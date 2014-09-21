import sys
import subprocess
import shutil
import os 
import re
import numpy as np
from scipy import linalg
amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
codon_dict  = {"AAA":"K", "AAC":"N", "AAG":"K", "AAT":"N", "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T", "AGA":"R", "AGC":"S", "AGG":"R", "AGT":"S", "ATA":"I", "ATC":"I", "ATG":"M", "ATT":"I", "CAA":"Q", "CAC":"H", "CAG":"Q", "CAT":"H", "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P", "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R", "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L", "GAA":"E", "GAC":"D", "GAG":"E", "GAT":"D", "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A", "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G", "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V", "TAC":"Y", "TAT":"Y", "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S", "TGC":"C", "TGG":"W", "TGT":"C", "TTA":"L", "TTC":"F", "TTG":"L", "TTT":"F"}
codons      = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
nucindex    = {'A':0, 'C':1, 'G':2, 'T':3}
purines     = ["A", "G"]
pyrims      = ["C", "T"]


def parse_input(arguments):
    usage_error = "\n\n Usage: python prefs_to_freqs.py <mu_scheme> <compute_preference_frequencies>.\n mu_scheme can either be np, yeast, or polio. The second argument is optional (provide 0 or 1 for false, true though.)"
    assert(len(arguments) == 2 or len(arguments) == 3), usage_error

    mu_type = arguments[1] # either np, yeast, or polio
    #### Mutation rates from Bloom (2014). "An Experimentally Determined Evolutionary Model Dramatically Improves Phylogenetic Fit." MBE. #### 
    if mu_type == 'np':
        mudict = {'AG':2.4e-5, 'TC':2.4e-5, 'GA':2.3e-5, 'CT':2.3e-5, 'AC':9.0e-6, 'TG':9.0e-6, 'CA':9.4e-6, 'GT':9.4e-6, 'AT':3.0e-6, 'TA':3.0e-6, 'GC':1.9e-6, 'CG':1.9e-6}
    
    #### Mutation rates from Zhu et al. (2014). "Precise estimates of mutation rate and spectrum in yeast." PNAS. ####
    elif mu_type == 'yeast':
        mu = 1.67e-10 # this is the mean per generation per nucleotide mutation rate. 
        mudict = {'AG':0.144/2*mu, 'TC':0.144/2*mu, 'GA':0.349/2*mu, 'CT':0.349/2*mu, 'AC':0.11/2*mu, 'TG':0.11/2*mu, 'CA':0.182/2*mu, 'GT':0.182/2*mu, 'AT':0.063/2*mu, 'TA':0.063/2*mu, 'GC':0.152/2*mu, 'CG':0.152/2*mu}

    #### Mutation rates from Acevedo, Brodskey, and Andino (2014). "Mutational and fitness landscapes of an RNA virus revealed through population sequencing." Nature. ####
    elif mu_type == 'polio':
        mudict = {'AG':2.495e-5, 'TC':6.886e-05, 'GA':1.259e-04, 'CT':2.602e-04, 'AC':1.721e-06, 'TG':1.177e-06, 'CA':9.072e-06, 'GT':1.472e-05, 'AT':3.812e-06, 'TA':3.981e-06, 'GC':6.301e-06, 'CG':1.633e-06}
    else:
        raise AssertionError(usage_error)
    
    try:
        compute_freqs_from_prefs = bool(int(arguments[2]))
    except:
        compute_freqs_from_prefs = True
    
    return mu_type, mudict, compute_freqs_from_prefs


def get_nuc_diff(source, target, grab_position = False):
    diff = ''
    position = 5
    for i in range(3):
        if source[i] != target[i]: 
            diff += source[i]+target[i]
            position = i
    if grab_position:
        return diff, position
    else:
        return diff

def build_metropolis_matrix(amino_prop_dict, mu_dict):
    ''' metropolis only, as this matrix definition more suits experimental propensities according to Bloom 2014. '''
    matrix = np.zeros([61,61])
    
    # off-diagonal entries
    for x in range(61):
        x_codon = codons[x]
        x_aa = codon_dict[x_codon]
        fx = amino_prop_dict[x_aa]
        for y in range(61):
            y_codon = codons[y]
            y_aa = codon_dict[y_codon]
            fy = amino_prop_dict[y_aa]
            diff = get_nuc_diff(x_codon, y_codon)
            if len(diff)==2:
                if x_aa == y_aa or fy >= fx:
                    matrix[x][y] =  mu_dict[diff]
                else:
                    matrix[x][y] = fy / fx  * mu_dict[diff]    
    # diagonal entries
    for i in range(61):
        matrix[i][i] = -1. * np.sum(matrix[i]) 
        assert( -1e-10 < np.sum(matrix[i]) < 1e-10 ), "diagonal fail"
    return matrix



def get_eq_from_eig(m):   
    ''' get the equilibrium frequencies from the matrix. the eq freqs are the *left* eigenvector corresponding to eigenvalue of 0. 
        Code here is largely taken from Bloom. See here - https://github.com/jbloom/phyloExpCM/blob/master/src/submatrix.py, specifically in the fxn StationaryStates
    '''
    (w, v) = linalg.eig(m, left=True, right=False)
    max_i = 0
    max_w = w[max_i]
    for i in range(1, len(w)):
        if w[i] > max_w:
            max_w = w[i]
            max_i = i
    assert( abs(max_w) < 1e-10 ), "Maximum eigenvalue is not close to zero."
    max_v = v[:,max_i]
    max_v /= np.sum(max_v)
    max_v = max_v.real # these are the stationary frequencies
    assert( abs(np.sum(max_v) - 1.) < 1e-8), "Eigenvector of equilibrium frequencies doesn't sum to 1."

    # SOME SANITY CHECKS
    assert np.allclose(np.zeros(61), np.dot(max_v, m)) # should be true since eigenvalue of zero
    pi_inv = np.diag(1.0 / max_v)
    s = np.dot(m, pi_inv)
    assert np.allclose(m, np.dot(s, np.diag(max_v)), atol=1e-10, rtol=1e-5), "exchangeability and equilibrium does not recover matrix"
    # additional overkill check
    for i in range(61):
        pi_i = max_v[i]
        for j in range(61):
            pi_j = max_v[j]
            forward  = pi_i * m[i][j] 
            backward = pi_j * m[j][i]
            assert(abs(forward - backward) < 1e-6), "Detailed balance violated."  # note that we need to use high tolerance here because the propensities have very few digits so we encounter more FLOP problems.
    return max_v


def get_eq_freqs(amino_prefs, mu_dict):
    amino_prefs_dict = dict(zip(amino_acids, amino_prefs)) 
    m = build_metropolis_matrix(amino_prefs_dict, mu_dict)
    cf = get_eq_from_eig(m) 
    assert( abs(np.sum(cf) - 1.) < 1e-8 ), "codon frequencies do not sum to 1" 
    return cf


def codon_to_nuc(codon_freqs):
    ''' Calculate nucleotide and positional nucleotide frequencies from codon frequencies. '''
    pos_nuc_freqs = np.zeros([4,3])
    for i in range(3):
        count = 0
        for codon in codons:
            pos_nuc_freqs[nucindex[codon[i]]][i] += codon_freqs[count]
            count += 1
    assert( np.allclose(np.sum(pos_nuc_freqs, axis=0),np.ones(3))), "bad positional nucleotide frequencies"
    return np.mean(pos_nuc_freqs, axis=1), pos_nuc_freqs



def calc_f3x4_freqs(pos_nuc_freqs):
    ''' Compute F3x4 codon frequencies from positional nucleotide frequencies. '''
    f3x4 = np.ones(61)
    pi_stop = pos_nuc_freqs[nucindex['T']][0]*pos_nuc_freqs[nucindex['A']][1]*pos_nuc_freqs[nucindex['G']][2]  +  pos_nuc_freqs[nucindex['T']][0]*pos_nuc_freqs[nucindex['G']][1]*pos_nuc_freqs[nucindex['A']][2]  +  pos_nuc_freqs[nucindex['T']][0]*pos_nuc_freqs[nucindex['A']][1]*pos_nuc_freqs[nucindex['A']][2] 

    for i in range(61):
        codon = codons[i]
        for j in range(3):
            f3x4[i] *= pos_nuc_freqs[ nucindex[codon[j]] ][j]
    f3x4 /= (1. - pi_stop)
    assert( abs(np.sum(f3x4) - 1.) < 1e-8), "Could not properly calculate (C)F3x4 frequencies."
    return f3x4   




def calc_f1x4_freqs(nuc_freqs):
    ''' Compute F1x4 codon frequencies from nucleotide frequencies. '''
    f1x4 = np.ones(61)
    pi_stop = nuc_freqs[nucindex['T']]*nuc_freqs[nucindex['A']]*nuc_freqs[nucindex['G']] + nuc_freqs[nucindex['T']]*nuc_freqs[nucindex['G']]*nuc_freqs[nucindex['A']] + nuc_freqs[nucindex['T']]*nuc_freqs[nucindex['A']]*nuc_freqs[nucindex['A']]
    for i in range(61):
        codon = codons[i]
        for j in range(3):
            f1x4[i] *= nuc_freqs[nucindex[codon[j]]]
    f1x4 /= (1. - pi_stop)
    assert( abs(np.sum(f1x4) - 1.) < 1e-8), "Could not properly caluclate F1x4 frequencies."
    return f1x4
        

def is_TI(source, target):
    check1 = source in purines and target in purines
    check2 = source in pyrims and target in pyrims
    if check1 or check2:
        return True
    else:
        return False
        


def build_mg_matrices(nuc_freqs, pos_nuc_freqs, outfile):
    ''' Create MG94-style matrices (use target nucleotide frequencies).  ''' 

    matrix_mg1  = 'MG1 = {61, 61, \n' # MG94
    matrix_mg3  = 'MG3 = {61, 61, \n' # MG94, positional
    
    for i in range(61):
        source = codons[i]
        for j in range(61):
            target = codons[j]
                
            diff, x = get_nuc_diff(source, target, grab_position=True)
            if len(diff) == 2:
                assert(len(str(x)) == 1), "Problem with determining nucleotide difference between codons"
                target_index = nucindex[diff[1]]   
                mg1  = str(nuc_freqs[nucindex[diff[1]]]) 
                mg3  = str(pos_nuc_freqs[nucindex[diff[1]]][x])
  
                # Create string for matrix element
                element = '{' + str(i) + ',' + str(j) + ',t'  
                if is_TI(diff[0], diff[1]):
                    element += '*k'
                if codon_dict[source] != codon_dict[target]:
                    element += '*w' 
                matrix_mg1  += element + '*' + mg1 + '}\n'
                matrix_mg3  += element + '*' + mg3 + '}\n'


    # And save to file.
    with open(outfile, 'w') as outf:
        outf.write(matrix_mg1 + '};\n\n\n' + matrix_mg3 + '};\n\n\n')



def array_to_hyphy_freq(f):
    ''' Convert array of codon frequencies to a hyphy frequency string to directly write into a batchfile. '''
    hyphy_f = "{"
    for freq in f:
        hyphy_f += "{" + str(freq) + "},"
    hyphy_f = hyphy_f[:-1]
    hyphy_f += "};"
    return hyphy_f

def array_to_hyphy_posfreq(f):
    ''' Convert array of positional nucleotide frequencies to a hyphy frequency string to directly write into a batchfile. '''
    hyphy_f = "{"
    for row in f:
        hyphy_f += "{"
        for entry in row:
            hyphy_f += str(entry) + ","
        hyphy_f = hyphy_f[:-1] + "},"
    hyphy_f = hyphy_f[:-1] + "};"
    return hyphy_f

    
    

def create_batchfile(basefile, outfile, pos_freqs, f61, f1x4, f3x4):
    ''' sed in the frequency specifications to create an output batchfile from the base/raw batchfile framework.'''
    cp_batch = subprocess.call("cp " + basefile + " " + outfile, shell=True)
    assert(cp_batch == 0), "couldn't copy batchfile"
    shutil.copy(basefile, outfile)
    flist = ['pos_freqs', 'f61', 'f1x4', 'f3x4']
    for i in range( len(flist) ):
        hyf = eval(flist[i])
        insert = flist[i].upper()
        setuphyphyf = "sed -i 's/INSERT_"+insert+"/"+hyf+"/g' " + outfile
        setupf = subprocess.call(setuphyphyf, shell = True)
        assert(setupf == 0), "couldn't properly add in frequencies"

           
    

def main():
    # First, we determine the equilibrium frequencies of the system on a per site basis. As we use amino acid preference data, we assign all synonymous codons the same fitness.
    # Second, we find the F61, F1x4, F3x4 frequencies, as well as the MG1 and MG3 matrices. These use the average dataset frequencies, analogous to taking global alignment frequencies.
    # Third, we set up the hyphy batch file which makes use of these frequencies. Note that the Fnuc matrices are saved to a separate matrix file, and are not placed directly into the hyphy batchfiles.
    
    # Parse input arguments and set up input/outfile files accordingly
    dataset, mudict, compute = parse_input(sys.argv)    
    data_dir      = "../experimental_data/"
    codon_freqs_outfile    = data_dir + dataset + "_codon_eqfreqs.txt"
    raw_batchfile = "globalDNDS_raw_exp.bf"
    batch_outfile = '../../hyphy_files/globalDNDS_' + dataset + '.bf'
    mg_outfile  = '../../hyphy_files/MG_' + dataset + '.mdl'
    
    # Compute codon equilibrium frequencies from amino acid preference data, if compute is True. Else, load the files of eq freqs which have already been calculated.
    # np_prefs are those taken from Bloom 2014 MBE paper. The np_prefs are directly from the paper's Supplementary_file_1.xls and refer to equilbrium amino acid propenisties. The best interpretation of these experimental propensities is metropolis.
    if compute:
        np_prefs = np.loadtxt(data_dir + 'nucleoprotein_amino_preferences.txt')
        nsites = len(np_prefs)
    
        # Site-wise equilibrium frequencies
        print "Calculating and saving site-wise equilibrium frequencies"
        final_codon_freqs = np.zeros([nsites, 61])
        for i in range(nsites):
            final_codon_freqs[i] = get_eq_freqs(np_prefs[i], mudict)
        np.savetxt(codon_freqs_outfile, final_codon_freqs)
    else:
        print "Loading equilibrium frequency data"
        final_codon_freqs = np.loadtxt(codon_freqs_outfile)


    # Calculate frequency parameterizations
    print "Calculating frequency parameterizations F61, F1x4, F3x4."
    f61_freqs = np.mean(final_codon_freqs, axis=0)
    nuc_freqs, pos_nuc_freqs = codon_to_nuc(f61_freqs)
    f1x4_freqs = calc_f1x4_freqs(nuc_freqs)
    f3x4_freqs = calc_f3x4_freqs(pos_nuc_freqs)

    # Convert frequency arrays to hyphy strings
    f61 = array_to_hyphy_freq(f61_freqs)
    f1x4 = array_to_hyphy_freq(f1x4_freqs)
    f3x4 = array_to_hyphy_freq(f3x4_freqs)
    pos_freqs_hyf = array_to_hyphy_posfreq(pos_nuc_freqs)
 
 
    # Create the hyphy batchfile to include the frequencies calculated here. Note that we need to do this since no actual alignment exists which includes all protein positions, so cannot be read in by hyphy file.       
    create_batchfile(raw_batchfile, batch_outfile, pos_freqs_hyf, f61, f1x4, f3x4)
          
          
    # Use nucleotide and positional nucleotide frequencies to construct MG-style matrices
    print "Building and saving the MG-style matrices"
    build_mg_matrices(nuc_freqs, pos_nuc_freqs, mg_outfile)

    
main() 

