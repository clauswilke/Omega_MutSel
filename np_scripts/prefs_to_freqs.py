import sys
import subprocess
import shutil
import os 
import numpy as np
from scipy import linalg

# Append path for code to calc GC content
initial_path = os.environ['HOME']
if initial_path == '/Users/sjspielman':
    initial_path += '/Research'
sys.path.append(initial_path + '/MutSel/Simulator/src')
from stateFreqs import *
codons=["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]



usage_error = "\n\n Usage: python prefs_to_freqs.py <mu_scheme>.\n mu_scheme can either be np, yeast, or polio."

assert(len(sys.argv) == 2), usage_error

mu_type = sys.argv[1] # either np, yeast, or polio
if mu_type == 'np':
    #### Mutation rates from Bloom (2014). "An Experimentally Determined Evolutionary Model Dramatically Improves Phylogenetic Fit." MBE. #### 
    mudict = {'AG':2.4e-5, 'TC':2.4e-5, 'GA':2.3e-5, 'CT':2.3e-5, 'AC':9.0e-6, 'TG':9.0e-6, 'CA':9.4e-6, 'GT':9.4e-6, 'AT':3.0e-6, 'TA':3.0e-6, 'GC':1.9e-6, 'CG':1.9e-6}
      
elif mu_type == 'yeast':
    #### Mutation rates from Zhu et al. (2014). "Precise estimates of mutation rate and spectrum in yeast." PNAS. ####
    mu = 1.67e-10 # this is the mean per generation per nucleotide mutation rate. 
    mudict = {'AG':0.144/2*mu, 'TC':0.144/2*mu, 'GA':0.349/2*mu, 'CT':0.349/2*mu, 'AC':0.11/2*mu, 'TG':0.11/2*mu, 'CA':0.182/2*mu, 'GT':0.182/2*mu, 'AT':0.063/2*mu, 'TA':0.063/2*mu, 'GC':0.152/2*mu, 'CG':0.152/2*mu}

elif mu_type == 'polio':
    #### Mutation rates from Acevedo, Brodskey, and Andino (2014). "Mutational and fitness landscapes of an RNA virus revealed through population sequencing." Nature. ####
    mudict = {'AG':2.495e-5, 'TC':6.886e-05, 'GA':1.259e-04, 'CT':2.602e-04, 'AC':1.721e-06, 'TG':1.177e-06, 'CA':9.072e-06, 'GT':1.472e-05, 'AT':3.812e-06, 'TA':3.981e-06, 'GC':6.301e-06, 'CG':1.633e-06}
else:
    raise AssertionError(usage_error)
    

# File names
data_dir      = "../experimental_data/"
cf_outfile    = data_dir + mu_type + "_codon_eqfreqs.txt"
raw_batchfile = "globalDNDS_raw_exp.bf"
batch_outfile = '../SelectionInference/hyphy/globalDNDS_' + mu_type + '.bf'


# np_prefs are those taken from Bloom 2014 paper (same as mu's above!). The np_prefs are directly from the paper's Supplementary_file_1.xls and refer to equilbrium amino acid propenisties. The best interpretation of these experimental propensities is metropolis.
np_prefs = np.loadtxt(data_dir + 'nucleoprotein_amino_preferences.txt')
nsites = len(np_prefs)

amino_acids  = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
codon_dict   = {"AAA":"K", "AAC":"N", "AAG":"K", "AAT":"N", "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T", "AGA":"R", "AGC":"S", "AGG":"R", "AGT":"S", "ATA":"I", "ATC":"I", "ATG":"M", "ATT":"I", "CAA":"Q", "CAC":"H", "CAG":"Q", "CAT":"H", "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P", "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R", "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L", "GAA":"E", "GAC":"D", "GAG":"E", "GAT":"D", "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A", "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G", "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V", "TAC":"Y", "TAT":"Y", "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S", "TGC":"C", "TGG":"W", "TGT":"C", "TTA":"L", "TTC":"F", "TTG":"L", "TTT":"F"}
codons       = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
nucleotides  = ["A", "C", "G", "T"]


def codon_to_posnuc(codon_freqs):
    ''' Calculate positional nucleotide frequencies from codon frequencies. '''
    pos_nuc_freqs = np.zeros([3,4])
    for i in range(3):
        count = 0
        for codon in codons:
            pos_nuc_freqs[i][nucleotides.index(codon[i])] += codon_freqs[count]
            count += 1
        print np.sum(pos_nuc_freqs[i])
        assert( abs(np.sum(pos_nuc_freqs[i]) - 1.) < 1e-8 ), "bad positional nucleotide frequencies"
    return pos_nuc_freqs
    

def get_nuc_diff(source, target):
    diff = ''
    for i in range(3):
        if source[i] != target[i]: 
            diff += source[i]+target[i]
    return diff


def build_matrix_nuc(nuc_freqs):
    kappa = 4.0
    omega = 0.5

    matrix = np.zeros([61,61])
    
    # off diagonal
    for x in range(61):
        source = codons[x]
        for y in range(61):
            target = codons[y]
            diff = get_nuc_diff(source, target)
            
            if len(diff) == 2:
                matrix[x][y] = nuc_freqs[diff[1]]
                if diff == 'AG' or diff == 'GA' or diff == 'CT' or diff == 'TC':
                    matrix[x][y] *= kappa
                if codon_dict[source] != codon_dict[target]:
                    matrix[x][y] *= omega
    # diagonal
    for i in range(61):
        matrix[i][i] = -1. * np.sum(matrix[i]) 
        assert( -1e-10 < np.sum(matrix[i]) < 1e-10 ), "diagonal fail"
    return matrix


def build_matrix(amino_prop_dict, mu_dict):
    ''' metropolis only, as this matrix definition more suits experimental propensities according to Bloom 2014. '''
    matrix = np.ones([61,61])
    
    # off-diagonal entries
    for x in range(61):
        source = codons[x]
        fx = amino_prop_dict[codon_dict[source]]
        for y in range(61):
            target = codons[y]
            fy = amino_prop_dict[codon_dict[target]]
            diff = get_nuc_diff(source, target)
            if len(diff)==2:
                if fx < fy:
                    matrix[x][y] = fy / fx  
                matrix[x][y] *= mu_dict[diff]
            else:
                matrix[x][y] = 0.      
    # diagonal entries
    for i in range(61):
        matrix[i][i] = -1. * np.sum(matrix[i]) 
        assert( -1e-10 < np.sum(matrix[i]) < 1e-10 ), "diagonal fail"
    return matrix



def get_eq_from_eig(m):   
    ''' get the equilibrium frequencies from the matrix. the eq freqs are the left eigenvector corresponding to eigenvalue of 0. 
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
            print abs(forward - backward)
            assert(-1e-5 < abs(forward - backward) < 1e-5), "Detailed balance violated." 
    return max_v


def get_eq_freqs(amino_prefs, mu_dict):
    amino_prefs_dict = dict(zip(amino_acids, amino_prefs)) 
    m = build_matrix(amino_prefs_dict, mu_dict)
    cf = get_eq_from_eig(m) 
    assert( -1e-10 < abs(np.sum(cf)) - 1. < 1e-10 ), "codon frequencies do not sum to 1" 
    return cf


def create_temp_alignment(codon_freqs, seqfile):
    ''' We just need to quick set up sequences so that HyPhy can read in an alignment to compute F3x4, CF3x4.
        No actual simulation needed, just a large dataset with codons in proportion to global frequencies.
    '''
    size = 1e6
    seq = ''
    for i in range(61):
        seq += codons[i] * int(round(codon_freqs[i] * size)) 
    sfile = open(seqfile, 'w')
    sfile.write('>taxon\n'+seq)
    sfile.close()

def array_to_hyphy_freq(f):
    ''' Convert a python/numpy list/array of codon frequencies to a hyphy frequency string to directly write into a batchfile. '''
    hyphy_f = "{"
    for freq in f:
        hyphy_f += "{" + str(freq) + "},"
    hyphy_f = hyphy_f[:-1]
    hyphy_f += "};"
    return hyphy_f
    
    
def hyout_to_hyin(type):
    ''' convert hyphy frequency output a string to give back to hyphy. hackish but ok...'''
    f = open(type + '.txt', 'r')
    raw = f.read()
    f.close()
    clean = raw.replace('{','').replace('}','')
    return array_to_hyphy_freq( np.fromstring(clean, sep='\n')  )
    

def create_batchfile(basefile, outfile, f61, fnull, f3x4, cf3x4):
    ''' sed in the frequency specifications to create an output batchfile from the base/raw batchfile framework.'''
    cp_batch = subprocess.call("cp " + basefile + " " + outfile, shell=True)
    assert(cp_batch == 0), "couldn't copy batchfile"
    shutil.copy(basefile, outfile)
    flist = ['f61', 'fnull', 'f3x4', 'cf3x4']
    sedlist = ['F61', 'NULL', 'F3x4', 'CF3x4']
    for i in range(4):
        hyf = eval(flist[i])
        setuphyphyf = "sed -i 's/INSERT"+sedlist[i]+"/"+hyf+"/g' " + outfile
        setupf = subprocess.call(setuphyphyf, shell = True)
        assert(setupf == 0), "couldn't properly add in frequencies"
    

    
    
    

def main():
    # First, we determine the equilibrium frequencies of the system on a per site basis.
    # Second, we find the global frequencies. We additionally use HyPhy to find the global F3x4 and CF3x4 frequencies. These will be manually hard-coded into the hyphy files in SelectionInference/hyphy.
    # Third, we determine the "null" equilibrium frequencies. These are the codon frequencies which would be expected in the ABSENCE of seletion.
    # Finally, we set up the hyphy batch file which makes use of much of these frequencies.
    
    # Site-wise equilibrium frequencies
    print "Determine and save site-wise equilibrium frequencies"
    final_codon_freqs = np.zeros([nsites, 61])
    for i in range(nsites):
        print i
        final_codon_freqs[i] = get_eq_freqs(np_prefs[i], mudict)
    np.savetxt(cf_outfile, final_codon_freqs)
    
    # Determine global f61, f3x4, cf3x4
    print "Determining global F61"
    global_freqs = np.mean(final_codon_freqs, axis=0)
    print "Creating temporary alignment for hyphy frequency estimation"
    create_temp_alignment(global_freqs, "temp.fasta")
    print "Determining global F3x4 and CF3x4 frequencies"
    run_hyphy = subprocess.call("HYPHYMP global_freqs.bf", shell = True)
    assert(run_hyphy == 0), "Hyphy fail"

    # Determine freqs in absence of selection. Additionally print out their GC contents. This will reveal the extent of compositional bias generated by these experimental mutation rates.
    print "Determining null frequencies"
    null_freqs = get_eq_freqs(np.ones(20) * 0.05 , mudict)
    fobj = UserFreqs(by = 'codon', freqs = dict(zip(codons, null_freqs)))
    nuc_freq_dict = dict(zip(['A', 'C', 'G', 'T'], fobj.calcFreqs(type = 'nuc')))    

##  these lines served to confirm that a matrix built using nucleotide frequencies instead of codon frequencies is indeed reversible.
#   matrix = build_matrix_nuc(nuc_freq_dict)
#   null_nuc_eqfreqs = get_eq_from_eig(matrix)


    # Create the hyphy batch file which will use much of this information
    # Note that first we must make the strings to put into this batchfile
    print "Creating HyPhy batchfiles"
    f61   = array_to_hyphy_freq(global_freqs)
    fnull  = array_to_hyphy_freq(null_freqs)
    f3x4  = hyout_to_hyin("f3x4")
    cf3x4 = hyout_to_hyin("cf3x4")
    create_batchfile(raw_batchfile, batch_outfile, f61, fnull, f3x4, cf3x4)
   
    print "Here's a dictionary of final nucleotide compositions for",mu_type,":"
    print nuc_freq

    

   
# Run and then cleanup.
main() 
os.remove("f3x4.txt")
os.remove("cf3x4.txt")
os.remove("temp.fasta")
os.remove("messages.log")
try:
    os.remove("errors.log") # this file is not always created.
except:
    pass







