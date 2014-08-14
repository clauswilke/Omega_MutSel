import sys
import subprocess
import shutil
import os 
import numpy as np
from scipy import linalg

usage_error = "\n\n Usage: python prefs_to_freqs.py <mu_scheme>.\n mu_scheme can either be np or yeast."

assert(len(sys.argv) == 2), usage_error

mu_type = sys.argv[1] # either np or yeast
if mu_type == 'np':
    ### Mutation rates from Bloom (2014). "An Experimentally Determined Evolutionary Model Dramatically Improves Phylogenetic Fit." MBE. #### 
    mudict = {'AG':2.4e-5, 'TC':2.4e-5, 'GA':2.3e-5, 'CT':2.3e-5, 'AC':9.0e-6, 'TG':9.0e-6, 'CA':9.4e-6, 'GT':9.4e-6, 'AT':3.0e-6, 'TA':3.0e-6, 'GC':1.9e-6, 'CG':1.9e-6}
      
elif mu_type == 'yeast':
    ### Mutation rates from Zhu et al. (2014). "Precise estimates of mutation rate and spectrum in yeast." PNAS. ####
    mu = 1.67e-10 # this is the mean per generation per nucleotide mutation rate. 
    mudict = {'AG':0.144/2*mu, 'TC':0.144/2*mu, 'GA':0.349/2*mu, 'CT':0.349/2*mu, 'AC':0.11/2*mu, 'TG':0.11/2*mu, 'CA':0.182/2*mu, 'GT':0.182/2*mu, 'AT':0.063/2*mu, 'TA':0.063/2*mu, 'GC':0.152/2*mu, 'CG':0.152/2*mu}

else:
    raise AssertionError(usage_error)
    

# Outfile names
path = 'data/'
cf_outfile      = path + mu_type + '_site_eqfreqs.txt'
mean_cf_outfile = path + mu_type + '_mean_eqfreqs.txt'
null_cf_outfile = path + mu_type + '_null_eqfreqs.txt'
f3x4_outfile    = path + mu_type + '_f3x4.txt' 
cf3x4_outfile   = path + mu_type + '_cf3x4.txt' 

# np_prefs are those taken from Bloom 2014 paper (same as mu's above!). The np_prefs are directly from the paper's Supplementary_file_1.xls and refer to equilbrium amino acid propenisties. The best interpretation of these experimental propensities is metropolis.
np_prefs = np.loadtxt(path + 'np_prefs.txt')
nsites = len(np_prefs)

amino_acids  = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
codon_dict   = {"AAA":"K", "AAC":"N", "AAG":"K", "AAT":"N", "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T", "AGA":"R", "AGC":"S", "AGG":"R", "AGT":"S", "ATA":"I", "ATC":"I", "ATG":"M", "ATT":"I", "CAA":"Q", "CAC":"H", "CAG":"Q", "CAT":"H", "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P", "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R", "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L", "GAA":"E", "GAC":"D", "GAG":"E", "GAT":"D", "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A", "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G", "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V", "TAC":"Y", "TAT":"Y", "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S", "TGC":"C", "TGG":"W", "TGT":"C", "TTA":"L", "TTC":"F", "TTG":"L", "TTT":"F"}
codons       = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]




def get_nuc_diff(source, target):
    diff = ''
    for i in range(3):
        if source[i] != target[i]: 
            diff += source[i]+target[i]
    return diff



def build_matrix(amino_prop_dict, mu_dict):
    ''' metropolis only, as this matrix definition more suits experimental propensities. '''
    matrix = np.zeros([61,61])
    
    # non-diagonal entries
    for x in range(61):
        source = codons[x]
        fx = amino_prop_dict[codon_dict[source]]
        for y in range(61):
            target = codons[y]
            fy = amino_prop_dict[codon_dict[target]]
            diff = get_nuc_diff(source, target)
            if len(diff) == 2 and source != target:
                if fy >= fx or codon_dict[source] == codon_dict[target]:
                    matrix[x][y] = mu_dict[diff] 
                else:
                    matrix[x][y] = fy / fx * mu_dict[diff]        
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
    return max_v


def get_eq_freqs(amino_prefs, mu_dict):
    amino_prefs_dict = dict(zip(amino_acids, amino_prefs)) 
    m = build_matrix(amino_prefs_dict, mu_dict)
    cf = get_eq_from_eig(m) 
    assert( -1e-8 < abs(np.sum(cf)) - 1. < 1e-8 ), "codon frequencies do not sum to 1" 
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


def main():
    # First, we determine the equilibrium frequencies of the system on a per site basis.
    # Second, we find the global frequencies. We additionally use HyPhy to find the global F3x4 and CF3x4 frequencies. These will be manually hard-coded into the hyphy files in SelectionInference/hyphy.
    # Third, we determine the "null" equilibrium frequencies. These are the codon frequencies which would be expected in the ABSENCE of seletion.
    
    # Site-wise equilibrium frequencies
    print "Determine site-wise equilibrium frequencies"
    final_codon_freqs = np.zeros([nsites, 61])
    for i in range(nsites):
        print i
        final_codon_freqs[i] = get_eq_freqs(np_prefs[i], mudict)
    
    # Determine global equilibrium frequencies. Also call HyPhy to retrieve F3x4, CF3x4 global frequencies.
    print "Determining global"
    global_freqs = np.mean(final_codon_freqs, axis=0)
    print "Simulating"
    create_temp_alignment(global_freqs, "temp.fasta")
    print "Determining F3x4, CF3x4 frequencies"
    run_hyphy = subprocess.call("HYPHYMP global_freqs.bf", shell = True)
    assert(run_hyphy == 0), "Hyphy fail"
    os.remove("temp.fasta")
    os.remove("messages.log")
    
    
    # Determine freqs in absence of selection.
    print "Determining null frequencies"
    null_freqs = get_eq_freqs(np.ones(20) * 0.05 , mudict)
   
    
    # Save all to files
    print "Saving to files"
    shutil.move("f3x4.txt", f3x4_outfile)
    shutil.move("cf3x4.txt", cf3x4_outfile)
    np.savetxt(cf_outfile, final_codon_freqs)
    np.savetxt(null_cf_outfile, null_freqs)
    np.savetxt(mean_cf_outfile, np.mean(final_codon_freqs, axis=0))
   
   
main() 







