# SJS. stephanie.spielman@gmail.com
# Code for simulating NYP sequences. 

######## Input parameters ########
import sys
if (len(sys.argv) != 4):
    print "\n\nUsage: python run_np.py <rep> <simdir> <dataset> \n."
    sys.exit()
rep = sys.argv[1]         # which rep we're on, for saving files. needs to run from 1-498, since 498 sites.
simdir = sys.argv[2]      # directory of simulation library
dataset = sys.argv[3]     # either np, yeast, or polio. determines the mutation scheme and eq freqs
from functions_omega_mutsel import *

seqfile   = "seqs"+str(rep)+".fasta"
write_treefile(treefile)
eq_codon_freqs = np.loadtxt(dataset + "_codon_eqfreqs.txt")
codon_freqs = eq_codon_freqs[ int(rep) - 1 ]
simulate(codon_freqs, seqfile, 'tree.tre', set_mu_dict(dataset), 500000, simdir)

