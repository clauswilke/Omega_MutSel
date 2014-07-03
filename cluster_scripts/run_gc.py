# SJS.
# Code specifically for generating data to demonstrate that trend holds across gc content levels.
# Parameters: fix kappa to 1.0, beta to 2.0, and seqlength to 1e6 in the qsub file.
# Conduct 100 reps for each of the following: 0.2-0.3, 0.3-0.4, 0.4-0.5, 0.5-0.6, 0.6-0.7, 0.7-0.8

import sys
from functions_simandinf import *

sys.path.append('src/')
from misc import *
from newick import *
from stateFreqs import *
from matrixBuilder import *
from evolver import *


# Input parameters and global stuff
if (len(sys.argv) != 9):
    print "\n\nUsage: python run_convergence.py <rep> <cpu> <beta> <mu> <kappa> <bl> <seqlen> <gc_start> \n."
    sys.exit()
rep = sys.argv[1]
cpu = sys.argv[2]
beta = float(sys.argv[3])
mu = float(sys.argv[4])
kappa = float(sys.argv[5])
bl = sys.argv[6]
seqlength = int(sys.argv[7])
gc_start = float(sys.argv[8])*0.1
gc_end = gc_start + 0.05

# Write tree given bl specifications
treefile = "tree.tre"
treef = open(treefile, 'w')
treef.write("(t1:" + str(bl) + ", t2:" + str(bl) + ");")
treef.close()

# set up output sequence and parameter files
freqfile = "codonFreqs" + str(rep) + ".txt"
outfile = "params" + str(rep) + ".txt"
seqfile = "seqs" + str(rep) + ".fasta"


# Simulate  
print "simulating"
f, num_pref_aa, gc_content = setFreqs(freqfile, beta, gc_start, gc_end) # last 2 args are gc min, gc max
simulate(f, seqfile, treefile, mu, kappa, seqlength, None) # omega is last argument. when None, sim via mutsel
    
# Derive omega
mu_dict = {'AT':mu, 'AC':mu, 'AG':mu, 'CG':mu, 'CT':mu, 'GT':mu}
mu_dict['AG'] = mu_dict['AG'] * kappa
mu_dict['CT'] = mu_dict['CT'] * kappa
derived_w = deriveOmega(f, mu_dict)

    
# HyPhy omega
print "ML"
gy94_w = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, kappa)
    
# Calculate relative error from derived omega. Not using abs() since plot nicer that way.
err = ( derived_w - gy94_w )/derived_w

# Save
outf = open(outfile, 'w')
outf.write(rep + '\t' + str(seqlength) + '\t' + str(bl) + '\t' + str(mu) + '\t' + str(kappa) + '\t' + str(num_pref_aa) + '\t' + str(round(gc_content, 6)) + '\t' + str(round(derived_w, 6)) + '\t' + str(round(gy94_w, 6)) + '\t' + str(round(err, 6)) + '\n')
outf.close()





