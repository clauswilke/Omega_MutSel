# SJS.
# Code specifically for generating data to demonstrate that trend holds across gc content levels.
# Parameters: fix kappa to 1.0, beta to 2.0, and seqlength to 1e6 in the qsub file.
# Conduct 100 reps for each of the following: 0.2-0.3, 0.3-0.4, 0.4-0.5, 0.5-0.6, 0.6-0.7, 0.7-0.8

import sys
from functions_simandinf import *
import random as rn
sys.path.append('src/')
from misc import *
from newick import *
from stateFreqs import *
from matrixBuilder import *
from evolver import *


# Input parameters and global stuff
if (len(sys.argv) != 7):
    print "\n\nUsage: python run_gc.py <rep> <mu> <kappa> <bl> <seqlen> <gc_start> \n."
    sys.exit()
rep = sys.argv[1]
mu = float(sys.argv[2])
kappa = float(sys.argv[3])
bl = sys.argv[4]
seqlength = int(sys.argv[5])
gc_raw = sys.argv[6]

gc_start = float(gc_raw)*0.1
gc_end = gc_start + 0.1

# Set up beta for getting amino acid frequencies
beta = rn.uniform(0.5,3.5)

# Write tree given bl specifications
treefile = "tree.tre"
treef = open(treefile, 'w')
treef.write("(t1:" + str(bl) + ", t2:" + str(bl) + ");")
treef.close()

# set up output sequence and parameter files
freqfile = "codonFreqs" + str(rep) + ".txt"
outfile = "params" + str(rep) + "_" + gc_raw + ".txt"
seqfile = "seqs" + str(rep) + ".fasta"


# Simulate  
print "simulating"
f, num_pref_aa, gc_content = setFreqs(freqfile, beta, gc_start, gc_end) # last 2 args are gc min, gc max
simulate(f, seqfile, treefile, mu, kappa, seqlength, None) # omega (well, beta, as in dN) is last argument. when None, sim via mutsel
    
# Derive omega
mu_dict = {'AT':mu, 'AC':mu, 'AG':mu, 'CG':mu, 'CT':mu, 'GT':mu}
mu_dict['AG'] = mu_dict['AG'] * kappa
mu_dict['CT'] = mu_dict['CT'] * kappa
derived_w = deriveOmega(f, mu_dict)

    
# HyPhy omega
print "ML"
gy94_w = float(runpaml(seqfile))
    
# Calculate relative error from derived omega. Not using abs() since plot nicer that way.
err = ( derived_w - gy94_w )/derived_w

# Save
outf = open(outfile, 'w')
outf.write(rep + '\t' + str(seqlength) + '\t' + str(bl) + '\t' + str(mu) + '\t' + str(kappa) + '\t' + str(beta) + '\t' + str(gc_content) + '\t' + str(derived_w) + '\t' + str(gy94_w) + '\t' + str(err) + '\n')
outf.close()

