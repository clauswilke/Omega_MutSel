# SJS.
# Code specifically for generating data to demonstrate omega convergence.
# Parameters: fix kappa to 1.0 and beta to 2.5 in the qsub file.

import sys
from random import uniform
from functions_simandinf import *

sys.path.append('src/')
from misc import *
from newick import *
from stateFreqs import *
from matrixBuilder import *
from evolver import *


# Input parameters and global stuff
if (len(sys.argv) != 7):
    print "\n\nUsage: python run_convergence.py <rep> <cpu> <beta> <mu> <kappa> <bl> \n."
    sys.exit()
rep = sys.argv[1]
cpu = sys.argv[2]
beta = float(sys.argv[3])
mu = float(sys.argv[4])
kappa = float(sys.argv[5])
bl = sys.argv[6]
seqlength = int(uniform(500,1e7)) # random sequence length between 1e2 and 1e7 

# Write tree given bl specifications
treefile = "tree.tre"
treef = open(treefile, 'w')
treef.write("(t1:" + str(bl) + ", t2:" + str(bl) + ");")
treef.close()

# set up output sequence and parameter files
freqfile = "codonFreqs" + str(rep)+".txt"
seqfile = "seqs"+str(rep)+".fasta"
outfile = "params"+str(rep)+".txt"
outf = open(outfile,'w')


# Now, simulate sequences and infer ML omegas
# Simulate
print "simulating"
f, num_pref_aa, gc_content = setFreqs(freqfile, beta, 0.0, 1.0) # last 2 args are gc min, gc max
simulate(f, seqfile, treefile, mu, kappa, seqlength, None) # omega is last argument. when None, sim via mutsel

# Derive
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
outf.write(rep + '\t' + str(seqlength) + '\t' + str(bl) + '\t' + str(mu) + '\t' + str(kappa) + '\t' + str(num_pref_aa) + '\t' + str(round(gc_content, 6)) + '\t' + str(round(derived_w, 6)) + '\t' + str(round(gy94_w, 6)) + '\t' + str(round(err, 6)) + '\n')
outf.close()





