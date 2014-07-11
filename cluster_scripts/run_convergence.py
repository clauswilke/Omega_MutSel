# SJS.
# Code specifically for generating data to demonstrate omega convergence.
# Parameters: fix kappa to 1.0 and beta to 2.5 in the qsub file.

import sys
from random import randint
from functions_simandinf import *

sys.path.append('src/')
from misc import *
from newick import *
from stateFreqs import *
from matrixBuilder import *
from evolver import *


# Input parameters and global stuff
if (len(sys.argv) != 6):
    print "\n\nUsage: python run_convergence.py <rep> <mu> <kappa> <bl> <expon> \n."
    sys.exit()
rep = sys.argv[1]
mu = float(sys.argv[2])
kappa = float(sys.argv[3])
bl = sys.argv[4]
expon = int(sys.argv[5])

# Set up beta for getting amino acid frequencies
beta = rn.uniform(0.5,3.5)

# Random sequence length between 5e2 and 1e6
if expon == 2:
    times = randint(5,10)
else:
    times = randint(1,10)
seqlength = int( times * 10**expon )

# Write tree given bl specifications
treefile = "tree.tre"
treef = open(treefile, 'w')
treef.write("(t1:" + str(bl) + ", t2:" + str(bl) + ");")
treef.close()

# set up output sequence and parameter files
freqfile = "codonFreqs" + str(rep) + "_" + str(seqlength) + ".txt"
seqfile = "seqs" + str(rep) + "_" + str(seqlength) + ".fasta"
outfile = "params" + str(rep) + "_" + str(seqlength) + ".txt"


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
gy94_w = float(runpaml(seqfile))

# Calculate relative error from derived omega. Not using abs() since plot nicer that way.
err = ( derived_w - gy94_w )/derived_w

# Save
outf = open(outfile,'w')
outf.write(rep + '\t' + str(seqlength) + '\t' + str(bl) + '\t' + str(mu) + '\t' + str(kappa) + '\t' + str(beta) + '\t' + str(gc_content) + '\t' + str(derived_w) + '\t' + str(gy94_w) + '\t' + str(err) + '\n')
outf.close()





