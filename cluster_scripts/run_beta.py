# SJS.
# Code specifically for generating data to demonstrate that trend holds across amino acid constraint, represented by beta, levels.
# Run 500 reps. k=1, bl=0.005, mu=1e-6, seqlen=1e6

import sys
from functions_simandinf import *
from random import uniform
sys.path.append('src/')
from misc import *
from newick import *
from stateFreqs import *
from matrixBuilder import *
from evolver import *


# Input parameters and global stuff
if (len(sys.argv) != 6):
    print "\n\nUsage: python run_beta.py <rep> <mu> <kappa> <bl> <seqlen> \n."
    sys.exit()
rep = sys.argv[1]
mu = float(sys.argv[2])
kappa = float(sys.argv[3])
bl = sys.argv[4]
seqlength = int(sys.argv[5])

beta = rn.uniform(0.5,3.5)

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
f, num_pref_aa, gc_content = setFreqs(freqfile, beta)
simulate(f, seqfile, treefile, mu, kappa, seqlength, None) # omega is last argument. when None, sim via mutsel
entropy = calcCodonEntropy(f)    

# Derive omega
print "deriving"
mu_dict = {'AT':mu, 'AC':mu, 'AG':mu, 'CG':mu, 'CT':mu, 'GT':mu}
mu_dict['AG'] = mu_dict['AG'] * kappa
mu_dict['CT'] = mu_dict['CT'] * kappa
derivedw = deriveOmega(f, mu_dict)

# HyPhy omega
print "ML"
mlw0 = float(runpaml(seqfile, codonFreq = "0", initw = 0.4))
mlw2   = float(runpaml(seqfile, codonFreq = "2", initw = 0.4))
mlw3  = float(runpaml(seqfile, codonFreq = "3", initw = 0.4))

# Calculate relative error from derived omega.
err0 = abs( derivedw - mlw0 )/derivedw
err2 = abs( derivedw - mlw2 )/derivedw
err3 = abs( derivedw - mlw3 )/derivedw

# Save
outf = open(outfile, 'w')
outf.write(rep + '\t' + str(seqlength) + '\t' + str(bl) + '\t' + str(mu) + '\t' + str(kappa) + '\t' + str(beta) + '\t' + str(entropy) + '\t' + str(gc_content) + '\t' + str(derivedw) + '\t' + str(mlw0) + '\t' + str(mlw2) + '\t' + str(mlw3) + '\t' + str(err0) + '\t' + str(err2) + '\t' + str(err3) + '\n')
outf.close()





