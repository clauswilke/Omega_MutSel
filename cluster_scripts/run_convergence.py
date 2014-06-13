# SJS
# Code to demonstrate convergence. Fix to mu = 1e-6, bl = 0.1. 
# Vary seqlength from 100->1000000.

import sys
import numpy as np
from functions_simandinf import *

sys.path.append('src/')
from misc import *
from newick import *
from stateFreqs import *
from matrixBuilder import *
from evolver import *

# Input parameters and global stuff
if (len(sys.argv) != 4):
    print "\n\nUsage: python run_convergence.py <rep> <cpu> <numaa>\n."
    sys.exit()
rep = sys.argv[1]
cpu = sys.argv[2]
numaa = int(sys.argv[3])
aadist = "exp"
mu = 1e-6
bl = 0.1
seqlengths = [1e2, 1e3, 1e4, 1e5, 1e6]

# Write tree given bl specifications
treefile = "tree.tre"
treef = open(treefile, 'w')
treef.write("(t1:" + str(bl) + ", t2:" + str(bl) + ");")
treef.close()

# Set frequencies 
f, mean_vol = setFreqs(aadist, numaa)
derived_w = None

outfile = "params"+str(rep)+".txt"
outf = open(outfile,'w') #outf.write('rep\tnumaa\taadist\tmu\tbl\tseqlength\tderived_w\tml_w\tmean_vol\n')

for seqlength in seqlengths:
    seqfile = "seqs"+str(rep)+"_"+str(seqlength)+".fasta"

    # Simulate
    print "simulating"
    simulate(f, seqfile, treefile, mu, int(seqlength), None) # omega is last arguement. when None, sim via mutsel

    # Derive omega
    if derived_w is None:
        print "deriving"
        derived_w, num_codons = deriveOmega(f)

    # HyPhy/PAML omega
    print "ML"
    ml_w = runhyphy("globalGY94.bf", "GY94_fixedkappa", seqfile, treefile, cpu)
    #ml_w = runpaml(seqfile, "0")

    # Save
    outf.write(rep + '\t' + str(numaa) + '\t' + str(aadist) + '\t' + str(mu) + '\t' + str(bl) + '\t' + str(seqlength) + '\t' + str(derived_w) + '\t' + str(ml_w) + '\t' + str(mean_vol) + '\n')

outf.close()


