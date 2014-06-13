# SJS
# Generic code for simulating and deriving omega via math, hyphy ML.
# Be sure to cp src/ directory (simulator), hyphy files, and the functions_simandinf.py script into working directory
# 6/12/14 -> best combo for 100000 sites where kappa=1 is mu=1e-6,bl=0.001

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
if (len(sys.argv) != 8):
    print "\n\nUsage: python run_siminf.py <rep> <cpu> <numaa> <aadist> <mu> <bl> <seqlength>\n."
    sys.exit()
rep = sys.argv[1]
cpu = sys.argv[2]
numaa = int(sys.argv[3])
aadist = sys.argv[4] # Either "exp", "equal", "random"
mu = float(sys.argv[5])
bl = sys.argv[6]
seqlength = int(sys.argv[7])

outfile = "params"+str(rep)+".txt"
outf = open(outfile,'w')
#outf.write('rep\tnumaa\taadist\tmu\tbl\tseqlength\tderived_w\tml_w\n')


seqfile = "seq"+str(rep)+".fasta"

# Write tree given bl specifications
treefile = "tree.tre"
treef = open(treefile, 'w')
treef.write("(t1:" + str(bl) + ", t2:" + str(bl) + ");")
treef.close()

# Simulate
print "simulating"
f, aminos_used = setFreqs(aadist, numaa)
simulate(f, seqfile, treefile, mu, seqlength, None) # omega is last arguement. when None, sim via mutsel

# Derive omega
print "deriving"
derived_w, num_codons = deriveOmega(f)

# HyPhy/PAML omega
print "ML"
ml_w = runhyphy("globalGY94.bf", "GY94_fixedkappa", seqfile, treefile, cpu)
#ml_w = runpaml(seqfile, "0")

# Save
outf.write(rep + '\t' + str(numaa) + '\t' + str(aadist) + '\t' + str(mu) + '\t' + str(bl) + '\t' + str(seqlength) + '\t' + str(derived_w) + '\t' + str(ml_w) + '\n')
outf.close()


