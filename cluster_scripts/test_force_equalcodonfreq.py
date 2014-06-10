## SJS.
## Script to see behavior when we force equal codon frequencies (1/61 for all, REGARDLESS of true dist). 
## Use PAML for this since easier!
## Can do equal or exponential amino acid distributions.


import os
import re
import sys
import subprocess
import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import *


sys.path.append('src/')
from misc import *
from newick import *
from stateFreqs import *
from matrixBuilder import *
from evolver import *
from mutation_counter import *
from site_counter import *
from predict_omega import *
from functions_simandinf import *


# Input parameters and global stuff
cpu = sys.argv[1]
rep = sys.argv[2]
numaa = int(sys.argv[3])
seqlength = 10000
treebl =  0.1
mu = 0.001

outf = open("tempout.txt", 'w')

# Write tree given bl specifications
treefile = "tree.tre"
treef = open(treefile, 'w')
treef.write("(t1:" + str(treebl) + ", t2:" + str(treebl) + ");")
treef.close()

seqfile = 'seqs'+str(rep)+'.fasta'

# Simulate
print "simulating"
f, aminos_used = setFreqs("user", numaa)
simulate(f, seqfile, treefile, mu, seqlength, None) # omega is last arguement. when None, sim via mutsel

# Use math to derive an omega for the simulated sequences. Also returns the number of codons theoretically sampled.
print "deriving"
derived_w, num_codons = deriveOmega(f)
    
# Nei-Gojobori Method
print "nei-gojobori"
nei_w, ns_mut, s_mut = run_neigojo(seqfile)

# HyPhy
print "ML"
paml_w = runpaml(seqfile, "1") # codonFreq=1 sets 1/61, regardless of reality.
outf.write(rep + '\t' + str(numaa) + '\t' + str(derived_w) + '\t' + str(nei_w) + '\t' + str(paml_w) + '\t' + aminos_used + '\n')

outf.close()

