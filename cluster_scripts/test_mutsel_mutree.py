## SJS.
## Script to see behavior when simulate according to mutsel model with various numaa, tree bl, and mu
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
length = 10000
branches =  [0.01, 0.05, 0.1, 0.5]
mutations = [0.0005, 0.001, 0.005, 0.01]

outf = open("tempout.txt", 'w')
for treebl in branches:

    # Write tree given bl specifications
    treefile = "tree.tre"
    treef = open(treefile, 'w')
    treef.write("(t1:" + str(treebl) + ", t2:" + str(treebl) + ");")
    treef.close()
    
    for mu in mutations:
        seqfile = 'seqs'+str(rep)+'_'+str(branches.index(treebl))+'_'+str(mutations.index(mu)) + '.fasta'

        # Simulate
        print "simulating"
        f, aminos_used = setFreqs("user", numaa)
        simulate(f, seqfile, treefile, mu, length) # do not specify last argument omega -> sim via mutsel
        
        # Use math to derive an omega for the simulated sequences. Also returns the number of codons theoretically sampled.
        print "deriving"
        derived_w, num_codons = deriveOmega(f)
    
        # Nei-Gojobori Method
        print "nei-gojobori"
        nei_w, ns_mut, s_mut = run_neigojo(seqfile)

        # HyPhy
        print "ML"
        hyphy_w = runhyphy("globalGY94.bf", "GY94_fixedkappa", seqfile, treefile, cpu, f)
        outf.write(rep + '\t' + str(numaa) + '\t' + str(treebl) + '\t' + str(mu) + '\t' + str(derived_w) + '\t' + str(nei_w) + '\t' + str(hyphy_w) + '\t' + str(ns_mut) + '\t' + str(s_mut) + '\t' + aminos_used + '\n')

outf.close()

