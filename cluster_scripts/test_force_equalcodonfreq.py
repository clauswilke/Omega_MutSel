## SJS.
## Script to see behavior when we force equal codon frequencies (1/61 for all, REGARDLESS of true dist). 
## Use PAML for this since easier!
## Runs with both equal and exponential amino acid distributions.


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

freqType = {"user": "expAA", "equal": "equalAA"}

# Write tree given bl specifications
treefile = "tree.tre"
treef = open(treefile, 'w')
treef.write("(t1:" + str(treebl) + ", t2:" + str(treebl) + ");")
treef.close()



outf = open("tempout.txt", 'w')


for ftype in freqType.keys():

    seqfile = 'seqs_'+freqType[ftype]+str(rep)+'.fasta'

    # Simulate
    print "simulating"
    f, aminos_used = setFreqs(ftype, numaa)
    simulate(f, seqfile, treefile, mu, seqlength, None) # omega is last arguement. when None, sim via mutsel
    
    # Use math to derive an omega for the simulated sequences. Also returns the number of codons theoretically sampled.
    print "deriving"
    derived_w, num_codons = deriveOmega(f)
       
    # Nei-Gojobori Method
    print "nei-gojobori"
    nei_w, ns_mut, s_mut = run_neigojo(seqfile)

    # PAML
    print "ML"
    pamlw = []
    for cf in range(4): # 1/61, f1x4, f3x4, data
        pamlw.append( runpaml(seqfile, str(cf)) )
    
    # Save. CURRECTLY HARDCODED THAT PAMLW HAS LENGTH 4. 
    outf.write(rep + '\t' + str(numaa) + '\t' + freqType[ftype] + '\t' + str(derived_w) + '\t' + str(nei_w) + '\t' + str(pamlw[0]) + '\t' + str(pamlw[1]) + '\t' + str(pamlw[2]) + '\t' + str(pamlw[3]) + '\n')
    
outf.close()

