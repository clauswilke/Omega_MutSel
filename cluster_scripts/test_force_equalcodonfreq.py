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

pamlFreq = {"0": "1/61", "2": "F3X4", "3": "empirical"}
freqType = {"user": "expAA", "equal": "equalAA"}

# Write tree given bl specifications
treefile = "tree.tre"
treef = open(treefile, 'w')
treef.write("(t1:" + str(treebl) + ", t2:" + str(treebl) + ");")
treef.close()



outf = open("tempout.txt", 'w')


for type in freqType.keys():

    seqfile = 'seqs_'+freqType+str(rep)+'.fasta'

    # Simulate
    print "simulating"
    f, aminos_used = setFreqs(freqType, numaa)
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
    for cf in pamlFreq.keys():
        pamlw.append( runpaml(seqfile, cf) )
    
    # Save. CURRECTLY HARDCODED THAT PAMLW HAS LENGTH 3. 
    outf.write(rep + '\t' + str(numaa) + '\t' + freqType[type] + '\t' + str(derived_w) + '\t' + str(nei_w) + '\t' + str(pamlw[0]) + '\t' + str(pamlw[1]) + '\t' + str(pamlw[2]) + '\n')
    
outf.close()

