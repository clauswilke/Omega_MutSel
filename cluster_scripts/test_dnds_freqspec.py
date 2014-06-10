## SJS.
## Script to test whether dnds is similarly susceptible to frequency specification
## PIPELINE:
#### Run 3-15 amino acids, 100 reps each, for omegas = [0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 2.5].


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


# Input parameters
cpu = sys.argv[1]
rep = sys.argv[2]
numaa = int(sys.argv[3])
seqlength = 10000
treebl =  0.1
mu = 0.001

# Global stuff
seqfile = "rep"+str(rep)+'.fasta'
mu = 0.001
length = 10000
omegas = [0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 2.5]
freqType = {"user": "expAA", "equal": "equalAA"}


# Write treestring to file
treefile = "tree.tre"
treef = open(treefile, 'w')
treef.write("(t1:0.01, t2:0.01);")
treef.close()

outf = open("tempout.txt", 'w')

for ftype in freqType.keys():

    # Frequencies for simulation. Same for all omegas in a given rep
    f, aminos_used = setFreqs(ftype, numaa)
    
    for omega in omegas:
   
        seqfile = 'seqs_'+freqType[ftype]+str(rep)+'_'+str(omega)+'.fasta'

        # Simulate
        print "simulating"
        simulate(f, seqfile, treefile, mu, length, omega)
    
        # Nei-Gojobori Method
        print "nei-gojobori"
        nei_w, ns_mut, s_mut = run_neigojo(seqfile)

        # PAML
        print "ML"
        pamlw = []
        for cf in range(4): # 1/61, f1x4, f3x4, data
            pamlw.append( runpaml(seqfile, str(cf)) )
    

        # Save. CURRECTLY HARDCODED THAT PAMLW HAS LENGTH 4. 
        outf.write(rep + '\t' + str(numaa) + '\t' + freqType[ftype] + '\t' + str(omega) + '\t' + str(nei_w) + '\t' + str(pamlw[0]) + '\t' + str(pamlw[1]) + '\t' + str(pamlw[2]) + '\t' + str(pamlw[3]) + '\n')
    
outf.close()




