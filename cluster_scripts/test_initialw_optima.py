## SJS.
## Script to test whether the initial omega guess causes local optima problems.
## PIPELINE:
#### Use prespecified amino acid list for [3,7] amino acids (set this up in optima.qsub).
#### Test the following 4 intial omega guesses - 0.01, 0.1, 0.5, 1.0
#### Examines if a local optima messes things up.



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
rdir = sys.argv[2]
if rdir[-1] != '/':
	rdir += '/'
final_outfile = rdir + "../" + sys.argv[3]
rep = sys.argv[4]
numaa = int(sys.argv[5])


# Global stuff
seqfile = "rep"+str(rep)+'.fasta'
mu = 0.001
length = 10000
intial_omega = [0.01, 0.1, 0.5, 1.0]


# Write treestring to file
treefile = "tree.tre"
treef = open(treefile, 'w')
treef.write("(t1:0.5, t2:0.5);")
treef.close()


# Simulate
print "simulating"
f, aminos_used = simulate(seqfile, numaa, "user", "amino", treefile, mu, length, True) # last argument = don't use prespecified amino acid choices.

# Use math to derive an omega for the simulated sequences
print "deriving"
derived_w = deriveOmega(f)

# Nei-Gojobori Method
print "nei-gojobori"
nei_w = run_neigojo(seqfile)


# PAML and HyPhy with several initial omega values. Save to file simultaneously (yes, derived and nei get written several times. no big.)
outf = open("tempout.txt", 'w')
for initw in initial_omega:
    paml_w = runpaml(seqfile, initw)
    hyphy_w_kappafixed = runhyphy("globalGY94.bf", "GY94_fixedkappa", seqfile, treefile, cpu, f, initw)
    outf.write(rep + '\t' + str(numaa) + '\t' + str(derived_w) + '\t' + str(nei_w) + '\t' + str(initw) + '\t' + str(paml_w) + '\t' + str(hyphy_w_kappafixed) + '\n')
outf.close()
		

# And now send to the final outfile
save = "cat tempout.txt >> "+final_outfile
saverun = subprocess.call(save, shell=True)
assert(saverun == 0), "couldn't save final file"
save2 = "cp "+seqfile+" "+rdir
saverun2 = subprocess.call(save2, shell=True)
assert(saverun2 == 0), "couldn't save seqfile"



