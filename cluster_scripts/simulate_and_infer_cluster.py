## SJS.
## Script to run mutsel simulations and infer dNdS using either Nei-Gojobori or HyPhy/PAML ML.
## NOTE: to use simulation library must copy it 

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


cpu = sys.argv[1]
rdir = sys.argv[2]
if rdir[-1] != '/':
	rdir += '/'
final_outfile = rdir + sys.argv[3]
rep = sys.argv[4]
numaa = int(sys.argv[5])

outfile = rdir+"mutsel_"+str(rep)+".txt"
seqfile = "rep"+str(rep)+'.fasta'
mu = 0.001
length = 10000

# Write treestring to file
treefile = "tree.tre"
treef = open(treefile, 'w')
treef.write("(t1:0.5, t2:0.5);")
treef.close()

# Simulate
print "simulating"
f = simulate(seqfile, numaa, "user", "amino", treefile, mu, length)

# Use math to derive an omega for the simulated sequences
print "deriving"
derived_w = deriveOmega(f)

# Nei-Gojobori Method
print "nei-gojobori"
nei_w = run_neigojo(seqfile)

# Send to PAML
print "paml"
paml_w = runpaml(seqfile)

# Send to HyPhy
print "hyphy"
hyphy_w_kappafree = runhyphy("globalGY94.bf", "GY94", seqfile, treefile, cpu, f)
hyphy_w_kappafixed = runhyphy("globalGY94.bf", "GY94_fixedkappa", seqfile, treefile, cpu, f)
		
# Now save everything to file
outf = open(outfile, 'w')
outf.write(str(numaa) + '\t' + str(derived_w) + '\t' + str(nei_w) + '\t' + str(paml_w) + '\t' + str(hyphy_w_kappafree) + '\t' + str(hyphy_w_kappafixed)+ '\n')
outf.close()

# And now send to the final outfile
save = "cat "+outfile+" >> "+final_outfile
saverun = subprocess.call(save, shell=True)
assert(saverun == 0), "couldn't save final file"
save2 = "cp "+seqfile+" "+rdir
saverun2 = subprocess.call(save2, shell=True)
assert(saverun2 == 0), "couldn't save seqfile"



