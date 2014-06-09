## SJS.
## Script to test whether, if I simulate according to a GY94 model but with similar frequencies to how we simulate MutSel sequences
## PIPELINE:
#### Run 3-15 amino acids, 100 reps each, for omegas = [0.05, 0.1, 0.25, 0.5, 0.75, 1.0].


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
mu = 0.005
length = 10000
omegas = [0.05, 0.1, 0.25, 0.5, 0.75, 1.0]


# Write treestring to file
treefile = "tree.tre"
treef = open(treefile, 'w')
treef.write("(t1:0.5, t2:0.5);")
treef.close()

outf = open("tempout.txt", 'w')

for omega in omegas:
    # Simulate
    print "simulating"
    f, aminos_used = setFreqs("user", "amino", numaa)
    simulate(f, seqfile, treefile, mu, length, omega)
    
    
    # Use math to derive an omega for the simulated sequences. Also returns the number of codons theoretically sampled.
    print "deriving"
    #derived_w, num_codons = deriveOmega(f)

    # Nei-Gojobori Method
    print "nei-gojobori"
    nei_w, ns_mut, s_mut = run_neigojo(seqfile)

    # PAML and HyPhy
    paml_w = runpaml(seqfile)
    hyphy_w_kappafixed = runhyphy("globalGY94.bf", "GY94_fixedkappa", seqfile, treefile, cpu, f)
    outf.write(rep + '\t' + str(numaa) + '\t' + str(omega) + '\t' + str(derived_w) + '\t' + str(nei_w) +  '\t' + str(paml_w) + '\t' + str(hyphy_w_kappafixed) + '\t' + aminos_used + '\t' + str(num_codons) + '\n')

outf.close()
            

# And now send to the final outfile
save = "cat tempout.txt >> "+final_outfile
saverun = subprocess.call(save, shell=True)
assert(saverun == 0), "couldn't save final file"
save2 = "cp "+seqfile+" "+rdir
saverun2 = subprocess.call(save2, shell=True)
assert(saverun2 == 0), "couldn't save seqfile"



