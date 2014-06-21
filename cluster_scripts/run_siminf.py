# SJS.
# Generic code for simulating and deriving omega via math, hyphy ML.
# Be sure to cp src/ directory (simulator), hyphy files, and the functions_simandinf.py script into working directory
# NOTE: very little to no sanity checking for input args


import sys
from functions_simandinf import *

sys.path.append('src/')
from misc import *
from newick import *
from stateFreqs import *
from matrixBuilder import *
from evolver import *


# Input parameters and global stuff
if (len(sys.argv) != 8):
    print "\n\nUsage: python run_siminf.py <rep> <cpu> <numaa> <mu> <kappa> <bl> <seqlength>\n."
    sys.exit()
rep = sys.argv[1]
cpu = sys.argv[2]
numaa = int(sys.argv[3])
mu = float(sys.argv[4])
kappa = float(sys.argv[5])
bl = sys.argv[6]
seqlength = int(sys.argv[7])


# Write tree given bl specifications
treefile = "tree.tre"
treef = open(treefile, 'w')
treef.write("(t1:" + str(bl) + ", t2:" + str(bl) + ");")
treef.close()

# set up output sequence and parameter files
seqfile = "seqs"+str(rep)+".fasta"
freqfile = "codonFreqs" + str(rep)+".txt"
outfile = "params"+str(rep)+".txt"

# Simulate
print "simulating"
f, gc_content, aminos_used = setFreqs(numaa, freqfile, 0., 1.) # last 2 args are gc min, gc max
simulate(f, seqfile, treefile, mu, kappa, seqlength, None) # omega is last argument. when None, sim via mutsel
    
# Derive omega
print "deriving"
mu_dict = {'AT':mu, 'AC':mu, 'AG':mu, 'CG':mu, 'CT':mu, 'GT':mu}
mu_dict['AG'] = mu_dict['AG'] * kappa
mu_dict['CT'] = mu_dict['CT'] * kappa
derived_w = deriveOmega(f, mu_dict)

# HyPhy omega
print "ML"
gy94_w = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, kappa)

# Save
outf = open(outfile,'w')
outf.write(rep + '\t' + str(numaa) + '\t' + aminos_used + '\t' + str(seqlength) + '\t' + str(bl) + '\t' + str(mu) + '\t' + str(kappa) + '\t' + str(round(gc_content, 5)) + '\t' + str(round(derived_w, 5)) + '\t' + str(round(gy94_w, 5)) + '\n')
outf.close()





