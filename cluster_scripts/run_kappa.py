# SJS.
# code to examine kappa influence 

import sys
sys.path.append('src/')
from misc import *
from newick import *
from stateFreqs import *
from matrixBuilder import *
from evolver import *
from functions_simandinf import *
import random as rn

# Input parameters and global stuff
if (len(sys.argv) != 5):
    print "\n\nUsage: python run_kappa.py <rep> <mu> <bl> <seqlength>\n."
    sys.exit()
rep = sys.argv[1]
mu = float(sys.argv[2])
bl = sys.argv[3]
seqlength = int(sys.argv[4])

# get kappa, beta
kappa = rn.uniform(1.0, 5.0)
# Set up beta for getting amino acid frequencies
beta = rn.uniform(0.5,3.5)


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
f, num_pref_aa, gc_content = setFreqs(freqfile, beta, 0.0, 1.0) # last 2 args are gc min, gc max
simulate(f, seqfile, treefile, mu, kappa, seqlength, None) # omega is last argument. when None, sim via mutsel
    
# Derive omega
print "deriving"
mu_dict = {'AT':mu, 'AC':mu, 'AG':mu, 'CG':mu, 'CT':mu, 'GT':mu}
mu_dict['AG'] = mu_dict['AG'] * kappa
mu_dict['CT'] = mu_dict['CT'] * kappa
derived_w = deriveOmega(f, mu_dict)

# ML omega
gy94_w = float(runpaml(seqfile))

err = (derived_w - gy94_w) / derived_w

# Save
outf = open(outfile,'w')
outf.write(rep + '\t' + str(seqlength) + '\t' + str(bl) + '\t' + str(mu) + '\t' + str(kappa) + '\t' + str(beta) + '\t' + str(derived_w) + '\t' + str(gy94_w) + '\t' + str(err) + '\n')
outf.close()





