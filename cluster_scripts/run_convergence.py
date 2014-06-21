# SJS.
# Code specifically for generating data to demonstrate omega convergence.
# Parameters: fix kappa to 1.0 and beta to 3.0 in the qsub file.
# Vary seq length from 10^2 - 10^6, for a single frequency distribution. Hence, for-loop within.

import sys
from functions_simandinf import *

sys.path.append('src/')
from misc import *
from newick import *
from stateFreqs import *
from matrixBuilder import *
from evolver import *


# Input parameters and global stuff
if (len(sys.argv) != 7):
    print "\n\nUsage: python run_convergence.py <rep> <cpu> <beta> <mu> <kappa> <bl> \n."
    sys.exit()
rep = sys.argv[1]
cpu = sys.argv[2]
beta = float(sys.argv[3])
mu = float(sys.argv[4])
kappa = float(sys.argv[5])
bl = sys.argv[6]
seqlens = [500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000, 10000000]

# Write tree given bl specifications
treefile = "tree.tre"
treef = open(treefile, 'w')
treef.write("(t1:" + str(bl) + ", t2:" + str(bl) + ");")
treef.close()

# set up output sequence and parameter files
freqfile = "codonFreqs" + str(rep)+".txt"
outfile = "params"+str(rep)+".txt"
outf = open(outfile,'w')

# First, determine frequencies and generate the derived omega from these
derived_w = 0.
while derived_w == 0.: # Must ensure that derived_w is greater than 0, because evolution is good and dividing by zero to determine relative error is bad.
    f, num_pref_aa, gc_content = setFreqs(freqfile, beta, 0.0, 1.0) # last 2 args are gc min, gc max
    mu_dict = {'AT':mu, 'AC':mu, 'AG':mu, 'CG':mu, 'CT':mu, 'GT':mu}
    mu_dict['AG'] = mu_dict['AG'] * kappa
    mu_dict['CT'] = mu_dict['CT'] * kappa
    derived_w = deriveOmega(f, mu_dict)

# Now, simulate sequences and infer ML omegas
for seqlength in seqlens:
    # Simulate
    seqfile = "seqs"+str(rep)+"_"+str(seqlens.index(seqlength))+".fasta"
    print "simulating"
    simulate(f, seqfile, treefile, mu, kappa, int(seqlength), None) # omega is last argument. when None, sim via mutsel
    
    # HyPhy omega
    print "ML"
    gy94_w = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, kappa)
    
    # Calculate relative error from derived omega. Not using abs() since plot nicer that way.
    err = ( derived_w - gy94_w )/derived_w

    # Save
    outf.write(rep + '\t' + str(seqlength) + '\t' + str(bl) + '\t' + str(mu) + '\t' + str(kappa) + '\t' + str(num_pref_aa) + '\t' + str(round(gc_content, 6)) + '\t' + str(round(derived_w, 6)) + '\t' + str(round(gy94_w, 6)) + '\t' + str(round(rel_err, 6)) + '\n')
outf.close()





