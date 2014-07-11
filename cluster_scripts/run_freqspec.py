# SJS.
# Script to simulate and demonstrate that the frequency specification really matters.
# Be sure to cp src/ directory (simulator), **paml** files, and the functions_simandinf.py script into working directory


import sys
sys.path.append('src/')
from misc import *
from newick import *
from stateFreqs import *
from matrixBuilder import *
from evolver import *
from functions_simandinf import *


# Input parameters and global stuff
if (len(sys.argv) != 6):
    print "\n\nUsage: python run_freqspec.py <rep> <mu> <bl> <seqlength> <cpu>\n."
    sys.exit()
rep = sys.argv[1]
mu = float(sys.argv[2])
bl = sys.argv[3]
seqlength = int(sys.argv[4])
cpu = sys.argv[5]


beta = rn.uniform(0.5,3.5)
kappa = rn.uniform(1.0, 5.0)

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

# Sequence entropy based on codon frequencies
entropy = calcCodonEntropy(f)


# Derive omega
print "deriving"
mu_dict = {'AT':mu, 'AC':mu, 'AG':mu, 'CG':mu, 'CT':mu, 'GT':mu}
mu_dict['AG'] = mu_dict['AG'] * kappa
mu_dict['CT'] = mu_dict['CT'] * kappa
derived_w = deriveOmega(f, mu_dict)


# PAML omegas. Using paml since makes frequency specification a lot easier.
print "ML"
gy94_w_equal = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, kappa, f, "equal")
gy94_w_3x4   = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, kappa, f, "f3x4")
gy94_w_data  = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, kappa, f, "data")

err_equal = (derived_w - gy94_w_equal) / derived_w
err_3x4 = (derived_w - gy94_w_3x4) / derived_w
err_data = (derived_w - gy94_w_data) / derived_w

# Save
outf = open(outfile,'w')
outf.write(rep + '\t' + str(seqlength) + '\t' + str(bl) + '\t' + str(mu) + '\t' + str(kappa) + '\t' + str(entropy) + '\t' + str(gc_content) + '\t' + str(derived_w) + '\t' + str(gy94_w_equal) + '\t' + str(gy94_w_3x4) + '\t' + str(gy94_w_data) + '\t' + str(err_equal) + '\t' + str(err_3x4) + '\t' + str(err_data) + '\n')
outf.close()





