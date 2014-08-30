# SJS. Code specifically for generating data to demonstrate omega convergence.


import sys
# Input parameters and global stuff
if (len(sys.argv) != 4):
    print "\n\nUsage: python run_convergence.py <rep> <simdir> <cpu> \n."
    sys.exit()
rep = sys.argv[1]
simdir = sys.argv[2]
cpu = sys.argv[3]
sys.path.append(simdir)
from functions_simandinf import *

# Set up output files and parameters
seqfile   = "seqs"+str(rep)+".fasta"
freqfile  = "codonFreqs" + str(rep)+".txt"
paramfile = "params"+str(rep)+".txt"
treefile = 'tree.tre'

mu = 1e-6
kappa = rn.uniform(1.0, 6.0)
sd = rn.uniform(0., 4.)
mu_dict = {'AT': mu, 'TA':mu, 'CG': mu, 'GC':mu, 'AC': mu, 'CA':mu, 'GT':mu, 'TG':mu, 'AG': kappa*mu, 'GA':kappa*mu, 'CT':kappa*mu, 'TC':kappa*mu}

# To test convergence, select random sequence length between 5e2 and 1e6
expon = rn.randint(2,5)
if expon == 2:
    times = rn.randint(5,10)
else:
    times = randint(1,10)
seqlength = int( times * 10**expon )


# Create treefile
write_treefile(treefile)

# Set up steady-state codon frequencies based on selection coefficients
print "Deriving equilibrium codon frequencies"
codon_freqs, codon_freqs_dict, gc_content, entropy = set_codon_freqs(sd, freqfile, 0.)


# Simulate according to MutSel model along phylogeny
print "Simulating"
simulate(codon_freqs, seqfile, treefile, mu_dict, seqlength)

# Derive dN/dS from selection coefficients (well, frequencies, but same deal)
print "Deriving dN/dS from selection coefficients"
dnds = derive_dnds(codon_freqs_dict, mu_dict)

# ML
print "Conducting ML inference with HyPhy"
mlw, k = run_hyphy_fequal(seqfile, treefile, cpu, kappa)
err = (dnds - mlw) / dnds

# Save
outf = open(paramfile,'w')
outf.write(rep + '\t' + str(seqlength) + '\t' + str(mu) + '\t' + str(kappa) + '\t' + str(sd) + '\t' + str(dnds) + '\t' + str(mlw) + '\t' + str(err) + '\n')
outf.close()





