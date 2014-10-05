# SJS. stephanie.spielman@gmail.com
# Code for using Jesse Bloom's NP amino acid preferene data, with mutation rates either from NP (Bloom 2014) or yeast.

######## Input parameters ########
import sys
if (len(sys.argv) != 4):
    print "\n\nUsage: python run_np.py <rep> <cpu> <dataset> \n."
    sys.exit()
rep = sys.argv[1]         # which rep we're on, for saving files. needs to run from 1-498, since 498 sites.
cpu = sys.argv[2]         # hyphy can use
dataset = sys.argv[3]     # either np, yeast, or polio. determines the mutation scheme and eq freqs

from functions_omega_mutsel import *

batchfile = 'batchfile.bf'
seqfile   = "seqs"+str(rep)+".fasta"
paramfile = "params"+str(rep)+".txt"
write_treefile(treefile)

# Read in equilibrium frequencies, determine entropy, and set up mu_dict
mu_dict = set_mu_dict(dataset)
eq_codon_freqs = np.loadtxt(dataset + "_codon_eqfreqs.txt")
codon_freqs = eq_codon_freqs[ int(rep) - 1 ]
codon_freqs_dict = dict(zip(codons, codon_freqs))
entropy = calc_entropy(codon_freqs)


# Derive omega from eq freqs
print "Deriving omega from equilibrium codon frequencies"
dnds = derive_dnds(codon_freqs_dict, mu_dict)


# Maximum likelihood omega inference across a variety of frequency specifications.
print "Conducting ML inference with HyPhy"

fspecs = ['f61',  'f1x4',  'f3x4',  'cf3x4',  'fnuc1',  'fnuc3']
k = 3. # all models have 3 free parameters (w,k,t). note that t is ok to be global, because simulation trees all have same branch lengths, so in fact t should be a global (not local/branch-specific) parameter.
lnliks, omegas, kappas = run_hyphy_nyp(batchfile, seqfile, treefile, cpu, fspecs)  
AICs = 2*(k - lnliks)
omega_errors = (dnds - omegas) / dnds


# Finally, save results
outstring_params = rep + '\t' + str(entropy) + '\t' + str(dnds)
with open(paramfile, 'w') as outf:
    for x in range(len(fspecs)):
        outf.write( outstring_params + '\t' + fspecs[x] + '\t' +  str(omegas[x]) + '\t' + str(omega_errors[x]) + '\t' + str(kappas[x]) + '\t' + str(lnliks[x]) + '\t' + str(int(k)) + '\t' + str(AICs[x]) + '\n')
