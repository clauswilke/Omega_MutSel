# SJS. stephanie.spielman@gmail.com
# Generic code for simulating and deriving dN/dS via selection coeffcients and hyphy ML.
# Note that we only run equal frequenies and true kappa here.

######## Input parameters ########
import sys
if (len(sys.argv) != 6):
    print "\n\nUsage: python run_siminf.py <rep> <treefile> <simdir> <cpu> <bias>\n."
    sys.exit()
rep = sys.argv[1]         # which rep we're on, for saving files
treefile = sys.argv[2]    # tree for simulation
simdir = sys.argv[3]      # directory of simulation library
cpu = sys.argv[4]         # hyphy can use
bias = float(sys.argv[5]) # codon bias factor
sys.path.append(simdir)
from functions_simandinf import *


# Set up output files and parameters
seqfile       = "seqs"+str(rep)+".fasta"
freqfile      = "codonFreqs" + str(rep)+".txt"
amino_sscfile = "aminoCoeffs" + str(rep)+".txt"
codon_sscfile = "codonCoeffs" + str(rep)+".txt"
paramfile     = "params"+str(rep)+".txt"


seqlength = 500000
mu = 1e-6
kappa = rn.uniform(1.0, 6.0)
sd = rn.uniform(0., 4.)
mu_dict = {'AT': mu, 'TA':mu, 'CG': mu, 'GC':mu, 'AC': mu, 'CA':mu, 'GT':mu, 'TG':mu, 'AG': kappa*mu, 'GA':kappa*mu, 'CT':kappa*mu, 'TC':kappa*mu}


# Set up steady-state codon frequencies based on selection coefficients
print "Deriving equilibrium codon frequencies"
codon_freqs_true, codon_freqs_true_dict, gc_content, entropy = set_codon_freqs(sd, freqfile, amino_sscfile, codon_sscfile, bias)


# Simulate according to MutSel model along phylogeny
print "Simulating"
simulate(codon_freqs_true, seqfile, treefile, mu_dict, seqlength)


# Derive omega from selection coefficients (well, frequencies, but same deal)
print "Deriving omega from selection coefficients"
derivedw = derive_omega(codon_freqs_true_dict, mu_dict, bias!=0.) # last argument as bool for function to know whether to compute dS.


# ML, with true kappa and equal frequencies
print "Conducting ML inference with HyPhy"
mlw = run_hyphy_convergence(seqfile, treefile, cpu, kappa)
werr = (derivedw - mlw) / derivedw

# Finally, save results
outstring_params = rep + '\t' + str(seqlength) + '\t' + str(mu) + '\t' + str(kappa) + '\t' + str(sd) + '\t' + str(bias) + '\t' + str(gc_content) + '\t' + str(entropy) + '\t' + str(derivedw)
outf = open(paramfile, 'w')
outf.write( outstring_params + '\t' + str(mlw) + '\t' + str(werr) + '\n')
outf.close()   







