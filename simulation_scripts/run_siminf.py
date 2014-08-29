# SJS. stephanie.spielman@gmail.com
# Generic code for simulating and deriving dN/dS via selection coeffcients and hyphy ML.
# Be sure to cp src/ directory (simulator), hyphy files, and the functions_simandinf.py script into working directory
# NOTE: very little (ok, none) sanity checking for input args..

######## Input parameters ########
import sys
if (len(sys.argv) != 6):
    print "\n\nUsage: python run_siminf.py <rep> <treefile> <simdir> <cpu> <bias>\n."
    sys.exit()
rep = sys.argv[1]         # which rep we're on, for saving files
treefile = sys.argv[2]    # tree for simulation
simdir = sys.argv[3]      # directory of simulation library
cpu = sys.argv[4]         # hyphy can use
bias = float(sys.argv[5]) # codon bias, yes or no?
sys.path.append(simdir)
from functions_simandinf import *


# Set up output files and parameters
seqfile       = "seqs"+str(rep)+".fasta"
freqfile      = "codonFreqs" + str(rep)+".txt"
paramfile     = "params"+str(rep)+".txt"


seqlength = 500000
if bias != 0.:
    bias = rn.uniform(ZERO,2) # 2 is a good top threshold. dN/dS typically <=2 this way, otherwise it gets absurd. 
mu = 1e-6
kappa = rn.uniform(1.0, 6.0)
sd = rn.uniform(0., 4.)
mu_dict = {'AT': mu, 'TA':mu, 'CG': mu, 'GC':mu, 'AC': mu, 'CA':mu, 'GT':mu, 'TG':mu, 'AG': kappa*mu, 'GA':kappa*mu, 'CT':kappa*mu, 'TC':kappa*mu}


# Derive equilibrium codon frequencies 
print "Deriving equilibrium codon frequencies"
codon_freqs_f61, codon_freqs_dict, gc_content, entropy = set_codon_freqs(sd, freqfile, bias)


# Simulate according to HB98 MutSel model
print "Simulating"
simulate(codon_freqs_f61, seqfile, treefile, mu_dict, seqlength)


# Derive dN/dS from equilibrium codon frequencies
print "Deriving dN/dS from equilibrium codon frequencies"
dnds = derive_dnds(codon_freqs_dict, mu_dict)


# Maximum likelihood omega inference for kappa={true,free} and Fequal only.
print "Conducting ML inference with HyPhy"

omegas = np.zeros(2)
kappas = np.zeros(2)

omegas[0], kappas[0] = run_hyphy_fequal(seqfile, treefile, cpu, kappa)  
omegas[1], kappas[1] = run_hyphy_fequal(seqfile, treefile, cpu, "free")  
omega_errors = (dnds - omegas)/dnds


# Finally, save results
outstring_params = rep + '\t' + str(seqlength) + '\t' + str(mu) + '\t' + str(kappa) + '\t' + str(sd) + '\t' + str(bias) + '\t' + str(gc_content) + '\t' + str(entropy) + '\t' + str(dnds)
outf = open(paramfile, 'w')
for k in range(2):
    x = kspecs.index(k)
    outf.write( outstring_params + '\t' + f + '\t' + k + '\t' + str(omegas[x]) + '\t' + str(omega_errors[x]) + '\t' + str(kappas[x]) + '\n')
outf.close()   








