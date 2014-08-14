# SJS. stephanie.spielman@gmail.com
# Code for using Jesse Bloom's NP data.

######## Input parameters ########
import sys
if (len(sys.argv) != 5):
    print "\n\nUsage: python run_siminf.py <rep> <treefile> <simdir> <cpu>\n."
    sys.exit()
rep = sys.argv[1]         # which rep we're on, for saving files. needs to run from 1-498, since 498 sites.
treefile = sys.argv[2]    # tree for simulation
simdir = sys.argv[3]      # directory of simulation library
cpu = sys.argv[4]         # hyphy can use
sys.path.append(simdir)
from functions_simandinf import *


# Set up output files and parameters
seqfile       = "seqs"+str(rep)+".fasta"
paramfile     = "params"+str(rep)+".txt"


seqlength = 500000
# from the NP paper, Bloom 2014 MBE.
mu_dict = {'AG':2.4e-5, 'TC':2.4e-5, 'GA':2.3e-5, 'CT':2.3e-5, 'AC':9.0e-6, 'TG':9.0e-6, 'CA':9.4e-6, 'GT':9.4e-6, 'AT':3.0e-6, 'TA':3.0e-6, 'GC':1.9e-6, 'CG':1.9e-6}
np_freqs = np.loadtxt("np_codon_eq_freqs.txt") # these have been calc'd in Omega_MutSel/np_scripts .

# Read in equilibrium frequencies and determine entropy
codon_freqs_true = np_freqs[ int(rep) - 1 ]
codon_freqs_true_dict = dict(zip(codons, codon_freqs_true))
entropy = calc_entropy(codon_freqs_true)


# Simulate according to MutSel model along phylogeny
print "Simulating"
simulate(codon_freqs_true, seqfile, treefile, mu_dict, seqlength)


# Derive omega from selection coefficients (well, frequencies, but same deal)
print "Deriving omega from selection coefficients"
derivedw = derive_omega(codon_freqs_true_dict, mu_dict, bias=False)


# Maximum likelihood omega inference across a variety of frequency, kappa specifications
print "Conducting ML inference with HyPhy"


# Lists for storing values and printing strings
krun = [1.0, 'free']
kspecs = ['one', 'free'] # no kappa true here, since there's no kappa.
fspecs = ['equal', 'f61_site', 'f61_global', 'f3x4_site', 'f3x4_global', 'cf3x4_site', 'cf3x4_global'] # DO NOT CHANGE THIS LIST !!!!
omegas = np.zeros([2,7])
kappas = np.zeros([2,7])
omega_errors = np.ones([2,7])


# First, set up F61 (data) frequency vector in the hyphy batchfile as this applies to all hyphy runs.
hyf = array_to_hyphy_freq(codon_freqs_true)
setuphyphyf = "sed -i 's/DATAFREQS/"+hyf+"/g' globalDNDS_np.bf"
setupf = subprocess.call(setuphyphyf, shell = True)
assert(setupf == 0), "couldn't properly add in F61 (data) frequencies"


# Run hyphy and save omegas, kappas (only sometimes returned, note), and omega errors along the way
kcount = 0
for kap in krun:
    wtemp, ktemp = run_hyphy_np(seqfile, treefile, cpu, kap, fspecs)  
    kappas[kcount] = ktemp
    omegas[kcount] = wtemp
    omega_errors[kcount] = (derivedw - wtemp) / derivedw
    kcount += 1


# Finally, save results
outstring_params = rep + '\t' + str(seqlength) + '\t' + str(entropy) + '\t' + str(derivedw)
outf = open(paramfile, 'w')
for f in fspecs:
    y =  fspecs.index(f)
    for k in kspecs:
        x = kspecs.index(k)
        outf.write( outstring_params + '\t' + f + '\t' + k + '\t' + str(omegas[x,y]) + '\t' + str(omega_errors[x,y]) + '\t' + str(kappas[x,y]) + '\n')
outf.close()   
