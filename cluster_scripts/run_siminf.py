# SJS.
# Generic code for simulating and deriving omega via math, hyphy ML.
# Be sure to cp src/ directory (simulator), hyphy files, and the functions_simandinf.py script into working directory
# NOTE: very little to no sanity checking for input args

######## Input parameters ########
import sys
if (len(sys.argv) != 7):
    print "\n\nUsage: python run_siminf.py <rep> <cpu> <mu> <bl> <seqlength> <simdir> \n."
    sys.exit()
rep = sys.argv[1]
cpu = sys.argv[2]
mu = float(sys.argv[3])
bl = sys.argv[4]
seqlength = int(sys.argv[5])
simdir = sys.argv[6]
sys.path.append(simdir)
from functions_simandinf import *
kappa = rn.uniform(1.0, 5.0)
beta = rn.uniform(0.5, 3.0)



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
entropy = calcCodonEntropy(f)

# Derive omega
print "deriving"
mu_dict = {'AT':mu, 'AC':mu, 'AG':mu*kappa, 'CG':mu, 'CT':mu*kappa, 'GT':mu}
derivedw = deriveOmega(f, mu_dict)


# HyPhy omega and kappa, forever and ever
print "ML"
ml_estimates = np.zeros(9)
kappas = np.zeros(3)
count = 0

# kappa constrained to 1
ml_estimates[0], noneK = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, 1., f, "equal")
ml_estimates[1], noneK  = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, 1., f, "f3x4")
ml_estimates[2], noneK  = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, 1., f, "data")


# kappa fixed as true value
ml_estimates[3], noneK = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, kappa , f, "equal")
ml_estimates[4], noneK  = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, kappa , f, "f3x4")
ml_estimates[5], noneK  = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, kappa , f, "f3x4")

# kappa as a free parameter
ml_estimates[6], kappas[0]  = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, 'free' , f, "equal")
ml_estimates[7], kappas[1]  = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, 'free' , f, "f3x4")
ml_estimates[8], kappas[2]  = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, 'free', f, "data")



out_string = rep + '\t' + str(seqlength) + '\t' + str(bl) + '\t' + str(mu) + '\t' + str(kappa) + '\t' + str(gc_content) + '\t' + str(beta) + '\t' + str(entropy) + '\t' + str(derivedw)

#mlw and error
for mlw in ml_estimates:
    err = (derivedw - mlw) / derivedw
    temp = '\t' + str(mlw) + '\t' + str(err)
    out_string += temp

#kappaw and error
for kap in kappas:
    err = (kappa - kap) / kappa
    temp = '\t' + str(kap) + '\t' + str(err)
    out_string += temp
    
out_string += '\n'

# Save
outf = open(outfile,'w')
outf.write(out_string)
outf.close()





