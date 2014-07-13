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
fspecs = ['equal', 'f3x4', 'data']
kspecs = {1.:'one', kappa:'true', 'free':'free'}

common_out_string = rep + '\t' + str(seqlength) + '\t' + str(bl) + '\t' + str(mu) + '\t' + str(kappa) + '\t' + str(gc_content) + '\t' + str(beta) + '\t' + str(entropy) + '\t' + str(derivedw)

outf = open(outfile, 'w')
for kapspec in kspecs:
    for freqspec in fspecs:
        mlw, mlk = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, kapspec, f, freqspec)
        w_err = (derivedw - mlw) / derivedw
        if mlk is not None:
            k_err = (kappa - mlk) / kappa
        else:
            mlk = "None"
            k_err = "None"   
        outf.write(common_out_string + '\t' + str(freqspec) + '\t' + str(kspecs[kapspec]) + '\t' + str(mlw) + '\t' + str(w_err) + '\t' + str(mlk) + '\t' + str(k_err) + '\n')
outf.close()





