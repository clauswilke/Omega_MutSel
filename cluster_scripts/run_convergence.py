# SJS. Code specifically for generating data to demonstrate omega convergence.




import sys
# Input parameters and global stuff
if (len(sys.argv) != 5):
    print "\n\nUsage: python run_convergence.py <rep> <treefile> <simdir> <cpu> \n."
    sys.exit()
rep = sys.argv[1]
treefile = sys.argv[2]
simdir = sys.argv[3]
cpu = sys.argv[4]
sys.path.append(simdir)
from functions_simandinf import *

mu = 1e-5
kappa = rn.uniform(1.0, 5.0)
lambda_ = rn.uniform(0.5, 2.0)

# Random sequence length between 5e2 and 1e6
expon = rn.randint(2,5)
if expon == 2:
    times = rn.randint(5,10)
else:
    times = randint(1,10)
seqlength = int( times * 10**expon )

# set up output sequence and parameter files
freqfile = "codonFreqs" + str(rep) + "_" + str(seqlength) + ".txt"
seqfile = "seqs" + str(rep) + "_" + str(seqlength) + ".fasta"
outfile = "params" + str(rep) + "_" + str(seqlength) + ".txt"


# Now, simulate sequences and infer ML omegas
# Simulate
print "simulating"
f, gc_content = setFreqs(freqfile, lambda_, 0.0, 1.0) # last 2 args are gc min, gc max
simulate(f, seqfile, treefile, mu, kappa, seqlength, None) # omega is last argument. when None, sim via mutsel

# Derive
mu_dict = {'AT':mu, 'AC':mu, 'AG':mu*kappa, 'CG':mu, 'CT':mu*kappa, 'GT':mu}
derivedw = deriveOmega(f, mu_dict)

# ML
f_equal = np.zeros(61)
f_equal[f_equal == 0.] = 1./61.
mlw, mlk = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, kappa, f_equal)
err = (derivedw - mlw) / derivedw

# Save
outf = open(outfile,'w')
outf.write(rep + '\t' + str(seqlength) + '\t' + str(kappa) + '\t' + str(lambda_) + '\t' + str(derivedw) + '\t' + str(mlw) + '\t' + str(err) + '\n')
outf.close()





