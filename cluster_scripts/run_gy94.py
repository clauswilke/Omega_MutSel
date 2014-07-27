import sys
if (len(sys.argv) != 5):
    print "\n\nUsage: python run_siminf.py <rep> <treefile> <simdir> <cpu> \n."
    sys.exit()
rep = sys.argv[1]
treefile = sys.argv[2]
simdir = sys.argv[3]
cpu = sys.argv[4]
sys.path.append(simdir)
from functions_simandinf import *

# set up output sequence and parameter files
seqfile = "seqs"+str(rep)+".fasta"
freqfile = "codonFreqs" + str(rep)+".txt"
outfile = "params"+str(rep)+".txt"

# More important parameters
seqlength = 50000
omega = rn.uniform(0.01, 0.99) # omega
kappa = rn.uniform(1.0, 5.0) # kappa
#lambda_ = rn.uniform(0.5, 3.5) # sets strength of selection, effectively. This parameter will be the stddev for the normal distribution from which we draw scaled selection coefficients. Larger stddev = larger fitness differences among amino acids.


# Simulate
print "simulating"
mu=1. # there isn't such a parameter for gy94 so set to 1.
#f_data, gc_content = setFreqs(freqfile, lambda_, 0.0, 1.0) # last 2 args are gc min, gc max
f_equal = np.zeros(61)
f_equal[f_equal == 0.] = 1./61.

simulate(f_equal, seqfile, treefile, mu, kappa, seqlength, omega) # omega is last argument. when None, sim via mutsel


mu_dict = {'AT':mu, 'AC':mu, 'AG':mu*kappa, 'CG':mu, 'CT':mu*kappa, 'GT':mu}
derivedw = deriveOmega(f_equal, mu_dict)

mlw, k = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, kappa, f_equal)

outf = open(outfile, 'w')
outf.write(str(rep) + '\t' + str(kappa) + '\t' + str(derivedw) + '\t' + str(omega) + '\t' + str(mlw) +  '\n')
outf.close()



