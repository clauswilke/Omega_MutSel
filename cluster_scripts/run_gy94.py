# SJS.
# Generic code for simulating and deriving omega via math, hyphy ML.
# Be sure to cp src/ directory (simulator), hyphy files, and the functions_simandinf.py script into working directory
# NOTE: very little to no sanity checking for input args

######## Input parameters ########
import sys
if (len(sys.argv) != 7):
    print "\n\nUsage: python run_siminf.py <rep> <cpu> <simdir> \n."
    sys.exit()
rep = sys.argv[1]
cpu = sys.argv[2]
simdir = sys.argv[3]
sys.path.append(simdir)
from functions_simandinf import *
seqlength = 10000
kappa = rn.uniform(1.0, 5.0)
lambda_ = rn.uniform(0.5, 2.0)
omega = rn.uniform(0.01, 0.99)



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
f_data, num_pref_aa, gc_content = setFreqs(freqfile, lambda_, 0.0, 1.0) # last 2 args are gc min, gc max
simulate(f_data, seqfile, treefile, mu, kappa, seqlength, omega) # omega is last argument. when None, sim via mutsel

f_equal = np.zeros(61)
f_equal[f_equal == 0.] = 1./61.

mlw1, mlk1 = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, 1., f_data)
mlw2, mlk2 = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, 'free', f_data)
mlw3, mlk3 = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, 1., f_equal)
mlw4, mlk4 = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, 'free', f_equal)


outf = open(outfile, 'w')
outf.write(str(rep) + '\t' + str(lambda_) + '\t' + str(kappa) + '\t' + str(omega) + '\t' + str(mlw1) + '\t' + str(mlw2) + '\t' + str(mlw3) + '\t' + str(mlw4) + '\t' + str(mlk2) + '\t' + str(mlk4) + '\n')
outf.close()



