# SJS.
# Generic code for simulating and deriving omega via math, hyphy ML.
# Be sure to cp src/ directory (simulator), hyphy files, and the functions_simandinf.py script into working directory
# NOTE: very little to no sanity checking for input args

######## Input parameters ########
import sys
if (len(sys.argv) != 4):
    print "\n\nUsage: python run_siminf.py <rep> <bl> <simdir> \n."
    sys.exit()
rep = sys.argv[1]
bl = float(sys.argv[2])
simdir = sys.argv[3]
sys.path.append(simdir)
from functions_simandinf import *
seqlength = 100000
#if int(rep) % 2 == 1:
#    kappa = rn.uniform(1.0, 5.0)
#else:
#    kappa = 1.

lambda_ = rn.uniform(0.5, 2.0)
omega = rn.uniform(0.01, 0.99)

# to get a params file named according to bl without decimals in it.
hi = str(str(bl).count('0')) # only works if the different branchlengths used are of different orders of magnitude.


# Write tree given bl specifications
treefile = "tree.tre"
treef = open(treefile, 'w')
treef.write("(t1:" + str(bl) + ", t2:" + str(bl) + ");")
treef.close()

# set up output sequence and parameter files
seqfile = "seqs"+str(rep)+".fasta"
freqfile = "codonFreqs" + str(rep)+".txt"
outfile = "params"+str(rep)+"_"+hi+".txt"

# Simulate
print "simulating"
mu=1. # there isn't such a parameter for gy94 so set to 1.
f_data, num_pref_aa, gc_content = setFreqs(freqfile, lambda_, 0.0, 1.0) # last 2 args are gc min, gc max
simulate(f_data, seqfile, treefile, mu, kappa, seqlength, omega) # omega is last argument. when None, sim via mutsel

f_equal = np.zeros(61)
f_equal[f_equal == 0.] = 1./61.

cpu="1"
mlw1, mlk1 = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, 1., f_data)
mlw2, mlk2 = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, 'free', f_data)
mlw3, mlk3 = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, 1., f_equal)
mlw4, mlk4 = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, 'free', f_equal)

neiw = run_neigojo(seqfile)

outf = open(outfile, 'w')
outf.write(str(rep) + '\t' + str(bl) + '\t' + str(lambda_) + '\t' + str(kappa) + '\t' + str(omega) + '\t' + str(neiw) + '\t' + str(mlw1) + '\t' + str(mlw2) + '\t' + str(mlw3) + '\t' + str(mlw4) + '\t' + str(mlk2) + '\t' + str(mlk4) + '\n')
outf.close()



