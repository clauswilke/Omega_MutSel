# codonfreqs to 1/61 an kappa to 1. simulate gy94. how does counting now work? does it agree with ml?

######## Input parameters ########
import sys
if (len(sys.argv) != 3):
    print "\n\nUsage: python run_siminf.py <rep> <simdir> \n."
    sys.exit()
rep = sys.argv[1]
simdir = sys.argv[2]
sys.path.append(simdir)
from functions_simandinf import *


# Write tree given bl specifications
treefile = "tree.tre"
treef = open(treefile, 'w')
treef.write("(t1:0.05, t2:0.05);")
treef.close()

# set up output sequence and parameter files
seqfile = "seqs"+str(rep)+".fasta"
freqfile = "codonFreqs" + str(rep)+".txt"
outfile = "params"+str(rep)+".txt"

# Simulate
print "simulating"
omega = rn.uniform(0.01, 0.99)
f_equal = np.zeros(61)
f_equal[f_equal == 0.] = 1./61.
simulate(f_equal, seqfile, treefile, 1., 1., 100000, omega) # omega is last argument. when None, sim via mutsel

# counting method omega
neiw = run_neigojo(seqfile)

# HyPhy omega and kappa, forever and ever
mlw, mlk = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, "1", 'free', f_data)

outf = open(outfile, 'w')
outf.write(str(rep) + '\t' + str(omega) + '\t' + str(neiw) + '\t' + str(mlk) + '\n')
outf.close()


