# SJS. Code specifically for generating data to demonstrate omega convergence.




import sys
# Input parameters and global stuff
if (len(sys.argv) != 5):
    print "\n\nUsage: python run_convergence.py <rep> <mu> <bl> <simdir> \n."
    sys.exit()
rep = sys.argv[1]
mu = float(sys.argv[2])
bl = sys.argv[3]
simdir = sys.argv[4]
sys.path.append(simdir)
from functions_simandinf import *
from random import randint

kappa = rn.uniform(1.0, 5.0)
lambda_ = rn.uniform(0.5, 2.0)

# Random sequence length between 5e2 and 1e6
expon = randint(2,5)
if expon == 2:
    times = randint(5,10)
else:
    times = randint(1,10)
seqlength = int( times * 10**expon )

# Write tree given bl specifications
treefile = "tree.tre"
treef = open(treefile, 'w')
treef.write("(t1:" + str(bl) + ", t2:" + str(bl) + ");")
treef.close()

# set up output sequence and parameter files
freqfile = "codonFreqs" + str(rep) + "_" + str(seqlength) + ".txt"
seqfile = "seqs" + str(rep) + "_" + str(seqlength) + ".fasta"
outfile = "params" + str(rep) + "_" + str(seqlength) + ".txt"


# Now, simulate sequences and infer ML omegas
# Simulate
print "simulating"
f, num_pref_aa, gc_content = setFreqs(freqfile, lambda_, 0.0, 1.0) # last 2 args are gc min, gc max
simulate(f, seqfile, treefile, mu, kappa, seqlength, None) # omega is last argument. when None, sim via mutsel

# Derive
mu_dict = {'AT':mu, 'AC':mu, 'AG':mu, 'CG':mu, 'CT':mu, 'GT':mu}
mu_dict['AG'] = mu_dict['AG'] * kappa
mu_dict['CT'] = mu_dict['CT'] * kappa
derivedw = deriveOmega(f, mu_dict)

# ng86 omega
neiw = runpaml_yn00(seqfile)

# Calculate relative error from derived omega. Not using abs() since plot nicer that way.
err = ( derived_w - neiw )/derived_w

# Save
outf = open(outfile,'w')
outf.write(rep + '\t' + str(seqlength) + '\t' + str(bl) + '\t' + str(mu) + '\t' + str(kappa) + '\t' + str(lambda_) + '\t' + str(derived_w) + '\t' + str(neiw) + '\t' + str(err) + '\n')
outf.close()





