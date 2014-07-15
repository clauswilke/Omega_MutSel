# SJS. Something funky is going on with kappa and it does not behave properly.
# plan: get state frequencies. one simulation with kappa=1, one simulation with kappa=[1,5].
# infer dnds with hyphy (k=1, k=true) and with counting method = 3 hyphy inferences and 2 counting inferences.

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




# Write tree given bl specifications
treefile = "tree.tre"
treef = open(treefile, 'w')
treef.write("(t1:" + str(bl) + ", t2:" + str(bl) + ");")
treef.close()

# set up output sequence and parameter files
seqfile = "seqs"+str(rep)+".fasta"
freqfile = "codonFreqs" + str(rep)+".txt"
outfile = "params"+str(rep)+".txt"

############ Set up state frequencies and also frequencies for hyphy #############
lambda_ = rn.uniform(1.0, 3.0)
f_data, num_pref_aa, gc_content = setFreqs(freqfile, lambda_, 0.0, 1.0) # last 2 args are gc min, gc max
f_equal = np.zeros(61)
f_equal[f_equal == 0.] = 1./61.


################ FIRST WITH TRUE KAPPA = 1.0 ##################

########## Simulate
simulate(f_data, seqfile, treefile, mu, 1.0, seqlength, None) # omega is last argument. when None, sim via mutsel

######### Derive omega with math
derivedw_1 = deriveOmega(f_data, {'AT':mu, 'AC':mu, 'AG':mu, 'CG':mu, 'CT':mu, 'GT':mu})

######### Infer omega with nei-gojobori
neiw_1 = run_neigojo(seqfile)

######### Infer hyphy
mlw_1, mlk = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, 1., f_equal)


################### NOW WITH TRUE KAPPA AS [2,5] #####################


########## Simulate
kappa = rn.uniform(2.0, 5.0)
simulate(f_data, seqfile, treefile, mu, kappa, seqlength, None) # omega is last argument. when None, sim via mutsel

######### Derive omega with math
derivedw_2 = deriveOmega(f_data, {'AT':mu, 'AC':mu, 'AG':mu*kappa, 'CG':mu, 'CT':mu*kappa, 'GT':mu})

######### Infer omega with nei-gojobori
neiw_2 = run_neigojo(seqfile)

######### Infer hyphy
mlw_2_k1, mlk = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, 1., f_equal)
mlw_2_ktrue, mlk = runhyphy("globalDNDS.bf", "GY94", seqfile, treefile, cpu, 'true', f_equal)


outf = open(outfile, 'w')
outf.write(rep + '\t' + str(kappa) + '\t' + str(derivedw_1) + '\t' + str(derivedw_2) + '\t' + str(neiw_1) + '\t' + str(neiw_2) + '\t' + str(mlw_1) + '\t' + str(mlw_2_k1) + '\t' + str(mlw_2_ktrue) + '\n')
outf.close()




