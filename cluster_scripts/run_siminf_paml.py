# SJS.
# Generic code for simulating and deriving omega via math, hyphy ML.
# Be sure to cp src/ directory (simulator), paml files, and the functions_simandinf.py script into working directory
# NOTE: very little to no sanity checking for input args

######## Input parameters ########
import sys
if (len(sys.argv) != 6):
    print "\n\nUsage: python run_siminf_paml.py <rep> <mu> <bl> <seqlength> <simdir> \n."
    sys.exit()
rep = sys.argv[1]
mu = float(sys.argv[2])
bl = sys.argv[3]
seqlength = int(sys.argv[4])
simdir = sys.argv[5]
sys.path.append(simdir)
from functions_simandinf import *
if int(rep)%2 == 1:
    kappa = rn.uniform(1.0, 5.0)
else:
    kappa = 1.0
lambda_ = rn.uniform(0.5, 3.0)



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
simulate(f_data, seqfile, treefile, mu, kappa, seqlength, None) # omega is last argument. when None, sim via mutsel

# codon entropies
entropy_data = calcCodonEntropy(f_data)
entropy_equal = 4.11087386417 # known, no need to calculate
entropy_f3x4 = calcCodonEntropy( calc_f3x4(f_data) )

# Our omega
print "deriving"
mu_dict = {'AT':mu, 'AC':mu, 'AG':mu*kappa, 'CG':mu, 'CT':mu*kappa, 'GT':mu}
derivedw = deriveOmega(f_data, mu_dict)

# NG86 omega
write_yn00 = '\t'.join( runpaml_yn00(seqfile) )


common_out_string = rep + '\t' + str(seqlength) + '\t' + str(bl) + '\t' + str(mu) + '\t' + str(kappa) + '\t' + str(gc_content) + '\t' + str(lambda_) + '\t' + str(entropy_data) + '\t' + str(entropy_f3x4) + '\t' + str(derivedw) + '\t' + write_yn00


# ML omegas
print "ML"
ml_estimates = np.zeros(9)
kappas = np.zeros(3)
count = 0
fspecs = {'0':"equal", '2':"f3x4", '3':"data"} # 1/61, f3x4, data
kspecs = {'0':'free', '1':'fixed'} #free, fixed
outf = open(outfile, 'w')
for freqspec in fspecs:
    for kapspec in kspecs:
        mlw = runpaml_codeml(seqfile, freqspec, kapspec, kappa)
        w_err = (derivedw - mlw) / derivedw
        outf.write(common_out_string + '\t' + str(fspecs[freqspec]) + '\t' + str(kspecs[kapspec]) + '\t' + str(mlw) + '\t' + str(w_err) + '\n')
outf.close()




