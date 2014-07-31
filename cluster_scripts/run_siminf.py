# SJS.
# Generic code for simulating and deriving omega via math, hyphy ML.
# Be sure to cp src/ directory (simulator), hyphy files, and the functions_simandinf.py script into working directory
# NOTE: very little to no sanity checking for input args

######## Input parameters ########
import sys
if (len(sys.argv) != 8):
    print "\n\nUsage: python run_siminf.py <rep> <treefile> <simdir> <cpu> <kappa> <gcbias> <selection>\n."
    sys.exit()
rep = sys.argv[1]
treefile = sys.argv[2]
simdir = sys.argv[3]
cpu = sys.argv[4]
kappa_runif = bool(int(sys.argv[5]))
asym = bool(int(sys.argv[6])) # 0 or 1
selection = bool(int(sys.argv[7])) # 0 or 1. 0 = no selection, so equal codon freqs. 1 = selection, so boltzmann amino acid, then convert to codon.
sys.path.append(simdir)
from functions_simandinf import *



# set up output sequence and parameter files
seqfile = "seqs_"+str(rep)+".fasta"
freqfile = "codonFreqs" + str(rep)+".txt"
outfile = "params"+str(rep)+".txt"



# Parameters
f_equal = np.zeros(61)
f_equal[f_equal == 0.] = 1./61.
seqlength = 500000
lambda_ = 1. #rn.uniform(0.5, 3.5) # sets strength of selection, effectively. This parameter will be the stddev for the normal distribution from which we draw scaled selection coefficients. Larger stddev = larger fitness differences among amino acids.
# Default mutation rates
mu = 1e-5
c = 1.0
kappa = 1.0

# Change defaults as needed
if kappa_runif:
    kappa = rn.uniform(1.0, 5.0)
if asym:
    c = rn.uniform(1,50)
    if rn.randint(0,1): #false (producing c>1) -> AT bias. true (producing c<1) -> GC bias.
        c =1/c
mu_dict = {'AT': mu, 'TA':mu, 'CG': mu, 'GC':mu, 'AC': mu, 'TG':mu, 'CA':c*mu, 'GT':c*mu, 'AG': kappa*mu, 'TC':kappa*mu, 'GA':c*kappa*mu, 'CT':c*kappa*mu}




# Simulate
print "freqs"
if selection and asym:
    f_true = setFreqsAsym(lambda_, mu_dict)
else:
    f_true = f_equal
    gc_true = 0.513661202
print "simulating"
simulate(f_true, seqfile, treefile, mu_dict, seqlength, None) # omega is last argument. when None, sim via mutsel



# Get empirical codon frequency and gc content, as well as codon entropies
f_data_obj = ReadFreqs(file = seqfile, by = 'codon')
f_data = f_data_obj.calcFreqs(type = 'codon')
nuc_data = f_data_obj.calcFreqs(type = 'nuc')
gc_data = nuc_data[1] + nuc_data[2]
entropy_true = calcCodonEntropy(f_true)
entropy_data = calcCodonEntropy(f_data)



# Derive omega
print "deriving"
derivedw_true = deriveOmega(f_true, mu_dict)
derivedw_data  = deriveOmega(f_data, mu_dict)


# ML
print "ML"
fspecs = {'equal':'globalDNDS_inputf.bf'} #, 'true':'globalDNDS_inputf.bf',  'data': 'globalDNDS_inputf.bf', 'f3x4':'globalDNDS_f3x4.bf', 'cf3x4':'globalDNDS_cf3x4.bf'}
kspecs = {kappa:'kappa_true'}#, 1.0:'kappa_one', 'free':'kappa_free'}

#common_out_string = rep + '\t' + str(seqlength) + '\t' + str(mu) + '\t' + str(kappa) + '\t' + str(bias) + '\t' + str(lambda_) + '\t' + str(gc_true) + '\t' + str(gc_data) + '\t' + str(entropy_true) + '\t' + str(entropy_data) + '\t' + str(derivedw_true) + '\t' + str(derivedw_data)

outf = open(outfile, 'w')
for freqspec in fspecs:
    try:
        hyfreq = eval('f_'+freqspec)
    except:
        hyfreq = None
    for kapspec in kspecs:
        mlw, mlk = runhyphy(fspecs[freqspec], "GY94", seqfile, treefile, cpu, kapspec, hyfreq)
        #w_err = (derivedw_true - mlw) / derivedw_true
        outf.write(str(derivedw_true) + '\t' + str(derivedw_data) + '\t' + str(mlw) + '\t' + str(c) + '\t' + str(kappa) + '\n')
        #outf.write(common_out_string + '\t' + 'freq_'+str(freqspec) + '\t' + str(kspecs[kapspec]) + '\t' + str(mlw) + '\t' + str(w_err) + '\t' + str(mlk) + '\n')
outf.close()



