# SJS.
# Generic code for simulating and deriving omega via math, hyphy ML.
# Be sure to cp src/ directory (simulator), hyphy files, and the functions_simandinf.py script into working directory
# NOTE: very little to no sanity checking for input args

######## Input parameters ########
import sys
if (len(sys.argv) != 6):
    print "\n\nUsage: python run_siminf.py <rep> <treefile> <simdir> <cpu> <gcbias>\n."
    sys.exit()
rep = sys.argv[1]
treefile = sys.argv[2]
simdir = sys.argv[3]
cpu = sys.argv[4]
gc_bias = bool(int(sys.argv[5])) # 0 or 1
sys.path.append(simdir)
from functions_simandinf import *




# set up output sequence and parameter files
seqfile = "seqs"+str(rep)+".fasta"
freqfile = "codonFreqs" + str(rep)+".txt"
outfile = "params"+str(rep)+".txt"

# Parameters
seqlength = 500000
lambda_ = rn.uniform(0.5, 3.5) # sets strength of selection, effectively. This parameter will be the stddev for the normal distribution from which we draw scaled selection coefficients. Larger stddev = larger fitness differences among amino acids.
mu = 1e-5
kappa = rn.uniform(1.0, 5.0)
if gc_bias:
    bias = rn.uniform(2., 5.)
else:
    bias = 1.0
mu_dict = {'AC': mu*bias, 'CA':mu, 'AG': mu*kappa*bias, 'GA':mu*kappa, 'AT': mu, 'TA':mu, 'CG': mu*bias, 'GC':mu*bias, 'CT': mu*kappa, 'TC':mu*kappa*bias, 'GT': mu, 'TG':mu*bias}
sym_mu_dict = {'AC': mu, 'AG': mu*kappa, 'AT': mu, 'CG': mu, 'CT': mu*kappa,  'GT': mu}




# Simulate
print "simulating"
f_true, gc_content = setFreqs(freqfile, lambda_, 0., 1.) # last 2 args are gc min, gc max
simulate(f_true, seqfile, treefile, mu_dict, kappa, seqlength, None) # omega is last argument. when None, sim via mutsel

f_data_obj = ReadFreqs(file = seqfile, by = 'codon')
f_data = f_data_obj.calcFreqs(type = 'codon')

# Derive omega
print "deriving"
derivedw_true = deriveOmega(f_true, mu_dict, False) #false = not symmetric
derivedw_true_sym = deriveOmega(f_true, sym_mu_dict, True)
derivedw_data  = deriveOmega(f_data, mu_dict, False) #false = not symmetric
derivedw_data_sym = deriveOmega(f_data, sym_mu_dict, True)

# Calculate entropy and setup the f_equal variable
f_equal = np.zeros(61)
f_equal[f_equal == 0.] = 1./61.
entropy_true = calcCodonEntropy(f_true)
entropy_data = calcCodonEntropy(f_data)

# ML
print "ML"
fspecs = {'equal':'globalDNDS_inputf.bf', 'true':'globalDNDS_inputf.bf',  'data': 'globalDNDS_inputf.bf', 'f3x4':'globalDNDS_f3x4.bf', 'cf3x4':'globalDNDS_cf3x4.bf'}
kspecs = {1.0: 'kappa_one', kappa:'kappa_true', 'free':'kappa_free'}

common_out_string = rep + '\t' + str(seqlength) + '\t' + str(mu) + '\t' + str(kappa) + '\t' + str(bias) + '\t' + str(lambda_) + '\t' + str(gc_content) + '\t' + str(entropy_data) + '\t' + str(entropy_true) + '\t' + str(derivedw_true) + '\t' + str(derivedw_true_sym) + '\t' + str(derivedw_data) + '\t' + str(derivedw_data_sym)


outf = open(outfile, 'w')
for freqspec in fspecs:
    try:
        hyfreq = eval('f_'+freqspec)
    except:
        hyfreq = None
    for kapspec in kspecs:
        mlw, mlk = runhyphy(fspecs[freqspec], "GY94", seqfile, treefile, cpu, kapspec, hyfreq)
        #w_err = (derivedw - mlw) / derivedw 
        outf.write(common_out_string + '\t' + 'freq_'+str(freqspec) + '\t' + str(kspecs[kapspec]) + '\t' + str(mlw) + '\n' ) #+ '\t' + str(w_err) + '\t' + str(mlk) + '\n')
outf.close()




