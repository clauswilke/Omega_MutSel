# SJS. stephanie.spielman@gmail.com
# Code for using Jesse Bloom's NP amino acid preferene data, with mutation rates either from NP (Bloom 2014) or yeast.

######## Input parameters ########
import sys
if (len(sys.argv) != 6):
    print "\n\nUsage: python run_np.py <rep> <simdir> <cpu> <dataset> <batchfile>\n."
    sys.exit()
rep = sys.argv[1]         # which rep we're on, for saving files. needs to run from 1-498, since 498 sites.
simdir = sys.argv[2]      # directory of simulation library
cpu = sys.argv[3]         # hyphy can use
dataset = sys.argv[4]     # either np, yeast, or polio. determines the mutation scheme and eq freqs
batchfile = sys.argv[5]   # hyphy batchfile name

sys.path.append(simdir)
from functions_simandinf import *

seqlength = 500
treefile = 'tree.tre'


# output files
seqfile   = "seqs"+str(rep)+".fasta"
paramfile = "params"+str(rep)+".txt"

# Set up mutation rates, frequencies, hyphy batchfile name based on the dataset specified
if dataset == 'np':
    mu_dict = {'AG':2.4e-5, 'TC':2.4e-5, 'GA':2.3e-5, 'CT':2.3e-5, 'AC':9.0e-6, 'TG':9.0e-6, 'CA':9.4e-6, 'GT':9.4e-6, 'AT':3.0e-6, 'TA':3.0e-6, 'GC':1.9e-6, 'CG':1.9e-6}
    nuc_freqs = {'pi_a': '0.238390045684', 'pi_c': '0.264432813113', 'pi_g': '0.253719016736', 'pi_t': '0.243458124467'}

elif dataset == 'yeast':
    mu = 1.67e-10 # this is the mean per generation per nucleotide mutation rate. 
    mu_dict = {'AG':0.144/2*mu, 'TC':0.144/2*mu, 'GA':0.349/2*mu, 'CT':0.349/2*mu, 'AC':0.11/2*mu, 'TG':0.11/2*mu, 'CA':0.182/2*mu, 'GT':0.182/2*mu, 'AT':0.063/2*mu, 'TA':0.063/2*mu, 'GC':0.152/2*mu, 'CG':0.152/2*mu}
    nuc_freqs = {'pi_a': '0.324036570484', 'pi_c': '0.174651194384', 'pi_g': '0.161406715511', 'pi_t': '0.33990551962'}

elif dataset == 'polio':
    mu_dict = {'AG':2.495e-5, 'TC':6.886e-05, 'GA':1.259e-04, 'CT':2.602e-04, 'AC':1.721e-06, 'TG':1.177e-06, 'CA':9.072e-06, 'GT':1.472e-05, 'AT':3.812e-06, 'TA':3.981e-06, 'GC':6.301e-06, 'CG':1.633e-06}
    nuc_freqs = {'pi_a': '0.356355793562', 'pi_c': '0.12632952606', 'pi_g': '0.0663827833931', 'pi_t': '0.450931896985'}

else:
    raise AssertionError("Dataset has to be np, yeast, or polio.")


# Create treefile
write_treefile(treefile)

# Read in equilibrium frequencies and determine entropy
eq_codon_freqs = np.loadtxt(dataset + "_codon_eqfreqs.txt")
codon_freqs_true = eq_codon_freqs[ int(rep) - 1 ]
codon_freqs_true_dict = dict(zip(codons, codon_freqs_true))
entropy = calc_entropy(codon_freqs_true)


# Simulate according to MutSel model along phylogeny
print "Simulating"
simulate(codon_freqs_true, seqfile, treefile, mu_dict, seqlength)


# Derive omega from eq freqs
print "Deriving omega from equilibrium codon frequencies"
dnds = derive_dnds(codon_freqs_true_dict, mu_dict)


# Run hyphy and save omega, kappa, omega error. Note we only run free kappa, so just a single call to hyphy.
for pi in nuc_freqs:
    sed_pi = "sed -i 's/"+pi+"/"+nuc_freqs[pi]+"/g' matrices_raw.mdl"
    assert(sed_pi == 0), "sed_pi didn't work"
shutil.copy('matrices_raw.mdl', 'matrices.mdl')
shutil.copy(seqfile, "temp.fasta")
setup_tree = subprocess.call("cat "+treefile+" >> temp.fasta", shell = True)
assert(setup_tree == 0), "couldn't add tree to hyphy infile"
runhyphy = subprocess.call( "./HYPHYMP globalDNDS_test.bf CPU="+cpu+" > hyout.txt", shell = True)
assert (runhyphy == 0), "hyphy fail"
w, k = parse_output_GY94("hyout.txt")
error = (dnds - w)/dnds

# Finally, save results
outstring_params = rep + '\t' + str(entropy) + '\t' + str(dnds)
outf = open(paramfile, 'w')
for f in fspecs:
    y =  fspecs.index(f)
    outf.write( outstring_params + '\t' + f + '\t' +  str(w) + '\t' + str(error) + '\t' + str(kappa) + '\n')
outf.close()   





