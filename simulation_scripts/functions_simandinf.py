## SJS. Functions that accompany run_siminf.py, etc.
# NOTE: to use simulation library, must cp the src/ directory (*not* contents, the whole directory!) into wdir.

import os
import re
import sys
import shutil
import subprocess
import numpy as np
from random import randint, shuffle
from scipy import linalg

# Simulation code
from misc import *
from newick import *
from stateFreqs import *
from matrixBuilder import *
from evolver import *

# Globals
ZERO = 1e-8
amino_acids  = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
codons=["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
codon_dict = {"AAA":"K", "AAC":"N", "AAG":"K", "AAT":"N", "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T", "AGA":"R", "AGC":"S", "AGG":"R", "AGT":"S", "ATA":"I", "ATC":"I", "ATG":"M", "ATT":"I", "CAA":"Q", "CAC":"H", "CAG":"Q", "CAT":"H", "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P", "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R", "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L", "GAA":"E", "GAC":"D", "GAG":"E", "GAT":"D", "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A", "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G", "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V", "TAC":"Y", "TAT":"Y", "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S", "TGC":"C", "TGG":"W", "TGT":"C", "TTA":"L", "TTC":"F", "TTG":"L", "TTT":"F"}
genetic_code = [["GCA", "GCC", "GCG", "GCT"], ["TGC","TGT"], ["GAC", "GAT"], ["GAA", "GAG"], ["TTC", "TTT"], ["GGA", "GGC", "GGG", "GGT"], ["CAC", "CAT"], ["ATA", "ATC", "ATT"], ["AAA", "AAG"], ["CTA", "CTC", "CTG", "CTT", "TTA", "TTG"], ["ATG"], ["AAC", "AAT"], ["CCA", "CCC", "CCG", "CCT"], ["CAA", "CAG"], ["AGA", "AGG", "CGA", "CGC", "CGG", "CGT"] , ["AGC", "AGT", "TCA", "TCC", "TCG", "TCT"], ["ACA", "ACC", "ACG", "ACT"], ["GTA", "GTC", "GTG", "GTT"], ["TGG"], ["TAC", "TAT"]]
family_size = [4., 2., 2., 2., 2., 4., 2., 3., 2., 6., 1., 2., 4., 2., 6., 6., 4., 4., 1., 2.] # alphabetical according to amino acids.


######################################################################################################################################

########################################################## SIMULATION ################################################################
def simulate(f, seqfile, tree, mu_dict, length):
    ''' Simulate single partition according homogeneous mutation-selection model.
    '''
    try:
        my_tree = readTree(file = tree, flags = False)
    except:
        my_tree = readTree(tree = tree, flags = False) 
          
    model = Model()
    params = {'stateFreqs':f, 'alpha':1.0, 'beta':1.0, 'mu': mu_dict}
    model.params = params
    mat = mutSel_MatrixBuilder(model)
    model.Q = mat.buildQ()
    partitions = [(length, {"rootModel":model})]        
    myEvolver = Evolver(partitions, "rootModel" )
    myEvolver.simulate(my_tree)
    myEvolver.writeSequences(outfile = seqfile)
######################################################################################################################################


################################### FUNCTIONS TO SET UP SCALED SEL COEFFS, CODON FREQUENCIES #########################################
def set_codon_freqs(sd, freqfile, bias):
    ''' Returns equilibrium codon frequencies, entropy, and gc content. Also saves codon frequencies to file. 
        We simulate codon frequencies according a Boltzmann distribution (eg, SellaHirsh 2005, eq 7). Thus, we must simulate values for the exponent. We draw them from a normal distribution.
        IMPORTANTLY, that expression (eq 7) for equilibrium frequencies applies only when the mutation matrix is symmetric.
        
        We draw 20 exponents to determine equilibrium codon frequencies. Each codon in an amino acid family is assigned the same exponent.
        These equilibrium frequencies are those which would exist in the absense of bias. We implement codon bias as follows - 
            Assume eq. freq of a codon in a given amino acid codon family is P, in the absense of codon bias.
            We randomly select one codon to be preferred, and the rest are non-preferred. The preferred codon has a frequency
            of P*(1+(k-1)\lambda) and the nonpreferred codons each have a frequency of P*(1-\lambda), where k=family size. 
            \lambda=0 means no codon bias and \lambda=1 means complete codon bias (there exists only one preferred codon).
    '''

    # Determine equilibrium codon frequencies in the absense of codon bias
    codon_freqs, codon_freqs_dict = calc_codon_freqs_nobias( np.random.normal(loc=0., scale=sd, size=20) )
    
    # If needed, implement codon bias
    if bias != 0.:
        codon_freqs, codon_freqs_dict = implement_bias(codon_freqs_dict, bias)
    
    # Save codon equilibrium frequencies to file  
    np.savetxt(freqfile, codon_freqs)
    
    # Determine gc content
    fobj = UserFreqs(by = 'codon', freqs = codon_freqs_dict)
    nuc_freq = fobj.calcFreqs(type = 'nuc')
    gc = nuc_freq[1] + nuc_freq[2]
    
    # Determine entropy
    entropy = calc_entropy(codon_freqs)

    return codon_freqs, codon_freqs_dict, gc, entropy
    
    
   

def calc_codon_freqs_nobias(exps):  
    ''' Calculate the equilibrium (steady-state) codon frequencies in the absence of codon bias. In other words, all synonymous codons will have the same frequency. '''
    codon_freqs = np.ones(61)
    for i in range(20):   
        for syn in genetic_code[i]:
            codon_index = codons.index(syn)
            codon_freqs[codon_index] = np.exp( exps[i] )
    codon_freqs /= np.sum(codon_freqs)
    assert(-1.*ZERO < np.sum(codon_freqs) - 1. < ZERO), "Bad no-bias codon frequencies."
    return codon_freqs, dict(zip(codons, codon_freqs))
    
         
    
def implement_bias(codon_freqs_dict, bias):
    ''' Implement codon bias. Argument raw_freqs_dict corresponds to what codon frequencies would be in *absence* of bias.
        The bias term ranges from [0,1], where 0=no bias and 1=complete bias.
    '''
    
    for i in range(20):       
        syn_codons = genetic_code[i]
        shuffle(syn_codons) # randomize this, otherwise the preferred codon will always be the first one alphabetically
        k = family_size[i]
        first=True
        for syn in syn_codons:
            if first:
                codon_freqs_dict[syn] *= (1. + (k-1.)*float(bias))
                first=False
            else:
                codon_freqs_dict[syn] *= (1. - float(bias))
    
    # Now need to convert this dictionary back to a list, properly ordered.
    codon_freqs = np.zeros(61)
    for i in range(61):
        codon_freqs[i] = codon_freqs_dict[codons[i]]
    assert(-1.*ZERO <  np.sum(codon_freqs) - 1. < ZERO), "Bad codon frequencies"  
    return codon_freqs, codon_freqs_dict



    
def calc_entropy(f):
    return -1. * np.sum ( f[f > ZERO] * np.log(f[f > ZERO]) )    
######################################################################################################################################





    

################################################# DN/DS DERIVATION FUNCTIONS #########################################################

def derive_dnds(codon_freqs_dict, mu_dict):
    ''' By default, calculate dS. If no bias and symmetric mutation rates, it will be 1 anyways at virtually no computational cost... '''
    
    numer_dn = 0.; denom_dn = 0.;
    numer_ds = 0.; denom_ds = 0.;

    for codon in codon_freqs_dict:
        if codon_freqs_dict[codon] > ZERO:  
        
            rate, sites = calc_nonsyn_paths(codon, codon_freqs_dict, mu_dict)
            numer_dn += rate
            denom_dn += sites
    
            rate, sites = calc_syn_paths(codon, codon_freqs_dict, mu_dict)
            numer_ds += rate
            denom_ds += sites
    
    assert( denom_dn != 0. and denom_ds != 0.), "Omega derivation, with bias, indicates no evolution, maybe?"
    return (numer_dn/denom_dn)/(numer_ds/denom_ds)
    


     
def calc_nonsyn_paths(source, cfreqs, mu_dict):
    rate = 0.
    sites = 0.
    source_freq = cfreqs[source]
    for target in codons:
        diff = get_nuc_diff(source, target) # only consider single nucleotide differences since are calculating instantaneous.
        if codon_dict[source] != codon_dict[target] and cfreqs[target] > ZERO and len(diff) == 2:
            rate  += calc_fixation_prob( source_freq, cfreqs[target] ) * mu_dict[diff]
            sites += mu_dict[diff]
    rate  *= source_freq
    sites *= source_freq
    return rate, sites
    

def calc_syn_paths(source, cfreqs, mu_dict):
    rate = 0.
    sites = 0.
    source_freq = cfreqs[source]
    for target in codons:
        diff = get_nuc_diff(source, target) # only consider single nucleotide differences since are calculating instantaneous.
        if codon_dict[source] == codon_dict[target] and cfreqs[target] > ZERO and len(diff) == 2:
            rate  += calc_fixation_prob( source_freq, cfreqs[target] ) * mu_dict[diff]
            sites += mu_dict[diff]
    rate  *= source_freq
    sites *= source_freq
    return rate, sites


def get_nuc_diff(source, target):
    diff = ''
    for i in range(3):
        if source[i] != target[i]: 
            diff += source[i]+target[i]
    return diff

    
def calc_fixation_prob(fi, fj):
    ''' fi is source and fj is target '''
    if fi == fj:
        return 1.
    elif fi == 0. or fj == 0.:
        return 0.
    else:
        return (np.log(fj) - np.log(fi)) / (1 - fi/fj)    
######################################################################################################################################




#################################################### HYPHY FUNCTIONS #################################################################
def run_hyphy_convergence(seqfile, treefile, cpu, kappa):
    ''' Run hyphy with kappa as true value and equal frequencies, to demonstrate convergence. '''
    
    # Set up sequence file with tree
    shutil.copy(seqfile, "temp.fasta")
    setup_tree = subprocess.call("cat "+treefile+" >> temp.fasta", shell = True)
    assert(setup_tree == 0), "couldn't add tree to hyphy infile"
            
    # Set up kappa in the matrices file
    runsedkappa = subprocess.call("sed 's/k/"+str(kappa)+"/g' matrices_raw.mdl > matrices.mdl", shell=True)
    assert(runsedkappa == 0), "couldn't set up kappa"
   
    # Run hyphy.
    runhyphy = subprocess.call( "./HYPHYMP globalDNDS_equalfreq_truekappa.bf CPU="+cpu+" > hyout.txt", shell = True)
    assert (runhyphy == 0), "hyphy fail"
    
    w, k = parse_output_GY94("hyout.txt")
    return w
    
    
    
    
    
def run_hyphy(seqfile, treefile, cpu, kappa, fspecs):
    ''' Run global omega inference according to GY94. The M0 model. 
        By default, conducts inferences for 4 sets of frequencies (equal, F61, F3x4, CF3x4) across a single kappa specification.
        DO NOT CHANGE FILE NAMES. THEY ARE HARDCODED HERE AND IN THE HYPHY BATCHFILE.
    '''
    
  
    # Set up sequence file with tree
    shutil.copy(seqfile, "temp.fasta")
    setup_tree = subprocess.call("cat "+treefile+" >> temp.fasta", shell = True)
    assert(setup_tree == 0), "couldn't add tree to hyphy infile"
            
        
    # Set up kappa in the matrices file
    if kappa != 'free':
        sedkappa = "sed 's/k/"+str(kappa)+"/g' matrices_raw.mdl > matrices.mdl"
        runsedkappa = subprocess.call(sedkappa, shell=True)
        assert(runsedkappa == 0), "couldn't set up kappa"
    else:
        shutil.copy('matrices_raw.mdl', 'matrices.mdl')

    # Run hyphy.
    runhyphy = subprocess.call("./HYPHYMP globalDNDS.bf CPU="+cpu, shell = True)
    assert (runhyphy == 0), "hyphy fail"
    
    # Retrieve omega, kappa MLEs from the hyout files Produces 4 output files, named suffix = {equal_, f61_, f3x4_, cf3x4_} + hyout.txt 
    omegas = np.zeros(4)
    kappas = np.zeros(4)
    count = 0
    for suffix in fspecs:
        file = suffix + '_hyout.txt'  
        mlw, mlk = parse_output_GY94(file)
	omegas[count] = mlw
        if mlk is None:
            kappas[count] = kappa
	else:
            kappas[count] = mlk
        count += 1
    return omegas, kappas



    
def run_hyphy_np(batchfile, seqfile, treefile, cpu, kappa, fspecs):
    ''' Run global omega inference according to GY94. The M0 model. FOR THE NUCLEOPROTEIN FREQUENCIES ONLY!!
        By default, conducts inferences for 8 sets of frequencies (equal, null, F61 site, F61 global, F3x4 site, F3x4 global, CF3x4 site, CF3x4 global) across a single kappa specification.
        DO NOT CHANGE FILE NAMES. THEY ARE HARDCODED HERE AND IN THE HYPHY BATCHFILE.
    '''
    
  
    # Set up sequence file with tree
    shutil.copy(seqfile, "temp.fasta")
    setup_tree = subprocess.call("cat "+treefile+" >> temp.fasta", shell = True)
    assert(setup_tree == 0), "couldn't add tree to hyphy infile"
            
        
    # Set up kappa in the matrices file
    if kappa != 'free':
        sedkappa = "sed 's/k/"+str(kappa)+"/g' matrices_raw.mdl > matrices.mdl"
        runsedkappa = subprocess.call(sedkappa, shell=True)
        assert(runsedkappa == 0), "couldn't set up kappa"
    else:
        shutil.copy('matrices_raw.mdl', 'matrices.mdl')

    # Run hyphy.
    runhyphy = subprocess.call("./HYPHYMP " + batchfile + " CPU="+cpu, shell = True)
    assert (runhyphy == 0), "hyphy fail"
    
    # Retrieve omega, kappa MLEs from the hyout files Produces 8 output files, names of which are hardcoded!!
    omegas = np.zeros(8)
    kappas = np.zeros(8)
    count = 0
    for suffix in fspecs:
        file = suffix + '_hyout.txt'  
        mlw, mlk = parse_output_GY94(file)
        omegas[count] = mlw
        if mlk is None:
            kappas[count] = kappa
	else:
            kappas[count] = mlk
        count += 1
    return omegas, kappas
     
    
def parse_output_GY94(file):
    hyout = open(file, 'r')
    hylines = hyout.readlines()
    hyout.close()
    hyphy_w = None
    hyphy_k = None
    for line in hylines:
        findw = re.search("^w=(\d+\.*\d*)", line)
        if findw:
            hyphy_w = float(findw.group(1))
        findk = re.search("^k=(\d+\.*\d*)", line)
        if findk:
            hyphy_k = float(findk.group(1))
    assert(hyphy_w is not None)
    return hyphy_w, hyphy_k


def array_to_hyphy_freq(f):
    ''' Convert codon frequencies to a hyphy frequency string. '''
    hyphy_f = "{"
    for freq in f:
        hyphy_f += "{"
        hyphy_f += str(freq)
        hyphy_f += "},"
    hyphy_f = hyphy_f[:-1]
    hyphy_f += "};"
    return hyphy_f



######################################################################################################################################
