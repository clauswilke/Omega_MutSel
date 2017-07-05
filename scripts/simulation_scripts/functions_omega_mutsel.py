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

import pyvolve


# Globals
ZERO = 1e-8
amino_acids  = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
codons=["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
codon_dict = {"AAA":"K", "AAC":"N", "AAG":"K", "AAT":"N", "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T", "AGA":"R", "AGC":"S", "AGG":"R", "AGT":"S", "ATA":"I", "ATC":"I", "ATG":"M", "ATT":"I", "CAA":"Q", "CAC":"H", "CAG":"Q", "CAT":"H", "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P", "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R", "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L", "GAA":"E", "GAC":"D", "GAG":"E", "GAT":"D", "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A", "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G", "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V", "TAC":"Y", "TAT":"Y", "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S", "TGC":"C", "TGG":"W", "TGT":"C", "TTA":"L", "TTC":"F", "TTG":"L", "TTT":"F"}
genetic_code = [["GCA", "GCC", "GCG", "GCT"], ["TGC","TGT"], ["GAC", "GAT"], ["GAA", "GAG"], ["TTC", "TTT"], ["GGA", "GGC", "GGG", "GGT"], ["CAC", "CAT"], ["ATA", "ATC", "ATT"], ["AAA", "AAG"], ["CTA", "CTC", "CTG", "CTT", "TTA", "TTG"], ["ATG"], ["AAC", "AAT"], ["CCA", "CCC", "CCG", "CCT"], ["CAA", "CAG"], ["AGA", "AGG", "CGA", "CGC", "CGG", "CGT"] , ["AGC", "AGT", "TCA", "TCC", "TCG", "TCT"], ["ACA", "ACC", "ACG", "ACT"], ["GTA", "GTC", "GTG", "GTT"], ["TGG"], ["TAC", "TAT"]]
family_size = [4., 2., 2., 2., 2., 4., 2., 3., 2., 6., 1., 2., 4., 2., 6., 6., 4., 4., 1., 2.] # alphabetical according to amino acids.



def write_treefile(filename):
    ''' write file containing 4-taxon tree'''
    treef = open(filename, 'w')
    treef.write("((t4:0.01,t1:0.01):0.01,(t3:0.01,t2:0.01):0.01);\n")
    treef.close()
    
def set_mu_dict(dataset):
    ''' Set up mutation rate dictionary for NP, yeast, or polio runs.'''
    if dataset == 'np':
        mu_dict = {'AG':2.4e-5, 'TC':2.4e-5, 'GA':2.3e-5, 'CT':2.3e-5, 'AC':9.0e-6, 'TG':9.0e-6, 'CA':9.4e-6, 'GT':9.4e-6, 'AT':3.0e-6, 'TA':3.0e-6, 'GC':1.9e-6, 'CG':1.9e-6}
    elif dataset == 'yeast':
        mu = 1.67e-10 # this is the mean per generation per nucleotide mutation rate. 
        mu_dict = {'AG':0.144/2*mu, 'TC':0.144/2*mu, 'GA':0.349/2*mu, 'CT':0.349/2*mu, 'AC':0.11/2*mu, 'TG':0.11/2*mu, 'CA':0.182/2*mu, 'GT':0.182/2*mu, 'AT':0.063/2*mu, 'TA':0.063/2*mu, 'GC':0.152/2*mu, 'CG':0.152/2*mu}
    elif dataset == 'polio':
        mu_dict = {'AG':2.495e-5, 'TC':6.886e-05, 'GA':1.259e-04, 'CT':2.602e-04, 'AC':1.721e-06, 'TG':1.177e-06, 'CA':9.072e-06, 'GT':1.472e-05, 'AT':3.812e-06, 'TA':3.981e-06, 'GC':6.301e-06, 'CG':1.633e-06}
    else:
        raise AssertionError("Dataset has to be np, yeast, or polio.")
    return mu_dict
    
    
###################################### SIMULATION CODE ##############################################
def simulate(f, seqfile, tree, mu_dict, length):
    ''' Simulate single partition according homogeneous mutation-selection model.
    '''
    
    try:
        my_tree = pyvolve.read_tree(file = tree)
    except:
        my_tree = pyvolve.read_tree(tree = tree) 

    model = pyvolve.Model("MutSel", {'state_freqs':f, 'mu': mu_dict})

    part = pyvolve.Partition(size = length, models = model)    
    e = pyvolve.Evolver(partitions = part, tree = my_tree)
    e(seqfile = seqfile, ratefile = None, infofile = None)
    

################################### FUNCTIONS TO SET UP SCALED SEL fitness, CODON FREQUENCIES #########################################
def set_codon_freqs(sd, freqfile, bias):
    ''' Returns equilibrium codon frequencies, entropy, and gc content. Also saves codon frequencies to file. 
        We simulate values for the amino acid scaled selection coefficients by drawing from N(0,x), where x~U(0,4). Note that x=0 means neutral evolution.
        We convert these values to codon frequencies via SellaHirsh 2005, eq 7 (Boltzmann).
        IMPORTANTLY, that expression (eq 7) for equilibrium frequencies applies only when the mutation matrix is symmetric.
        
        We implement codon bias as follows - 
            Assume ssc of a given amino acid codon family is P, in the absense of codon bias.
            We randomly select one codon to be preferred, and the rest are non-preferred. 
            The preferred codon will be assigned an ssc = P*(1+(k-1)\lambda).
            All nonpreferred codons will be assigned an ssc = P*(1-\lambda), where k=family size. 
            \lambda=0 means no codon bias and \lambda=1 means complete codon bias (there exists only one preferred codon).
    '''

    # Draw amino acid ssc values and assign randomly to amino acids.
    aminos = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    shuffle(aminos)  # To randomly assign coefficients, shuffle aminos acids.
    aa_fitness = dict(zip(aminos, np.random.normal(loc = 0, scale = sd, size = 20)))

    # Convert amino acid coefficients to codon coefficients
    codon_fitness = aa_to_codon_fitness(aa_fitness, bias)

    # Convert codon coefficients to steady-state frequencies
    codon_freqs = codon_fitness_to_freqs(codon_fitness)
    codon_freqs_dict = dict(zip(codons, codon_freqs))
        
    # Save codon equilibrium frequencies to file  
    np.savetxt(freqfile, codon_freqs)
    
    # Determine gc content
    f = pyvolve.CustomFrequencies("codon", freq_dict = codon_freqs_dict)
    nuc_freq = f.compute_frequencies(type = "nucleotide")
    gc = nuc_freq[1] + nuc_freq[2]
    
    # Determine entropy
    entropy = calc_entropy(codon_freqs)

    return codon_freqs, codon_freqs_dict, gc, entropy
    

    
def aa_to_codon_fitness(aa_fitness, lambda_):
    ''' Assign amino acid selection coefficients to codon ssc values. lambda_ is the bias term. '''
    codon_fitness = {}
    for aa in aa_fitness:
        syn_codons = genetic_code[ amino_acids.index(aa) ]
        shuffle(syn_codons) # randomize otherwise the preferred will be the first one alphabetically
        k = float(len(syn_codons) - 1.)
        first=True
        for syn in syn_codons:
            if first:
                codon_fitness[syn] = aa_fitness[aa] + lambda_
                first=False
            else:
                codon_fitness[syn] = aa_fitness[aa] - lambda_
    return codon_fitness         
    
def codon_fitness_to_freqs(codon_fitness):
    codon_freqs = np.zeros(61)
    count = 0
    for codon in codons:
        codon_freqs[count] = np.exp( codon_fitness[codon] )
        count += 1
    codon_freqs /= np.sum(codon_freqs)                   
    assert(-1*ZERO < np.sum(codon_freqs) - 1.0 < ZERO), "codon_freq doesn't sum to 1 in codon_fitness_to_freqs"
    return codon_freqs    


    
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
        
            rate, sites = calc_paths(codon, codon_freqs_dict, mu_dict, 'nonsyn')
            numer_dn += rate
            denom_dn += sites
    
            rate, sites = calc_paths(codon, codon_freqs_dict, mu_dict, 'syn')
            numer_ds += rate
            denom_ds += sites
    
    assert( denom_dn != 0. and denom_ds != 0.), "Omega derivation, with bias, indicates no evolution, maybe?"
    return (numer_dn/denom_dn)/(numer_ds/denom_ds)
    


     
def calc_paths(source, cfreqs, mu_dict, type):
    ''' type is either syn or nonsyn '''
    rate = 0.
    sites = 0.
    source_freq = cfreqs[source]
    for target in codons:
        diff = get_nuc_diff(source, target) # only consider single nucleotide differences since are calculating instantaneous.
        if (type == 'nonsyn' and codon_dict[source] != codon_dict[target]) or (type == 'syn' and codon_dict[source] == codon_dict[target]):
            if cfreqs[target] > ZERO and len(diff) == 2:
                rate  += calc_subst_prob( source_freq, cfreqs[target], mu_dict[diff], mu_dict[diff[1]+diff[0]] )
                sites += mu_dict[diff]
        else:
            continue
    rate  *= source_freq
    sites *= source_freq
    return rate, sites



def get_nuc_diff(source, target):
    diff = ''
    for i in range(3):
        if source[i] != target[i]: 
            diff += source[i]+target[i]
    return diff

    
def calc_subst_prob(pi, pj, mu_ij, mu_ji):
    if pi == 0. or pj == 0.:
        return 0.
    else:
        ratio = (pj*mu_ji)/(pi*mu_ij)
        if ratio == 1.:
            return mu_ij
        else:
            return np.log(ratio)/(1. - 1./ratio) * mu_ij
######################################################################################################################################




#################################################### HYPHY FUNCTIONS #################################################################
def run_hyphy_fequal(seqfile, treefile, cpu, kappa):
    ''' Run hyphy with kappa as true value and equal frequencies, to demonstrate convergence. '''
    
    # Set up sequence file with tree
    shutil.copy(seqfile, "temp.fasta")
    setup_tree = subprocess.call("cat "+treefile+" >> temp.fasta", shell = True)
    assert(setup_tree == 0), "couldn't add tree to hyphy infile"
            
    # Set up kappa in the matrices file
    if kappa != 'free':
        sedkappa = "sed 's/k/"+str(kappa)+"/g' GY94_raw.mdl > GY94_raw.mdl"
        runsedkappa = subprocess.call(sedkappa, shell=True)
        assert(runsedkappa == 0), "couldn't set up kappa"
    else:
        shutil.copy('GY94_raw.mdl', 'GY94.mdl')
   
    # Run hyphy.
    runhyphy = subprocess.call( "HYPHYMP globalDNDS_fequal.bf CPU="+cpu+" > hyout.txt", shell = True)
    assert (runhyphy == 0), "hyphy fail"
    
    lk, w, k = parse_output_GY94("hyout.txt")
    if k is None:
        k = kappa
    return w, k
    




    
def run_hyphy_nyp(batchfile, seqfile, treefile, cpu, fspecs):
    ''' Run global omega inference according to GY94 and GY94_Fnuc. To be used with experimental datasets that use NP, yeast, and polio (hence, nyp) mutation rates.'''

    # Set up sequence file with tree
    shutil.copy(seqfile, "temp.fasta")
    setup_tree = subprocess.call("cat "+treefile+" >> temp.fasta", shell = True)
    assert(setup_tree == 0), "couldn't add tree to hyphy infile"

    # Run hyphy.
    runhyphy = subprocess.call("./HYPHYMP " + batchfile + " CPU="+cpu, shell = True)
    assert (runhyphy == 0), "hyphy fail"
    
    # Retrieve likelihood, omega, kappa from the hyout files, names of which are hardcoded!!
    omegas = np.zeros(len(fspecs))
    kappas = np.zeros(len(fspecs))
    lnliks = np.zeros(len(fspecs)) # log likelihood values
    count = 0
    for suffix in fspecs:
        file = suffix + '_hyout.txt'  
        lnliks[count], omegas[count], kappas[count] = parse_output_GY94(file)
        count += 1
    return lnliks, omegas, kappas
     
    
def parse_output_GY94(file):
    with open(file, 'r') as hyout:
        hylines = hyout.readlines()
    lnlik = None; hyphy_w = None; hyphy_k = None;
    for line in hylines:
        findlk = re.search("^Likelihood Function's Current Value\s+=\s+(-\d+\.*\d*)", line)
        if findlk:
            lnlik = float(findlk.group(1))
        findw = re.search("^w=(\d+\.*\d*)", line)
        if findw:
            hyphy_w = float(findw.group(1))
        findk = re.search("^k=(\d+\.*\d*)", line)
        if findk:
            hyphy_k = float(findk.group(1))
    assert(lnlik is not None),   "Couldn't retrieve log likelihood from hyphy output file."
    assert(hyphy_w is not None), "Couldn't retrieve omega from hyphy output file."
    assert(hyphy_k is not None), "Couldn't retrieve kappa from hyphy output file."
    return lnlik, hyphy_w, hyphy_k



def calc_AIC(k, lnlk):
    return 2.*(float(k) - float(lnlk))

######################################################################################################################################
