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
def set_codon_freqs(sd, freqfile, aafile, codonfile, bias):
    ''' Returns equilibrium codon frequencies, entropy, and gc content. Also saves codon frequencies, aa coefficients, and codon coefficients to file. '''
        
    # Draw amino acid ssc values and assign randomly to amino acids. NOTE: the sd argument is either the sd or it is the pre-determined amino acid selection coefficients.
    if type(sd) is float:
        aminos = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
        shuffle(aminos)  # To randomly assign coefficients, shuffle aminos acids.
        aa_coeffs = dict(zip(aminos, draw_amino_coeffs(sd)))
    elif type(sd) is dict:
        aa_coeffs = sd

    # Convert amino acid coefficients to codon coefficients
    codon_coeffs = aa_to_codon_coeffs(aa_coeffs, bias)


    # Convert codon coefficients to steady-state frequencies
    codon_freqs = codon_coeffs_to_freqs(codon_coeffs)
    codon_freqs_dict = dict(zip(codons, codon_freqs))
            
     
    # Save frequencies and selection coeffs to file  
    np.savetxt(freqfile, codon_freqs)
    np.savetxt(aafile, [aa_coeffs[key] for key in sorted(aa_coeffs)])
    np.savetxt(codonfile, [codon_coeffs[key] for key in sorted(codon_coeffs)])
    
    # Determine gc content
    fobj = UserFreqs(by = 'codon', freqs = codon_freqs_dict)
    nuc_freq = fobj.calcFreqs(type = 'nuc')
    gc = nuc_freq[1] + nuc_freq[2]
    
    # Determine entropy
    entropy = calc_entropy(codon_freqs)

    return codon_freqs, codon_freqs_dict, gc, entropy
    

def draw_amino_coeffs(sd):
    ssc_values = np.random.normal(loc=0., scale=sd, size=20) # note: default np fxn is mean=0, sd=1.
    ssc_values[0] = 0.
    return ssc_values


def aa_to_codon_coeffs(aa_coeffs, bias):
    ''' Assign amino acid selection coefficients to codons. Note that if bias=0, then there's no codon bias. '''
    codon_coeffs = {}
    for aa in aa_coeffs:
        syn_codons = genetic_code[ amino_acids.index(aa) ]
        shuffle(syn_codons) # randomize otherwise the preferred will be the first one alphabetically
        k = float(len(syn_codons) - 1.)
        first=True
        for syn in syn_codons:
            if first:
                codon_coeffs[syn] = aa_coeffs[aa] + bias
                first=False
            else:
                codon_coeffs[syn] = aa_coeffs[aa] - bias/k
    return codon_coeffs



def codon_coeffs_to_freqs(codon_coeffs):
    codon_freqs = np.zeros(61)
    count = 0
    for codon in codons:
        codon_freqs[count] = np.exp( codon_coeffs[codon] )
        count += 1
    codon_freqs /= np.sum(codon_freqs)                   
    assert(np.sum(codon_freqs) - 1.0 < ZERO), "codon_freq doesn't sum to 1 in codon_coeffs_to_freqs"
    return codon_freqs
    
def calc_entropy(f):
    return -1. * np.sum ( f[f > ZERO] * np.log(f[f > ZERO]) )    
######################################################################################################################################





    

################################################# OMEGA DERIVATION FUNCTIONS #########################################################

def derive_omega(codon_freqs_dict, mu_dict, bias):
    ''' NOTE: bias=0 assumes dS = 1, which follows from a symmetric mutation rate and no codon bias (synonymous have same fitness). 
        When bias is 0, then dS is also computed. Mutation matrix still symmetric though.
    '''
    
    numer_dn = 0.; denom_dn = 0.;
    numer_ds = 0.; denom_ds = 0.;

    for codon in codon_freqs_dict:
        if codon_freqs_dict[codon] > ZERO:  
        
            rate, sites = calc_nonsyn_paths(codon, codon_freqs_dict, mu_dict)
            numer_dn += rate
            denom_dn += sites
    
            if bias:
		rate, sites = calc_syn_paths(codon, codon_freqs_dict, mu_dict)
                numer_ds += rate
                denom_ds += sites
    
    if bias is False:
        assert( denom_dn != 0. ), "Omega derivation indicates no evolution, maybe?"
        return numer_dn/denom_dn
        
    else:
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
        kappas[count] = mlk
        count += 1
    return omegas, kappas



    
def run_hyphy_np(seqfile, treefile, cpu, kappa, fspecs):
    ''' Run global omega inference according to GY94. The M0 model. FOR THE NUCLEOPROTEIN FREQUENCIES ONLY!!
        By default, conducts inferences for 7 sets of frequencies (equal, F61 site, F61 global, F3x4 site, F3x4 global, CF3x4 site, CF3x4 global) across a single kappa specification.
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
    runhyphy = subprocess.call("./HYPHYMP globalDNDS_np.bf CPU="+cpu, shell = True)
    assert (runhyphy == 0), "hyphy fail"
    
    # Retrieve omega, kappa MLEs from the hyout files Produces 7 output files, names of which are hardcoded!!
    omegas = np.zeros(7)
    kappas = np.zeros(7)
    count = 0
    for suffix in fspecs:
        file = suffix + '_hyout.txt'  
        mlw, mlk = parse_output_GY94(file)
        omegas[count] = mlw
        kappas[count] = mlk
        count += 1
    return omegas, kappas
     
    
def parse_output_GY94(file):
    hyout = open(file, 'r')
    hylines = hyout.readlines()
    hyout.close()
    hyphy_w = None
    hyphy_k = 1000 # so R will read column as numeric
    for line in hylines:
        findw = re.search("^w=(\d+\.*\d*)", line)
        if findw:
            hyphy_w = float(findw.group(1))
        findk = re.search("^k=(\d+\.*\d*)", line)
        if findk:
            hyphy_k = float(findk.group(1))
    assert(hyphy_w is not None)
    return hyphy_w, hyphy_k
######################################################################################################################################

