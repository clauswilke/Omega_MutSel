## SJS. Functions that accompany simulate_and_infer.py and all derivatives.
# NOTE: to use simulation library, must cp the src/ directory (*not* contents, the whole directory!) into wdir.

import os
import re
import sys
import subprocess
import numpy as np
from random import randint


# Simulation code
sys.path.append('src/')
from misc import *
from newick import *
from stateFreqs import *
from matrixBuilder import *
from evolver import *

# Nei-Gojobori code
from mutation_counter import *
from site_counter import *

# Globals
zero = 1e-8
amino_acids  = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
codons=["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
nslist = [['CAA', 'GAA', 'ACA', 'ATA', 'AGA', 'AAC', 'AAT'], ['CAC', 'TAC', 'GAC', 'ACC', 'ATC', 'AGC', 'AAA', 'AAG'], ['CAG', 'GAG', 'ACG', 'ATG', 'AGG', 'AAC', 'AAT'], ['CAT', 'TAT', 'GAT', 'ACT', 'ATT', 'AGT', 'AAA', 'AAG'], ['CCA', 'TCA', 'GCA', 'AAA', 'ATA', 'AGA'], ['CCC', 'TCC', 'GCC', 'AAC', 'ATC', 'AGC'], ['CCG', 'TCG', 'GCG', 'AAG', 'ATG', 'AGG'], ['CCT', 'TCT', 'GCT', 'AAT', 'ATT', 'AGT'], ['GGA', 'AAA', 'ACA', 'ATA', 'AGC', 'AGT'], ['CGC', 'TGC', 'GGC', 'AAC', 'ACC', 'ATC', 'AGA', 'AGG'], ['TGG', 'GGG', 'AAG', 'ACG', 'ATG', 'AGC', 'AGT'], ['CGT', 'TGT', 'GGT', 'AAT', 'ACT', 'ATT', 'AGA', 'AGG'], ['CTA', 'TTA', 'GTA', 'AAA', 'ACA', 'AGA', 'ATG'], ['CTC', 'TTC', 'GTC', 'AAC', 'ACC', 'AGC', 'ATG'], ['CTG', 'TTG', 'GTG', 'AAG', 'ACG', 'AGG', 'ATA', 'ATC', 'ATT'], ['CTT', 'TTT', 'GTT', 'AAT', 'ACT', 'AGT', 'ATG'], ['AAA', 'GAA', 'CCA', 'CTA', 'CGA', 'CAC', 'CAT'], ['AAC', 'TAC', 'GAC', 'CCC', 'CTC', 'CGC', 'CAA', 'CAG'], ['AAG', 'GAG', 'CCG', 'CTG', 'CGG', 'CAC', 'CAT'], ['AAT', 'TAT', 'GAT', 'CCT', 'CTT', 'CGT', 'CAA', 'CAG'], ['ACA', 'TCA', 'GCA', 'CAA', 'CTA', 'CGA'], ['ACC', 'TCC', 'GCC', 'CAC', 'CTC', 'CGC'], ['ACG', 'TCG', 'GCG', 'CAG', 'CTG', 'CGG'], ['ACT', 'TCT', 'GCT', 'CAT', 'CTT', 'CGT'], ['GGA', 'CAA', 'CCA', 'CTA'], ['AGC', 'TGC', 'GGC', 'CAC', 'CCC', 'CTC'], ['TGG', 'GGG', 'CAG', 'CCG', 'CTG'], ['AGT', 'TGT', 'GGT', 'CAT', 'CCT', 'CTT'], ['ATA', 'GTA', 'CAA', 'CCA', 'CGA'], ['ATC', 'TTC', 'GTC', 'CAC', 'CCC', 'CGC'], ['ATG', 'GTG', 'CAG', 'CCG', 'CGG'], ['ATT', 'TTT', 'GTT', 'CAT', 'CCT', 'CGT'], ['AAA', 'CAA', 'GCA', 'GTA', 'GGA', 'GAC', 'GAT'], ['AAC', 'CAC', 'TAC', 'GCC', 'GTC', 'GGC', 'GAA', 'GAG'], ['AAG', 'CAG', 'GCG', 'GTG', 'GGG', 'GAC', 'GAT'], ['AAT', 'CAT', 'TAT', 'GCT', 'GTT', 'GGT', 'GAA', 'GAG'], ['ACA', 'CCA', 'TCA', 'GAA', 'GTA', 'GGA'], ['ACC', 'CCC', 'TCC', 'GAC', 'GTC', 'GGC'], ['ACG', 'CCG', 'TCG', 'GAG', 'GTG', 'GGG'], ['ACT', 'CCT', 'TCT', 'GAT', 'GTT', 'GGT'], ['AGA', 'CGA', 'GAA', 'GCA', 'GTA'], ['AGC', 'CGC', 'TGC', 'GAC', 'GCC', 'GTC'], ['AGG', 'CGG', 'TGG', 'GAG', 'GCG', 'GTG'], ['AGT', 'CGT', 'TGT', 'GAT', 'GCT', 'GTT'], ['ATA', 'CTA', 'TTA', 'GAA', 'GCA', 'GGA'], ['ATC', 'CTC', 'TTC', 'GAC', 'GCC', 'GGC'], ['ATG', 'CTG', 'TTG', 'GAG', 'GCG', 'GGG'], ['ATT', 'CTT', 'TTT', 'GAT', 'GCT', 'GGT'], ['AAC', 'CAC', 'GAC', 'TCC', 'TTC', 'TGC'], ['AAT', 'CAT', 'GAT', 'TCT', 'TTT', 'TGT'], ['ACA', 'CCA', 'GCA', 'TTA'], ['ACC', 'CCC', 'GCC', 'TAC', 'TTC', 'TGC'], ['ACG', 'CCG', 'GCG', 'TTG', 'TGG'], ['ACT', 'CCT', 'GCT', 'TAT', 'TTT', 'TGT'], ['AGC', 'CGC', 'GGC', 'TAC', 'TCC', 'TTC', 'TGG'], ['AGG', 'CGG', 'GGG', 'TCG', 'TTG', 'TGC', 'TGT'], ['AGT', 'CGT', 'GGT', 'TAT', 'TCT', 'TTT', 'TGG'], ['ATA', 'GTA', 'TCA', 'TTC', 'TTT'], ['ATC', 'CTC', 'GTC', 'TAC', 'TCC', 'TGC', 'TTA', 'TTG'], ['ATG', 'GTG', 'TCG', 'TGG', 'TTC', 'TTT'], ['ATT', 'CTT', 'GTT', 'TAT', 'TCT', 'TGT', 'TTA', 'TTG']]

genetic_code = [["GCA", "GCC", "GCG", "GCT"], ["TGC","TGT"], ["GAC", "GAT"], ["GAA", "GAG"], ["TTC", "TTT"], ["GGA", "GGC", "GGG", "GGT"], ["CAC", "CAT"], ["ATA", "ATC", "ATT"], ["AAA", "AAG"], ["CTA", "CTC", "CTG", "CTT", "TTA", "TTG"], ["ATG"], ["AAC", "AAT"], ["CCA", "CCC", "CCG", "CCT"], ["CAA", "CAG"], ["AGA", "AGG", "CGA", "CGC", "CGG", "CGT"] , ["AGC", "AGT", "TCA", "TCC", "TCG", "TCT"], ["ACA", "ACC", "ACG", "ACT"], ["GTA", "GTC", "GTG", "GTT"], ["TGG"], ["TAC", "TAT"]]





        

############################# SIMULATION FUNCTIONS #######################################

def simulate(seqfile, numaa, freqClass, freqBy, tree, mu, length, prespec = False):
    ''' Simulate single partition according to mutsel model.
        Uses equal mutation rates.
    '''
    try:
        my_tree = readTree(file = tree)
    except:
        my_tree = readTree(tree = tree)
    
    # Equal frequencies
    if freqClass == 'equal':
        fobj = EqualFreqs(by = freqBy, type = 'codon')
        aminos_used = ''
    
    # Random frequencies
    elif freqClass == 'random':
        fobj = RandFreqs(by = freqBy, type = 'codon')
        aminos_used = ''
        
    # User frequencies
    elif freqClass == 'user':
        userFreq, aminos_used = generateExpFreqDict(numaa, prespec)
        fobj = UserFreqs(by = freqBy, type = 'codon', freqs = userFreq)

    else:
        raise AssertionError("Bad freqClass specification. Byebye.")
   
    f = fobj.calcFreqs() 
    
    model = Model()
    params = {'alpha':1.0, 'beta':1.0, 'mu': {'AC': mu, 'CA':mu, 'AG': mu, 'GA':mu, 'AT': mu, 'TA':mu, 'CG': mu, 'GC':mu, 'CT': mu, 'TC':mu, 'GT': mu, 'TG':mu}}
    params['stateFreqs'] = f
    model.params = params
    mat = mutSel_MatrixBuilder(model)
    model.Q = mat.buildQ()
    partitions = [(length, model)]        
    
    myEvolver = StaticEvolver(partitions = partitions, tree = my_tree, outfile = seqfile)
    myEvolver.sim_sub_tree(my_tree)
    myEvolver.writeSequences()
    
    return f, aminos_used


def generateExpFreqDict(size, specified=True):
    ''' Generate a dictionary of exponentially distributed amino acid frequencies.
        size = number of amino acids
        specified = use the prespecified amino acids for a given size (True) or get random amino acids (False)
        NOTE: if size==1 or size>6, we will just get random amino acids, NOT SPECIFIED. When that many, probably physiochemical properties don't mean much anymore.
        If size==1, need to make sure that the amino picked has MULTIPLE synonymous codons, otherwise evolution is broken. = NOT AMINO ACID 10,18 (M and W)
    '''

    # Create the amino acid frequency distribution
    final_dict = {}
    raw = np.random.exponential(size=size)
    final = raw/np.sum(raw)
    
    # Create a dictionary of frequencies using "final"
    if specified and 2<=size<=6:    
        prespec_aa = { 2: ['I', 'V'], 3: ['H', 'K', 'R'], 4: ['A', 'S', 'C', 'V'], 5: ['A', 'S', 'T', 'G', 'V'], 6: ['A', 'S', 'T', 'G', 'V', 'C']}
        count = 0
        for aa in prespec_aa[size]:
            final_dict[aa] = final[count]
            count +=1
    else:
        assert (1 <= size <= 20), "that's a silly size."
        aminos = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
        for i in range(size):
            n = randint(0,len(aminos)-1)
            # make sure evolution is allowed. no met (10) or trp (18) allowed.
            if size == 1:
                while n == 10 or n == 18:
                    n = randint(0,len(aminos)-1)   
            final_dict[aminos[n]] = final[i]
            aminos.pop(n)        
    
    return final_dict, "".join(final_dict.keys()) 



############################ HYPHY-RELATED FUNCTIONS #####################################
def runhyphy(batchfile, matrix_name, seqfile, treefile, cpu, codonfreq, initw=0.5):
    ''' pretty specific function.'''
    setuphyphy1 = "cp "+seqfile+" temp.fasta"
    setup1 = subprocess.call(setuphyphy1, shell = True)
    assert(setup1 == 0), "couldn't create temp.fasta"
    
    setuphyphy2 = "cat "+treefile+" >> temp.fasta"
    setup2 = subprocess.call(setuphyphy2, shell = True)
    assert(setup2 == 0), "couldn't add tree to hyphy infile"
    
    hyf = freq2Hyphy(codonfreq)
    setuphyphy3 = "sed 's/MYFREQUENCIES/"+hyf+"/g' "+batchfile+" > run.bf"
    setup3 = subprocess.call(setuphyphy3, shell = True)
    assert(setup3 == 0), "couldn't properly add in frequencies"
    
    setuphyphy4 = "sed -i 's/MYMATRIX/"+matrix_name+"/g' run.bf"
    setup4 = subprocess.call(setuphyphy4, shell = True)
    assert(setup4 == 0), "couldn't properly define matrix"
    
    setuphyphy5 = "sed -i 's/MYINITIALW/"+str(initw)+"/g' run.bf"
    setup5 = subprocess.call(setuphyphy5, shell = True)
    assert(setup5 == 0), "couldn't properly define intial omega guess"
    
    
    hyphy = "./HYPHYMP run.bf CPU="+cpu+" > hyout.txt"
    runhyphy = subprocess.call(hyphy, shell = True)
    assert (runhyphy == 0), "hyphy fail"
    
    # grab hyphy output
    hyout = open('hyout.txt', 'r')
    hylines = hyout.readlines()
    hyout.close()
    for line in hylines:
        findw = re.search("^w=(\d+\.*\d*)", line)
        if findw:
            hyphy_w = findw.group(1)
            break
    return hyphy_w
    
def freq2Hyphy(f):
    ''' Convert codon frequencies to a form hyphy can use. '''
    hyphy_f = "{"
    for freq in f:
        hyphy_f += "{"
        hyphy_f += str(freq)
        hyphy_f += "},"
    hyphy_f = hyphy_f[:-1]
    hyphy_f += "}"
    return hyphy_f




############################ PAML-RELATED FUNCTIONS ###############################
def runpaml(seqfile, initw=0.5):
    setuppaml1 = "cp "+seqfile+" temp.fasta"
    setup1 = subprocess.call(setuppaml1, shell = True)
    assert(setup1 == 0), "couldn't create temp.fasta"
    
    setuppaml2 = 'sed "s/MYINITIALW/'+str(initw)+'/g" codeml_raw.txt > codeml.ctl' 
    setup2 = subprocess.call(setuppaml2, shell = True)
    assert(setup2 == 0), "couldn't set paml initial w"
    
    runpaml = subprocess.call("./codeml", shell=True)
    assert (runpaml == 0), "paml fail"

    # Grab paml output
    paml_w = parsePAML("pamloutfile")
    return paml_w
    
def parsePAML(pamlfile):
    ''' get the omega from a paml file. model run is single omega for an entire alignment. '''
    paml = open(pamlfile, 'rU')
    pamlines = paml.readlines()
    paml.close()
    omega = None
    for line in pamlines:
        findw = re.search("^omega \(dN\/dS\)\s*=\s*(\d+\.*\d*)", line)
        if findw:
            omega = findw.group(1)
            break
    assert (omega is not None), "couldn't get omega from paml file"
    return omega




############################# NEI-GOJOBORI FUNCTIONS ##################################
def run_neigojo(seqfile):
    ''' Get omega using counting method '''
    M = MutationCounter()
    S = SiteCounter()
    records = list(SeqIO.parse(seqfile, 'fasta'))
    s1 = records[0].seq
    s2 = records[1].seq
    ( ns_mut, s_mut ) = M.countMutations( s1, s2 )
    ( ns_sites1, s_sites1 ) = S.countSites( s1 )
    ( ns_sites2, s_sites2 ) = S.countSites( s2 )
    dS = 2*sum( s_mut )/(sum( s_sites1 ) + sum( s_sites2 ))
    dN = 2*sum( ns_mut )/(sum( ns_sites2 ) + sum( ns_sites2 ))
    return dN/dS



############################# OMEGA DERIVATION FUNCTIONS ##############################

def deriveOmega(codonFreq):
    ''' Derive an omega using codon frequencies. ''' 
    nonZero = getNonZeroFreqs(codonFreq) # get indices which aren't zero.
    
    kN=0. #dN numerator
    nN=0. #dN denominator. NOTE: Does not correct for consider number of nonsyn options

    # Calculations
    for i in nonZero:
        fix_sum=0.
        
        ### Nonsynonymous.
        for nscodon in nslist[i]:
            nscodon_freq = codonFreq[codons.index(nscodon)]
            fix_sum += fix(float(codonFreq[i]), float(nscodon_freq))                    
            nN += codonFreq[i]
        kN += fix_sum*codonFreq[i]

    # Final dN/dS
    if kN < zero:
        return 0., len(nonZero)
    else:
        return kN/nN, len(nonZero)
        
        
        

def getNonZeroFreqs(freq):
    ''' Return indices whose frequencies are not 0.'''
    nonZero = []
    for i in range(len(freq)):
        if freq[i] > zero:
            nonZero.append(i)
    return nonZero




def fix(fi, fj):
    if fi == fj:
        return 1.
    elif fi == 0.  or fj == 0.:
        return 0.
    else:
        return (np.log(fj) - np.log(fi)) / (1 - fi/fj)
#########################################################################################
