## SJS. Functions that accompany run_siminf.py, etc.
# NOTE: to use simulation library, must cp the src/ directory (*not* contents, the whole directory!) into wdir.

import os
import re
import sys
import shutil
import subprocess
import numpy as np
from random import randint


# Simulation code
#sys.path.append('src/')
from misc import *
from newick import *
from stateFreqs import *
from matrixBuilder import *
from evolver import *

# Globals
zero = 1e-8
amino_acids  = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
codons=["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
codon_dict = {"AAA":"K", "AAC":"N", "AAG":"K", "AAT":"N", "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T", "AGA":"R", "AGC":"S", "AGG":"R", "AGT":"S", "ATA":"I", "ATC":"I", "ATG":"M", "ATT":"I", "CAA":"Q", "CAC":"H", "CAG":"Q", "CAT":"H", "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P", "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R", "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L", "GAA":"E", "GAC":"D", "GAG":"E", "GAT":"D", "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A", "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G", "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V", "TAC":"Y", "TAT":"Y", "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S", "TGC":"C", "TGG":"W", "TGT":"C", "TTA":"L", "TTC":"F", "TTG":"L", "TTT":"F"}
genetic_code = [["GCA", "GCC", "GCG", "GCT"], ["TGC","TGT"], ["GAC", "GAT"], ["GAA", "GAG"], ["TTC", "TTT"], ["GGA", "GGC", "GGG", "GGT"], ["CAC", "CAT"], ["ATA", "ATC", "ATT"], ["AAA", "AAG"], ["CTA", "CTC", "CTG", "CTT", "TTA", "TTG"], ["ATG"], ["AAC", "AAT"], ["CCA", "CCC", "CCG", "CCT"], ["CAA", "CAG"], ["AGA", "AGG", "CGA", "CGC", "CGG", "CGT"] , ["AGC", "AGT", "TCA", "TCC", "TCG", "TCT"], ["ACA", "ACC", "ACG", "ACT"], ["GTA", "GTC", "GTG", "GTT"], ["TGG"], ["TAC", "TAT"]]
grantham = {'AA':0, 'AC':195, 'AD':126, 'AE':107, 'AF':113, 'AG':60, 'AH':86, 'AI':94, 'AK':106, 'AL':96, 'AM':84, 'AN':111, 'AP':27, 'AQ':91, 'AR':112, 'AS':99, 'AT':58, 'AV':64, 'AW':148, 'AY':112, 'CA':195, 'CC':0, 'CD':154, 'CE':170, 'CF':205, 'CG':159, 'CH':174, 'CI':198, 'CK':202, 'CL':198, 'CM':196, 'CN':139, 'CP':169, 'CQ':154, 'CR':180, 'CS':112, 'CT':149, 'CV':192, 'CW':215, 'CY':194, 'DA':126, 'DC':154, 'DD':0, 'DE':45, 'DF':177, 'DG':94, 'DH':81, 'DI':168, 'DK':101, 'DL':172, 'DM':160, 'DN':23, 'DP':108, 'DQ':61, 'DR':96, 'DS':65, 'DT':85, 'DV':152, 'DW':181, 'DY':160, 'EA':107, 'EC':170, 'ED':45, 'EE':0, 'EF':140, 'EG':98, 'EH':40, 'EI':134, 'EK':56, 'EL':138, 'EM':126, 'EN':42, 'EP':93, 'EQ':29, 'ER':54, 'ES':80, 'ET':65, 'EV':121, 'EW':152, 'EY':122, 'FA':113, 'FC':205, 'FD':177, 'FE':140, 'FF':0, 'FG':153, 'FH':100, 'FI':21, 'FK':102, 'FL':22, 'FM':28, 'FN':158, 'FP':114, 'FQ':116, 'FR':97, 'FS':155, 'FT':103, 'FV':50, 'FW':40, 'FY':22, 'GA':60, 'GC':159, 'GD':94, 'GE':98, 'GF':153, 'GG':0, 'GH':98, 'GI':135, 'GK':127, 'GL':138, 'GM':127, 'GN':80, 'GP':42, 'GQ':87, 'GR':125, 'GS':56, 'GT':59, 'GV':109, 'GW':184, 'GY':147, 'HA':86, 'HC':174, 'HD':81, 'HE':40, 'HF':100, 'HG':98, 'HH':0, 'HI':94, 'HK':32, 'HL':99, 'HM':87, 'HN':68, 'HP':77, 'HQ':24, 'HR':29, 'HS':89, 'HT':47, 'HV':84, 'HW':115, 'HY':83, 'IA':94, 'IC':198, 'ID':168, 'IE':134, 'IF':21, 'IG':135, 'IH':94, 'II':0, 'IK':102, 'IL':5, 'IM':10, 'IN':149, 'IP':95, 'IQ':109, 'IR':97, 'IS':142, 'IT':89, 'IV':29, 'IW':61, 'IY':33, 'KA':106, 'KC':202, 'KD':101, 'KE':56, 'KF':102, 'KG':127, 'KH':32, 'KI':102, 'KK':0, 'KL':107, 'KM':95, 'KN':94, 'KP':103, 'KQ':53, 'KR':26, 'KS':121, 'KT':78, 'KV':97, 'KW':110, 'KY':85, 'LA':96, 'LC':198, 'LD':172, 'LE':138, 'LF':22, 'LG':138, 'LH':99, 'LI':5, 'LK':107, 'LL':0, 'LM':15, 'LN':153, 'LP':98, 'LQ':113, 'LR':102, 'LS':145, 'LT':92, 'LV':32, 'LW':61, 'LY':36, 'MA':84, 'MC':196, 'MD':160, 'ME':126, 'MF':28, 'MG':127, 'MH':87, 'MI':10, 'MK':95, 'ML':15, 'MM':0, 'MN':142, 'MP':87, 'MQ':101, 'MR':91, 'MS':135, 'MT':81, 'MV':21, 'MW':67, 'MY':36, 'NA':111, 'NC':139, 'ND':23, 'NE':42, 'NF':158, 'NG':80, 'NH':68, 'NI':149, 'NK':94, 'NL':153, 'NM':142, 'NN':0, 'NP':91, 'NQ':46, 'NR':86, 'NS':46, 'NT':65, 'NV':133, 'NW':174, 'NY':143, 'PA':27, 'PC':169, 'PD':108, 'PE':93, 'PF':114, 'PG':42, 'PH':77, 'PI':95, 'PK':103, 'PL':98, 'PM':87, 'PN':91, 'PP':0, 'PQ':76, 'PR':103, 'PS':74, 'PT':38, 'PV':68, 'PW':147, 'PY':110, 'QA':91, 'QC':154, 'QD':61, 'QE':29, 'QF':116, 'QG':87, 'QH':24, 'QI':109, 'QK':53, 'QL':113, 'QM':101, 'QN':46, 'QP':76, 'QQ':0, 'QR':43, 'QS':68, 'QT':42, 'QV':96, 'QW':130, 'QY':99, 'RA':112, 'RC':180, 'RD':96, 'RE':54, 'RF':97, 'RG':125, 'RH':29, 'RI':97, 'RK':26, 'RL':102, 'RM':91, 'RN':86, 'RP':103, 'RQ':43, 'RR':0, 'RS':110, 'RT':71, 'RV':96, 'RW':101, 'RY':77, 'SA':99, 'SC':112, 'SD':65, 'SE':80, 'SF':155, 'SG':56, 'SH':89, 'SI':142, 'SK':121, 'SL':145, 'SM':135, 'SN':46, 'SP':74, 'SQ':68, 'SR':110, 'SS':0, 'ST':58, 'SV':124, 'SW':177, 'SY':144, 'TA':58, 'TC':149, 'TD':85, 'TE':65, 'TF':103, 'TG':59, 'TH':47, 'TI':89, 'TK':78, 'TL':92, 'TM':81, 'TN':65, 'TP':38, 'TQ':42, 'TR':71, 'TS':58, 'TT':0, 'TV':69, 'TW':128, 'TY':92, 'VA':64, 'VC':192, 'VD':152, 'VE':121, 'VF':50, 'VG':109, 'VH':84, 'VI':29, 'VK':97, 'VL':32, 'VM':21, 'VN':133, 'VP':68, 'VQ':96, 'VR':96, 'VS':124, 'VT':69, 'VV':0, 'VW':88, 'VY':55, 'WA':148, 'WC':215, 'WD':181, 'WE':152, 'WF':40, 'WG':184, 'WH':115, 'WI':61, 'WK':110, 'WL':61, 'WM':67, 'WN':174, 'WP':147, 'WQ':130, 'WR':101, 'WS':177, 'WT':128, 'WV':88, 'WW':0, 'WY':37, 'YA':112, 'YC':194, 'YD':160, 'YE':122, 'YF':22, 'YG':147, 'YH':83, 'YI':33, 'YK':85, 'YL':36, 'YM':36, 'YN':143, 'YP':110, 'YQ':99, 'YR':77, 'YS':144, 'YT':92, 'YV':55, 'YW':37, 'YY':0}



        

################################################ SIMULATION FUNCTIONS #####################################################

def simulate(f, seqfile, tree, mu_dict, length, beta=None):
    ''' Simulate single partition according to either codon or mutsel model (check beta (dN) value for which model).
    '''
    try:
        my_tree = readTree(file = tree, flags = False)
    except:
        my_tree = readTree(tree = tree, flags = False) 
          
    model = Model()
    if beta:
    	params = {'stateFreqs':f, 'alpha':1.0, 'beta':float(beta), 'mu': mu_dict}
        model.params = params
        mat = mechCodon_MatrixBuilder(model)
    else:
        print "Setting up MutSel model"
        params = {'stateFreqs':f, 'alpha':1.0, 'beta':1.0, 'mu': mu_dict}
        model.params = params
        mat = mutSel_MatrixBuilder(model)
    model.Q = mat.buildQ()
    partitions = [(length, {"rootModel":model})]        
    myEvolver = Evolver(partitions, "rootModel" )
    myEvolver.simulate(my_tree)
    myEvolver.writeSequences(outfile = seqfile)



def setFreqs(freqfile, lambda_):
    ''' Returns codon frequencies and gc content. '''
 
    redo = True
    while redo:
        # Frequencies based on boltzmann dist, such that amino acid frequencies distributed exponentially.
        raw_aafreqs = setBoltzFreqs(lambda_) # gets frequencies 
        assert(np.sum(raw_aafreqs) - 1.0 < zero), "bad amino freq calculation"
        aaFreq = dict(zip(amino_acids, raw_aafreqs))
    
        # Calculate codon state frequencies given amino acid frequencies, above.
        fobj = UserFreqs(by = 'amino', freqs = aaFreq)
        codonFreq = fobj.calcFreqs(type = 'codon', savefile = freqfile)
        
        # Should I redo based on excessive codon freq stringency?
        redo = np.any(codonFreq >= 0.985)
        
        # Get gc content
        nucFreq = fobj.calcFreqs(type = 'nuc')
        gc = nucFreq[1] + nucFreq[2]
    return codonFreq, gc


def setBoltzFreqs(lambda_):
    ''' Use Boltzmann distribution to get amino acid frequencies for a certain number of amino acids.'''
    # lambda_ basically determines the strength of selection. We are now using it as the stddev of the normal distribution from which ssc's are drawn.
    ssc_values = np.random.normal(loc=0., scale=lambda_, size = 20) #ssc = scaled selection coefficient.
    ssc_values[0] = 0. # set one value to zero to make these values *scaled* selection coeffs 
    numer_list = np.zeros(20)
    for ssc in range(20):
        val = np.exp(-1. * ssc_values[ssc])
        numer_list[ssc] = val
    return numer_list/np.sum(numer_list)    

def calcCodonEntropy(f):
    sum = 0.
    for entry in f:
        if entry > 1e-8:
            sum += entry*np.log(entry)
    return -1. * sum
    

####################################### OMEGA DERIVATION FUNCTIONS ######################################

def deriveOmega(codonFreqs, mu_dict):
    ''' NOTE: assumes dS = 1. If codon bias is introduced, this will be violated. '''
    
    numer = 0.; denom = 0.;
    cfreqs = dict(zip(codons, codonFreqs))

    for codon in cfreqs:
        if cfreqs[codon] > zero:  
            rate, sites = calcNonSyn(codon, cfreqs, mu_dict)
            numer += rate
            denom += sites
    assert( denom != 0. ), "Omega derivation indicates no evolution, maybe?"
    return numer/denom


     
def calcNonSyn(source, cfreqs, mu_dict):
    rate = 0.
    sites = 0.
    sourceFreq = cfreqs[source]
    for target in codons:
        diff = getNucleotideDiff(source, target) # only consider single nucleotide differences since are calculating instantaneous.
        if codon_dict[source] != codon_dict[target] and cfreqs[target] > zero and len(diff) == 2:
            rate  += calcFix( sourceFreq, cfreqs[target] ) * mu_dict[diff]
            sites += mu_dict[diff]
    rate  *= sourceFreq
    sites *= sourceFreq
    return rate, sites
    

def getNucleotideDiff(source, target):
    diff = ''
    for i in range(3):
        if source[i] != target[i]: 
            diff += "".join( source[i]+target[i] )
    return diff

    
def calcFix(fi, fj):
    if fi == fj:
        return 1.
    elif fi == 0. or fj == 0.:
        return 0.
    else:
        return (np.log(fj) - np.log(fi)) / (1 - fi/fj)    
#########################################################################################




################################################# HYPHY-RELATED FUNCTIONS ############################################################
def runhyphy(batchfile, matrix_name, seqfile, treefile, cpu, kappa, codonFreqs):
    ''' pretty specific function.
        codonFreqs = what was simulated. note that empirical are not signif different from true, so this is fine.
    '''
    
  
    # Set up sequence file with tree
    setuphyphy1 = "cp "+seqfile+" temp.fasta"
    setup1 = subprocess.call(setuphyphy1, shell = True)
    assert(setup1 == 0), "couldn't create temp.fasta"
    setuphyphy2 = "cat "+treefile+" >> temp.fasta"
    setup2 = subprocess.call(setuphyphy2, shell = True)
    assert(setup2 == 0), "couldn't add tree to hyphy infile"
    
    # Set up matrix (GY94 or MG94), within run.bf
    setuphyphy3 = "sed 's/MYMATRIX/"+matrix_name+"/g' " + batchfile + " > run.bf"
    setup3 = subprocess.call(setuphyphy3, shell = True)
    assert(setup3 == 0), "couldn't properly define matrix"
    
    # Set up frequencies.
    if codonFreqs is not None:
        hyf = freq2hyphy(codonFreqs)
        setuphyphyf = "sed -i 's/MYFREQUENCIES/"+hyf+"/g' run.bf"
        setupf = subprocess.call(setuphyphyf, shell = True)
        assert(setupf == 0), "couldn't properly add in frequencies"
    
    # Set up kappa
    if kappa != 'free':
        sedkappa = "sed 's/k/"+str(kappa)+"/g' matrices_raw.mdl > matrices.mdl"
        runsedkappa = subprocess.call(sedkappa, shell=True)
        assert(runsedkappa == 0), "couldn't set up kappa"
    else:
        shutil.copy('matrices_raw.mdl', 'matrices.mdl')

    # Run hyphy
    hyphy = "./HYPHYMP run.bf CPU="+cpu+" > hyout.txt"
    runhyphy = subprocess.call(hyphy, shell = True)
    assert (runhyphy == 0), "hyphy fail"
    
    # grab hyphy output
    return parseHyphyGY94('hyout.txt')



def freq2hyphy(f):
    ''' Convert codon frequencies to a form hyphy can use. '''
    hyphy_f = "{"
    for freq in f:
        hyphy_f += "{"
        hyphy_f += str(freq)
        hyphy_f += "},"
    hyphy_f = hyphy_f[:-1]
    hyphy_f += "};"
    return hyphy_f
    
    
def parseHyphyGY94(file):
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

