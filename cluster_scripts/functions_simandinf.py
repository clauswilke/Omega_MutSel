## SJS. Functions that accompany run_siminf.py, etc.
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

# Globals
zero = 1e-8
amino_acids  = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
codons=["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
nslist = [['CAA', 'GAA', 'ACA', 'ATA', 'AGA', 'AAC', 'AAT'], ['CAC', 'TAC', 'GAC', 'ACC', 'ATC', 'AGC', 'AAA', 'AAG'], ['CAG', 'GAG', 'ACG', 'ATG', 'AGG', 'AAC', 'AAT'], ['CAT', 'TAT', 'GAT', 'ACT', 'ATT', 'AGT', 'AAA', 'AAG'], ['CCA', 'TCA', 'GCA', 'AAA', 'ATA', 'AGA'], ['CCC', 'TCC', 'GCC', 'AAC', 'ATC', 'AGC'], ['CCG', 'TCG', 'GCG', 'AAG', 'ATG', 'AGG'], ['CCT', 'TCT', 'GCT', 'AAT', 'ATT', 'AGT'], ['GGA', 'AAA', 'ACA', 'ATA', 'AGC', 'AGT'], ['CGC', 'TGC', 'GGC', 'AAC', 'ACC', 'ATC', 'AGA', 'AGG'], ['TGG', 'GGG', 'AAG', 'ACG', 'ATG', 'AGC', 'AGT'], ['CGT', 'TGT', 'GGT', 'AAT', 'ACT', 'ATT', 'AGA', 'AGG'], ['CTA', 'TTA', 'GTA', 'AAA', 'ACA', 'AGA', 'ATG'], ['CTC', 'TTC', 'GTC', 'AAC', 'ACC', 'AGC', 'ATG'], ['CTG', 'TTG', 'GTG', 'AAG', 'ACG', 'AGG', 'ATA', 'ATC', 'ATT'], ['CTT', 'TTT', 'GTT', 'AAT', 'ACT', 'AGT', 'ATG'], ['AAA', 'GAA', 'CCA', 'CTA', 'CGA', 'CAC', 'CAT'], ['AAC', 'TAC', 'GAC', 'CCC', 'CTC', 'CGC', 'CAA', 'CAG'], ['AAG', 'GAG', 'CCG', 'CTG', 'CGG', 'CAC', 'CAT'], ['AAT', 'TAT', 'GAT', 'CCT', 'CTT', 'CGT', 'CAA', 'CAG'], ['ACA', 'TCA', 'GCA', 'CAA', 'CTA', 'CGA'], ['ACC', 'TCC', 'GCC', 'CAC', 'CTC', 'CGC'], ['ACG', 'TCG', 'GCG', 'CAG', 'CTG', 'CGG'], ['ACT', 'TCT', 'GCT', 'CAT', 'CTT', 'CGT'], ['GGA', 'CAA', 'CCA', 'CTA'], ['AGC', 'TGC', 'GGC', 'CAC', 'CCC', 'CTC'], ['TGG', 'GGG', 'CAG', 'CCG', 'CTG'], ['AGT', 'TGT', 'GGT', 'CAT', 'CCT', 'CTT'], ['ATA', 'GTA', 'CAA', 'CCA', 'CGA'], ['ATC', 'TTC', 'GTC', 'CAC', 'CCC', 'CGC'], ['ATG', 'GTG', 'CAG', 'CCG', 'CGG'], ['ATT', 'TTT', 'GTT', 'CAT', 'CCT', 'CGT'], ['AAA', 'CAA', 'GCA', 'GTA', 'GGA', 'GAC', 'GAT'], ['AAC', 'CAC', 'TAC', 'GCC', 'GTC', 'GGC', 'GAA', 'GAG'], ['AAG', 'CAG', 'GCG', 'GTG', 'GGG', 'GAC', 'GAT'], ['AAT', 'CAT', 'TAT', 'GCT', 'GTT', 'GGT', 'GAA', 'GAG'], ['ACA', 'CCA', 'TCA', 'GAA', 'GTA', 'GGA'], ['ACC', 'CCC', 'TCC', 'GAC', 'GTC', 'GGC'], ['ACG', 'CCG', 'TCG', 'GAG', 'GTG', 'GGG'], ['ACT', 'CCT', 'TCT', 'GAT', 'GTT', 'GGT'], ['AGA', 'CGA', 'GAA', 'GCA', 'GTA'], ['AGC', 'CGC', 'TGC', 'GAC', 'GCC', 'GTC'], ['AGG', 'CGG', 'TGG', 'GAG', 'GCG', 'GTG'], ['AGT', 'CGT', 'TGT', 'GAT', 'GCT', 'GTT'], ['ATA', 'CTA', 'TTA', 'GAA', 'GCA', 'GGA'], ['ATC', 'CTC', 'TTC', 'GAC', 'GCC', 'GGC'], ['ATG', 'CTG', 'TTG', 'GAG', 'GCG', 'GGG'], ['ATT', 'CTT', 'TTT', 'GAT', 'GCT', 'GGT'], ['AAC', 'CAC', 'GAC', 'TCC', 'TTC', 'TGC'], ['AAT', 'CAT', 'GAT', 'TCT', 'TTT', 'TGT'], ['ACA', 'CCA', 'GCA', 'TTA'], ['ACC', 'CCC', 'GCC', 'TAC', 'TTC', 'TGC'], ['ACG', 'CCG', 'GCG', 'TTG', 'TGG'], ['ACT', 'CCT', 'GCT', 'TAT', 'TTT', 'TGT'], ['AGC', 'CGC', 'GGC', 'TAC', 'TCC', 'TTC', 'TGG'], ['AGG', 'CGG', 'GGG', 'TCG', 'TTG', 'TGC', 'TGT'], ['AGT', 'CGT', 'GGT', 'TAT', 'TCT', 'TTT', 'TGG'], ['ATA', 'GTA', 'TCA', 'TTC', 'TTT'], ['ATC', 'CTC', 'GTC', 'TAC', 'TCC', 'TGC', 'TTA', 'TTG'], ['ATG', 'GTG', 'TCG', 'TGG', 'TTC', 'TTT'], ['ATT', 'CTT', 'GTT', 'TAT', 'TCT', 'TGT', 'TTA', 'TTG']]
slist=[['AAG'], ['AAT'], ['AAA'], ['AAC'], ['ACC', 'ACT', 'ACG'], ['ACA', 'ACT', 'ACG'], ['ACA', 'ACC', 'ACT'], ['ACA', 'ACC', 'ACG'], ['CGA', 'AGG'], ['AGT'], ['CGG', 'AGA'], ['AGC'], ['ATC', 'ATT'], ['ATA', 'ATT'], [], ['ATA', 'ATC'], ['CAG'], ['CAT'], ['CAA'], ['CAC'], ['CCC', 'CCT', 'CCG'], ['CCA', 'CCT', 'CCG'], ['CCA', 'CCC', 'CCT'], ['CCA', 'CCC', 'CCG'], ['AGA', 'CGC', 'CGT', 'CGG'], ['CGA', 'CGT', 'CGG'], ['AGG', 'CGA', 'CGC', 'CGT'], ['CGA', 'CGC', 'CGG'], ['TTA', 'CTC', 'CTT', 'CTG'], ['CTA', 'CTT', 'CTG'], ['TTG', 'CTA', 'CTC', 'CTT'], ['CTA', 'CTC', 'CTG'], ['GAG'], ['GAT'], ['GAA'], ['GAC'], ['GCC', 'GCT', 'GCG'], ['GCA', 'GCT', 'GCG'], ['GCA', 'GCC', 'GCT'], ['GCA', 'GCC', 'GCG'], ['GGC', 'GGT', 'GGG'], ['GGA', 'GGT', 'GGG'], ['GGA', 'GGC', 'GGT'], ['GGA', 'GGC', 'GGG'], ['GTC', 'GTT', 'GTG'], ['GTA', 'GTT', 'GTG'], ['GTA', 'GTC', 'GTT'], ['GTA', 'GTC', 'GTG'], ['TAT'], ['TAC'], ['TCC', 'TCT', 'TCG'], ['TCA', 'TCT', 'TCG'], ['TCA', 'TCC', 'TCT'], ['TCA', 'TCC', 'TCG'], ['TGT'], [], ['TGC'], ['CTA', 'TTG'], ['TTT'], ['CTG', 'TTA'], ['TTC']]
genetic_code = [["GCA", "GCC", "GCG", "GCT"], ["TGC","TGT"], ["GAC", "GAT"], ["GAA", "GAG"], ["TTC", "TTT"], ["GGA", "GGC", "GGG", "GGT"], ["CAC", "CAT"], ["ATA", "ATC", "ATT"], ["AAA", "AAG"], ["CTA", "CTC", "CTG", "CTT", "TTA", "TTG"], ["ATG"], ["AAC", "AAT"], ["CCA", "CCC", "CCG", "CCT"], ["CAA", "CAG"], ["AGA", "AGG", "CGA", "CGC", "CGG", "CGT"] , ["AGC", "AGT", "TCA", "TCC", "TCG", "TCT"], ["ACA", "ACC", "ACG", "ACT"], ["GTA", "GTC", "GTG", "GTT"], ["TGG"], ["TAC", "TAT"]]
grantham = {'AA':0, 'AC':195, 'AD':126, 'AE':107, 'AF':113, 'AG':60, 'AH':86, 'AI':94, 'AK':106, 'AL':96, 'AM':84, 'AN':111, 'AP':27, 'AQ':91, 'AR':112, 'AS':99, 'AT':58, 'AV':64, 'AW':148, 'AY':112, 'CA':195, 'CC':0, 'CD':154, 'CE':170, 'CF':205, 'CG':159, 'CH':174, 'CI':198, 'CK':202, 'CL':198, 'CM':196, 'CN':139, 'CP':169, 'CQ':154, 'CR':180, 'CS':112, 'CT':149, 'CV':192, 'CW':215, 'CY':194, 'DA':126, 'DC':154, 'DD':0, 'DE':45, 'DF':177, 'DG':94, 'DH':81, 'DI':168, 'DK':101, 'DL':172, 'DM':160, 'DN':23, 'DP':108, 'DQ':61, 'DR':96, 'DS':65, 'DT':85, 'DV':152, 'DW':181, 'DY':160, 'EA':107, 'EC':170, 'ED':45, 'EE':0, 'EF':140, 'EG':98, 'EH':40, 'EI':134, 'EK':56, 'EL':138, 'EM':126, 'EN':42, 'EP':93, 'EQ':29, 'ER':54, 'ES':80, 'ET':65, 'EV':121, 'EW':152, 'EY':122, 'FA':113, 'FC':205, 'FD':177, 'FE':140, 'FF':0, 'FG':153, 'FH':100, 'FI':21, 'FK':102, 'FL':22, 'FM':28, 'FN':158, 'FP':114, 'FQ':116, 'FR':97, 'FS':155, 'FT':103, 'FV':50, 'FW':40, 'FY':22, 'GA':60, 'GC':159, 'GD':94, 'GE':98, 'GF':153, 'GG':0, 'GH':98, 'GI':135, 'GK':127, 'GL':138, 'GM':127, 'GN':80, 'GP':42, 'GQ':87, 'GR':125, 'GS':56, 'GT':59, 'GV':109, 'GW':184, 'GY':147, 'HA':86, 'HC':174, 'HD':81, 'HE':40, 'HF':100, 'HG':98, 'HH':0, 'HI':94, 'HK':32, 'HL':99, 'HM':87, 'HN':68, 'HP':77, 'HQ':24, 'HR':29, 'HS':89, 'HT':47, 'HV':84, 'HW':115, 'HY':83, 'IA':94, 'IC':198, 'ID':168, 'IE':134, 'IF':21, 'IG':135, 'IH':94, 'II':0, 'IK':102, 'IL':5, 'IM':10, 'IN':149, 'IP':95, 'IQ':109, 'IR':97, 'IS':142, 'IT':89, 'IV':29, 'IW':61, 'IY':33, 'KA':106, 'KC':202, 'KD':101, 'KE':56, 'KF':102, 'KG':127, 'KH':32, 'KI':102, 'KK':0, 'KL':107, 'KM':95, 'KN':94, 'KP':103, 'KQ':53, 'KR':26, 'KS':121, 'KT':78, 'KV':97, 'KW':110, 'KY':85, 'LA':96, 'LC':198, 'LD':172, 'LE':138, 'LF':22, 'LG':138, 'LH':99, 'LI':5, 'LK':107, 'LL':0, 'LM':15, 'LN':153, 'LP':98, 'LQ':113, 'LR':102, 'LS':145, 'LT':92, 'LV':32, 'LW':61, 'LY':36, 'MA':84, 'MC':196, 'MD':160, 'ME':126, 'MF':28, 'MG':127, 'MH':87, 'MI':10, 'MK':95, 'ML':15, 'MM':0, 'MN':142, 'MP':87, 'MQ':101, 'MR':91, 'MS':135, 'MT':81, 'MV':21, 'MW':67, 'MY':36, 'NA':111, 'NC':139, 'ND':23, 'NE':42, 'NF':158, 'NG':80, 'NH':68, 'NI':149, 'NK':94, 'NL':153, 'NM':142, 'NN':0, 'NP':91, 'NQ':46, 'NR':86, 'NS':46, 'NT':65, 'NV':133, 'NW':174, 'NY':143, 'PA':27, 'PC':169, 'PD':108, 'PE':93, 'PF':114, 'PG':42, 'PH':77, 'PI':95, 'PK':103, 'PL':98, 'PM':87, 'PN':91, 'PP':0, 'PQ':76, 'PR':103, 'PS':74, 'PT':38, 'PV':68, 'PW':147, 'PY':110, 'QA':91, 'QC':154, 'QD':61, 'QE':29, 'QF':116, 'QG':87, 'QH':24, 'QI':109, 'QK':53, 'QL':113, 'QM':101, 'QN':46, 'QP':76, 'QQ':0, 'QR':43, 'QS':68, 'QT':42, 'QV':96, 'QW':130, 'QY':99, 'RA':112, 'RC':180, 'RD':96, 'RE':54, 'RF':97, 'RG':125, 'RH':29, 'RI':97, 'RK':26, 'RL':102, 'RM':91, 'RN':86, 'RP':103, 'RQ':43, 'RR':0, 'RS':110, 'RT':71, 'RV':96, 'RW':101, 'RY':77, 'SA':99, 'SC':112, 'SD':65, 'SE':80, 'SF':155, 'SG':56, 'SH':89, 'SI':142, 'SK':121, 'SL':145, 'SM':135, 'SN':46, 'SP':74, 'SQ':68, 'SR':110, 'SS':0, 'ST':58, 'SV':124, 'SW':177, 'SY':144, 'TA':58, 'TC':149, 'TD':85, 'TE':65, 'TF':103, 'TG':59, 'TH':47, 'TI':89, 'TK':78, 'TL':92, 'TM':81, 'TN':65, 'TP':38, 'TQ':42, 'TR':71, 'TS':58, 'TT':0, 'TV':69, 'TW':128, 'TY':92, 'VA':64, 'VC':192, 'VD':152, 'VE':121, 'VF':50, 'VG':109, 'VH':84, 'VI':29, 'VK':97, 'VL':32, 'VM':21, 'VN':133, 'VP':68, 'VQ':96, 'VR':96, 'VS':124, 'VT':69, 'VV':0, 'VW':88, 'VY':55, 'WA':148, 'WC':215, 'WD':181, 'WE':152, 'WF':40, 'WG':184, 'WH':115, 'WI':61, 'WK':110, 'WL':61, 'WM':67, 'WN':174, 'WP':147, 'WQ':130, 'WR':101, 'WS':177, 'WT':128, 'WV':88, 'WW':0, 'WY':37, 'YA':112, 'YC':194, 'YD':160, 'YE':122, 'YF':22, 'YG':147, 'YH':83, 'YI':33, 'YK':85, 'YL':36, 'YM':36, 'YN':143, 'YP':110, 'YQ':99, 'YR':77, 'YS':144, 'YT':92, 'YV':55, 'YW':37, 'YY':0}



        

################################################ SIMULATION FUNCTIONS #####################################################

def simulate(f, seqfile, tree, mu, kappa, length, beta=None):
    ''' Simulate single partition according to either codon or mutsel model (check beta value for which model).
    '''
    try:
        my_tree = readTree(file = tree)
    except:
        my_tree = readTree(tree = tree) 
        
    model = Model()
    if beta:
        params = {'alpha':1.0, 'beta':float(beta), 'mu': {'AC': mu, 'AG': mu*kappa, 'AT': mu, 'CG': mu, 'CT': mu*kappa, 'GT': mu}}
        params['stateFreqs'] = f
        model.params = params
        mat = mechCodon_MatrixBuilder(model)
    else:
        params = {'alpha':1.0, 'beta':1.0, 'mu': {'AC': mu, 'CA':mu, 'AG': mu*kappa, 'GA':mu*kappa, 'AT': mu, 'TA':mu, 'CG': mu, 'GC':mu, 'CT': mu*kappa, 'TC':mu*kappa, 'GT': mu, 'TG':mu}}
        params['stateFreqs'] = f
        model.params = params
        mat = mutSel_MatrixBuilder(model)
    
    model.Q = mat.buildQ()
    partitions = [(length, model)]        
    myEvolver = StaticEvolver(partitions = partitions, tree = my_tree, outfile = seqfile)
    myEvolver.sim_sub_tree(my_tree)
    myEvolver.writeSequences()

def setFreqs(freqClass, numaa, bias):
    ''' Returns codon frequencies and gc content '''
    
    assert(numaa != 1), "Can't only have a single amino acid! omega will always be zero. we're not interested in that.'''
    
    if freqClass == 'exp':
        userFreq = setFreqDict(numaa)
        fobj = UserFreqs(by = 'amino', type = 'codon', freqs = userFreq, bias = bias)
    else:
        aalist = generateAAlist(numaa) 
        if freqClass == 'equal':
            fobj = EqualFreqs(by = 'amino', type = 'codon', restrict = aalist, bias = bias)
        elif freqClass == 'random':
            fobj = RandFreqs(by = 'amino', type = 'codon', restrict = aalist, bias = bias)              
        else:
            raise AssertionError("Bad freqClass specification. Byebye.")
    codonFreq = fobj.calcFreqs()
    nucFreq = fobj.codon2nuc()
    return codonFreq, nucFreq[1] + nucFreq[2], "".join(userFreq.keys())


def setFreqDict(size):
    ''' Generate a dictionary of exponentially distributed amino acid frequencies.
        size = number of amino acids. please never use 1.
    '''
    rawfreqs = setBoltzFreqs(size)
    aalist = generateAAlist(size)
    return dict(zip(aalist, rawfreqs))
    

def setBoltzFreqs(size, beta = 1.5):
    ''' Use Boltzmann distribution to get amino acid frequencies for a certain number of amino acids'''
    ddg_values = np.random.normal(size = size) 
    numer_list = np.zeros(size)
    denom = 0.
    for d in range(size):
        val = np.exp(-1. * beta * ddg_values[d])
        denom += val
        numer_list[d] = val
    return numer_list/(np.sum(numer_list))     


def checkGrantham(aalist, size, cutoff):
    ''' Given a list of amino acids, ensure that they can reasonably co-occur based on Grantham indices.
        The mean similarity score should be <= 100. Based on mean distance for non-ID amino acids is 99.97. Good enough for now.
    '''
    scores = 0.
    total = float(size*(size-1)/2)
    for i in range(len(aalist)):
        for j in range(i, len(aalist)):
            if i != j:
                key = "".join(sorted(aalist[i] + aalist[j]))
                scores += grantham[key]
    if scores/total <= float(cutoff):
        return True
    else:
        return False
 
 
def generateAAlist(size):
    ''' Generate a list of size of reasonable co-occuring amino acids.
        Ensure acceptable choices by making mean pairwise Grantham <=100.
        If size>=10, no need to check Grantham since so many are allowed, properties probably don't matter so much.
    '''
    list_is_ok = False
    while not list_is_ok:
        aalist = []
        amino = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
        if size == 20:
            return amino
        else:         
            for i in range(size):
                n = randint(0,len(amino)-1)
                aalist.append(amino[n])
                amino.pop(n)
            if size >= 10:
                list_is_ok = True
            else:
                list_is_ok = checkGrantham(aalist, 100.)  
    return aalist

################################################# HYPHY-RELATED FUNCTIONS ############################################################
def runhyphy(batchfile, matrix_name, seqfile, treefile, cpu, kappa = 1.0):
    ''' pretty specific function.'''
    
  
    # Set up sequence file with tree
    setuphyphy1 = "cp "+seqfile+" temp.fasta"
    setup1 = subprocess.call(setuphyphy1, shell = True)
    assert(setup1 == 0), "couldn't create temp.fasta"
    setuphyphy2 = "cat "+treefile+" >> temp.fasta"
    setup2 = subprocess.call(setuphyphy2, shell = True)
    assert(setup2 == 0), "couldn't add tree to hyphy infile"
    

    hyf = freq2hyphy(eqf)
    setuphyphy3 = "sed 's/MYFREQUENCIES/"+hyf+"/g' "+batchfile+" > run.bf"
    setup3 = subprocess.call(setuphyphy3, shell = True)
    assert(setup3 == 0), "couldn't properly add in frequencies"
    
    # Set up kappa
    if kappa != 'free':
        sedkappa = "sed -i 's/k/"+str(kappa)+"/g' matrices.mdl"
        runsedkappa = subprocess.call(sedkappa, shell=True)
        assert(runsedkappa == 0), "couldn't set up kappa"
    
    # Set up matrix, within run.bf
    setuphyphy4 = "sed -i 's/MYMATRIX/"+matrix_name+"/g' run.bf"
    setup4 = subprocess.call(setuphyphy4, shell = True)
    assert(setup4 == 0), "couldn't properly define matrix"
    

    # Run hyphy
    hyphy = "./HYPHYMP run.bf CPU="+cpu+" > hyout.txt"
    runhyphy = subprocess.call(hyphy, shell = True)
    assert (runhyphy == 0), "hyphy fail"
    
    # grab hyphy output
    if matrix_name == 'GY94':
        return parseHyphyGY94('hyout.txt')
    else:
        return parseHyphyMG94('hyout.txt')


def parseHyphyGY94(file):
    hyout = open(file, 'r')
    hylines = hyout.readlines()
    hyout.close()
    for line in hylines:
        findw = re.search("^w=(\d+\.*\d*)", line)
        if findw:
            hyphy_w = findw.group(1)
    return float(hyphy_w)

def parseHyphyMG94(file):
    hyout = open(file, 'r')
    hylines = hyout.readlines()
    hyout.close()
    for line in hylines:
        finda = re.search("^a=(\d+\.*\d*)", line)
        findb = re.search("^b=(\d+\.*\d*)", line)
        if finda:
             hyphy_alpha = finda.group(1)
        if findb:
            hyphy_beta = findb.group(1)
    return float(hyphy_alpha), float(hyphy_beta)

def freq2hyphy(f):
    ''' Convert codon frequencies to a form hyphy can use. '''
    hyphy_f = "{"
    for freq in f:
        hyphy_f += "{"
        hyphy_f += str(freq)
        hyphy_f += "},"
    hyphy_f = hyphy_f[:-1]
    hyphy_f += "}"
    return hyphy_f

def setFreqsForHyphy(codonFreq, aminos_used, bias):
    ''' Build input hyphy frequencies.
        In a **highly** hackish way, will also manipulate the equal freqs by any bias factors.
    '''
    eqf = np.zeros(61)
    for i in range(61):
        eqf[i] = 1./61.
    if bias is not None:
        for amino in aminos_used:
            aa_ind = amino_acids.index(amino)
            syn_codons = genetic_code[aa_ind]
            cfreqs = []
            cinds=[]
            for syn in syn_codons:
                cinds.append(codons.index(syn))
                cfreqs.append(codonFreq[ codons.index(syn) ])

            pref_factor = (1. + bias)
            unpref_factor = 1. - bias/(len(syn_codons)-1.)
            maxFreq = max(cfreqs)
            for c in range(len(syn_codons)):
                if cfreqs[c] < maxFreq:
                    eqf[cinds[c]] = 1./61. * unpref_factor
                elif cfreqs[c] == maxFreq:
                    eqf[cinds[c]] = 1./61. * pref_factor   
    return freq2hyphy(eqf)

####################################### OMEGA DERIVATION FUNCTIONS ######################################

def deriveOmega(codonFreq, mu_dict = {'AT':1.0, 'AC':1.0, 'AG':1.0, 'CG':1.0, 'CT':1.0, 'GT':1.0}):
    ''' Derive an omega using codon frequencies. Single calculation. Get numerator, get denominator, divide once at end.
        Default mutational scheme is all mu's are equal.
        Requires symmetric mutational scheme.
    ''' 
    
    kN = 0.; nN = 0.; kS = 0.; nS = 0.
    
    # codon indices which do not have zero frequency
    nonZero = getNonZeroFreqs(codonFreq)
    
    # Calculations
    for i in nonZero:
        codon = codons[i]
        # Nonsynonymous calculation
        num, den = calcNS(codon, codonFreq, i, nslist, mu_dict)
        kN += num
        nN += den
        
        # Synonymous calculation
        num, den = calcNS(codon, codonFreq, i, slist, mu_dict)
        kS += num
        nS += den
   
    assert( nS != 0. and nN != 0.), "Omega derivation indicates no evolution, maybe?"
    dN = kN/nN
    dS = kS/nS
    return dN, dS, dN/dS



def getNonZeroFreqs(freq):
    ''' Return indices whose frequencies are not 0.'''
    nonZero = []
    for i in range(len(freq)):
        if freq[i] > zero:
            nonZero.append(i)
    return nonZero

def getNucleotideDiff(source, target):
    for i in range(3):
        if source[i] != target[i]:    
            return "".join(sorted(source[i]+target[i]))
    
def calcFix(fi, fj):
    if fi == fj:
        return 1.
    elif fi == 0.  or fj == 0.:
        return 0.
    else:
        return (np.log(fj) - np.log(fi)) / (1 - fi/fj)

def calcNS(codon, codonFreq, i, list, mu_dict):
    numer = 0.
    denom = 0.
    fix_sum=0.
    freq = codonFreq[i]
    for other_codon in list[i]:
        if codon != other_codon:
            diff = getNucleotideDiff(codon,other_codon)
            tempfreq = codonFreq[codons.index(other_codon)]
            fix_sum += calcFix( float(freq), float(tempfreq) ) * mu_dict[diff]                  
            denom += codonFreq[i]
    numer += fix_sum*codonFreq[i]
    return numer, denom
#########################################################################################









 

############################# NEI-GOJOBORI FUNCTIONS ##################################
def run_neigojo(seqfile):
    ''' Get omega using counting method '''
    from mutation_counter import *
    from site_counter import *
    
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
    return dN/dS, np.mean(ns_mut), np.mean(s_mut)



############################ PAML-RELATED FUNCTIONS ###############################
def runpaml(seqfile, codonFreq = "3", initw = 0.4):
    
    # Set up sequence file
    setuppaml1 = "cp "+seqfile+" temp.fasta"
    setup1 = subprocess.call(setuppaml1, shell = True)
    assert(setup1 == 0), "couldn't create temp.fasta"
    
    # Set up initial omega
    setuppaml2 = 'sed "s/MYINITIALW/'+str(initw)+'/g" codeml_raw.txt > codeml.ctl' 
    setup2 = subprocess.call(setuppaml2, shell = True)
    assert(setup2 == 0), "couldn't set paml initial w"
    
    # Set up codon frequency specification NOTE: 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
    setuppaml3 = 'sed -i "s/MYCODONFREQ/'+str(codonFreq)+'/g" codeml.ctl' 
    setup3 = subprocess.call(setuppaml3, shell = True)
    assert(setup3 == 0), "couldn't set paml codon frequencies"
    
    # Run paml
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

