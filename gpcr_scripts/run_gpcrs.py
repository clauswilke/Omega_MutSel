import numpy as np
from Bio import SeqIO
import os
import sys
import subprocess
path_to_mutsel = "/Users/sjspielman/Research/MutSel/Simulator/src/"
sys.path.append(path_to_mutsel)
from stateFreqs import *

codons      = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
nucindex    = {'A':0, 'C':1, 'G':2, 'T':3}
purines     = ["A", "G"]
pyrims      = ["C", "T"]


def calc_nuc_freqs(alnfile):
    ''' From an alignment file, return the empirical nucleotide frequencies. Note that as the input alignment file is actually a hyphy input file (has tree), we have to remove it. '''
    
    nuc_freqs = np.zeros(4)
    
    with open(alnfile, 'r') as infile:
        for line in infile:
            if "(" not in line and ">" not in line:
                    nuc_freqs[0] += float(line.count('A'))
                    nuc_freqs[1] += float(line.count('C'))
                    nuc_freqs[2] += float(line.count('G'))
                    nuc_freqs[3] += float(line.count('T'))
    nuc_freqs /= np.sum(nuc_freqs)
    return nuc_freqs
    
    
    

def codon_to_nuc(codon_freqs):
    ''' Calculate positional nucleotide frequencies from codon frequencies. '''
    pos_nuc_freqs = np.zeros( [3, 4] ) # row is position, column in nucleotide
    for i in range(3):
        count = 0
        for codon in codons:
            pos_nuc_freqs[i][nucindex[codon[i]]] += codon_freqs[count]
            count += 1
        assert( abs(np.sum(pos_nuc_freqs[i]) - 1.) < 1e-8 ), "bad positional nucleotide frequencies"
    return pos_nuc_freqs


def calc_f3x4_freqs(pos_nuc_freqs):
    ''' Compute F3x4 codon frequencies from positional nucleotide frequencies. '''
    
    f3x4 = np.ones(61)
    pi_stop = pos_nuc_freqs[0][nucindex['T']]*pos_nuc_freqs[1][nucindex['A']]*pos_nuc_freqs[2][nucindex['G']]  +  pos_nuc_freqs[0][nucindex['T']]*pos_nuc_freqs[1][nucindex['G']]*pos_nuc_freqs[2][nucindex['A']]  +  pos_nuc_freqs[0][nucindex['T']]*pos_nuc_freqs[1][nucindex['A']]*pos_nuc_freqs[2][nucindex['A']] 

    for i in range(61):
        codon = codons[i]
        for j in range(3):
            f3x4[i] *= pos_nuc_freqs[j][ nucindex[codon[j]] ]
    f3x4 /= (1. - pi_stop)
    assert( abs(np.sum(f3x4) - 1.) < 1e-8), "Could not properly caluclate (C)F3x4 frequencies."
    return f3x4   



def calc_f1x4_freqs(nuc_freqs):
    ''' Compute F1x4 codon frequencies from nucleotide frequencies. '''
    f1x4 = np.ones(61)
    pi_stop = nuc_freqs[nucindex['T']]*nuc_freqs[nucindex['A']]*nuc_freqs[nucindex['G']] + nuc_freqs[nucindex['T']]*nuc_freqs[nucindex['G']]*nuc_freqs[nucindex['A']] + nuc_freqs[nucindex['T']]*nuc_freqs[nucindex['A']]*nuc_freqs[nucindex['A']]
    for i in range(61):
        codon = codons[i]
        for j in range(3):
            f1x4[i] *= nuc_freqs[nucindex[codon[j]]]
    f1x4 /= (1. - pi_stop)
    assert( abs(np.sum(f1x4) - 1.) < 1e-8), "Could not properly caluclate F1x4 frequencies."
    return f1x4
    

def calc_cf3x4_freqs(pos_nuc_freqs):
    assert( pos_nuc_freqs.shape == (4,3)), "You need to provide hyphy with the transpose of your pos_nuc_freqs matrix!!"
    
    # Create positional nucleotide matrix string to sed into hyphy cf3x4 batchfile
    posnuc = '{'
    for row in pos_nuc_freqs:
        p = str(row).replace("[","").replace("]","").replace("","").strip()
        p = re.sub("\s+", ",", p)
        posnuc += '{' + p + '},'
    posnuc = posnuc[:-1]+'};'

    # sed into hyphy
    sed_cf3x4 = subprocess.call("sed 's/INSERT_POS_FREQS/"+posnuc+"/g' cf3x4_raw.bf > cf3x4.bf", shell=True)
    assert(sed_cf3x4 == 0), "Couldn't sed positional nucleotide frequencies into cf3x4.bf"
    
    # run hyphy 
    run_hyphy = subprocess.call("HYPHYMP cf3x4.bf > cf3x4.out", shell=True)
    assert(run_hyphy == 0), "Couldn't get hyphy to run to compute cf3x4"

    # parse hyphy output file and save new positional frequencies
    cf3x4_pos_freqs = np.zeros([3,4])
    with open("cf3x4.out", "r") as file:
        hystring = file.readlines()
    for i in range(4):
        line = str(hystring[i+1]).replace("{","").replace("}","").rstrip()
        freqs = line.split(',')
        for j in range(3):
            cf3x4_pos_freqs[j][i] = float(freqs[j])
    assert( np.allclose( np.sum(cf3x4_pos_freqs, axis=1), np.ones(3)) ), "Bad CF3x4 positional frequencies."

    # Finally convert these cf3x4 positional frequencies to codon frequencies
    return calc_f3x4_freqs(cf3x4_pos_freqs)





def get_nuc_diff(source, target):
    diff = ''
    for i in range(3):
        if source[i] != target[i]: 
            diff += source[i]+target[i]
    return diff



def is_TI(source, target):
    check1 = source in purines and target in purines
    check2 = source in pyrims and target in pyrims
    if check1 or check2:
        return True
    else:
        return False
    
        
def build_fnuc(nuc_freqs):
    ''' Create hyphy matrix which use target nucleotide frequencies (not codon frequencies). ''' 

    matrix = 'Fnuc = {61, 61, \n'    
    for i in range(61):
        source = codons[i]
        for j in range(61):
            target = codons[j]
                
            diff = get_nuc_diff(source, target)
            if len(diff) == 2:
                target_nuc_index = nucindex[diff[1]]   
                    
                # Create string for matrix element
                element = '{' + str(i) + ',' + str(j) + ',t*'                    
                if is_TI(diff[0], diff[1]):
                    element += 'k*'
                if codon_dict[source] != codon_dict[target]:
                    element += 'w*'
                matrix_data += element + str(nuc_freqs[target_nuc_index]) + '}\n'

    # And save to file.
    with open('fnuc.mdl', 'w') as outf:
        outf.write(matrix+ '};\n\n\n')
        
        
        
def array_to_hyphy_freq(f):
    ''' Convert a python/numpy list/array of codon frequencies to a hyphy frequency string to directly write into a batchfile. '''
    hyphy_f = "{"
    for freq in f:
        hyphy_f += "{" + str(freq) + "},"
    hyphy_f = hyphy_f[:-1]
    hyphy_f += "};"
    return hyphy_f
    
    

def create_batchfile(basefile, outfile, f61, f1x4, f3x4, cf3x4):
    ''' sed in the frequency specifications to create an output batchfile from the base/raw batchfile framework.'''
    cp_batch = subprocess.call("cp " + basefile + " " + outfile, shell=True)
    assert(cp_batch == 0), "couldn't copy batchfile"
    shutil.copy(basefile, outfile)
    flist = ['f61', 'f1x4', 'f3x4', 'cf3x4']
    for i in range( len(flist) ):
        hyf = eval(flist[i])
        insert = flist[i].upper()
        setuphyphyf = "sed -i 's/INSERT_"+insert+"/"+hyf+"/g' " + outfile
        setupf = subprocess.call(setuphyphyf, shell = True)
        assert(setupf == 0), "couldn't properly add in frequencies"
        


def run_hyphy(batchfile, seqfile, fspecs):

    # Set up sequence file with tree
    shutil.copy(seqfile, "temp.fasta")

    # Run hyphy.
    runhyphy = subprocess.call("./HYPHYMP " + batchfile, shell = True)
    assert (runhyphy == 0), "hyphy fail"
    
    lnliks = np.zeros(len(fspecs)) # log likelihood values
    count = 0
    for suffix in fspecs:
        file = suffix + '_hyout.txt'  
        lnliks[count] = parse_output_GY94(file)
        count += 1
    return lnliks
     
    
def parse_output_GY94(file):
    with open(file, 'r') as hyout:
        hylines = hyout.readlines()
    lnlik = None;
    for line in hylines:
        findlk = re.search("^Likelihood Function's Current Value = (-\d+\.*\d*)", line)
        if findlk:
            lnlik = float(findlk.group(1))
    assert(lnlik is not None),   "Couldn't retrieve log likelihood from hyphy output file."
    return lnlik
        
        

def main():
    rep = sys.argv[1]
    seqfile = "gpcr"+rep+".fasta"
    outfile = "output"+rep+".txt"
    
    # Use MutSel simulator code to get f61_freqs, nuc_freqs
    # First need to have a temporary file w/out the tree at the bottom
    rm_tree = subprocess.call("sed '$d' " + seqfile + " > notree.fasta", shell=True)
    assert(rm_tree == 0), "Could not remove tree from fasta file."
    freqObject = ReadFreqs(by ='codon', file = "notree.fasta")
    f61_freqs = freqObject.calcFreqs(type = 'codon')
    f61 = array_to_hyphy_freq(f61_freqs)
    nuc_freqs = freqObject.calcFreqs(type = 'nuc')
    
    # F1x4, F3x4, CF3x4
    f1x4 = array_to_hyphy(calc_f1x4_freqs(nuc_freqs))
    pos_nuc_freqs = codon_to_nuc(f61_freqs)
    f3x4 = array_to_hyphy(calc_f3x4_freqs(pos_nuc_freqs))
    cf3x4 = array_to_hyphy_freq( calc_cf3x4_freqs(pos_nuc_freqs.T) ) 
    
    # Save batchfile
    create_batchfile("raw_batchfile.bf", "batchfile.bf", f61, f1x4, f3x4, cf3x4)
    
    # Fnuc 
    build_fnuc_matrices(nuc_freqs)
    
    # Call hyphy
    numparams = np.array([62., 11., 11., 5., 5.])
    lnliks = run_hyphy("batchfile.bf", seqfile, ['f61', 'f3x4', 'cf3x4', 'f1x4', 'fnuc'])
    AIC = 2*(numparms - lnliks)
    
    with open(outfile, 'w') as outf:
        outf.write(rep + '\t' + AIC[0] + '\t' + AIC[1] + '\t' + AIC[2] + '\t' + AIC[3] + '\t' + AIC[4] + '\t' + AIC[5] + '\n')
        
        
        

main()    
    
    




















        
        
