import numpy as np
from Bio import AlignIO
import sys
import re
import shutil
import subprocess

codon_dict  = {"AAA":"K", "AAC":"N", "AAG":"K", "AAT":"N", "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T", "AGA":"R", "AGC":"S", "AGG":"R", "AGT":"S", "ATA":"I", "ATC":"I", "ATG":"M", "ATT":"I", "CAA":"Q", "CAC":"H", "CAG":"Q", "CAT":"H", "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P", "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R", "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L", "GAA":"E", "GAC":"D", "GAG":"E", "GAT":"D", "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A", "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G", "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V", "TAC":"Y", "TAT":"Y", "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S", "TGC":"C", "TGG":"W", "TGT":"C", "TTA":"L", "TTC":"F", "TTG":"L", "TTT":"F"}
codons      = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
nucindex    = {'A':0, 'C':1, 'G':2, 'T':3}
purines     = ["A", "G"]
pyrims      = ["C", "T"]
  
    
    
    
def calc_nuc_freqs(infile):
    pos_nuc_freqs = np.zeros([4,3])
    
    full_seq = ''
    aln = AlignIO.read(infile, "fasta")
    for record in aln:
        full_seq += str(record.seq)  
        
    # Calc positional nucleotide frequencies      
    for i in range(0,len(full_seq),3):
        codon = full_seq[i:i+3]
        if codon in codons:
            for x in range(3):
                pos_nuc_freqs[nucindex[codon[x]]][x] += 1.
        
    pos_nuc_freqs /= np.sum(pos_nuc_freqs, axis=0)
    nuc_freqs = np.mean(pos_nuc_freqs, axis=1)
    
    assert( abs(np.sum(nuc_freqs) - 1.) < 1e-8), "bad global nucleotide frequencies"
    assert( np.allclose(np.sum(pos_nuc_freqs, axis=0),np.ones(3))), "bad positional nucleotide frequencies"
    return nuc_freqs, pos_nuc_freqs


        
def calc_f3x4_freqs(pos_nuc_freqs):
    ''' Compute F3x4 codon frequencies from positional nucleotide frequencies. '''
    f3x4 = np.ones(61)
    pi_stop = pos_nuc_freqs[nucindex['T']][0]*pos_nuc_freqs[nucindex['A']][1]*pos_nuc_freqs[nucindex['G']][2]  +  pos_nuc_freqs[nucindex['T']][0]*pos_nuc_freqs[nucindex['G']][1]*pos_nuc_freqs[nucindex['A']][2]  +  pos_nuc_freqs[nucindex['T']][0]*pos_nuc_freqs[nucindex['A']][1]*pos_nuc_freqs[nucindex['A']][2] 

    for i in range(61):
        codon = codons[i]
        for j in range(3):
            f3x4[i] *= pos_nuc_freqs[ nucindex[codon[j]] ][j]
    f3x4 /= (1. - pi_stop)
    assert( abs(np.sum(f3x4) - 1.) < 1e-8), "Could not properly calculate F3x4 frequencies."
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
        

def is_TI(source, target):
    check1 = source in purines and target in purines
    check2 = source in pyrims and target in pyrims
    if check1 or check2:
        return True
    else:
        return False

def get_nuc_diff(source, target, grab_position = False):
    diff = ''
    position = 5
    for i in range(3):
        if source[i] != target[i]: 
            diff += source[i]+target[i]
            position = i
    if grab_position:
        return diff, position
    else:
        return diff
        
def build_fnuc_matrices(nuc_freqs, pos_nuc_freqs, outfile):
    ''' Create matrices which use target nucleotide frequencies (not codon frequencies). 
        Note: pos means "positional", we use 12 positional nucleotide frequencies. Goes with F3x4 stationary distribution.
              glob means "global", we use 4 global (position-free) nucleotide frequencies. Goes with F1x4 stationary distribution.
    ''' 

    matrix_pos  = 'Fnuc_pos = {61, 61, \n'    
    matrix_glob = 'Fnuc_glob = {61, 61, \n'

    for i in range(61):
        source = codons[i]
        for j in range(61):
            target = codons[j]
                
            diff, x = get_nuc_diff(source, target, grab_position=True)
            if len(diff) == 2:
                assert(len(str(x)) == 1), "Problem with determining nucleotide difference between codons"
                target_index = nucindex[diff[1]]   
                glob = str(nuc_freqs[nucindex[diff[1]]])
                pos  = str(pos_nuc_freqs[nucindex[diff[1]]][x])
    
                # Create string for matrix element
                element = '{' + str(i) + ',' + str(j) + ',t'  

                if is_TI(diff[0], diff[1]):
                    element += '*k'
                if codon_dict[source] != codon_dict[target]:
                    element += '*w'
                
                matrix_pos  += element + '*' + pos + '}\n'
                matrix_glob += element + '*' + glob + '}\n'

    # And save to file.
    with open(outfile, 'w') as outf:
        outf.write(matrix_pos + '};\n\n\n' + matrix_glob + '};\n\n\n')


def run_hyphy(batchfile, cpu, seqfile, fspecs):

    # Set up sequence file with tree
    shutil.copy(seqfile, "temp.fasta")

    # Run hyphy.
    runhyphy = subprocess.call("./HYPHYMP CPU="cpu + " " + batchfile, shell = True)
    assert (runhyphy == 0), "hyphy fail"
    
    lnliks = np.zeros(len(fspecs)) # log likelihood values
    numparams = np.zeros(len(fspecs)) # free parameters in model
    
    count = 0
    for suffix in fspecs:
        file = suffix + '_hyout.txt'  
        numparams[count], lnliks[count] = parse_output_GY94(file)
        count += 1
    return 2*(numparams - lnliks)
     
    
def parse_output_GY94(file):
    with open(file, 'r') as hyout:
        hylines = hyout.readlines()
    numparams = float(len(hylines)-1)
    lnlik = None;
    for line in hylines:
        findlk = re.search("^Likelihood Function's Current Value = (-\d+\.*\d*)", line)
        if findlk:
            lnlik = float(findlk.group(1))
    assert(lnlik is not None),   "Couldn't retrieve log likelihood from hyphy output file."
    return numparams, lnlik

        

def main():
    
    rep = sys.argv[1]
    cpu = sys.argv[2]
    seqfile = "gpcr"+rep+".fasta"
    outfile = "output"+rep+".txt"
    
    # Build the Fnuc matrices
    rm_tree = subprocess.call("sed '$d' " + seqfile + " > notree.fasta", shell=True)
    assert(rm_tree == 0), "Could not remove tree from fasta file."
    nuc_freqs, pos_nuc_freqs = calc_nuc_freqs('notree.fasta')
    f1x4 = calc_f1x4_freqs(nuc_freqs)
    f3x4 = calc_f3x4_freqs(pos_nuc_freqs)
    build_fnuc_matrices(nuc_freqs, pos_nuc_freqs, "fnuc.mdl")
        
    # Call hyphy
    AICs = run_hyphy("globalDNDS_gpcr.bf", cpu, seqfile, ['f61', 'f3x4', 'cf3x4', 'f1x4', 'fnuc_pos', 'fnuc_glob'])
    
    outstring = rep
    for aic in AICs:
        outstring += '\t' + str(aic)
    with open(outfile, 'w') as outf:
        outf.write(outstring + '\n')
        
        
        

main()    
    
    




















        
        
