#! /usr/bin/env python
"""
@author: beamimc

"""
# global variables

#list of aminoacids. X is STOP 
aas = ['I','M','T','N','K','S','R','L','P','H','Q','V',
       'A','D','E','G','F','Y','C','W','X']
#diccionary with codon (just with ATCG) and the aminoacid they code
codon_aa={
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'X', 'TAG':'X',
        'TGC':'C', 'TGT':'C', 'TGA':'X', 'TGG':'W'
        }
#diccionary with degenerate nucleotides and the list of nucleotides they work as 
deg_nucl={
         'A':'A',
         'T':'T',
         'C':'C',
         'G':'G',
         'B':'GTC',
         'D':'GAT',
         'H':'ATC',
         'K':'GT',
         'M':'AC',
         'N':'ACGT',
         'R':'AG',
         'S':'GC',
         'V':'ACG',
         'W':'AT',
         'Y':'CT'
        }
#relative paths to csv 
#these csv wewre created with the script create_deg_codons_DB.py
CODON_DB ='../datasets/deg_codons_DB.csv'
CLEAN_CODON_DB = '../datasets/deg_codons_DB_clean.csv'
