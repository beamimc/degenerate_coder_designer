#! /usr/bin/env python

"""
@author: beamimc

"""
import pandas as pd
import math 
import time
import os 

from dcg import gcm_multiple
# global variables
aas = ['I','M','T','N','K','S','R','L','P','H','Q','V',
       'A','D','E','G','F','Y','C','W','STOP']
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
        'TAC':'Y', 'TAT':'Y', 'TAA':'STOP', 'TAG':'STOP',
        'TGC':'C', 'TGT':'C', 'TGA':'STOP', 'TGG':'W'
        }
deg_nucl={
         'A':['A'],
         'T':['T'],
         'C':['C'],
         'G':['G'],
         'B':['G','T','C'],
         'D':['G','A','T'],
         'H':['A','T','C'],
         'K':['G','T'],
         'M':['A','C'],
         'N':['A','C','G','T'],
         'R':['A','G'],
         'S':['G','C'],
         'V':['A','C','G'],
         'W':['A','T'],
         'Y':['C','T']
        }

def deg_codon_to_aa(deg_codon):
    def get_all_real_codons(deg_codon):
         posibilities = []
         for n1 in deg_nucl[deg_codon[0]]:
             for n2 in deg_nucl[deg_codon[1]]:
                 for n3 in deg_nucl[deg_codon[2]]:
                     posibilities.append(n1+n2+n3)
         return posibilities
    posibilities = get_all_real_codons(deg_codon)
    df_aas = pd.DataFrame(columns = ['count'], index =aas).fillna(0)   
    for codon in posibilities:    
        aa = codon_aa[codon]
        df_aas.loc[aa,'count'] +=1

    ##normalize proportions 
    prop = (df_aas[df_aas['count']>0]['count'].tolist())
    gcm = gcm_multiple(prop)
    df_aas.loc[:,'count'] = df_aas.loc[:,'count'] / gcm 

    dc = df_aas['count'].to_dict()
    dc['deg_codon'] = deg_codon
    return dc

def create_db_codon():
    def get_degenerated_codons():
        ##15x15x15 different posible deg_codons
        all_deg_codons=[]
        bases = list(deg_nucl.keys())
        for n1 in bases:
             for n2 in bases:
                 for n3 in bases:
                     all_deg_codons.append(n1+n2+n3)
        return all_deg_codons
    
    #get all posible combinations into df
    all_deg_codons = get_degenerated_codons()
    columns = ['deg_codon']
    columns.extend(aas)
    df = pd.DataFrame(columns=columns)
    #generate df of codon and aa they code proportionally
    for deg_codon in all_deg_codons:
        dc_aa = deg_codon_to_aa(deg_codon)
        df = df.append(dc_aa,ignore_index=True)

    return df

if __name__ == '__main__':
    startTime = time.time()
    current_path = os.getcwd()
    print('Creating DB')  
    df = create_db_codon()
    df.to_csv(current_path+'\\deg_codons_DB.csv',index = False)
    print(f'File saved at {current_path}\\deg_codons_DB.csv')
    executionTime = (time.time() - startTime)