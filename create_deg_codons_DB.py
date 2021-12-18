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
       'A','D','E','G','F','Y','C','W','X']
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
         # print(posibilities)
         return posibilities
    posibilities = get_all_real_codons(deg_codon)
    df_aas = pd.DataFrame(columns = ['count'], index =aas).fillna(0)  
    
    ## remove redundant deg_codons ej TTR codes L 2 times, so we only leave the 
    ## redundant normal codons TTA, TTG, 
    ok = True
    if len(posibilities)>1:
        coded = []
        for codon in posibilities:    
            coded.append(codon_aa[codon])
        # print(coded)
        if len(set(coded)) == 1 and deg_codon not in set(codon_aa.keys()):
            ok = False
    ## remove redundant deg_codons ej TTA and TTG both code L 1 times
    ## we only keep 1 option
    deg_codon not in set(codon_aa.keys())
    
    for codon in posibilities:    
        aa = codon_aa[codon]
        df_aas.loc[aa,'count'] +=1

    ##normalize proportions 
    prop = (df_aas[df_aas['count']>0]['count'].tolist())
    gcm = gcm_multiple(prop)
    df_aas.loc[:,'count'] = df_aas.loc[:,'count'] / gcm 

    dc = df_aas['count'].to_dict()
    dc['deg_codon'] = deg_codon
    return dc, ok 


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
    aux = {}
    for aa in aas:
        aux[aa]=0
    #get all posible combinations into df
    all_deg_codons = get_degenerated_codons()
    columns = ['deg_codon']
    columns.extend(aas)
    df = pd.DataFrame(columns=columns)
    df_clean = pd.DataFrame(columns=columns)
    #generate df of codon and aa they code proportionally
    for deg_codon in all_deg_codons:
        dc_aa, ok = deg_codon_to_aa(deg_codon)
        df = df.append(dc_aa,ignore_index=True)
        if ok:
            df_clean = df_clean.append(dc_aa,ignore_index=True)
 
    ### remove codons that do the same 
    df_clean = df_clean.drop_duplicates(
          subset = aas,
          keep = 'first').reset_index(drop = True)
    return df, df_clean

if __name__ == '__main__':
    current_path = os.getcwd()
    directory = 'datasets'
    
    if not os.path.exists(directory):
        os.mkdir(directory)
        
    print('Creating DB and cleaned DB')  
    filename = 'deg_codons_DB.csv'
    filename_clean = 'deg_codons_DB_clean.csv'
    
    df, df_clean = create_db_codon()
    
    df.to_csv(current_path+'\\'+directory+'\\'+filename,index = False)
    print(f'File saved at {current_path}\\{directory}\\{filename}')
    
    df_clean.to_csv(current_path+'\\'+directory+'\\'+filename_clean,index = False)
    print(f'File saved at {current_path}\\{directory}\\{filename_clean}')

