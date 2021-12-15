# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 16:12:00 2021

@author: es0055
"""
import pandas as pd
import math 
import itertools as it

__all__ = [ "generateCodon" ]

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

def get_degenerated_codons():
    ##15x15x15 different posible deg_codons
    all_deg_codons=[]
    bases = list(deg_nucl.keys())
    for n1 in bases:
         for n2 in bases:
             for n3 in bases:
                 all_deg_codons.append(n1+n2+n3)
    return all_deg_codons

def get_all_real_codons(deg_codon):
     posibilities = []
     for n1 in deg_nucl[deg_codon[0]]:
         for n2 in deg_nucl[deg_codon[1]]:
             for n3 in deg_nucl[deg_codon[2]]:
                 posibilities.append(n1+n2+n3)
     return posibilities

def deg_codon_to_aa(deg_codon):
    posibilities = get_all_real_codons(deg_codon)
    df_aas = pd.DataFrame(columns = ['count','probab'], index =aas).fillna(0)   
    # print(df_aas)
    for codon in posibilities:    
        aa = codon_aa[codon]
        df_aas.loc[aa,'count'] +=1
    df_aas.loc[:,'probab'] = df_aas.loc[:,'count']/df_aas.loc[:,'count'].sum()
    coded_aas = df_aas[df_aas['count']>0].index.tolist()
    return df_aas, coded_aas

def check_restriction_1(list_deg_codons):
    ## codons cannot repeat AAs -- 1aa - 1 codon -- 
    all_=[]
    for deg_codon in list_deg_codons:
        _,coded_aas= deg_codon_to_aa(deg_codon)
        # print(deg_codon)
        # print(coded_aas,'\n')
        all_.extend(coded_aas)
    # print(all_)
    if len(all_)>len(set(all_)):
        return False
    else:
        return True    
    
def check_restriction_3(deg_codon):
    ## check if deg codong codes all aa with the same probabilities
    df_aas,coded_aas= deg_codon_to_aa(deg_codon)
    unique_probabilities = set(df_aas['probab'][df_aas['probab']>0])
    # print(deg_codon)
    # print(unique_probabilities)
    if len(unique_probabilities)>1:
        return False
    else:
        return True    
         
def check_restriction_2(deg_codon, AAset):
    _, coded_aas= deg_codon_to_aa(deg_codon)    
    ##restriction 2: deg_codon must not encoded aa not in AAset
    ## this also retricts those deg_codon that code a STOP
    for a in coded_aas:
        if a not in AAset:
            return False
    return True

def gcm_multiple(list_numbers):
    if 1 in list_numbers:
        return 1
    elif len(list_numbers)==1:
        return list_numbers[0]
    elif len(list_numbers) == 2:
        return math.gcd(list_numbers[-1],list_numbers[-2]) 
    else:
        new_list = list_numbers[:-2]
        last_gcm = math.gcd(list_numbers[-1],list_numbers[-2])
        new_list.append(last_gcm)
        return gcm_multiple(new_list)
    
def same_probability(list_deg_codons):
    len_list =[]
    for deg_codon in list_deg_codons:
        _,coded_aas= deg_codon_to_aa(deg_codon)
        len_list.append(len(coded_aas))
    print(len_list)
    gcm = gcm_multiple(len_list)
    proportions = []
    for ele in len_list:
        proportions.append(ele / gcm)
    print(list_deg_codons)
    print(proportions)
    return proportions


def generateCodon(AAset):
    all_deg_codons = get_degenerated_codons()
    ### first search all the deg_codons that code aa inside of the AAset 
    ### and do not code any aa out of this set
    i = 0
    pull = []
    for deg_codon in all_deg_codons:
        i +=1
        ##check if they do not code AAs not in the AAset
        if check_restriction_2(deg_codon,AAset):
            _, aas = deg_codon_to_aa(deg_codon)
            # print( deg_codon, aas)
            #check if the code all the AAs with same probability
            if check_restriction_3(deg_codon):
                # print('yes\n')
                pull.append(deg_codon)
            # else:
            #     print('no\n')
    print(len(pull))
    ## at most, the max num of deg_codon needed woul be the nº of AA in AAset
    max_codons = len(AAset)
    
    
    #####solve by brute forze --> works but too much time or OutOfMemoryError
    #########################
    for i in range(2,max_codons+1):
        combinations = list(it.combinations(pull, i))
        for combi in combinations:
            #check if combination meet restriction 1: 1aa encodeb by just 1 codon
            if check_restriction_1(combi):
                encoded_aas_by_combination = []
                for deg_codon in combi:
                    _, aas = deg_codon_to_aa(deg_codon)
                    encoded_aas_by_combination.extend(aas)
                    print(encoded_aas_by_combination)
                    if len(encoded_aas_by_combination) != len(set(encoded_aas_by_combination)):
                        print('RESTR 1 gone wrong!')
                        return
                    else:
                        print('check combi')
                        valid_combi = True
                        # print(AAset)
                        for AA in AAset:
                            # print(AA)
                            if AA not in encoded_aas_by_combination:
                                print('not valid')
                                valid_combi = False
                if valid_combi:
                    ## if found a valid combi we stop bc it is the smallest
                    ## since started checking combination 1 to 1, 2 to 2, etc
                    return combi


if __name__ == '__main__':
    deg_codon = 'ANT'
    # # print(codon_aa('TAC'))
    # print(get_all_real_codons(deg_codon))
    # print(deg_codon_to_aa(deg_codon))
    # df_aas= deg_codon_to_aa(deg_codon)
    
    # a = check_restriction_2(deg_codon,'KJKELK')
    # print(a)
    
    # a= check_restriction_1([deg_codon, 'VMA', 'ARA','TGG'])
    # print(a)
    # restriction_3(['NDT','VMA','ATG','TGG'])
    all_AA = 'SNIRHLGDVCYFKTQPEAMW'
    combi = generateCodon('KRED')