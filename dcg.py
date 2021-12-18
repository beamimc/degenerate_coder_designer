#! /usr/bin/env python

"""
@author: beamimc

"""
import pandas as pd
import math 
import time
import itertools as it
import random
import argparse
# from metaheuristic import SA
__all__ = [ "generateCodon" ]

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
    
def check_restriction_1(df_all,list_deg_codons):
    # ## codons cannot repeat AAs -- 1aa - 1 codon -- 
    #get only codons that are being cheked right now
    # df = df_all[df_all['deg_codon'].isin(list_deg_codons)]
    total_coded_aas = []
    for deg_codon in list_deg_codons:
        coded_aas=get_coded_aas(df_all, deg_codon)
        total_coded_aas.extend(coded_aas)
    if len(total_coded_aas)>len(set(total_coded_aas)):
        ##if any aa is coded more than one time -> more than 1 codong are 
        ## coding the same aa 
        return False,set(total_coded_aas)
    else:
        return True,set(total_coded_aas)

def restriction_3(df_all):
    ##restrict to only viable codons
    ##restriction 3: codon must encode all aas with same prob
    for aa in aas:
        df_all = df_all[df_all[aa]<=1]
    df_all = df_all.reset_index(drop=True)
    return df_all

def get_coded_aas(df_all, deg_codon):
    df = df_all[aas][df_all['deg_codon'] == deg_codon]
    df = df.reset_index(drop=True)
    coded_aas = []
    for aa in aas:
        if df.loc[0,aa]>0:
            coded_aas.append(aa)
    return coded_aas

def restriction_2(df_all,AAset):
    # coded_aas = get_coded_aas(df_all, deg_codon)
    ##restriction 2: deg_codon must not encoded aa not in AAset
    ## this also retricts those deg_codon that code a STOP
    to_drop = []
    for aa in aas:
        index_drop=[]
        if aa not in AAset:
            index_drop = df_all.index[df_all[aa] > 0].tolist()
        to_drop.extend(index_drop)
    df = df_all.drop(to_drop, axis =0)
    df = df.reset_index(drop=True)
    return df

def same_probability(df_all,list_deg_codons):
    len_list =[]
    for deg_codon in list_deg_codons:
        coded_aas= get_coded_aas(df_all, deg_codon)
        len_list.append(len(coded_aas))
    gcm = gcm_multiple(len_list)
    proportions = []
    cod_prop={}
    for i in range(len(len_list)):
        cod_prop[list_deg_codons[i]] = len_list[i] / gcm
        # proportions.append(ele/gcm)
    return cod_prop

def get_equivalent_codons(df_all, deg_codon):
    # newdf1 = df_all.groupby(aas)
    for _, df in  df_all.groupby(aas):
        codons = df['deg_codon'].tolist()
        if deg_codon in codons:
            return codons

def generateCodon(AAset):
    #get all posible deg_codons
    df_all = pd.read_csv('./datasets/deg_codons_DB_clean.csv')
    
    #filter by restriction 3: all aas coded by same probl
    df_all = restriction_3(df_all)
    
    #filter by restrition 2: get deg_codons that only code aa in the AAset 
    df_all = restriction_2(df_all,AAset)

    #get pull of deg_codons that can be used to create a combination
    pull = df_all['deg_codon'].tolist()
    print(len(pull))

    ###check if 1 deg_codon codes all AAset
    #####################################
    for combi in pull:
        #check if combination meet restriction 1: 1aa encodeb by just 1 codon
        rest_1, total_coded_aas = check_restriction_1(df_all,[combi])
        if rest_1:
            ##if valid combi, check if all aas in AAset are codeb by combi
            valid_combi = True
            for AA in AAset:
                if AA not in total_coded_aas:
                    # print('not valid')
                    valid_combi = False
            if valid_combi:
                # if found a valid combi we stop bc it is the smallest
                # since started checking combination 1 to 1, 2 to 2, etc
                combi_prop = same_probability(df_all,[combi])
                return combi_prop
                # return combi    
    if len(pull) <= 50 : 
        print('brute force')
        #####solve by brute force --> works but too much time or OutOfMemoryError
        ##### only do if nº of posible codons is <50 and didnt find only 1 codon
        ##### best solution guaranteed
        ############################################
        #at most, the max num of deg_codon needed would be the nº of AA in AAset
        max_codons = len(AAset)
        for i in range(2,max_codons+1):
            combinations = list(it.combinations(pull, i))
            for combi in combinations:
                #check if combination meet restriction 1: 1aa encodeb by just 1 codon
                rest_1, total_coded_aas = check_restriction_1(df_all,combi)
                if rest_1:
                    ##if valid combi, check if all aas in AAset are codeb by combi
                    valid_combi = True
                    for AA in AAset:
                        if AA not in total_coded_aas:
                            # print('not valid')
                            valid_combi = False
                    if valid_combi:
                        # if found a valid combi we stop bc it is the smallest
                        # since started checking combination 1 to 1, 2 to 2, etc
                        combi_prop = same_probability(df_all,combi)
                        return combi_prop
                        # return combi
    else:
        ######## greedy search--> fast
        ######## best solution not guaranteed
        #####################################
        AAset = set(AAset)
        print('greedy')
        stop = False
        best = pull
        best_n = len(pull)
        i = 0
        not_better = 0
        while not stop:
            ##inicialize on random codon, and then greddy search 
            combi = []
            AAset_ = set(AAset)
            initial_codon =random.choice(pull)
            combi.append(initial_codon)
            # print(initial_codon)
            already_coded = set(get_coded_aas(df_all, initial_codon))
            # print(already_coded)
            #new AAset with aas not yet coded
            AAset_ -=already_coded
            ##make meet restriction 2
            df = restriction_2(df_all,AAset_)
            while len(AAset_)>0:
                # first , get deg_codon that codes the most aas in AAset
                id_max = df[aas].sum(axis=1).idxmax() 
                max_codon = df.iloc[id_max,0]
                already_coded = set(get_coded_aas(df, max_codon))
                combi.extend([max_codon])
                # print(max_codon,already_coded)
                #new AAset with aas not yet coded
                AAset_ -=already_coded
                ##make meet restriction 2
                df = restriction_2(df,AAset_)
            # print(combi)
            if len(combi)<len(best):
                best = combi.copy()
                best_n = len(combi)
            else:
                not_better +=1
            i+=1
            if not_better >150:
                stop = True 
        combi = best.copy()
        ## get proportions of the best combination found
        combi_prop = same_probability(df_all,combi)
        return combi_prop

def main():
    all_AA = 'SNIRHLGDVCYFKTQPEAMW'
    charged = 'KRED'
    charged_ = 'KREDH'
    non_polar = 'GAVCPLIMWF'
    polar = 'STYNQ'
    # #get arguments from cmd    
    ok = True                           
    AAset = input('Introducde AAset: ').upper()
    if AAset == 'ALL':
        AAset =set(all_AA)
    elif AAset =='CHARGED':
        AAset = set(charged)
    elif AAset == 'POLAR':
        AAset = set(polar)
    elif AAset == 'NON POLAR':
        AAset = set(non_polar)
    else:
        ## check if input is correct
        AAset_ = set(AAset)
        for aa in AAset_:
            if aa not in aas:
                return print('wrong AAset')
    
    if ok:
        ###execute 
        startTime = time.time()
        combi_prop = generateCodon(AAset)
        print(combi_prop)
        sorted_combi_prop=sorted(combi_prop.items(), key=lambda kv: kv[1],reverse=True)
        combi_str = ":".join(str(elem[0]) for elem in sorted_combi_prop)
        prop_str = ":".join(str(int(elem[1])) for elem in sorted_combi_prop)

        print(f'Output: {combi_str}  ratio {prop_str}\n')
        executionTime = (time.time() - startTime)
        # print(executionTime,'\n')
        
        ## provide equivalent codons and aas coded by each codon 
        df_all = pd.read_csv('./datasets/deg_codons_DB.csv')
        eq_str = ''
        cd_str = ''
        for ele in sorted_combi_prop:
            codon = ele[0]
            equivalents = get_equivalent_codons(df_all, codon)
            equivalents.remove(codon)
            s =", ".join(str(elem) for elem in equivalents)
            if len(equivalents)>0:
                # print(f'Instead of {codon} you could equivalently use: {s}\n')
                eq_str+=(f'Instead of {codon} you could equivalently use: {s}\n')
            coded_aas = get_coded_aas(df_all, codon)
            s ="".join(str(elem) for elem in coded_aas)
            cd_str+=(f'Codon {codon} codes for: {s}\n')
        print(eq_str,'\n',cd_str)
    
    else:
        print('wrong input',AAset)
        return

if __name__ == '__main__':
    main()