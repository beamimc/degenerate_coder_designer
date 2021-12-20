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

from global_variables import *

        
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
    
def get_num_deg_nucl(list_deg_codons):
    count=0
    for deg_codon in list_deg_codons:
        for nucl in deg_codon:
            if nucl not in 'ATCG':
                count+=1
    return(count)

def generateCodon(AAset):
    #get all posible deg_codons without 'duplicated codons' i.e codons that 
    #code the same aminoacids with same probabilities
    df_all = pd.read_csv(CLEAN_CODON_DB)
    
    #filter by restriction 3: all aas coded by same probl
    df_all = restriction_3(df_all)
    
    #filter by restrition 2: get deg_codons that only code aa in the AAset 
    df_all = restriction_2(df_all,AAset)

    #get pull of deg_codons that can be used to create a combination
    pull = df_all['deg_codon'].tolist()
    print(len(pull))

    ###check if 1 deg_codon codes all AAset
    #######################################
    for combi in pull:
        #check if combi meets restriction 1: each aa coded by just 1 codon
        rest_1, total_coded_aas = check_restriction_1(df_all,[combi])
        if rest_1:
            ##if rest_1 is met, check if all aas in AAset are codeb by combi
            valid_combi = True
            for AA in AAset:
                if AA not in total_coded_aas:
                    # print('not valid')
                    valid_combi = False
            if valid_combi:
                # if found a valid combi we stop bc it is the smallest (1codon)
                combi_prop = same_probability(df_all,[combi])
                return combi_prop
    greedy = False
    if len(pull) < 65 : 
        print('brute force')
        #solve by brute force -> works but too much time or OutOfMemoryError
        #only do if nº of posible combinations is < 2000        
        #directly skip if num of codons in pull is > 65 bc combinations of 2 
        #out of 65 items is already >2000
        #best solution guaranteed
        ############################################
        #at most, the max num of deg_codon needed would be the nº of AA in AAset
        max_codons = len(AAset)
        for i in range(2,max_codons+1):
            print('combinations of, ',i)
            combinations = list(it.combinations(pull, i))
            print(len(combinations))
            if len(combinations) < 2000:
                for combi in combinations:
                    #check if combi meets restrict 1: each aa coded by just 1 codon
                    rest_1, total_coded_aas = check_restriction_1(df_all,combi)
                    if rest_1:
                        ##if rest_1 is met, check if all aas in AAset are coded
                        valid_combi = True
                        for AA in AAset:
                            if AA not in total_coded_aas:
                                # print('not valid')
                                valid_combi = False
                        if valid_combi:
                            #if found a valid combi, stop bc it is the smallest
                            #since started checking combinations of 1, then 2,etc
                            combi_prop = same_probability(df_all,combi)
                            return combi_prop
            else:
                greedy = True
                break
    else:
        greedy = True
    if greedy:
        #variation of greedy search until stop condition is met -> fast
        #best solution not guaranteed
        #####################################
        print('greedy')
        stop = False
        best = pull
        best_n = len(pull)
        i = 0
        not_better = 0
        while not stop:
            ##inicialize on random codon, then greddy search to get full combi
            combi = []
            AAset_ = AAset.copy()
            initial_codon =random.choice(pull)
            combi.append(initial_codon)
            
            #get aas coded by the initial codon
            already_coded = set(get_coded_aas(df_all, initial_codon))

            #new AAset with aas not yet coded
            AAset_ -=already_coded
            #meet restriction 2: 
            #only use deg_codons that do not code aas outside AAset 
            df = restriction_2(df_all,AAset_)
            
            #greddy search to finish the combination
            #stop when there are no more aas in AAset_: all aas coded by combi
            while len(AAset_)>0:
                #first, get deg_codon that codes the most aas in AAset
                id_max = df[aas].sum(axis=1).idxmax() 
                max_codon = df.iloc[id_max,0]
                already_coded = set(get_coded_aas(df, max_codon))
                #save codon into combi
                combi.extend([max_codon])
                #new AAset with aas not yet coded                
                already_coded = set(get_coded_aas(df, max_codon))
                AAset_ -=already_coded
                #meet restriction 2: 
                #only use deg_codons that do not code aas outside AAset 
                df = restriction_2(df,AAset_)
            
            
            if len(combi)<len(best):
                best = combi.copy()
                best_n = len(combi)
            else:
                not_better +=1
            i+=1
            if not_better > 200:
                stop = True 
        combi = best.copy()
        ## get proportions of the best combination found
        combi_prop = same_probability(df_all,combi)
        return combi_prop

def main():
    
    good_example = 'GAVCPLIMWF' #### best by brute force is 4 codons {'ATK': 1.0, 'TKT': 1.0, 'KGG': 1.0, 'SYA': 2.0}
    example = 'GAVCPLIMWFKRED' ## best is 4 codons
    exampl = 'STYNQKRED'
    
    all_AA = 'SNIRHLGDVCYFKTQPEAMW'
    charged = 'KREDH'
    polar = 'STYNQ'  
    non_polar = 'GAVCPLIMWF'
    
    ok = True       
    run = True                   
    
    #welcome message to user
    print('\nDEGENERATE CODON DESIGNER')
    print('Program to desgin the minimun combination of degenerate codons that \
code all the aminoacids given (AAset) in equal proportions\n')
    print('Please, provide de AAset you want to code. Aminoacids have to be \
introduced by their one letter representation')
    print('E.g: input KRED for coding Lys, Arg, Glu and Asp\n')
    print('Shortcut inputs:\n \
    -- all: SNIRHLGDVCYFKTQPEAMW \n \
    -- charged: KREDH \n \
    -- polar: STYNQ \n \
    -- non polar: GAVCPLIMWF')
    print('E.g: input "all" would be equal to input SNIRHLGDVCYFKTQPEAMW\n')
    #get input
    AAset = input('Introducde AAset: ').upper()
    # check input
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
        AAset = set(AAset)
        for aa in AAset:
            if aa not in aas:
                return print('The introduced input was not a valid AAset')
            
    #if input is ok run the program
    # startTime = time.time()
    #execution
    combi_prop = generateCodon(AAset)
    
    #sort result in order to provide ratios in descendent order
    sorted_combi_prop=sorted(combi_prop.items(), key=lambda kv: kv[1],reverse=True)
    combi_str = ":".join(str(elem[0]) for elem in sorted_combi_prop)
    prop_str = ":".join(str(int(elem[1])) for elem in sorted_combi_prop)
    
    #show result
    print(f'Output: {combi_str}  ratio {prop_str}\n')
    # executionTime = (time.time() - startTime)
    
    #provide equivalent codons and aas coded by each codon 
    df_all = pd.read_csv(CODON_DB)
    eq_str = ''
    cd_str = ''
    for ele in sorted_combi_prop:
        codon = ele[0]
        equivalents = get_equivalent_codons(df_all, codon)
        equivalents.remove(codon)
        s =", ".join(str(elem) for elem in equivalents)
        if len(equivalents)>0:
            eq_str+=(f'Instead of {codon} you could equivalently use: {s}\n')
        coded_aas = get_coded_aas(df_all, codon)
        s ="".join(str(elem) for elem in coded_aas)
        cd_str+=(f'Codon {codon} codes for: {s}\n')
    print(eq_str)
    print(cd_str)
    return

if __name__ == '__main__':
    main()