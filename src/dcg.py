#! /usr/bin/env python

"""
@author: beamimc

"""
import pandas as pd
import math 
import time
import itertools as it
import random

from .global_variables import *

def gcd_multiple(list_numbers):
    '''
    get the greatest common divisor of a list of numbers

    Parameters
    ----------
    list_numbers : list
            list of numbers

    Returns
    -------
    int
        greatest common divisor of the list of numbers.

    '''
    #if there is a 1 in the list, gcm is 1
    if 1 in list_numbers:
        return 1
    #if only one number in list, gcm is that number
    elif len(list_numbers)==1:
        return list_numbers[0]
    #if two numbers in list, calculate gcm with math.gcd
    elif len(list_numbers) == 2:
        return math.gcd(list_numbers[-1],list_numbers[-2]) 
    #if there are more than 2 numbers
    else:
        #remove last two numers
        new_list = list_numbers[:-2]
        #calculate gcd of those two numbers
        last_gcd = math.gcd(list_numbers[-1],list_numbers[-2])
        #append the gcd
        new_list.append(last_gcd)
        #recursively will call function until 1 in list or only 2 numbers left
        return gcd_multiple(new_list)
    
def check_restriction_1(df_all,list_deg_codons):
    '''
    check if codons from a list of codons repeat aminoacids and calculate
    the aminoacids coded by the whole set of codons
    restriction 1 : codons cannot repeat AAs -> each aa coded by just 1 codon 

    Parameters
    ----------
    df_all : pandas DataFrame
        df with codons as rows and the number of times they code each aa as columns
    list_deg_codons : list
        list of codon we want to chek if they code the same aas or not

    Returns
    -------
    bool
        True if codons do not repeat aminoacids. False if they do
    set
        Set of the total aminoacids coded by all the codons in list_deg_codons

    '''
    total_coded_aas = []
    for deg_codon in list_deg_codons:
        #get aas coded by each codon 
        coded_aas=get_coded_aas(df_all, deg_codon)
        total_coded_aas.extend(coded_aas)
    if len(total_coded_aas)>len(set(total_coded_aas)):
        #if any aa is coded more than one time -> more than 1 codong are 
        #coding the same aa 
        return False,set(total_coded_aas)
    else:
        return True,set(total_coded_aas)

def restriction_3(df_all):
    '''
    create df with only viable codons, i.e remove codons that code more than
    one aminoacid with different probabilities
    restriction 3: codons must encode all aas with same probability
    
    Parameters
    ----------
    df_all : pandas DataFrame
        original df of codons (rows) and their aminoacids (columns)

    Returns
    -------
    df_all : pandas DataFrame
        filtered df with only viable codons

    '''
    for aa in aas:
        df_all = df_all[df_all[aa]<=1]
    df_all = df_all.reset_index(drop=True)
    return df_all

def get_coded_aas(df_all, deg_codon):
    '''
    calculate aminoacids coded by a given degenerated codon

    Parameters
    ----------
    df_all : pandas DataFrame
        df of codons (rows) and their aminoacids (columns)
    deg_codon : string
        degenerated codon

    Returns
    -------
    coded_aas : list
        list of aminoacids coded by deg_codon

    '''
    df = df_all[aas][df_all['deg_codon'] == deg_codon]
    df = df.reset_index(drop=True)
    coded_aas = []
    for aa in aas:
        if df.loc[0,aa]>0:
            coded_aas.append(aa)
    return coded_aas

def restriction_2(df_all,AAset):
    '''
    create df with only codons that do not code aas outside AAset
    restriction 2: deg_codon must not encoded aa not in AAset

    Parameters
    ----------
    df_all : pandas DataFrame
        df of codons (rows) and their aminoacids (columns)
    AAset : set
        set of aminoacids that need to be coded

    Returns
    -------
    df : pandas DataFrame
        filtered df with only codons that do not code aminoacids outside AAset

    '''
    to_drop = []
    #check all aminoacids
    for aa in aas:
        index_drop=[]
        #if aminoacid not in AAset, remove codons that code them
        if aa not in AAset:
            #get index of rows (codons) that code that aa more than 0 times 
            index_drop = df_all.index[df_all[aa] > 0].tolist()
        to_drop.extend(index_drop)
    df = df_all.drop(to_drop, axis =0)
    df = df.reset_index(drop=True)
    return df

def same_probability(df_all,list_deg_codons):
    '''
    given a set of degenerated codons, calculate ratio, so that all aminoacids
    coded by all of them are coded with same probabilities

    Parameters
    ----------
    df_all : pandas DataFrame
        df of codons (rows) and their aminoacids (columns)
    list_deg_codons : list
        list of degenerated codons

    Returns
    -------
    cod_prop : dictionary
        dictionary with codons as keys and ratios as values.

    '''
    len_list =[]
    for deg_codon in list_deg_codons:
        coded_aas= get_coded_aas(df_all, deg_codon)
        len_list.append(len(coded_aas))
    gcm = gcd_multiple(len_list)
    cod_prop={}
    for i in range(len(len_list)):
        cod_prop[list_deg_codons[i]] = len_list[i] / gcm
    return cod_prop

def get_num_deg_nucl(list_deg_codons):
    '''
    calculate number of degenerated nucleotides used in a set of codons

    Parameters
    ----------
    list_deg_codons : list
        list of degenerated codons

    Returns
    -------
    int
        number of degenerated nucleotides in the list of degenerated codons

    '''
    count=0
    for deg_codon in list_deg_codons:
        for nucl in deg_codon:
            #if nucleotide is not A T C or G, it is a degenerated nucleotide
            if nucl not in 'ATCG':
                count+=1
    return(count)

def get_equivalent_codons(df_all, deg_codon):
    '''
    given a degenerated codon, find all equivalent codons i.e those that code
    the same aminoacids with same probabilities, and also have the same number
    of degenerated nucleotides

    Parameters
    ----------
    df_all : pandas DataFrame
        df of codons (rows) and their aminoacids (columns)
    deg_codon : string
        degenerated codon

    Returns
    -------
    equivalents : list
        list of equivalent codons

    '''
    n_deg_nucl = get_num_deg_nucl([deg_codon])
    equivalents = []
    for _, df in  df_all.groupby(aas):
        codons = df['deg_codon'].tolist()
        if deg_codon in codons:
            for codon in codons:
                if get_num_deg_nucl([codon]) == n_deg_nucl:
                    equivalents.append(codon)
            return equivalents

def generateCodon(AAset):
    '''
    calculate the minimun set of degenerated codons to code all the aminoacids
    in a given AAset 

    Parameters
    ----------
    AAset : set
        set of aminoacids that to be coded

    Returns
    -------
    combi_prop : dictionary
        dictionary with codons as keys, and their proportions/ratio needed so as
        all aminoacids are coded in the same probability as values.

    '''
    #get all posible deg_codons without 'duplicated codons' i.e codons that 
    #code the same aminoacids with same probabilities
    df_all = pd.read_csv(CLEAN_CODON_DB)
    
    #filter by restriction 3: all aas coded by same probl
    df_all = restriction_3(df_all)
    
    #filter by restrition 2: get deg_codons that only code aa in the AAset 
    df_all = restriction_2(df_all,AAset)

    #get pull of deg_codons that can be used to create a combination
    pull = df_all['deg_codon'].tolist()

    ###check if 1 deg_codon codes all AAset
    #######################################
    best_n_deg_nucl = 15
    best = pull[0]
    found_valid = False
    for combi in pull:
        #check if combi meets restriction 1: each aa coded by just 1 codon
        rest_1, total_coded_aas = check_restriction_1(df_all,[combi])
        if rest_1:
            ##if rest_1 is met, check if all aas in AAset are codeb by combi
            valid_combi = True
            for AA in AAset:
                if AA not in total_coded_aas:
                    valid_combi = False
                    break
            if valid_combi:
                found_valid = True
                n_deg_nucl = get_num_deg_nucl([combi])
                # print(n_deg_nucl)
                if n_deg_nucl == 0:
                    #if found a valid combi without deg nucleotides stop bc it 
                    #is already the best: the smallest (1codon) with 0 deg nucl
                    #(would only happen when trying to code 1 aminoacid)
                    combi_prop = same_probability(df_all,[combi])
                    return combi_prop
                elif n_deg_nucl < best_n_deg_nucl:
                    best = combi
                    best_n_deg_nucl = n_deg_nucl
    #after cheking all combinations of 1:
    #if found any validwe stop bc it is the smallest (1codon)
    if found_valid:
        #return the best one, i.e the one with less deg nucleotides
        # and 0 deg nucl
        combi_prop = same_probability(df_all,[best])
        return combi_prop
    
    ###brute force in case combinations of 1 codon were not found
    #############################################################
    #solve by brute force -> works but too much time or OutOfMemoryError
    #only do if nº of posible combinations is < 2000 bc checking all this
    #combinations takes a reasonable time 
    #directly skip if num of codons in pull is > 65 bc combinations of 2 
    #out of 65 items is already >2000
    #best solution guaranteed
    greedy = False
    if len(pull) < 65 : 
        print('trying brute force...')
        #at most, the max num of deg_codon needed would be num of AA in AAset
        max_codons = len(AAset)
        #start looking for combis of 2 bc combinations of 1 already checked
        for i in range(2,max_codons+1):
            combinations = list(it.combinations(pull, i))
            if len(combinations) < 2000:
                best_n_deg_nucl = 15
                best = pull[0]
                found_valid = False
                for combi in combinations:
                    #check if combi meets restr 1: each aa coded by just 1 codon
                    rest_1, total_coded_aas = check_restriction_1(df_all,combi)
                    if rest_1:
                        ##if rest_1 is met, check if all aas in AAset are coded
                        valid_combi = True
                        for AA in AAset:
                            if AA not in total_coded_aas:
                                # print('not valid')
                                valid_combi = False
                                break
                        if valid_combi:
                            found_valid = True
                            n_deg_nucl = get_num_deg_nucl(combi)
                            if n_deg_nucl < best_n_deg_nucl:
                                best = combi
                                best_n_deg_nucl = n_deg_nucl
                if found_valid:
                    #if found a valid combi, stop bc it is the smallest
                    #since started checking combinations of 1, then 2,etc
                    #return the best one found (min num deg nucl)
                    combi_prop = same_probability(df_all,best)
                    return combi_prop
            else:
                greedy = True
                break
    else:
        greedy = True
        
    #variation of greedy search until stop condition is met -> fast
    ###############################################################  
    #in cases combinations of 1 did not work and searching by brute force
    #takes too much computaional power
    #best solution not guaranteed
    if greedy:
        print('greedy search...')
        stop = False
        best = pull
        best_n = len(pull)
        best_n_deg_nucl=15
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
            #if new combi has less codons than the best so far, keep it
            if len(combi)<len(best):
                n_deg_nucl = get_num_deg_nucl(combi)
                best = combi.copy()
                best_n = len(combi)
                best_n_deg_nucl = n_deg_nucl
                #restart stop condition counter
                not_better = 0
            #if new combi has same number of codons than the best so far
            #keep it only if it has less deg nucleotides
            elif len(combi)==len(best):
                n_deg_nucl = get_num_deg_nucl(combi)
                if n_deg_nucl < best_n_deg_nucl:
                    best = combi.copy()
                    best_n = len(combi)
                    best_n_deg_nucl = n_deg_nucl
                    #restart stop condition counter
                    not_better = 0
                else:
                    #increment stop condition counter
                    not_better +=1
            else:
                #increment stop condition counter
                not_better +=1
            #stop condition: stop if not found a better solution after 300 iters
            if not_better > 300:
                stop = True 
        # print(best_n, best_n_deg_nucl)
        
        combi = best.copy()
        #get proportions (ratio) of the best combination found
        combi_prop = same_probability(df_all,combi)
        return combi_prop

def run_generateCodon(AAset):
    '''
    runs program to calculate minimun degenerated codons given an AAset,
    calculates which aminoacids are coded by each codon, and calculates 
    equivalent codons of each.

    Parameters
    ----------
    AAset : set
        set of aminoacids that to be coded

    Returns
    -------
    sorted_combi_prop : dictionary
        dictionary with codons as keys, and their proportions/ratio as values
        sorted by values in descendent way
    eq_str : string
        string to print user the equivalent codons of each codon.
    cd_str : string
        string to print user the aminoacids coded by each codon.

    '''
    #execution
    combi_prop = generateCodon(AAset)
    
    #sort result to provide ratios in descendent order
    sorted_combi_prop=sorted(combi_prop.items(), 
                             key=lambda kv: kv[1],
                             reverse=True)
    
    #provide equivalent codons and aas coded by each codon 
    df_all = pd.read_csv(CODON_DB)
    eq_str = ''
    cd_str = ''
    for ele in sorted_combi_prop:
        codon = ele[0]
        equivalents = get_equivalent_codons(df_all, codon)
        equivalents.remove(codon)
        #create string to show user
        s =", ".join(str(elem) for elem in equivalents)
        if len(equivalents)>0:
            eq_str+=(f'Instead of {codon} you could equivalently use: {s}\n')
        coded_aas = get_coded_aas(df_all, codon)
        #create string to show user
        s ="".join(str(elem) for elem in coded_aas)
        cd_str+=(f'Codon {codon} codes for: {s}\n')

    return sorted_combi_prop, eq_str, cd_str
