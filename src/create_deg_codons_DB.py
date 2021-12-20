#! /usr/bin/env python
"""
@author: beamimc

"""
import pandas as pd
import os 

from global_variables import *
from dcg import gcm_multiple

def get_all_real_codons(deg_codon):
    '''
    given a degenerated codon, get all (non degenerated) codons they work as
    e.g given AAY, calculte that it works as AAT and AAC
    if given a non deg codon, it would also work, returning that one
    e.g given AAT, it would retutn AAT
    
    Parameters
    ----------
    deg_codon: string
               codon that we want to know non degenerated codons they work as
            
    Returns    
    -------
    posibilities: list
                  list of possible codons            
    '''
    posibilities = []
    for n1 in deg_nucl[deg_codon[0]]:
        for n2 in deg_nucl[deg_codon[1]]:
            for n3 in deg_nucl[deg_codon[2]]:
                posibilities.append(n1+n2+n3)
    return posibilities
def deg_codon_to_aa(deg_codon):
    '''
    given a codon, get the numer of times it codes each aminoacid 
    
    Parameters
    ----------
    deg_codon : string 
                codon that we want to know the aminoacids it codes 

    Returns
    -------
    dc : dict
        dict with aas as keys and the numer of times it is coded as value.

    '''
    posibilities = get_all_real_codons(deg_codon)
    df_aas = pd.DataFrame(columns = ['count'], index =aas).fillna(0)  
    
    for codon in posibilities:    
        aa = codon_aa[codon]
        df_aas.loc[aa,'count'] +=1

    ##normalize proportions 
    # e.g AAY codes N 2 times and no other AA, so it would be equivalent to put
    # it codes it 1 times, since the prob of coding N is 1
    # e.g AAN codes K 2 times and N 2 times, so it would be equivalent to 
    # put K 1 time and N 1 time, since the prob of coding each is 0.5 
    prop = (df_aas[df_aas['count']>0]['count'].tolist())
    gcm = gcm_multiple(prop)
    df_aas.loc[:,'count'] = df_aas.loc[:,'count'] / gcm 
    
    #create dictionaty with result, so it is easier to append to a df
    dc = df_aas['count'].to_dict()
    dc['deg_codon'] = deg_codon
    return dc
def get_degenerated_codons():
    ''' 
    create all posible combinations (codons) of 3 degenerated nucleotides
    
    Returns    
    -------
    all_deg_codons: list
                    list with all posible degenerated codons
    '''
    #there are 15 degenerate nucleotides
    #15x15x15 different posible deg_codons
    all_deg_codons=[]
    nucls = list(deg_nucl.keys())
    for n1 in nucls:
         for n2 in nucls:
             for n3 in nucls:
                 all_deg_codons.append(n1+n2+n3)
    return all_deg_codons

def create_db_codon():
    '''
    create 
    each value represents how many times that aminoacid is coded by that codon
    e.g df.loc[6,'N'] = 2 would mean that codon 0 (i.e AAH) codes N two times
    
    Returns
    -------
    df : pandas DataFrame
        dataframe with all posible degenerated codons
    df_clean : pandas DataFrame
        dataframe after removing 'duplicates' of equivalent codons .

    '''
    #get all posible combinations into df
    all_deg_codons = get_degenerated_codons()
    columns = ['deg_codon']
    columns.extend(aas)
    df = pd.DataFrame(columns=columns)
    df_clean = pd.DataFrame(columns=columns)
    
    #generate df of codon and aa they code proportionally
    for deg_codon in all_deg_codons:
        dc_aa = deg_codon_to_aa(deg_codon)
        df = df.append(dc_aa,ignore_index=True)
        
    #remove codons that do the same i.e equivalent codons
    # e.g. AAT, AAC and AAY all code N with probability 1
    # so only AAT would be kept in the df_clean
    df_clean = df.copy()
    df_clean = df_clean.drop_duplicates(
          subset = aas,
          keep = 'first').reset_index(drop = True)
    return df, df_clean

def main():
    #get current path
    current_path = os.getcwd()
    directory = 'datasets'
    
    #create directory if it does not exist yet
    if not os.path.exists(directory):
        os.mkdir(directory)
        
    #new file names
    filename = 'deg_codons_DB.csv'
    filename_clean = 'deg_codons_DB_clean.csv'
    
    #execute 
    print('Creating DB and cleaned DB...') 
    df, df_clean = create_db_codon()
    
    #save results to csv
    df.to_csv(current_path+'\\'+directory+'\\'+filename,index = False)
    print(f'File saved at {current_path}\\{directory}\\{filename}')
    df_clean.to_csv(current_path+'\\'+directory+'\\'+filename_clean,index = False)
    print(f'File saved at {current_path}\\{directory}\\{filename_clean}')
    return
    
if __name__ == '__main__':
    main()
