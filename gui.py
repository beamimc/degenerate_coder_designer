#! /usr/bin/env python

"""
@author: beamimc

"""

import streamlit as st
import pandas as pd
import numpy as np
import time
from dcg import generateCodon, get_equivalent_codons, get_coded_aas

st.title('Degenerate codon designer')

options ={'all':'SNIRHLGDVCYFKTQPEAMW',
          'charged':'KRED',
          'charged with H':'KREDH',
          'non_polar':'GAVCPLIMWF',
          'polar':'STYNQ',
          'custom':''
            }
st.subheader('Select AA set')
option = st.radio(
    '',
     list(options.keys()))

if option == 'custom':
    AAset = st.text_input('Input custom AA set')
else:
    AAset = options[option]
    
st.write('You selected: ', AAset)

run = st.button('Design')
if run:
    with st.spinner('Designing degenerated codons...'):
        combi_prop = generateCodon(AAset)
    st.success('Designed!')
    # st.balloons()
    print(combi_prop)
    # print(combi_prop)
    sorted_combi_prop=sorted(combi_prop.items(), key=lambda kv: kv[1],reverse=True)
    combi_str = ":".join(str(elem[0]) for elem in sorted_combi_prop)
    prop_str = ":".join(str(int(elem[1])) for elem in sorted_combi_prop)
    
    st.write(f'Output: {combi_str}  ratio {prop_str}\n')

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
    st.markdown(eq_str)
    st.markdown(cd_str)