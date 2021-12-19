#! /usr/bin/env python

"""
@author: beamimc

"""

import pandas as pd
import streamlit as st
import plotly 

import plotly.graph_objects as go

from plotly.subplots import make_subplots

from dcg import generateCodon, get_equivalent_codons, get_coded_aas, check_restriction_1
# st.set_page_config(layout="wide")
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
    
    st.text(f'{combi_str}  ratio {prop_str}\n')

    cols = st.columns(7)
    for j in range(7):
        cols[j].metric('','','')

    for i in range(len(sorted_combi_prop)):
        j=i
        cols[j].metric("Codon",
                       str(sorted_combi_prop[i][0]), 
                       '')
                       # str(sorted_combi_prop[i][1]))
        cols[j].metric("Ratio",
                     str(sorted_combi_prop[i][1]), 
                     '')
                     # str(sorted_combi_prop[i][1]))

    ## provide equivalent codons and aas coded by each codon 
    df_all = pd.read_csv('./datasets/deg_codons_DB.csv')
    eq_str = ''
    cd_str = ''
    aux = {}
    for ele in sorted_combi_prop:
        codon = ele[0]
        equivalents = get_equivalent_codons(df_all, codon)
        equivalents.remove(codon)
        s =", ".join(str(elem) for elem in equivalents)
        if len(equivalents)>0:
            # print(f'Instead of {codon} you could equivalently use: {s}\n')
            eq_str+=(f'Instead of {codon} you could equivalently use: {s} <br />  ')
        coded_aas = get_coded_aas(df_all, codon)
        s ="".join(str(elem) for elem in coded_aas)
        cd_str+=(f'Codon {codon} codes for: {s}<br />  ')
        
        aux[codon] = len(coded_aas)
        # df = pd.DataFrame(index = list(combi_prop.keys()),columns = ['count'])
        # print(df)
        # df.append({'count':len(coded_aas)}, ignore_index = True)
    st.markdown(eq_str, True)
    st.markdown(cd_str, True)
    
    ###bar chart
    # df = pd.DataFrame.from_dict(aux,orient ='index')
    # st.bar_chart(df, width=0, height=0, use_container_width=True)
    
    
    
    _, total_coded_aas = check_restriction_1(df_all,list(aux.keys()))
    
    col1, col2, col3 = st.columns(3)
    
    fig = go.Figure(
        go.Pie(
        labels = list(aux.keys()),
        values = list(aux.values()),
        hoverinfo = "label+percent",
        textinfo = "value"
        ))
    fig.update_layout(
    autosize=False,
    width=300,
    height=300,)
    col1.plotly_chart(fig)
    
    fig = go.Figure(
        go.Pie(labels=['Coded AAs','Non coded AAs'], 
               values=[len(total_coded_aas), len(AAset)-len(total_coded_aas)],
            hoverinfo = "label+percent",
            textinfo = "value"
    ))
    fig.update( layout_showlegend=False)
    fig.update_layout(
    autosize=False,
    width=300,
    height=300,)
    col2.plotly_chart(fig)


    fig = go.Figure(
        go.Pie(labels=['Num of AAs outside AAset','Num inside'], 
                          values=[0,1],
            hoverinfo = "label+percent",
            textinfo = "value"
    ))
    fig.update(layout_showlegend=False)
    fig.update_layout(
    autosize=False,
    width=300,
    height=300,)
    col3.plotly_chart(fig)

    