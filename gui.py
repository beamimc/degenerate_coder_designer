#! /usr/bin/env python

"""
@author: beamimc

"""

import pandas as pd
import streamlit as st
from PIL import Image
import plotly.graph_objects as go

from plotly.subplots import make_subplots

from dcg import generateCodon, get_equivalent_codons, get_coded_aas, check_restriction_1

import base64
import textwrap

from global_variables import CODON_DB, aas

#render_svg function code from https://discuss.streamlit.io/t/is-there-a-way-to-show-svg-url-as-an-image-in-streamlit/12461
def render_svg(svg):
    """Renders the given svg string."""
    b64 = base64.b64encode(svg.encode('utf-8')).decode("utf-8")
    html = r'<img src="data:image/svg+xml;base64,%s"/>' % b64
    st.write(html, unsafe_allow_html=True)


st.set_page_config(
     page_title="Degenerate Codon Designer",
     page_icon="👓",
     )

with st.container():
    option = st.selectbox(
     '',
     ('Main page', 'See an example', ))



f = open("logo.svg","r")
lines = f.readlines()
line_string=''.join(lines)

render_svg(line_string)
f.close()

# st.title('Degenerate Codon Designer')

options ={'All (SNIRHLGDVCYFKTQPEAMW)':'SNIRHLGDVCYFKTQPEAMW',
          'Charged (KREDH)':'KREDH',
          'KRED':'KRED',
          'Non polar (GAVCPLIMWF) ':'GAVCPLIMWF',
          'Polar (STYNQ)':'STYNQ',
          'Custom':''
            }
st.subheader('Select AA set')
option = st.radio(
    '',
     list(options.keys()))

ok = True
if option == 'Custom':
    AAset = st.text_input('Input custom AA set').upper()
    AAset_ = set(AAset)
    for aa in AAset_:
        if aa not in aas:
            ok = False
else:
    AAset = options[option]

st.write('You selected: ', AAset)
AAset = set(AAset)

if ok:
    run = st.button('Design')
    if run:
        with st.spinner('Designing degenerated codons...'):
            combi_prop = generateCodon(AAset)
        sorted_combi_prop=sorted(combi_prop.items(), key=lambda kv: kv[1],reverse=True)
        combi_str = ",".join(str(elem[0]) for elem in sorted_combi_prop)
        prop_str = ":".join(str(int(elem[1])) for elem in sorted_combi_prop)
        
        st.success(f'{combi_str}  ratio {prop_str}\n')
    
        cols = st.columns(7)
        for j in range(7):
            cols[j].metric('','','')
    
        for i in range(len(sorted_combi_prop)):
            j=i
            cols[j].metric("Codon",
                           str(sorted_combi_prop[i][0]), 
                           '')
            cols[j].metric("Ratio",
                         str(int(sorted_combi_prop[i][1])), 
                         '')
    
        ## provide equivalent codons and aas coded by each codon 
        df_all = pd.read_csv(CODON_DB)
        aux = {}  
        mark_str = ''
        for ele in sorted_combi_prop:
            codon = ele[0]
            equivalents = get_equivalent_codons(df_all, codon)
            equivalents.remove(codon)
            s =", ".join(str(elem) for elem in equivalents)
            coded_aas = get_coded_aas(df_all, codon)
            ss ="".join(str(elem) for elem in coded_aas)
            mark_str+=(f'Codon <b>{codon}</b> codes for: {ss}<br />  ')
            if len(equivalents)>0:
                mark_str+=(f'<ul><li>You could equivalently use: {s} </li></ul><br />  ')
            else:
                mark_str+=('<ul><li>There are no equivalent codons </li></ul><br />  ')
     
            
            aux[codon] = len(coded_aas)
        
        with st.expander("See more information"):
            st.markdown(mark_str,True)
        
        
        
        _, total_coded_aas = check_restriction_1(df_all,list(aux.keys()))    
        specs = [[{'type':'domain'}, {'type':'domain'}, {'type':'domain'}]]
        in_colors = ['rgb(248,196,103)', 'rgb(252,237,209)']
        
        out_colors = ['rgb(238,86,86)', 'rgb(248,187,187)']
        fig = make_subplots(rows=1, cols=3, 
                        specs=specs,
                        subplot_titles=['Num. of AAs coded <br /> by each codon', 
                                        'Percent of <br /> AAset coded', 
                                        'Percent of AAs <br /> coded outside AAset'])
        fig.add_trace(go.Pie(
            labels = list(aux.keys()),
            values = list(aux.values()),
            hoverinfo = "label+percent",
            textinfo = "label+value"
            ),1,1)
    
        fig.add_trace(
            go.Pie(labels=['Coded AAs','Non coded AAs'], 
                   values=[len(total_coded_aas), len(AAset)-len(total_coded_aas)],
                hoverinfo = "label+value",
                textinfo = "label+percent",marker_colors=in_colors
        ),1,2)
    
        fig.add_trace(go.Pie(labels=['Inside AAset\n 100%','Outside AAset\n 0%'], 
                              values=[0,1],
                hoverinfo = "label",
                textinfo = "label",marker_colors=out_colors),1,3)
        fig.update( layout_showlegend=False)
        st.plotly_chart(fig)
else:
    st.error('Invalid input')

    