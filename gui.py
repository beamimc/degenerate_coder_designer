#! /usr/bin/env python

"""
@author: beamimc

"""

import pandas as pd
import streamlit as st
import plotly 
from PIL import Image
import plotly.graph_objects as go

from plotly.subplots import make_subplots

from dcg import generateCodon, get_equivalent_codons, get_coded_aas, check_restriction_1
# st.set_page_config(layout="wide")

import base64
import textwrap

def render_svg(svg):
    """Renders the given svg string."""
    b64 = base64.b64encode(svg.encode('utf-8')).decode("utf-8")
    html = r'<img src="data:image/svg+xml;base64,%s"/>' % b64
    st.write(html, unsafe_allow_html=True)

def render_svg_example():
    svg = """
        <svg xmlns="http://www.w3.org/2000/svg" width="100" height="100">
            <circle cx="50" cy="50" r="40" stroke="green" stroke-width="4" fill="yellow" />
        </svg>
    """
    st.write('## Rendering an SVG in Streamlit')

    st.write('### SVG Input')
    st.code(textwrap.dedent(svg), 'svg')

    st.write('### SVG Output')
    render_svg(svg)
    
def header(url):
        st.markdown(f'<p style="background-color:#ededed;color:#000000;font-size:18px;border-radius:2%;border-color:#000000">{url}</p>', unsafe_allow_html=True)
  
st.set_page_config(
     page_title="Degenerate Codon Designer",
     page_icon="🍩",
     # initial_sidebar_state="expanded",
       # layout="wide",
 )

with st.container():
    option = st.selectbox(
     '',
     ('Main page', 'See an example', ))



f = open("logoo.svg","r")
lines = f.readlines()
line_string=''.join(lines)

render_svg(line_string)
f.close()

# st.title('Degenerate codon designer')

options ={'all':'SNIRHLGDVCYFKTQPEAMW',
          'charged':'KRED',
          'charged with H':'KREDH',
          'non polar':'GAVCPLIMWF',
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
    # st.success('Designed!')
    # st.balloons()
    print(combi_prop)
    # print(combi_prop)
    sorted_combi_prop=sorted(combi_prop.items(), key=lambda kv: kv[1],reverse=True)
    combi_str = ":".join(str(elem[0]) for elem in sorted_combi_prop)
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
            # print(f'Instead of {codon} you could equivalently use: {s}\n')
            eq_str+=(f'  Instead of {codon} you could equivalently use: {s} <br />  ')
            mark_str+=(f'<ul><li>You could equivalently use: {s} </li></ul><br />  ')
        else:
            mark_str+=(f'<ul><li>There are no equivalent codons </li></ul><br />  ')
       
        cd_str+=(f'  Codon {codon} codes for: {s}<br />  ')
        
        aux[codon] = len(coded_aas)
        # df = pd.DataFrame(index = list(combi_prop.keys()),columns = ['count'])
        # print(df)
        # df.append({'count':len(coded_aas)}, ignore_index = True)
    
    # def header(url):
    #     st.markdown(f'<p style="background-color:#ededed;color:#000000;font-size:18px;border-radius:2%;border-color:#000000">{url}</p>', unsafe_allow_html=True)

    # col2.markdown(f'<p style="background-color:#0066cc;color:#33ff33;font-size:24px;border-radius:2%;">{eq_str}</p>')
    with st.expander("See more information"):
        st.markdown(mark_str,True)
    ###bar chart
    # df = pd.DataFrame.from_dict(aux,orient ='index')
    # st.bar_chart(df, width=0, height=0, use_container_width=True)
    
    
    
    _, total_coded_aas = check_restriction_1(df_all,list(aux.keys()))
    
    # col1, col2, col3 = st.columns(3)
    
    specs = [[{'type':'domain'}, {'type':'domain'}, {'type':'domain'}]]
    out_colors = ['rgb(248,196,103)', 'rgb(252,237,209)']
    fig = make_subplots(rows=1, cols=3, 
                    specs=specs,
                    subplot_titles=['Nº of AAs coded <br /> by each codon', 
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
            textinfo = "label+percent"
    ),1,2)

    fig.add_trace(go.Pie(labels=['AAs outside AAset','AAs inside AAset'], 
                          values=[0,1],
            hoverinfo = "label+percent",
            textinfo = "label+percent",marker_colors=out_colors),1,3)
    fig.update( layout_showlegend=False)
    st.plotly_chart(fig)

    