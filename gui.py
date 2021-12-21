#! /usr/bin/env python

"""
@author: beamimc

"""

import pandas as pd
import streamlit as st
import base64
import plotly.graph_objects as go

from plotly.subplots import make_subplots

from src.dcg import *
from src.global_variables import *

#variation of render_svg() function code from https://discuss.streamlit.io/t/is-there-a-way-to-show-svg-url-as-an-image-in-streamlit/12461
def render_svg():
    """Renders the given svg string."""
    f = open("./static/logo.svg","r")
    lines = f.readlines()
    svg=''.join(lines)
    f.close()

    b64 = base64.b64encode(svg.encode('utf-8')).decode("utf-8")
    html = r'<img src="data:image/svg+xml;base64,%s"/>' % b64
    return html

#set logo (page title) as global so its not calculated each time the program runs
LOGO_HTML = render_svg()

def main():
    ''' To run GUI of Degenerate Codon Designer on streamlit '''
    
    ##strings to show explanation of examples
    example_all='This is known as the ‘MAX’ degenerate codon set. This set could\
    be used when we want to mutate the a certain position of the protein with all\
    amino acids. That position will mutate to be any of the 20 natural AAs, \
    with equal probabilities.'
    example_KRED='This set could be used to mutate in a case we want to mutate\
    a charged residue of a protein (e.g a voltage indicator) to either \
    positively charged (KR) or negatively charged residues (ED) with equal \
    probabilities.'
    example_polar='This set could be used in order to mutate a known polar residue\
    in a protein, to see if any other polar group would enhance its action. E.g.\
    we have a certain protein in which a residue is known to be key in the binding \
    of a ligand. With this set, we would be able to generate a library of mutants \
    with all the polar aminoacids on that position, thus search for another one\
    that could make the protein-ligand interaction stronger.'
    example_polar='This set could be used in order to mutate a known polar residue\
    in a protein, to see if any other polar group would enhance its action. E.g.\
    we have a certain protein in which a residue is known to be key on the binding \
    of a ligand. With this set, we would be able to generate a library of mutants \
    with all the polar aminoacids on that position, thus search for another one\
    that could make the protein-ligand interaction stronger.'
    example_nonpolar='This set could be used in order to mutate a polar residue\
    known to be key in the binding of a ligand. By generating a library \
    of mutants with non polar aminoacids on that position, we could search for \
    the best non polar aminoacid that disrupts this interaction the most.'
    example_aromatic='This set could be used in a case we want to mutate an \
    aromatic residue of a protein to either of aromatic amino acids (FYW). \
    Within this mutants, we could search for new properties or stronger/weaker\
    interactions due to changes in Pi-Stacking.'
    examples = {'SNIRHLGDVCYFKTQPEAMW':example_all,
                'KRED':example_KRED,
                'STCNQ':example_polar,
                'GAVLMIP':example_nonpolar,
                'FYW':example_aromatic}
   
    #config page
    st.set_page_config(
         page_title="Degenerate Codon Designer",
         page_icon="🔀",
         )

    
    #show page title 
    #########
    st.write(LOGO_HTML, unsafe_allow_html=True)
    st.markdown('<b>Desgin the minimun set of degenerate codons that \
    encode all the amino acids given (AAset) in equal probabilities</b>'
                , True)
       
    #show main menu 
    ###############
    col1,col2 = st.columns(2)
    col1.subheader('Select AA set')
    #set options and get input AAset from user
    options ={'All (SNIRHLGDVCYFKTQPEAMW)':'SNIRHLGDVCYFKTQPEAMW',
              'Charged (KRED)':'KRED',
              'Non polar aliphatic (GAVLMIP) ':'GAVLMIP',
              'Polar (STCNQ)':'STCNQ',
              'Aromatic (FYW)':'FWY',
              'Custom':''
                }
    option = col1.radio('Select one or more amino acids\n', 
                  list(options.keys()))
    col2.image('./static/chart.png','Source: https://che.gg/3Hhut5V')
    
       
    #check if AAset is a valid input
    ok = True
    if option == 'Custom':
        AAset = st.text_input('Input custom AA set. Aminoacids have to be \
                      introduced by their one letter codification.\
                      E.g KR for Lys and Arg. X represents STOP codon.').upper()
        AAset_ = set(AAset)
        for aa in AAset_:
            if aa not in aas:
                ok = False
    else:
        AAset = options[option]
    AAset_ = AAset
    AAset = set(AAset)   
    
    #if valid iput, run program 
    if ok:
        info = st.info(f'Amino acids selected: {AAset_}')
        #run button
        ###########
        run = st.button('Design')
        if run:
            
            #show message while executing generateCodon(AAset)
            ############
            with st.spinner('Designing degenerate codons...'):
                combi_prop = generateCodon(AAset)
            info.empty()
            #show input
            ############
            st.info(f'Amino acids selected: {AAset_}')
                
            #sort result to provide ratios in descendent order
            sorted_combi_prop=sorted(combi_prop.items(), 
                                     key=lambda kv: kv[1],
                                     reverse=True)
            combi_str = ", ".join(str(ele[0]) for ele in sorted_combi_prop)
            prop_str = ":".join(str(int(ele[1])) for ele in sorted_combi_prop)
            
            #show result in green box
            ############
            st.success(f'{combi_str}  ratio {prop_str}\n')
            
            #show result cooler 
            ############
            #create 7 columns to get a better format
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
        
            #get equivalent codons and aas coded by each codon 
            df_all = pd.read_csv(CODON_DB)
            aux = {}  
            mark_str = ''
            for ele in sorted_combi_prop:
                codon = ele[0]
                equivalents = get_equivalent_codons(df_all, codon)
                #create strings to show user
                equivalents.remove(codon)
                s =", ".join(str(elem) for elem in equivalents)
                coded_aas = get_coded_aas(df_all, codon)
                ss ="".join(str(elem) for elem in coded_aas)
                mark_str+=(f'Codon <b>{codon}</b> codes for: {ss}<br/>')
                if len(equivalents)>0:
                    mark_str+=(f'<ul><li>You could equivalently use: {s}</li></ul><br/>')
                else:
                    mark_str+=('<ul><li>There are no equivalent codons</li></ul><br/>')
                aux[codon] = len(coded_aas)
                
            list_deg_codons= list(combi_prop.keys())
            n_deg_nucl = get_num_deg_nucl(list_deg_codons)
            #get total coded aas by codon set 
            _, total_coded_aas = check_restriction_1(df_all,list(aux.keys()))  
            #coutn how many aas are coded by codon set that are outside AAset
            not_in= 0
            for aa in total_coded_aas:
                if aa not in AAset:
                    not_in+=1
                    
            #show infomation in an expander container 
            ############
            with st.expander("See more information"):
                col1, col2 = st.columns(2)
                col1.markdown(mark_str,True)
                with col2:
                    st.metric('Number of degenerate codons',
                                f'{int(len(list_deg_codons))}',
                                '')
                    st.metric('Number of degenerate nucleotides',
                                f'{int(n_deg_nucl)}',
                                '')
                #show use-case explanation for all, KRED, aromatic, and polar
                ############
                explain=True
                if AAset == set('SNIRHLGDVCYFKTQPEAMW'):
                    str_example = examples['SNIRHLGDVCYFKTQPEAMW']
                elif AAset == set('KRED'):
                    str_example = examples['KRED']
                elif AAset == set('STCNQ'):
                    str_example = examples['STCNQ']
                elif AAset == set('GAVLMIP'):
                    str_example = examples['GAVLMIP']
                elif AAset == set('FYW'):
                    str_example = examples['FYW']
                else:
                    explain = False
                if explain:
                    st.markdown(f'<b>Real-world application</b>\
                                <br/>{str_example}',
                                True)
            #show legend in an expander container
            ############
            with st.expander("Legend"):
                col1,col2 = st.columns(2)
                #explain what is meant by equivalent codons
                col1.markdown(f'<br/><b>Equivalent codons:</b><ul><li> code \
                       the same aminoacids in same probabilities</li></ul>\
                       <ul><li>could be exchanged and the ratio would not \
                       change</li></ul><ul><li>have the same number of\
                       degenerated nucleotides</li></ul>', True)
                
                nucl_df = pd.DataFrame(list(deg_nucl.items()),
                                columns=(['Base IUB Code','Bases Represented']))
                #show table with deg nucleotides and bases they represent
                ############
                col2.table(nucl_df)
                
            
                    
            #create pie charts
            specs = [[{'type':'domain'}, {'type':'domain'}, {'type':'domain'}]]
            in_colors = ['rgb(248,196,103)', 'rgb(252,237,209)']
            out_colors = ['rgb(238,86,86)', 'rgb(248,187,187)']
            fig = make_subplots(rows=1, 
                        cols=3, 
                        specs=specs,
                        subplot_titles=[
                            '(A) Number of AAs <br /> coded by each codon', 
                            '(B) Percentage of <br /> AAset coded', 
                            '(C) Percentage of AAs <br /> coded outside AAset'
                            ])
            fig.add_trace(go.Pie(
                        labels = list(aux.keys()),
                        values = list(aux.values()),
                        hoverinfo = "label+percent",
                        textinfo = "label+value"
                        ),1,1)
        
            fig.add_trace(go.Pie(
                        labels=['Coded AAs','Non coded AAs'], 
                        values=[len(total_coded_aas) - not_in, 
                                len(AAset)-(len(total_coded_aas)-not_in)],
                        hoverinfo = "label+value",
                        textinfo = "percent+label",marker_colors=in_colors
                        ),1,2)
            
            fig.add_trace(go.Pie(
                        labels=[f'Inside AAset {int((1-not_in/len(AAset))*100)}%',
                                f'Outside AAset {int((not_in/len(AAset))*100)}%'], 
                        values=[not_in,
                                len(total_coded_aas) - not_in],
                        hoverinfo = "label",
                        textinfo = "label",
                        marker_colors=out_colors
                        ),1,3)
            fig.update( layout_showlegend=False)
            
            #show pie charts
            ############
            st.plotly_chart(fig)
            
            #show legend for pie charts
            ############
            with st.expander("See figure interpretion"):
                st.markdown("<b>Pie Chart A:</b> shows the number of \
                    aminoacids that each degenerated codon codes for.<br />\
                    <b>Pie Chart B:</b>shows the percentage of the aminoacids,\
                    from the AA set input, that are coded by the codon set \
                    designed.<br /><b>Pie Chart C:</b> shows the percentage of \
                    the aminoacids, coded by the codon set designed, that are \
                    not included in the input AA set.<br /> *'AAs' means 'amino acids'"
                    , True)
                    
    #show error if invalid input
    ############
    else:
        st.error('Invalid input')

if __name__ == '__main__':
    main()

    
