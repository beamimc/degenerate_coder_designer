#! /usr/bin/env python
"""
@author: beamimc

"""
from src.dcg import run_generateCodon
from src.global_variables import aas

def main():
    """ To run Degenerate Codon Designer from cmd """
    
    good_example = 'GAVCPLIMWF' #### best by brute force is 4 codons {'ATK': 1.0, 'TKT': 1.0, 'KGG': 1.0, 'SYA': 2.0}
    example = 'GAVCPLIMWFKRED' ## best is 4 codons
    exampl = 'STYNQKRED'
    
    all_AA = 'SNIRHLGDVCYFKTQPEAMW'
    charged = 'KREDH'
    polar = 'STCNQ'  
    non_polar = 'GAVLMIP'
    aromatic = 'FWY'
    
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
    -- polar: STCNQ \n \
    -- non polar: GAVLMIP (alifatic)\n \
    -- aromatic: FYW')
    print('E.g: input "all" would be equal to input SNIRHLGDVCYFKTQPEAMW\n')
    #get input
    AAset = input('Introduce AAset: ').upper()
    # check input
    if AAset == 'ALL':
        AAset =set(all_AA)
    elif AAset =='CHARGED':
        AAset = set(charged)
    elif AAset == 'POLAR':
        AAset = set(polar)
    elif AAset == 'NON POLAR':
        AAset = set(non_polar)
    elif AAset == 'AROMATIC':
        AAset = set(aromatic)
    else:
        ## check if input is correct
        AAset = set(AAset)
        for aa in AAset:
            if aa not in aas:
                return print('The introduced input was not a valid AAset')
            
    #if input is ok run the program
    # startTime = time.time()
    #execution
    sorted_combi_prop, eq_str, cd_str = run_generateCodon(AAset)
    
    #create string to show user
    combi_str = ",".join(str(elem[0]) for elem in sorted_combi_prop)
    prop_str = ":".join(str(int(elem[1])) for elem in sorted_combi_prop)
    
    #show result
    print(f'\nOutput: {combi_str}  ratio {prop_str}\n')
    # executionTime = (time.time() - startTime)
    
    #provide equivalent codons and aas coded by each codon 
    print(cd_str)
    print(eq_str)
    
    return
if __name__ == '__main__':
    main()
