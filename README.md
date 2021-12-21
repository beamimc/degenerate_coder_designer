# Degenerate Codon Designer
A codon is a sequence of three DNA/RNA nucleotides that corresponds with an specific aminoacid or stop signal during protein synthesis. Naturally, a codon codes just for one aminoacid. Degenerate nucleotides represent that a certain position can have multiple alternatives from the canonical bases (ATCG). A degenerate codon is a sequence of three degenerate nucleotides which, as a result, encodes for (one or) more than one amino acid.

**This program lets you calculate the minimun set of degenerate codons that would code a given set of amino acids with equal probabilities.**

## Installation
- Download the repository or clone it into prefered location

`git clone https://github.com/beamimc/degenerate_codon_designer.git`

### Requirements
- python version >= 3.8
- libraries needed: 
  ```
  pandas==1.3.3
  plotly==4.14.3
  streamlit==1.3.0
  ```
  - install each by pip install 
  - install all by: `pip install -r requirements.txt` 

## Description of files
 - requirements.txt to install all dependencies
 - run_dcg.py is for running program on python enviroment / terminal
 - gui.py is for running program on streamlit
 - create_deg_codons_DB.py is program made to create datasets folder
 ### /src (module)
- global_variables.py has global variables of the module
- dcg.py has the degenerateCodon function that calculates the minimum set of degenerated codons given an AAset
### /datasets 
- Has deg_codons_DB_clean.csv and deg_codons_DB.csv used during the program, created by create_deg_codons_DB.py
### /static
Has images used for the GUI version
- chart.png source is: https://che.gg/3Hhut5V
- logo.svg was made by me
### /.streamlit 
- config.toml for setting the main theme for the GUI version
# How to run 
A GUI version of this program was developed, creating an web app using Streamlit (https://streamlit.io/), an open-source Python library.

A Terminal version was also implemented, so as to run the program on the terminal or python interpreter.

## Run GUI version Locally (Recomended)
1. On Terminal, move into the folder of the repository downloaded/cloned
2. Type command `streamlit run gui.py` 
3. First time running streamlit a message like this will show. Press enter to skip.

<img width="445" alt="image" src="https://user-images.githubusercontent.com/59894638/146967097-c5956c7e-2a2a-45db-9c68-d31ff189adc4.png">

4. A message like the following should show up

<img width="308" alt="image" src="https://user-images.githubusercontent.com/59894638/146967373-06624fc1-ee49-495e-a146-097f1abcb84a.png">

5. A new tab in your default browser will automatically open with the app. You can also open it by copying the **Local URL** on your browser. 
Please, note that Safari is not completely compatible. Preferably use Chrome or Firefox.
6. Main page will show, where the aminoacids can be selected from a list, or from an input of the user (Custom option)

![image](https://user-images.githubusercontent.com/59894638/146930685-f2231532-e74b-44eb-96a5-6fca045a0f3f.png)
![image](https://user-images.githubusercontent.com/59894638/146930447-7f993bae-dd46-4a10-a41d-ead48ba991c3.png)

7. Click Design button to run the program and obtain the set of optimized degenerated codons 
8. Results will show up, along with some extra information
9. In case we select any of the provided options (i.e all/KRED/polar/non polar/aromatic), the result will come with a brief real-work application example

![image](https://user-images.githubusercontent.com/59894638/146986142-27a9fc0a-b9a1-4c9a-8997-8ce79817fdc7.png)

 ## Run Terminal version
 1. On Terminal, move into the folder of the repository downloaded/cloned
 
 2. Type command `python run_dcg.py` 

 3. Program will start and ask AAset to the user 

![image](https://user-images.githubusercontent.com/59894638/146986336-abb68a5d-4179-46e2-9602-555c6986bd90.png)

4. After input the AAset, program will run and show the results

![image](https://user-images.githubusercontent.com/59894638/146986492-060aade6-3451-4b40-9870-cd6dff7aac15.png)
