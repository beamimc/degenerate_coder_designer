# degenerate_codon_designer

## Instalation
- Download the repository or clone it into prefered location

`git clone https://github.com/beamimc/degenerate_codon_designer.git`

### Requirements
- python version >= 3.8.5
- libraries needed: 
  ```
  pandas==1.3.3
  plotly==4.14.3
  streamlit==1.3.0
  ```
  - install each by pip
  - install all by: `pip -r requirements.txt` 
 
 ## Run CMD version
 1. On cmd input command `python run_dcg.py` 
 2. Program will start and ask AAset to the user 

![image](https://user-images.githubusercontent.com/59894638/146851593-e025dad9-566c-421f-800d-445ce048feab.png)

3. After getting the AAset, program will run and show the results

![image](https://user-images.githubusercontent.com/59894638/146852342-924f5d8d-fabd-4e1a-8012-372ff9dfd437.png)

 ## Run GUI version
 
1. On cmd input command `streamlit run gui.py` 
2. A message like the following should show up

 ![image](https://user-images.githubusercontent.com/59894638/146851086-f0b54661-f8a2-4307-bec4-57109e26407e.png)
 
3. A new tab in your default browser will automatically open with the app. You can also open it by copying the Local URL on your browser. 
Please, note that Safari is not completely compatible. Preferably use Chrome or Firefox.
4. Main page will show, where the aminoacids can be selected from a list, or from an input of the user (Custom option)

![image](https://user-images.githubusercontent.com/59894638/146851926-da35c997-fe35-4564-b49c-32e4239287c0.png)
![image](https://user-images.githubusercontent.com/59894638/146852166-1638b2e6-9760-45cc-9129-bf2da1977199.png)

5. Click Design button to run the program and obtain the set of optimized degenerated codons 
6. Results will show up, along with some extra information

![image](https://user-images.githubusercontent.com/59894638/146852015-f289dbe8-49f8-4c5e-9e7e-223afcc90278.png)
![image](https://user-images.githubusercontent.com/59894638/146852046-423bc46a-da08-46e1-8c03-fe1cb444b149.png)


