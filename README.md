# degenerate_codon_designer

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
# Run program 
 ## Run GUI version running on server (Recomended)
 1. Open this url in your browser 
 
 ` `
 
## Run GUI version locally
1. On Terminal, move into the folder of the repository downloaded/cloned
2. Type command `streamlit run gui.py` 

3. A message like the following should show up![146851593-e025dad9-566c-421f-800d-445ce048feab](https://user-images.githubusercontent.com/59894638/146930763-a0a605a0-d346-4de0-a436-389cc4859176.png)


![image](https://user-images.githubusercontent.com/59894638/146930087-808efcd6-4a7c-4ac6-b16a-4201914d9988.png)
 
4. A new tab in your default browser will automatically open with the app. You can also open it by copying the **Local URL** on your browser. 
Please, note that Safari is not completely compatible. Preferably use Chrome or Firefox.
5. Main page will show, where the aminoacids can be selected from a list, or from an input of the user (Custom option)

![image](https://user-images.githubusercontent.com/59894638/146930685-f2231532-e74b-44eb-96a5-6fca045a0f3f.png)
![image](https://user-images.githubusercontent.com/59894638/146930447-7f993bae-dd46-4a10-a41d-ead48ba991c3.png)

6. Click Design button to run the program and obtain the set of optimized degenerated codons 
7. Results will show up, along with some extra information


 ## Run Terminal version
 1. On Terminal, move into the folder of the repository downloaded/cloned
 
 2. Type command `python run_dcg.py` 

 3. Program will start and ask AAset to the user 

![image](https://user-images.githubusercontent.com/59894638/146851593-e025dad9-566c-421f-800d-445ce048feab.png)

4. After getting the AAset, program will run and show the results

![image](https://user-images.githubusercontent.com/59894638/146852342-924f5d8d-fabd-4e1a-8012-372ff9dfd437.png)

