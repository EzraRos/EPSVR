from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import StaleElementReferenceException, TimeoutException
import pandas as pd
import time
import os
import shutil
import sys
import pickle
import os

def click_with_retry(driver, by, value, retries=3):
    for _ in range(retries):
        try:
            element = WebDriverWait(driver, 20).until(EC.visibility_of_element_located((by, value)))
            element.click()
            return
        except StaleElementReferenceException:
            print("Stale element. Retrying...")
            time.sleep(1)
        except TimeoutException as e:
            print(f"Timed out waiting for element: {e}")
            sys.exit("Script terminated due to timeout")


#Set directories
filepath= "/Users/ezrarosenbaum/Desktop/Research"

info_df=pd.read_csv(filepath+"/Test_Set/test_with_cutoff_number_annotated_and_size.csv")
annotation_folder=(filepath+"/unbound_annotations")# the folder that contains the experimental results(only residues at surface level noted)
pdb_folder=(filepath+"/unbound_pdbs")# the folder that contains the pdbs of the protein test sets
protein_outputfile=(filepath+"/Test_Set/EPSVR_with_annotated_tset.csv")
pickle_file=filepath+"/Test_Set/EPSVR_count_tset.pickle"

# Saving filenames to list
annotated_og=[filename for filename in os.listdir(annotation_folder)]# saving original filenames to list
annotation=[((filename.replace(".txt","")).replace("_annotated","")).replace("_",".") for filename in os.listdir(annotation_folder)]# formatting to just keep protein name
annotated_dict= {annotation[i]:annotated_og[i] for i in range(len(annotated_og))} # creates a dict so the new name(just protein name) can be later matched up with actual file name

pdb_og=[filename for filename in os.listdir(pdb_folder)]# saving original filenames to list
pdbs=[(((filename.replace("A_","")).replace("B_","")).replace(".pdb","")).replace("D_","") for filename in os.listdir(pdb_folder)]# formatting
pdb_dict= {pdbs[i]:pdb_og[i] for i in range(len(pdb_og))} # creates a dict so the new name(just protein name) can be later matched up with actual file name

count=0 # keep track of the proteins that are missing pdbs

for protein in info_df["Name"]:
    try:
        with open(pickle_file, "rb") as file:
            proteins_done= pickle.load(file)
    except FileNotFoundError:
        proteins_done= ["Let's","begin"]
    
    done=False
    for finished in proteins_done:
        if protein == finished:
            done= True
    if not done:
        print("Begin protein analysis...")
        position= info_df[info_df["Name"]== protein].index[0]# notes the position of the protein in info_df
        print(f"{position} {protein}")

        pdb_oi= [pdb_dict[pdb] for pdb in pdbs if protein in pdb]
        annotate_oi= [annotated_dict[annotate] for annotate in annotation if protein in annotate]

        #makes sure it's looking at the correct files
        print(protein)
        print(pdb_oi)
        print(annotate_oi)

        cutoff=info_df[info_df["Name"]==protein]["Cutoff"].values
        if len(pdb_oi)>0 and len(annotate_oi)>0 and len(cutoff)>0:
            with open(filepath+f"/unbound_annotations/{annotate_oi[0]}", "r") as file:
                lines=file.readlines()
                annotate_info=[(line.split())[0] for line in lines]
                print(annotate_info)


            # Set up webDriver
            driver = webdriver.Chrome(options=webdriver.ChromeOptions())#sets up the webDriver and gives it the library of capabilities
            start_time=time.time()#sets a timer

            # Navigate to Spatom site to send the pdbs for analysis
            driver.get("http://sysbio.unl.edu/EPSVR/")
            click_with_retry(driver, By.CSS_SELECTOR, '#fileupload') #switches from PDBID to upload file

            choose_file = driver.find_element(By.CSS_SELECTOR, "#filename")
            choose_file.send_keys(filepath+f'/unbound_pdbs/{pdb_oi[0]}')#sends the pdb of interest to input file category
            chain_names = driver.find_element(By.XPATH, '/html/body/center[1]/ul/form/table/tbody/tr[8]/td[2]/input')
            chain_names.send_keys("*")
            time.sleep(4)
            click_with_retry(driver, By.XPATH, '/html/body/center[1]/ul/form/input[1]')
            driver.maximize_window()
            
            start=time.time()
            dwnld_btn_presence= WebDriverWait(driver, 2400).until(EC.presence_of_element_located((By.XPATH, "/html/body/div[2]/div/div[2]/form/input[2]")))
            end=time.time()
            print(f"Button was found after {round((end-start)/60,2)} minutes")
            time.sleep(3)
            click_dwnld=click_with_retry(driver, By.XPATH, "/html/body/div[2]/div/div[2]/form/input[2]")

            #This then renames the download to include the protein name(differentiating it from future downloads)
            time.sleep(10)
            downloaded=False
            while not downloaded:
                try:
                    rename_dwnld=os.rename("/Users/ezrarosenbaum/Downloads/prediction.pdb.crdownload",f"/Users/ezrarosenbaum/Downloads/prediction_{protein}.pdb")
                except FileNotFoundError:
                    pass
                else:
                    downloaded=True

            print("Donwloaded and moved...")
            #moving it from downloads section to backup folder for EPSVR
            pred_src= f"/Users/ezrarosenbaum/Downloads/prediction_{protein}.pdb"
            pred_dst= "/Users/ezrarosenbaum/Desktop/Research/Test_Set/EPSVR_tset_backups"
            shutil.move(pred_src,pred_dst)

            new_protein=[protein]
            proteins_done.extend(new_protein)
            with open(pickle_file, "wb") as file:
                pickle.dump(proteins_done, file)
            print(f"Finished collection {protein}")
            with open(pickle_file, "rb") as file:
                check=pickle.load(file)
                print(check)
            driver.quit()
                        
    #         xp_dwnld = WebDriverWait(driver, 120).until(EC.visibility_of_element_located((By.CLASS_NAME,'info')))#sets a timer of 120s to wait, cause Spatom takes a bit of time to load

    #         # Retrieve data from each row and save it to a csv file
    #         data_rows = driver.find_elements(By.CLASS_NAME, 'info')

    #         with open(protein_tempfile, 'w') as inputfile:
    #             inputfile.write("Number,Residue,Site,Score\n")
    #             for row in data_rows:
    #                 # td is the tag name given to each value in a row 
    #                 td_elements = row.find_elements(By.TAG_NAME, 'td')
                                        
    #                 # Loops through each value in a given row and extracts it's text
    #                 data_values = [td.text for td in td_elements]
                                        
    #                 # Join the data values with a comma and write to the file
    #                 output_line = ','.join(data_values) + '\n'
    #                 inputfile.write(output_line)

    #         driver.quit()
                    
    #         # Scaling
    #         file=pd.read_csv(protein_tempfile)
    #         values=file["Score"].tolist()
    #         highest=float(max(values))

    #         # Creating lists of the columns of interest(just Number and Score) these will be renamed into Residue and Spatom respectively
    #         res_numb= [''.join([res for res in str(val) if res.isnumeric()]) for val in file["Number"]] #removes any unecessary prefixes such as letters or dashes, we just want the numbers(represent res. position)
    #         print(res_numb)
    #         residues=[(f"{protein}-{name}") for name in res_numb]      
    #         scaled_scores=[round(float(row)/highest, 3) for row in values]         
    #         #Looping through and annotating
    #         print(residues)
    #         annotate_score=[]
    #         for residue in residues:
    #             prot=residue.split("-",1)
    #             found=False
    #             for annotate in annotate_info:
    #                 if prot[1]==annotate:
    #                     annotate_score.append("1")
    #                     found=True
    #             if not found:
    #                 annotate_score.append("0")

    #         data={
    #             "Residue":residues,
    #             "Spatom":scaled_scores,
    #             "Annotation":annotate_score
    #         }
    #         df=pd.DataFrame(data)
    #         df.to_csv(protein_outputfile,mode="a",index=False,header=False)
    #         output= pd.read_csv(protein_outputfile)
    #     else:
    #         print("It's coming here")
    #     #Adds finished protein to pickle file to keep track of the one's done
    #     new_protein=[protein]
    #     proteins_done.extend(new_protein)
    #     with open(pickle_file, "wb") as file:
    #         pickle.dump(proteins_done, file)
    #     print(f"Finished collection {protein}")
    #     with open(pickle_file, "rb") as file:
    #         check=pickle.load(file)
    #         print(check)
    # else:
    #     print(f"Already collected {protein}")
    # count+=1
    # print(f"{count}/29\n")