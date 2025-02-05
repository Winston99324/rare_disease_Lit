# https://pmc.ncbi.nlm.nih.gov/articles/PMC8374293/


from Bio import Entrez
import time
import json
import os  
import traceback
from tqdm import tqdm 

Entrez.email = "yilianz@uci.edu"

# Rename disease names for recorded files 
def sanitize_filename(name):
    
    return "".join(c if c.isalnum() or c in ("-", "_") else "_" for c in name).strip()

# Separate list of diseases into n piese 
def split_list(lst, n):
    
    k, m = divmod(len(lst), n)
    return [lst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n)]

def error_log(disease, error_msg):
    error_file = "./pubmed_search/error_log.txt"
    with open(error_file, "a") as f:
        f.write(f"Error processing {disease}: {error_msg}\n")

def fetch_pubmed_data(webenv, query_key, total_count, batch_size=10000):
    all_data = ""
    for start in range(0, total_count, batch_size):
        try:
            with Entrez.efetch(
                db="pubmed",
                rettype="medline",
                retmode="text",
                retstart=start,
                retmax=batch_size,
                webenv=webenv,
                query_key=query_key
            ) as fetch_handle:
                all_data += fetch_handle.read()
            time.sleep(0.01)  
        except Exception as e:
            error_message = f"{e}\n{traceback.format_exc()}"
            error_log("Fetching batch", error_message)
            return None  
    return all_data

# This function will output 4 files: 1. Total research count and valid count in a JSON 2. abnormal finding log 3. each diseases' search separately txt file 4. 25 groups of combined txt file for uploading
def search_and_save_all_pubmed(disease, output_filename, json_filename, abnormal_file_log, start_year=2015):
    time.sleep(0.01)
    print(f"Processing: {disease}")
    try:
        # set up base search term
        base_search_term = (
            '("electronic health record" OR "electronic health records" OR "electronic medical record" OR "electronic medical records" OR "EHR" OR "EMR") '
            #'AND ("rare disease" OR "rare diseases" OR "undiagnosed disease" OR "undiagnosed diseases") '
            'AND ("{rare_disease}"[Title/Abstract])'
            'AND ( "language model" OR "large language model" OR llm OR llms OR "transformer model" OR "machine learning" OR nlp OR "natural language processing" OR "artificial intelligence" OR phenotyp*)'
            #"AND (large language model OR language model OR NLP OR natural language processing OR language processing OR machine learning OR artificial intelligence OR predictive modeling OR deep learning)"
        )
        search_term = base_search_term.format(rare_disease=f'{disease.lower()}')
        
        date_range = f"{start_year}/01/01:3000/12/31[dp]"  
        search_query = f"{search_term} AND {date_range}"

        
        with Entrez.esearch(db="pubmed", term=search_query, retmax=100000, usehistory="y") as search_handle:
            search_results = Entrez.read(search_handle)
        
        disease_total_count = int(search_results["Count"])


        # Append research reults number to Json 
        if os.path.exists(json_filename):
            with open(json_filename, "r", encoding="utf-8") as json_file:
                existing_data = json.load(json_file)
        else:
            existing_data = {"valid_count": 0,"total_count": 0, "results": []}

        existing_data["results"].append({"disease_name": disease,"search_query": search_query, "count": disease_total_count})

        # Get Abnormal count out of calculations, specified how many total research achived and how many real we need  
        if disease_total_count >= 350:
            existing_data["total_count"] += disease_total_count
        else:
            existing_data["valid_count"] += disease_total_count
            existing_data["total_count"] += disease_total_count
            
        # Record every search query result't coutn num 
        with open(json_filename, "w", encoding="utf-8") as json_file:
            json.dump(existing_data, json_file, indent=4)

        # Stop if none result found     
        if disease_total_count == 0:
            return
        
        # Check for abnormal & record 
        if disease_total_count >= 350:
                #print(f"Abnormal count detected for {disease}: {disease_total_count}.")
                with open(abnormal_file_log, "a") as ab_file:
                    ab_file.write(f"Abnormal count detected for {disease}: {disease_total_count}.\n")
                return
        
        # Use WebEnv and QueryKey for batch fetching
        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]

        all_data = fetch_pubmed_data(webenv, query_key, disease_total_count)
        if all_data is None:
            error_log(disease, "Fetching failed. Skipping.")
            return
        
                
        # Create separate file for recording 
        disease_dir = './pubmed_search/all_diseases'
        os.makedirs(disease_dir, exist_ok=True)

        disease_file = os.path.join(disease_dir, f"{sanitize_filename(disease)}.txt")
        with open(disease_file, "w") as file:
            file.write(all_data)
            
        with open(output_filename, "a") as all_file:
            all_file.write(all_data)

    except Exception as e:
        error_message = f"{e}\n{traceback.format_exc()}"
        error_log(disease, error_message)
        return
    
if __name__ == "__main__":

    path_file = "./disease_json/data.json"

    ab_dir = "./pubmed_search/abnormal_diseases"
    os.makedirs(ab_dir, exist_ok=True)

    abnormal_file_log = os.path.join(ab_dir, f"abnormal_log.txt")
    output_json_filename = f"./pubmed_search/pubmed_search_summary.json"

    # Read Orphanet Rare Disease names 
    with open(path_file, "r", encoding="utf-8") as file:
        disease_data = json.load(file)
    rare_diseases = [d['Name'] for d in disease_data['disease']]

    total_diseases = len(rare_diseases)
    start_time = time.time()

    upload_dir = "./pubmed_search/file_4_upload"
    os.makedirs(upload_dir, exist_ok=True)
    
    #split_disease_groups = split_list(rare_diseases, 25)
    #for idx, group in enumerate(split_disease_groups):
    
    output_filename = os.path.join(upload_dir, f"pubmed_output_all.txt")
  
    #for disease in rare_diseases:
    for disease in tqdm(rare_diseases, desc="Processing Diseases", unit="disease", ncols=80):
        search_and_save_all_pubmed(disease, output_filename, output_json_filename, abnormal_file_log)
        


    