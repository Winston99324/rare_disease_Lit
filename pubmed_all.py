# https://pmc.ncbi.nlm.nih.gov/articles/PMC8374293/


from Bio import Entrez
from datetime import datetime
import json
import os  

Entrez.email = "yilianz@uci.edu"





def search_and_save_all_pubmed(output_filename, json_filename,  start_year=2015):

    # set up base search term
    base_search_term = (
    '("electronic health record" OR "electronic health records" OR "electronic medical record" OR "electronic medical records" OR "EHR" OR "EMR") '
    'AND ("rare disease" OR "rare diseases" OR "undiagnosed disease" OR "undiagnosed diseases") '
    'AND ( "language model" OR "large language model" OR llm OR llms OR "transformer model" OR "machine learning" OR nlp OR "natural language processing" OR "artificial intelligence" OR phenotyp*)'
    #"AND (large language model OR language model OR NLP OR natural language processing OR language processing OR machine learning OR artificial intelligence OR predictive modeling OR deep learning)"
)
    search_term = base_search_term#.format(rare_disease=f'{disease.lower()}')
    
    date_range = f"{start_year}/01/01:3000/12/31[dp]"  
    search_query = f"{search_term} AND {date_range}"

    
    with Entrez.esearch(db="pubmed", term=search_query, retmax=100000, usehistory="y") as search_handle:
        search_results = Entrez.read(search_handle)
    
    total_count = int(search_results["Count"])
    print(total_count)
    #print(f"Search term: '{search_query}'")
    #print(f"Total results found: {total_count}")
    search_metadata = {
        "search_query": search_query,
        "total_results": total_count
    }

    with open(json_filename, "w") as json_file:
        json.dump(search_metadata, json_file, indent=4)

    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]

    all_data = ""
    batch_size = 10000
    
    count =0 
    for start in range(0, total_count, batch_size):
        
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
            count += 1 
    
   
    with open(output_filename, "w") as file:
        file.write(all_data)
    
    print(count)

   

if __name__ == "__main__":

    path_file = "./disease_json/data.json"
    #abnormal_file = "./pubmed_search/abnormal_diseases/abnormal_counts.txt"
    output__filename = f"./pubmed_search/pubmed_search_all.txt"
    json__filename = f"./pubmed_search/pubmed_search_all.json"
    search_and_save_all_pubmed(output__filename,json__filename)

    