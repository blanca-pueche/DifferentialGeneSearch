import streamlit as st
import pandas as pd
import requests
import gseapy as gp
from Bio import Entrez
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import urllib.parse
import plotly.express as px
import base64
import threading
import os
import urllib.parse
import gseapy as gp
import re
import io, zipfile, datetime
from collections import defaultdict
from streamlit_autorefresh import st_autorefresh
import re
import urllib.parse
import streamlit.components.v1 as components


GRAPHQL_URL = "https://dgidb.org/api/graphql"
ENSEMBL_LOOKUP_URL = "https://rest.ensembl.org/lookup/id/"


# Methods
def get_disease_name(mesh_id):
    """
    Returns the corresponding name of a disease given a MESH id provided by the user. 
    """
    try:
        search_handle = Entrez.esearch(db="mesh", term=mesh_id)
        search_record = Entrez.read(search_handle)
        search_handle.close()
        if not search_record['IdList']:
            return None
        uid = search_record['IdList'][0]
        summary_handle = Entrez.esummary(db="mesh", id=uid)
        summary_record = Entrez.read(summary_handle)
        summary_handle.close()
        disease_url = f"https://www.ncbi.nlm.nih.gov/mesh/{uid}"
        return summary_record[0]['DS_MeshTerms'][0], disease_url
    except:
        return None
    
def generate_expression_atlas_link(disease_name):
    """
    Generates link to Expression Atlas for a given disease to extract gene information. 
    """
    encoded_disease = urllib.parse.quote(disease_name)
    base_url = (
        "https://www.ebi.ac.uk/gxa/search?geneQuery=%5B%5D"
        "&species=Homo%20sapiens"
        f"&conditionQuery=[%7B%22value%22%3A%22{encoded_disease}%22%7D]"
        "&ds=%7B%22kingdom%22%3A%5B%22animals%22%5D%2C%22regulation%22%3A%5B%22UP%22%5D%7D"
        "&bs=%7B%22homo%20sapiens%22%3A%5B%22ORGANISM_PART%22%5D%7D"
        "#differential"
    )
    return base_url

def get_gene_name_from_ensembl(ensembl_id):
    """
    Given an ensembl id it extracts the corresponding gene name.
    """
    response = requests.get(f"{ENSEMBL_LOOKUP_URL}{ensembl_id}?content-type=application/json")
    if response.status_code == 200:
        data = response.json()
        return data.get("display_name", "Not Found")
    elif response.status_code == 429:
        st.warning("Error 429: Resource Exceeded. Please try again later.")
        st.stop()
    return "Not Found"

def fetch_gene_names(df):
    """
    Fetches gene names from ensembl and groups them according to log_2 fold change.
    """

    df["Gene Name"] = df["Gene"].apply(get_gene_name_from_ensembl)
    df_filtered = df[df["Gene Name"] != "Not Found"]

    # Sum log2 fold changes for repeated genes
    grouped = df_filtered.groupby(["Gene", "Gene Name"], as_index=False).agg({
        "log_2 fold change": "sum"
    })

    return grouped


def find_possible_target_of_drugs(ensembl_id):
    """
    Given an Ensembl Gene ID, it asks the OpenTragetS API to check if it's a drug target.
    Returns the gene symbol, whether it's a known drug target, and associated approved drugs.
    """
    #url = f"https://platform-api.opentargets.io/v3/platform/public/target/{ensembl_id}/associations"
    #url= f"https://api.platform.opentargets.org/api/v4/graphql/browser?query={ensembl_id}"
    


    # Build query string to get general information about AR and genetic constraint and tractability assessments 
    query_string = """
      query target($ensemblId: String!){
        target(ensemblId: $ensemblId){
          id
          approvedSymbol
          biotype
          geneticConstraint {
            constraintType
            exp
            obs
            score
            oe
            oeLower
            oeUpper
          }
          tractability {
            label
            modality
            value
          }
        }
      }"""

    # Set variables object of arguments to be passed to endpoint
    variables = {"ensemblId": ensembl_id}

    # Set base URL of GraphQL API endpoint
    base_url = "https://api.platform.opentargets.org/api/v4/graphql"

    # Perform POST request and check status code of response
    try:
        r = requests.post(base_url, json={"query": query_string, "variables": variables})
        if r.status_code != 200:
            return None
        data = r.json()['data']['target']
        return {
            "Gene Symbol": data.get("approvedSymbol", ""),
            "Ensembl ID": data.get("id", ""),
            "Name": data.get("approvedName", ""),
            "Biotype": data.get("biotype", ""),
            "Tractability": [
                t["label"] for t in data.get("tractability", []) if t["value"]
            ]
        }
    except:
        return None


def analyze_pathways(df, number):
    """
    Analyzes n (given number) pathways in which the the genes interact.
    """
    # Prepare gene list
    df_genes = df["Gene Name"].dropna().astype(str).str.strip().str.upper().unique().tolist()

    # Enrichment
    enr = gp.enrichr(gene_list=df_genes, gene_sets="Reactome_2022", organism="Human", outdir=None)
    if enr.results.empty:
        return None, None

    top_pathways = enr.results.copy()

    # Extract Reactome ID and create direct link
    def make_reactome_link(term):
        match = re.search(r'(R-HSA-\d+)', term)
        if match:
            rid = match.group(1)
            return f'<a href="https://reactome.org/PathwayBrowser/#/{rid}" target="_blank">{term}</a>'
        return term

    top_pathways["Reactome Link"] = top_pathways["Term"].apply(make_reactome_link)


    # Parse overlap info
    top_pathways[["Input Genes", "Pathway Genes"]] = top_pathways["Overlap"].str.split("/", expand=True).astype(int)
    top_pathways["Input %"] = top_pathways["Input Genes"] / top_pathways["Pathway Genes"] * 100
    top_pathways = top_pathways.sort_values("Adjusted P-value", ascending=True).head(number)

    # Sum log2fc for overlapping genes per pathway
    df["Gene Name"] = df["Gene Name"].astype(str).str.strip().str.upper()
    sum_fc = []
    for _, row in top_pathways.iterrows():
        genes = [g.strip().upper() for g in row["Genes"].split(";")]
        overlap_df = df[df["Gene Name"].isin(genes)]
        sum_fc.append(overlap_df["log_2 fold change"].sum())

    top_pathways["Sum log2fc"] = sum_fc

    return top_pathways

def get_overlapping_genes(df, selected_pathway_row):
    """
    Get overlapped genes as a way to identify most important genes in a given specific pathway
    """
    # Normalize input genes
    df["Gene Name"] = df["Gene Name"].astype(str).str.strip().str.upper()

    # Get genes from selected pathway
    pathway_genes = [g.strip().upper() for g in selected_pathway_row["Genes"].split(";")]

    # Get overlapping genes from input
    overlap_df = df[df["Gene Name"].isin(pathway_genes)].copy()
    overlap_df["abs_fc"] = overlap_df["log_2 fold change"].abs()
    overlap_df = overlap_df.sort_values("abs_fc", ascending=False)

    return overlap_df



def get_drug_targets_dgidb_graphql(gene_names):

    """
    Queries the DGIdb GraphQL API for drug–gene interaction data for a list of gene names.

    For each gene, the function retrieves associated drugs, interaction types, 
    directionality, interaction scores, sources, and PMIDs (if available), 
    and compiles the results into a pandas DataFrame.
    """

    all_results = []
    for gene_name in gene_names:
        graphql_query = f"""
        {{
          genes(names: ["{gene_name}"]) {{
            nodes {{
              interactions {{
                drug {{
                  name
                  conceptId
                }}
                interactionScore
                interactionTypes {{
                  type
                  directionality
                }}
                interactionAttributes {{
                  name
                  value
                }}
                publications {{
                  pmid
                }}
                sources {{
                  sourceDbName
                }}
              }}
            }}
          }}
        }}
        """
        response = requests.post(GRAPHQL_URL, json={"query": graphql_query})
        if response.status_code == 200:
            try:
                data = response.json()
                if 'data' in data and 'genes' in data['data'] and len(data['data']['genes']['nodes']) > 0:
                    interactions = data['data']['genes']['nodes'][0].get('interactions', [])
                    for interaction in interactions:
                        drug_name = interaction['drug'].get('name', 'Unknown')
                        score = interaction.get('interactionScore', 'N/A')
                        types = interaction.get('interactionTypes', [])
                        interaction_type = types[0].get('type', 'N/A') if types else 'N/A'
                        directionality = types[0].get('directionality', 'N/A') if types else 'N/A'
                        sources = interaction.get('sources', [])
                        source = sources[0]['sourceDbName'] if sources else 'N/A'
                        pmids = interaction.get('publications', [])
                        pmid = pmids[0]['pmid'] if pmids else 'N/A'

                        all_results.append({
                            'Gene': gene_name,
                            'Drug': drug_name,
                            'Interaction Type': interaction_type,
                            'Directionality': directionality,
                            'Source': source,
                            'PMID': pmid,
                            'Interaction Score': score
                        })
            except:
                continue
    
    return pd.DataFrame(all_results)

def drug_with_links(df):
    df_with_links = df.copy()
    df_with_links["Gene"] = df_with_links["Gene"].apply(
            lambda gene: f'<a href="https://dgidb.org/results?searchType=gene&searchTerms={urllib.parse.quote(gene)}" target="_blank">{gene}</a>'
            if gene else ""
    )
    df_with_links["Drug"] = df_with_links["Drug"].apply(
        lambda drug: f'<a href="https://dgidb.org/results?searchType=drug&searchTerms={urllib.parse.quote(drug)}" target="_blank">{drug}</a>'
    )
    df_with_links["PMID"] = df_with_links["PMID"].apply(
        lambda pmid: f'<a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}/" target="_blank">{pmid}</a>'
        if pmid else "NaN"
    )
    return df_with_links
    
    
# Estimate table height based on number of rows --> dynamic change for better presentation
def estimate_table_height(df, max_visible_rows=10, base_row_height=50, tall_row_height=70, overhead=70, min_height=350, max_height=1200): 
    visible_rows = min(len(df), max_visible_rows) 
    
    if df.shape[1] >= 2: 
        lengths = df.iloc[:visible_rows, :2].astype(str).map(len) 
    else: 
        lengths = df.iloc[:visible_rows].astype(str).map(len) 
    
    tall_rows = (lengths > 60).any(axis=1).sum() 
    normal_rows = visible_rows - tall_rows 
    row_height = (normal_rows * base_row_height) + (tall_rows * tall_row_height) 
    
    estimated_height = row_height + overhead 
    # Clamp height between min and max 
    return max(min(estimated_height, max_height), min_height)

def normalize_disease_name(name: str) -> str:
    # If there's a comma, move the part after the comma to the front
    if "," in name:
        parts = [part.strip() for part in name.split(",")]
        # e.g., ["Muscular Dystrophy", "Duchenne"] → "Duchenne Muscular Dystrophy"
        return f"{parts[1]} {parts[0]}"
    return name


def add_links_to_final_table(df):
    df = df.copy()
    
    # Gene → Ensembl
    if "Ensembl ID" in df.columns:
        df["Ensembl ID"] = df["Ensembl ID"].apply(
            lambda gene_id: f'<a href="https://www.ensembl.org/Multi/Search/Results?q={urllib.parse.quote(str(gene_id))}" target="_blank">{gene_id}</a>'
        )
    
    # Drug → DGIdb
    if "Drug" in df.columns:
        def drug_links(drugs_str):
            if pd.isna(drugs_str):
                return ""
            drugs = [d.strip() for d in drugs_str.split(";")]
            return "; ".join([f'<a href="https://dgidb.org/results?searchType=drug&searchTerms={urllib.parse.quote(d)}" target="_blank">{d}</a>' for d in drugs])
        df["Drug"] = df["Drug"].apply(drug_links)
    
    # PMID → PubMed
    if "PMID" in df.columns:
        def pmid_links(pmids_str):
            if pd.isna(pmids_str):
                return ""
            pmids = [p.strip() for p in str(pmids_str).split(";")]
            return "; ".join([f'<a href="https://pubmed.ncbi.nlm.nih.gov/{p}" target="_blank">{p}</a>' for p in pmids])
        df["PMID"] = df["PMID"].apply(pmid_links)
    
    # Pathways → Reactome
    if "Pathways" in df.columns:
        def pathway_links(pathways_str):
            if pd.isna(pathways_str):
                return ""
            pathways = [p.strip() for p in pathways_str.split(";")]
            linked_pathways = []
            for p in pathways:
                match = re.search(r"(R-HSA-\d+)", p)  # Extract Reactome ID
                if match:
                    reactome_id = match.group(1)
                    linked_pathways.append(f'<a href="https://reactome.org/content/detail/{reactome_id}" target="_blank">{p}</a>')
                else:
                    linked_pathways.append(p)
            return "; ".join(linked_pathways)
        df["Pathways"] = df["Pathways"].apply(pathway_links)
    
    return df

def save_pathway_csvs(df_selected, top_pathways, output_dir="all_pathway_genes_csvs"):
    os.makedirs(output_dir, exist_ok=True)
    
    for _, row in top_pathways.iterrows():
        pathway_name = row["Term"]
        pathway_genes = get_overlapping_genes(df_selected, row)
        pathway_genes = pathway_genes.drop(columns=["Gene Name_raw", "abs_fc"])
        pathway_genes_newNames = pathway_genes.rename(columns={"Gene":"Ensembl ID", "Gene Name":"Gene"})
        safe_name = pathway_name.replace("/", "_").replace(" ", "_")
        csv_path = os.path.join(output_dir, f"{safe_name}.csv")
        pathway_genes_newNames.to_csv(csv_path, index=False)

    # Signal completion
    with open(os.path.join(output_dir, ".done"), "w") as f:
        f.write("done")

def save_drug_csvs(df_selected, top_pathways, output_dir="all_drug_csvs"):
    os.makedirs(output_dir, exist_ok=True)
    
    for _, row in top_pathways.iterrows():
        pathway_name = row["Term"]
        pathway_genes = get_overlapping_genes(df_selected, row)
        drug_df = get_drug_targets_dgidb_graphql(pathway_genes["Gene Name"].tolist())
        if not drug_df.empty:
            drug_df = drug_with_links(drug_df)
            safe_name = pathway_name.replace("/", "_").replace(" ", "_")
            csv_path = os.path.join(output_dir, f"{safe_name}_drugs.csv")
            drug_df.to_csv(csv_path, index=False)

    # Signal completion
    with open(os.path.join(output_dir, ".done"), "w") as f:
        f.write("done")


def create_combined_zip(zip_path, folders=[], extra_files=[]):
    """
    Create a zip containing all CSVs from the given folders + any extra files.
    
    Parameters:
        zip_path (str): Path where the zip will be saved.
        folders (list of str): List of folder paths containing CSVs.
        extra_files (list of str): List of extra CSV file paths to include.
    """
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        # Add all CSVs from each folder
        for folder in folders:
            for root, _, files in os.walk(folder):
                for file in files:
                    if file.endswith(".csv"):
                        file_path = os.path.join(root, file)
                        zipf.write(file_path, arcname=file)
        
        # Add any extra individual files
        for file_path in extra_files:
            if os.path.exists(file_path):
                zipf.write(file_path, arcname=os.path.basename(file_path))


