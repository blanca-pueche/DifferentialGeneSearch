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

Entrez.email = "blanca.puechegranados@usp.ceu.es"
GRAPHQL_URL = "https://dgidb.org/api/graphql"
ENSEMBL_LOOKUP_URL = "https://rest.ensembl.org/lookup/id/"

st.markdown("""
    <style>
    .main {
        background: linear-gradient(to bottom right, #e3f2fd, #fce4ec);
        font-family: 'Arial', sans-serif;
        padding: 20px;
    }
    .stButton > button {
        background-color: #6a1b9a;
        color: white;
        font-weight: bold;
        border-radius: 10px;
    }
    .stDownloadButton > button {
        background-color: #2e7d32;
        color: white;
        font-weight: bold;
        border-radius: 10px;
    }
    </style>
""", unsafe_allow_html=True)

# Methods
def get_disease_name(mesh_id):
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
        return summary_record[0]['DS_MeshTerms'][0]
    except:
        return None
    
def generate_expression_atlas_link(disease_name):
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
    response = requests.get(f"{ENSEMBL_LOOKUP_URL}{ensembl_id}?content-type=application/json")
    if response.status_code == 200:
        data = response.json()
        return data.get("display_name", "Not Found")
    return "Not Found"

def fetch_gene_names(df):
    df["Gene Name"] = df["Gene"].apply(get_gene_name_from_ensembl)
    df_filtered = df[df["Gene Name"] != "Not Found"]

    # Sum log2 fold changes for repeated genes
    grouped = df_filtered.groupby(["Gene", "Gene Name"], as_index=False).agg({
        "log_2 fold change": "sum"
    })

    return grouped


    
    
    

st.title("🧬 Disease-Gene-Drug Analysis Dashboard") #TODO change name

# Step 1: Ask for e-mail
email = st.text_input("Enter user e-mail:")

#Step 2: Name from MeshID
mesh_id = st.text_input("🔍 Enter MeSH ID (e.g., D003920 for Diabetes Mellitus):")

if mesh_id:
    with st.spinner("Fetching disease information..."):
        disease = get_disease_name(mesh_id)
        
        if disease:
            st.success(f"🎯 Disease **{disease}** identified")
            
            # Step 3: File upload only after disease name is given
            url = generate_expression_atlas_link(disease_name=disease)
            uploaded_file = st.file_uploader(
                f"📁 Upload Differential Expression File (TSV from Expression Atlas: [link]({url}))",
                type=["tsv"]
            )
            
            if uploaded_file:
                df_raw = pd.read_csv(uploaded_file, sep="\t")
                st.write("📥 Downloaded correctly")
                
                with st.spinner("🔄 Mapping Ensembl IDs to gene names..."):
                    df_selected = fetch_gene_names(df_raw)
                # Convert the Ensembl gene IDs into HTML links
                df_selected_with_links = df_selected.copy()
                df_selected_with_links["Gene"] = df_selected_with_links["Gene"].apply(
                    lambda gene_id: f'<a href="https://www.ensembl.org/Multi/Search/Results?q={gene_id}" target="_blank">{gene_id}</a>'
                )
                df_selected_with_links = df_selected_with_links.sort_values(by="log_2 fold change", ascending=False)
                
                st.markdown("### 🧬 Gene Table with Links to Ensembl")

                # Crear la tabla HTML con escape=False para mantener los links
                html_table = df_selected_with_links.to_html(escape=False, index=False)

                # Envolver la tabla en un div con scroll y altura fija
                html_with_scroll = f"""
                    <div style="max-height: 400px; overflow-y: auto;">
                        {html_table}
                    </div>
                    """

                # Mostrar la tabla con scroll y permitiendo HTML
                st.markdown(html_with_scroll, unsafe_allow_html=True)

                # Botón para descargar CSV
                st.download_button(
                    "📥 Download Genes CSV",
                    df_selected_with_links.to_csv(index=False),
                    "genes.csv",
                    "text/csv"
                )

