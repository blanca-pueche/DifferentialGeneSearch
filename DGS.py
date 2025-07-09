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
    html, body, .stApp {
        width: 100%;
        margin: 0;
        padding: 0;
        overflow-x: hidden;
    }
    .main {
        background: linear-gradient(to bottom right, #e3f2fd, #fce4ec);
        font-family: 'Arial', sans-serif;
        padding: 20px;
    }
    .stButton > button {
        background-color: #6a1b9a;
        color: white;
        font-weight: bold;
        border-radius: 5px;
    }
    .stDownloadButton > button {
        background-color: #2e7d32;
        color: white;
        font-weight: bold;
        border-radius: 5px;
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
    gene_list = df["Gene Name"].dropna().unique().tolist()
    enr = gp.enrichr(gene_list=gene_list, gene_sets="Reactome_2022", organism="Human", outdir=None)
    if enr.results.empty:
        return None, None
    top_pathways = enr.results.sort_values("Adjusted P-value").head(number)

    top_pathways["Reactome Link"] = top_pathways["Term"].apply(
        lambda term: f'<a href="https://reactome.org/content/query?q={urllib.parse.quote(term)}" target="_blank">{term}</a>'
    )

    top_pathway_genes = top_pathways.iloc[0]["Genes"].split(";")
    important_genes = df[df["Gene Name"].isin(top_pathway_genes)].copy()
    important_genes["abs_fc"] = important_genes["log_2 fold change"].abs()
    important_genes = important_genes.sort_values("abs_fc", ascending=False)
    return top_pathways, important_genes
    
    

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
                
                # Crear una copia para conservar los datos originales si es necesario
                df_selected["Gene Name_raw"] = df_selected["Gene Name"]

                # Agrupar por "Gene Name" y sumar los valores de "log_2 fold change"
                gname_fc = df_selected.groupby("Gene Name_raw", as_index=False)["log_2 fold change"].sum()

                # Crear gráfico de barras con plotly
                fig = px.bar(
                    gname_fc,
                    x="Gene Name_raw",
                    y="log_2 fold change",
                    color="log_2 fold change",
                    color_continuous_scale="sunset",
                    title="Sum of log₂ fold change per gene",
                    labels={"Gene Name_raw": "Gene Name", "log_2 fold change": "Sum of log₂ fold change"},
                    height=500,
                    width=2000
                )

                st.plotly_chart(fig, use_container_width=True)

                         
                
                # Buscar genes con posible tractabilidad farmacológica
                with st.spinner("🔍 Checking Open Targets..."):
                    openTargets_results = df_selected["Gene"].apply(find_possible_target_of_drugs)
                    openTargets_df = pd.DataFrame([r for r in openTargets_results if r is not None])

                if not openTargets_df.empty:
                    # Eliminar la columna "Name" si existe
                    if "Name" in openTargets_df.columns:
                        openTargets_df = openTargets_df.drop(columns=["Name"])

                        openTargets_df = openTargets_df.sort_values(by="Gene Symbol")


                        # Formatear para mostrar en la tabla con píldoras
                        def format_tractability(tags):
                            if not tags:
                               return ""
                            return " ".join(
                               f'<span style="background-color:#d1c4e9; color:#4a148c; padding:4px 8px; border-radius:10px; margin:2px; display:inline-block;">{tag}</span>'
                               for tag in tags
                            )

                        openTargets_df["Tractability"] = openTargets_df["Tractability"].apply(format_tractability)
                        openTargets_df["Gene Symbol"] = openTargets_df["Gene Symbol"].apply(
                            lambda x: f'<span title="{x}">{x[:10]}...</span>' if len(x) > 10 else x
                        )


                        st.markdown("### 💊 Genes with Tractability (Open Targets)")

                        # Crear la tabla HTML bonita
                        html_ot_table = openTargets_df.to_html(escape=False, index=False)

                        html_ot_scroll = f"""
                            <div style="max-height: 500px; overflow-y: auto;"">
                                {html_ot_table}
                            </div>
                        """

                        st.markdown(html_ot_scroll, unsafe_allow_html=True)

                        # --- Aquí generamos el gráfico de barras de tractabilidad ---
                        # Guardar la lista original de tags en otra columna para análisis
                        openTargets_df["Tractability_raw"] = openTargets_df["Tractability"]
                        openTargets_df["Biotype_raw"] = openTargets_df["Biotype"]

                        # Explode para obtener cada tag en una fila
                        tract_tags = openTargets_df["Tractability_raw"].explode()
                        bio_tags = openTargets_df["Biotype_raw"].explode()

                        # Contar frecuencia de cada tipo de tractabilidad
                        tract_counts = tract_tags.value_counts().reset_index()
                        tract_counts.columns = ["Tractability", "Count"]
                        bio_counts = bio_tags.value_counts().reset_index()
                        bio_counts.columns = ["Biotype", "Count"]

                        # Crear gráfico de barras con plotly
                        fig = px.bar(
                            tract_counts,
                            x="Tractability",
                            y="Count",
                            color="Count",
                            color_continuous_scale="burg",
                            title="Count of genes by Tractability",
                            labels={"Tractability": "Type of tractability", "Count": "Number of genes"},
                            height=700,
                            width= 1000
                        )

                        st.plotly_chart(fig, use_container_width=True)
                        
                        fig = px.bar(
                            bio_counts,
                            x="Biotype",
                            y="Count",
                            color="Count",
                            color_continuous_scale="blues",
                            title="Count of genes by Biotype",
                            labels={"Biotype": "Biotype", "Count": "Number of genes"},
                            height=500,
                            width= 1000
                        )

                        st.plotly_chart(fig, use_container_width=True)


                        # Botón para descargar CSV
                        st.download_button(
                            "📥 Download results from Open Targets",
                            openTargets_df.to_csv(index=False),
                            "drug_target_genes.csv",
                            "text/csv"
                        )
                    else:
                        st.warning("⚠️ No results found in Open Targets.")

                with st.spinner("🧠 Performing pathway analysis..."):
                    number_pathways = st.text_input("🔍 Enter number of pathways to retrieve", value="10")
                    if number_pathways:
                        try:
                           number_pathways = int(number_pathways)
                           top_pathways, important_genes = analyze_pathways(df_selected, number_pathways)
                           top_pathways["-log10(Adj P)"] = -np.log10(top_pathways["Adjusted P-value"])
            
            
                           if top_pathways is not None:
                              st.subheader("📈 Top Enriched Reactome Pathways")
                              #st.write(top_pathways[["Reactome Link", "Adjusted P-value", "-log10(Adj P)"]].to_html(escape=False, index=False), unsafe_allow_html=True)


                              # Crear la tabla HTML bonita
                              html_pathway_table = top_pathways[["Reactome Link", "Adjusted P-value", "-log10(Adj P)"]].to_html(escape=False, index=False)

                              html_pathway_scroll = f"""
                                 <div style="max-height: 500px; overflow-y: auto;">
                                    {html_pathway_table}
                                 </div>
                              """
                              st.markdown(html_pathway_scroll, unsafe_allow_html=True)

                              selected_pathway = st.selectbox(
                                  "🔍 Select a pathway to see the top genes",
                                  top_pathways["Term"].tolist(),
                                  key="pathway"
                               )

                              pathway_genes = important_genes[important_genes["Gene Name"].isin(
                                  top_pathways[top_pathways["Term"] == selected_pathway]["Genes"].iloc[0].split(";")
                               )]
    
                              st.subheader(f"🔥 Important genes in pathway: {selected_pathway}")
                              pathway_genes["Gene"] = pathway_genes["Gene"].apply(
                                  lambda gene_id: f'<a href="https://www.ensembl.org/Multi/Search/Results?q={gene_id}" target="_blank">{gene_id}</a>'
                              )
                              #st.write(pathway_genes.to_html(escape=False, index=False), unsafe_allow_html=True)
                              html_genespathway_table = pathway_genes.to_html(escape=False, index=False)

                              html_genespathway_scroll = f"""
                                 <div style="max-height: 500px; overflow-y: auto;">
                                    {html_genespathway_table}
                                 </div>
                              """
                              st.markdown(html_genespathway_scroll, unsafe_allow_html=True)
                              st.download_button("📥 Download Important Genes CSV", pathway_genes.to_csv(index=False), "genes_pathway.csv", "text/csv")
                        except ValueError:
                           st.error("⚠️ Please enter a valid integer.")
                            
                    else:
                        st.warning("⚠️ No enriched pathways found.")

