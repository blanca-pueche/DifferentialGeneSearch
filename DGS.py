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

Entrez.email = "blanca.puechegranados@usp.ceu.es"
GRAPHQL_URL = "https://dgidb.org/api/graphql"
ENSEMBL_LOOKUP_URL = "https://rest.ensembl.org/lookup/id/"


st.markdown("""
    <style>
    .fixed-logo {
        position: fixed;
        top: 0;
        left: 0;
        padding: 5px;
        z-index: 9999;
        background: transparent;
    }
    html, body, .stApp {
        width: 100%;
        margin: 0;
        padding: 0;
        overflow-x: hidden;
        font-size: 18px;
    }
    .block-container {
        max-width: 80% !important;
        padding-left: 2rem;
        padding-right: 2rem;
    }
    .main {
        background: linear-gradient(to bottom right, #e3f2fd, #fce4ec);
        font-family: 'Arial', sans-serif;
        padding: 10px;
    }
    .stButton > button {
        background-color: #6a1b9a;
        color: white;
        font-weight: bold;
        border-radius: 2px;
    }
    .stDownloadButton > button {
        background-color: #2e7d32;
        color: white;
        font-weight: bold;
        border-radius: 2px;
    }
    </style>
""", unsafe_allow_html=True)

logo_placeholder = st.empty()
logo_placeholder.markdown(
    """
    <div class="fixed-logo"></div>
    """,
    unsafe_allow_html=True
)

# Pero la imagen la ponemos con st.image(), para que Streamlit la maneje y respete tamaño
st.image("https://upload.wikimedia.org/wikipedia/commons/7/7f/Logo_CNB.jpg", width=200)


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
    
import urllib.parse
import gseapy as gp

def analyze_pathways(df, number):
    # Prepare gene list
    df_genes = df["Gene Name"].dropna().astype(str).str.strip().str.upper().unique().tolist()

    # Enrichment
    enr = gp.enrichr(gene_list=df_genes, gene_sets="Reactome_2022", organism="Human", outdir=None)
    if enr.results.empty:
        return None, None

    top_pathways = enr.results.copy()

    # Add Reactome links
    top_pathways["Reactome Link"] = top_pathways["Term"].apply(
        lambda term: f'<a href="https://reactome.org/content/query?q={urllib.parse.quote(term)}" target="_blank">{term}</a>'
    )

    # Parse overlap info
    top_pathways[["Input Genes", "Pathway Genes"]] = top_pathways["Overlap"].str.split("/", expand=True).astype(int)
    top_pathways["Input %"] = top_pathways["Input Genes"] / top_pathways["Pathway Genes"] * 100
    top_pathways = top_pathways.sort_values("Input %", ascending=False).head(number)

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
    df_with_links["Drug"] = df_with_links["Drug"].apply(
        lambda drug: f'<a href="https://dgidb.org/results?searchType=drug&searchTerms={urllib.parse.quote(drug)}" target="_blank">{drug}</a>'
    )
    return df_with_links
    
    

st.title("Disease-Gene-Drug Analysis Dashboard")

import streamlit as st

# Estilo y layout del label con tooltip
st.markdown("""
<style>
.email-label-wrapper {
  display: flex;
  align-items: center;
  gap: 6px;
  font-weight: 500;
  font-size: 16px;
  margin-bottom: -10px;
}

.tooltip {
  position: relative;
  display: inline-block;
  font-size: 14px;
  vertical-align: middle;
  cursor: pointer;
  color: #1E90FF;
}

.tooltip .tooltiptext {
  visibility: hidden;
  width: 280px;
  background-color: #555;
  color: #fff;
  text-align: left;
  border-radius: 6px;
  padding: 6px;
  position: absolute;
  z-index: 1;
  bottom: 125%;  /* aparece encima del icono */
  left: 50%;
  margin-left: -140px;
  opacity: 0;
  transition: opacity 0.3s;
}

.tooltip .tooltiptext::after {
  content: "";
  position: absolute;
  top: 100%;
  left: 50%;
  margin-left: -5px;
  border-width: 5px;
  border-style: solid;
  border-color: #555 transparent transparent transparent;
}

.tooltip:hover .tooltiptext {
  visibility: visible;
  opacity: 1;
}
</style>

<div class="email-label-wrapper">
  <span>Enter user e-mail:</span>
  <div class="tooltip">ℹ
    <span class="tooltiptext">
      Se necesita el correo electrónico para hacer búsquedas en Entrez.
    </span>
  </div>
</div>
""", unsafe_allow_html=True)

# Campo de texto sin label visible, pero con accesibilidad
email = st.text_input("Correo electrónico", key="user_email", label_visibility="collapsed")


#Step 2: Name from MeshID
mesh_id = st.text_input("🔍 Enter MeSH ID (e.g., D003920 for Diabetes Mellitus):")

spinner = st.spinner

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
        st.write("Uploaded correctly")
        
        with st.spinner("Mapping Ensembl IDs to gene names..."):
            df_selected = fetch_gene_names(df_raw)
        
        # Step 4: Obtain Ensembl ID with links for the genes
        df_selected_with_links = df_selected.copy()
        df_selected_with_links["Gene"] = df_selected_with_links["Gene"].apply(
            lambda gene_id: f'<a href="https://www.ensembl.org/Multi/Search/Results?q={gene_id}" target="_blank">{gene_id}</a>'
        )
        df_selected_with_links = df_selected_with_links.sort_values(by="log_2 fold change", ascending=False)
        
        import streamlit.components.v1 as components

        # Título
        st.markdown("## Gene Table with Links to Ensembl")

        # Convertir tu dataframe a HTML (como antes)
        html_table = df_selected_with_links.to_html(escape=False, index=False, table_id="geneTable")

        # Código HTML con CSS fijo en blanco/negro
        html_code = f"""
        <link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
        <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
        <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>

        <style>
        /* Fuerza fondo blanco y texto negro */
        #geneTable, #geneTable thead, #geneTable tbody, #geneTable tr, #geneTable td, #geneTable th {{
            background-color: white !important;
            color: black !important;
            border-color: #ccc !important;
        }}

        #geneTable input {{
            background-color: white !important;
            color: black !important;
            border: 1px solid #ccc !important;
        }}

        .dataTables_wrapper {{
            color: black !important;
        }}

        div.dataTables_wrapper div.dataTables_paginate, 
        div.dataTables_wrapper div.dataTables_info {{
            color: black !important;
        }}
        </style>

        <script>
        $(document).ready(function() {{
            $('#geneTable thead tr').clone(true).appendTo('#geneTable thead');
            $('#geneTable thead tr:eq(1) th').each(function(i) {{
                var title = $(this).text();
                $(this).html('<input type="text" placeholder="Search ' + title + '" />');

                $('input', this).on('keyup change', function () {{
                    if (table.column(i).search() !== this.value) {{
                        table.column(i).search(this.value).draw();
                    }}
                }});
            }});

            var table = $('#geneTable').DataTable({{
                scrollY: '400px',
                scrollCollapse: true,
                paging: true,
                orderCellsTop: true,
                fixedHeader: true,
                searching: false
            }});
        }});
        </script>

        <div style="overflow-x:auto">
        {html_table}
        </div>
        """
        # Mostrar tabla con colores forzados
        components.html(html_code, height=650, scrolling=True)
        st.download_button(
            "📥 Download Genes CSV",
            df_selected_with_links.to_csv(index=False),
            "genes.csv",
            "text/csv"
        )
        
        df_selected["Gene Name_raw"] = df_selected["Gene Name"]

        gname_fc = df_selected.groupby("Gene Name_raw", as_index=False)["log_2 fold change"].sum()
        gname_fc = gname_fc.sort_values("log_2 fold change", ascending=False)

        # Graph for the genes and their log2 fold change
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

        
        # Step 5: Search biotype and tractability for the genes in Open Targets
        with st.spinner("Checking Open Targets..."):
            openTargets_results = df_selected["Gene"].apply(find_possible_target_of_drugs)
            openTargets_df = pd.DataFrame([r for r in openTargets_results if r is not None])

        if not openTargets_df.empty:
            if "Name" in openTargets_df.columns:
                openTargets_df = openTargets_df.drop(columns=["Name"])

                openTargets_df = openTargets_df.sort_values(by="Gene Symbol")

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


                st.markdown("# Genes with Tractability (Open Targets)")

                #scroll
                html_ot_table = openTargets_df.to_html(escape=False, index=False)

                html_ot_scroll = f"""
                    <div style="max-height: 500px; overflow-y: auto;"">
                        {html_ot_table}
                    </div>
                """

                st.markdown(html_ot_scroll, unsafe_allow_html=True)
                
                st.download_button(
                    "📥 Download results from Open Targets",
                    openTargets_df.to_csv(index=False),
                    "openTargets_genes.csv",
                    "text/csv"
                )

                openTargets_df["Tractability_raw"] = openTargets_df["Tractability"]
                openTargets_df["Biotype_raw"] = openTargets_df["Biotype"]

                tract_tags = openTargets_df["Tractability_raw"].explode()
                bio_tags = openTargets_df["Biotype_raw"].explode()

                tract_counts = tract_tags.value_counts().reset_index()
                tract_counts.columns = ["Tractability", "Count"]
                bio_counts = bio_tags.value_counts().reset_index()
                bio_counts.columns = ["Biotype", "Count"]

                # Tractability graph
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
                fig.update_xaxes(showticklabels=False)
                st.plotly_chart(fig, use_container_width=True)
                
                # Biotype graph
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

            else:
                st.warning("No results found in Open Targets.")

        
        with st.spinner("Performing pathway analysis..."):
            
            # Step 6: Search pathways for the genes
            st.markdown("# Top Enriched Reactome Pathways")
            number_pathways = st.text_input("🔍 Enter number of pathways to retrieve", value="10")
            if number_pathways:
                try:
                    number_pathways = int(number_pathways)
                    top_pathways = analyze_pathways(df_selected, number_pathways)
                    top_pathways["-log10(Adj P)"] = -np.log10(top_pathways["Adjusted P-value"])


                    if top_pathways is not None:
                        html_pathway_table = top_pathways[["Reactome Link", "Adjusted P-value", "-log10(Adj P)", "Overlap", "Input %", "Sum log2fc"]].to_html(escape=False, index=False)

                        html_pathway_scroll = f"""
                            <div style="max-height: 500px; overflow-y: auto;">
                            {html_pathway_table}
                            </div>
                        """
                        st.markdown(html_pathway_scroll, unsafe_allow_html=True)
                                    
                        st.markdown("# Important genes in pathway: ")
                        selected_pathway = st.selectbox(
                            "🔍 Select a pathway to see the top genes",
                            top_pathways["Term"].tolist(),
                            key="pathway"
                        )

                        # Get full row of selected pathway
                        selected_pathway_row = top_pathways[top_pathways["Term"] == selected_pathway].iloc[0]

                        # Get overlapping genes for that pathway
                        pathway_genes = get_overlapping_genes(df_selected, selected_pathway_row)

                        # Display pathway name
                        st.subheader(f"{selected_pathway}")

                        # Turn Gene IDs into Ensembl search links
                        pathway_genes["Gene Name"] = pathway_genes["Gene Name"].apply(
                            lambda gene_id: f'<a href="https://www.ensembl.org/Multi/Search/Results?q={gene_id}" target="_blank">{gene_id}</a>'
                        )

                        # Render HTML table with scroll
                        html_genespathway_table = pathway_genes.to_html(escape=False, index=False)

                        html_genespathway_scroll = f"""
                            <div style="max-height: 500px; overflow-y: auto;">
                            {html_genespathway_table}
                            </div>
                        """

                        st.markdown(html_genespathway_scroll, unsafe_allow_html=True)
                        st.download_button("📥 Download Important Genes CSV", pathway_genes.to_csv(index=False), "genes_pathway.csv", "text/csv")
                except ValueError:
                    st.error("Please enter a valid integer.")
                
                st.markdown("# Drug-Gene Interactions")
                with st.spinner("Searching drug-gene interactions..."):
                
                    selected_pathway = st.selectbox(
                        "🔍 Select a pathway to see drugs",
                    top_pathways["Term"].tolist(),
                    key="drugs"
                    )
                    selected_pathway_row = top_pathways[top_pathways["Term"] == selected_pathway].iloc[0]

                    pathway_genes = get_overlapping_genes(df_selected, selected_pathway_row)
                    print(pathway_genes["Gene Name"].to_list())
                    drug_df = get_drug_targets_dgidb_graphql(pathway_genes["Gene Name"].tolist())

                if not drug_df.empty:
                    
                    drug_df_with_links = drug_with_links(drug_df)
        
                    #st.write(drug_df_with_links.to_html(escape=False, index=False), unsafe_allow_html=True)
                    
                    # Crear la tabla HTML bonita
                    html_drug_table = drug_df_with_links.to_html(escape=False, index=False)

                    html_drug_scroll = f"""
                    <div style="max-height: 500px; overflow-y: auto;">
                        {html_drug_table}
                    </div>
                        """
                    st.markdown(html_drug_scroll, unsafe_allow_html=True)
        
                    st.download_button("📥 Download Drug Interactions CSV", drug_df.to_csv(index=False), "drug_interactions.csv", "text/csv")
                    
                    # Gene selection and plot — BELOW the table
                    unique_genes = drug_df_with_links['Gene'].unique()
                    selected_gene = st.selectbox("Select a gene to view its drug interaction scores", unique_genes)

                    # Filter by selected gene
                    gene_df = drug_df_with_links[drug_df['Gene'] == selected_gene]

                    # Plot
                    fig = px.bar(
                        gene_df,
                        x='Drug',
                        y='Interaction Score',
                        hover_data=['Interaction Type', 'Source', 'PMID'],
                        title=f"Drug Interactions for {selected_gene}",
                        labels={'Interaction Score': 'Interaction Score', 'Drug': 'Drug'},
                    )
                    fig.update_layout(xaxis_tickangle=-45)

                    st.plotly_chart(fig)
                    
                    if 'Interaction Type' in drug_df.columns:
                        st.markdown("# Distribution of Interaction Types")

                        # Count each type
                        interaction_counts = drug_df['Interaction Type'].value_counts().reset_index()
                        interaction_counts.columns = ['Interaction Type', 'Count']

                        # Pie chart
                        pie_fig = px.pie(
                            interaction_counts,
                            values='Count',
                            names='Interaction Type',
                            title='Percentage of Interaction Types',
                            hole=0.4
                        )

                        # Agrupar relaciones por tipo
                        interaction_grouped = (
                            drug_df.groupby('Interaction Type')
                            .apply(lambda df: [f"{g} → {d}" for g, d in zip(df['Gene'], df['Drug'])])
                        )

                        # Convertir a DataFrame con columnas por tipo de interacción
                        max_len = interaction_grouped.map(len).max()
                        interaction_wide = pd.DataFrame({
                            interaction_type: values + [""] * (max_len - len(values))  # rellenar con "" para igualar longitudes
                            for interaction_type, values in interaction_grouped.items()
                        })

                        # Layout: dos columnas
                        col1, col2 = st.columns([1.5, 1.8])

                        with col1:
                            st.plotly_chart(pie_fig)

                        with col2:
                            st.markdown("**Gene–Drug Interactions by Interaction Type**")
                            st.dataframe(interaction_wide, use_container_width=True)
                            
                        # -- Merge all data into a full report table --

                        # First, clean and prepare dataframes
                        merged = df_selected.copy()

                        # Merge with Open Targets
                        if not openTargets_df.empty:
                            ot_clean = openTargets_df.copy()
                            ot_clean["Tractability"] = ot_clean["Tractability"].str.replace(r'<.*?>', '', regex=True)
                            merged = pd.merge(merged, ot_clean.drop(columns=["Biotype_raw"]), how="left", left_on="Gene Name", right_on="Gene Symbol")

                        # Merge with pathway genes (optional: only if pathway is selected)
                        if not pathway_genes.empty:
                            merged = pd.merge(merged, pathway_genes[["Gene", "abs_fc"]], on="Gene", how="left")

                        # Merge with drug interactions
                        if not drug_df.empty:
                            drug_clean = drug_df.copy()
                            drug_summary = drug_clean.groupby("Gene").agg({
                                "Drug": lambda x: "; ".join(sorted(set(x))),
                                "Interaction Type": lambda x: "; ".join(sorted(set(x))),
                                "PMID": lambda x: "; ".join(sorted(set(str(p) for p in x if p != 'N/A'))),
                                "Interaction Score": "mean"
                            }).reset_index()
                            merged = pd.merge(merged, drug_summary, how="left", left_on="Gene Name", right_on="Gene")
                        
                        if not top_pathways.empty:
                            pathway_genes_unique = pathway_genes[["Gene Name"]].drop_duplicates()
                            pathway_genes_unique["Pathway Match"] = True
                            merged = pd.merge(merged, pathway_genes, how ="left", on="Gene Name")
                        
                        #Delete duplicated columns
                        merged = merged.drop(columns=["Gene_y", "Gene Symbol", "Tractability_raw", "Gene Name_raw_x", "Gene Name_raw_y", "abs_fc_y", "log_2 fold change_y"], errors="ignore")
                        merged = merged.rename(columns={"Gene_x": "Gene", "abs_fc_x" : "abs_fc"})
                        dupes = merged.columns[merged.columns.duplicated()].tolist()
                        merged = merged.loc[:, ~merged.columns.duplicated(keep="first")]

                        # Optional: reorder columns
                        #merged = merged.drop(columns=["Gene Symbol", "Gene_y", "Tractability_raw"], errors="ignore")
                        #merged = merged.rename(columns={"Gene_x": "Gene"})

                        # Show and offer download
                        st.markdown("## 📦 Download Full Results Table")
                        st.dataframe(merged, use_container_width=True)

                        st.download_button(
                            "📥 Download Full Combined CSV",
                            merged.to_csv(index=False),
                            "full_results_table.csv",
                            "text/csv"
                        )
                    else:
                        st.info("No 'Interaction Type' column found.")

                else:
                    st.warning("No drug-gene interactions found.")
                    
            else:
                st.warning("No enriched pathways found.")

