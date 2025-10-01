import streamlit as st
import pandas as pd
from Bio import Entrez
import plotly.express as px
import os
import re
import io, zipfile, datetime
import re
import streamlit.components.v1 as components
from pipeline import (estimate_table_height, normalize_disease_name, save_pathway_csvs,  save_drug_csvs)
import time
from pathlib import Path


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

# Image set with st.image() 
st.image("https://upload.wikimedia.org/wikipedia/commons/7/7f/Logo_CNB.jpg", width=200)


st.title("Tractome")
st.markdown(
    "<p style='font-size:18px; font-weight:bold;'>Integrative analysis of upregulated disease genes, pathways, and drug interactions.</p>",
    unsafe_allow_html=True
)


# Style and layout (tooltip)
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
  <div class="tooltip">ℹ️
    <span class="tooltiptext">
      Email adress is necessary for Entrez searches.
    </span>
  </div>
</div>
""", unsafe_allow_html=True)

#Step 1: Text input email empty
email = st.text_input("User e-mail:", key="user_email", label_visibility="collapsed", value="")
Entrez.email = email

#Step 2: Name from MeshID default D003110
mesh_id = st.text_input("🔍 Enter MeSH ID (e.g., D003920 for Diabetes Mellitus):", value="D003110")

spinner = st.spinner

if mesh_id:
    with st.spinner("Fetching disease information..."):
        disease = "Colonic Neoplasms"
        disease_url = "https://www.ncbi.nlm.nih.gov/mesh/68003110"
        
        if disease:
            normalized_disease = normalize_disease_name(disease)
            st.success(f"🎯 Disease **[{normalized_disease}]({disease_url})** identified")
            
            # Step 3: File upload only after disease name is given default link to Expression Atlas
            url = "https://www.ebi.ac.uk/gxa/search?geneQuery=%5B%5D&species=Homo%20sapiens&conditionQuery=%5B%7B%22value%22%3A%22Colonic%20Neoplasms%22%7D%5D&ds=%7B%22kingdom%22%3A%5B%22animals%22%5D%2C%22regulation%22%3A%5B%22UP%22%5D%7D&bs=%7B%22homo%20sapiens%22%3A%5B%22ORGANISM_PART%22%5D%7D#differential"
            uploaded_file = st.file_uploader(
                f"📁 Upload Differential Expression File (TSV from Expression Atlas: [link]({url}))",
                type=["tsv"]
            )

            if uploaded_file is None:
                #demo Expression Atlas file
                uploaded_file = "../demoData/colorectal.tsv"

            
            
            
    if uploaded_file:
        df_raw = pd.read_csv(uploaded_file, sep="\t")
        st.write("Uploaded correctly")
        
        with st.spinner("Mapping Ensembl IDs to gene names..."):
            time.sleep(1)
        
        # Step 4: Obtain Ensembl ID with links for the genes. Default csv file
        st.markdown("## Gene Table with Links to Ensembl")

        # Conversion of dataframe to HTML
        df_genes = pd.read_csv("../demoData/genes.csv", sep=",")
        html_table = df_genes.to_html(
            escape=False, index=False, table_id="geneTable"
        )

        # HTML code with CSS fixed in black and white
        html_code = f"""
        <link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
        <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
        <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>

        <style>
        /* Apply style to whole table */
        .dataTables_wrapper, .dataTables_wrapper * {{
            font-family: "Segoe UI", "Helvetica", "Arial", sans-serif !important;
            color: #31333f !important;
            font-size: 15px !important;
        }}

        /* Table */
        #geneTable, #geneTable thead, #geneTable tbody, #geneTable tr, #geneTable td, #geneTable th {{
            background-color: white !important;
            color: #31333f !important;
            border-color: #ccc !important;
        }}

        /* Headers */
        #geneTable thead th {{
            font-weight: 600 !important;
        }}

        /* Search inputs */
        #geneTable input {{
            background-color: white !important;
            color: #31333f !important;
            font-size: 14px !important;
            border: 1px solid #ccc !important;
        }}

        /* Entry selections */
        .dataTables_length select {{
            background-color: white !important;
            color: #31333f !important;
            border: 1px solid #ccc !important;
            font-size: 14px !important;
        }}

        /* Results info */
        .dataTables_info {{
            color: #31333f !important;
        }}

        /* Pages */
        .dataTables_paginate a {{
            color: #31333f !important;
            font-size: 14px !important;
            font-weight: normal !important;
            padding: 6px 12px;
            border-radius: 6px;
            text-decoration: none;
            margin: 0 2px;
        }}

        .dataTables_paginate a.current {{
            background-color: #f0f0f0 !important;
            font-weight: bold !important;
            border: 1px solid #aaa !important;
        }}

        .dataTables_paginate a:hover {{
            background-color: #e3e3e3 !important;
        }}
        </style>

        <script>
        $(document).ready(function() {{
            var table = $('#geneTable').DataTable({{
                scrollCollapse: true,
                paging: true,
                orderCellsTop: true,
                fixedHeader: true
            }});

            // Add search inputs to each column
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
        }});
        </script>

        <div style="overflow-x:auto; margin-bottom: 0px;"">
        {html_table}
        </div>
        """

        height = estimate_table_height(df_genes)
        
        # Show table
        components.html(html_code, height=height, scrolling=True)
        
        # Download button
        st.download_button(
            "📥 Download Genes CSV",
            df_genes.to_csv(index=False),
            "genes.csv",
            "text/csv")
        
        df_genes["Gene Name_raw"] = df_genes["Gene Name"]

        gname_fc = df_genes.groupby("Gene Name_raw", as_index=False)["log_2 fold change"].sum()
        gname_fc = gname_fc.sort_values("log_2 fold change", ascending=False)

        # Graph for genes and their log2 fold change
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
            time.sleep(1)
        #default file
        openTargets_df = pd.read_csv("../demoData/openTargets_genes.csv", sep=",")
        if True:
            if True:
                st.markdown("# Genes with Tractability (Open Targets)")
                #helpful link to the overview page
                st.markdown(
                    """
                    ℹ️ For details on how tractability is defined, see the 
                    [Open Targets Tractability Overview](https://platform-docs.opentargets.org/target/tractability).
                    """
                )
                html_ot_table = openTargets_df.to_html(escape=False, index=False, table_id="openTargetsTable")

                html_ot_scroll = f"""
                    <link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
                    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
                    <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>

                    <style>
                    /* Apply style to whole table */
                    .dataTables_wrapper, .dataTables_wrapper * {{
                        font-family: "Segoe UI", "Helvetica", "Arial", sans-serif !important;
                        color: #31333f !important;
                        font-size: 15px !important;
                    }}

                    /* Table */
                    #openTargetsTable, #openTargetsTable thead, #openTargetsTable tbody, #openTargetsTable tr, #openTargetsTable td, #openTargetsTable th {{
                        background-color: white !important;
                        color: #31333f !important;
                        border-color: #ccc !important;
                    }}

                    /* Headers */
                    #openTargetsTable thead th {{
                        font-weight: 600 !important;
                    }}

                    /* Search inputs */
                    #openTargetsTable input {{
                        background-color: white !important;
                        color: #31333f !important;
                        font-size: 14px !important;
                        border: 1px solid #ccc !important;
                    }}

                    /* Entry selections */
                    .dataTables_length select {{
                        background-color: white !important;
                        color: #31333f !important;
                        border: 1px solid #ccc !important;
                        font-size: 14px !important;
                    }}

                    /* Results info */
                    .dataTables_info {{
                        color: #31333f !important;
                    }}

                    /* Pages */
                    .dataTables_paginate a {{
                        color: #31333f !important;
                        font-size: 14px !important;
                        font-weight: normal !important;
                        padding: 6px 12px;
                        border-radius: 6px;
                        text-decoration: none;
                        margin: 0 2px;
                    }}

                    .dataTables_paginate a.current {{
                        background-color: #f0f0f0 !important;
                        font-weight: bold !important;
                        border: 1px solid #aaa !important;
                    }}

                    .dataTables_paginate a:hover {{
                        background-color: #e3e3e3 !important;
                    }}
                    </style>

                    <script>
                    $(document).ready(function() {{
                        var table = $('#openTargetsTable').DataTable({{
                            scrollCollapse: true,
                            paging: true,
                            orderCellsTop: true,
                            fixedHeader: true
                        }});

                        // Add search inputs to each column
                        $('#openTargetsTable thead tr').clone(true).appendTo('#openTargetsTable thead');
                        $('#openTargetsTable thead tr:eq(1) th').each(function(i) {{
                            var title = $(this).text();
                            $(this).html('<input type="text" placeholder="Search ' + title + '" />');

                            $('input', this).on('keyup change', function () {{
                                if (table.column(i).search() !== this.value) {{
                                    table.column(i).search(this.value).draw();
                                }}
                            }});
                        }});
                    }});
                    </script>
                    
                    <div style="overflow-x:auto">
                        {html_ot_table}
                    </div>
                """
                
                height = estimate_table_height(openTargets_df)
                components.html(html_ot_scroll, height=height+200, scrolling=True)
                
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

        
        with st.spinner("Performing pathway analysis..."):
            
            # Step 6: Search pathways for the genes
            st.markdown("# Top Enriched Reactome Pathways")
            # Add a helpful link to the overview page
            st.markdown(
                """
                ℹ️ For details on how the columns are defined, see the 
                [Enrichr Help Center](https://maayanlab.cloud/Enrichr/help#background).
                """
            )
            st.markdown(
                """
                ℹ️ **Custom column in this table**  
                - **Input %**: Percentage of the pathway’s genes present in your input list. 
                """
            )

            number_pathways = st.text_input("🔍 Enter number of pathways to retrieve", value="10")
            if number_pathways:
                try:
                    number_pathways = int(number_pathways)
                    #default file
                    top_pathways = pd.read_csv("../demoData/topPathways.csv", sep=",")

                    if top_pathways is not None:
                        html_pathway_table = top_pathways[["Reactome Link", "Adjusted P-value", "-log10(Adj P)", "Overlap", "Input %", "Sum log2fc"]].to_html(escape=False, index=False, table_id="topPathwayTable")
                        
                        html_code_top_pathways = f"""
                        <link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
                        <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
                        <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>

                        <style>
                            body, .dataTables_wrapper, .dataTables_wrapper * {{
                                font-family: "Segoe UI", "Helvetica", "Arial", sans-serif !important;
                                color: #31333f !important;
                                font-size: 15px !important;
                            }}

                            #topPathwayTable {{
                                width: 100% !important;
                            }}

                            #topPathwayTable, #topPathwayTable thead, #topPathwayTable tbody, #topPathwayTable tr, #topPathwayTable td, #topPathwayTable th {{
                                background-color: white !important;
                                color: #31333f !important;
                                border-color: #ccc !important;
                            }}

                            #topPathwayTable thead th {{
                                font-weight: 600 !important;
                            }}

                            #topPathwayTable input {{
                                width: 100%;
                                box-sizing: border-box;
                                background-color: white !important;
                                color: #31333f !important;
                                font-size: 14px !important;
                                border: 1px solid #ccc !important;
                            }}

                            .dataTables_length select {{
                                background-color: white !important;
                                color: #31333f !important;
                                border: 1px solid #ccc !important;
                                font-size: 14px !important;
                            }}

                            .dataTables_info {{
                                color: #31333f !important;
                            }}

                            .dataTables_paginate a {{
                                color: #31333f !important;
                                font-size: 14px !important;
                                font-weight: normal !important;
                                padding: 6px 12px;
                                border-radius: 6px;
                                text-decoration: none;
                                margin: 0 2px;
                            }}

                            .dataTables_paginate a.current {{
                                background-color: #f0f0f0 !important;
                                font-weight: bold !important;
                                border: 1px solid #aaa !important;
                            }}

                            .dataTables_paginate a:hover {{
                                background-color: #e3e3e3 !important;
                            }}
                        </style>

                        <script>
                            $(document).ready(function() {{
                                var table = $('#topPathwayTable').DataTable({{
                                    scrollCollapse: true,
                                    paging: true,
                                    orderCellsTop: true,
                                    fixedHeader: false,
                                    autoWidth: false
                                }});

                                $('#topPathwayTable thead tr').clone(true).appendTo('#topPathwayTable thead');
                                $('#topPathwayTable thead tr:eq(1) th').each(function(i) {{
                                    var title = $(this).text();
                                    $(this).html('<input type="text" placeholder="Search ' + title + '" />');
                                }});

                                $('#topPathwayTable thead').on('keyup change', 'input', function () {{
                                    let i = $(this).parent().index();
                                    table.column(i).search(this.value).draw();
                                }});
                            }});
                        </script>

                        <div style="overflow-x:auto">
                            {html_pathway_table}
                        </div>
                        """
                        height = estimate_table_height(top_pathways)
                        components.html(html_code_top_pathways, height=height, scrolling=True)
                        
                        st.download_button(
                            "📥 Download Top N pathways CSV",
                            top_pathways.to_csv(index=False),
                            "topPathways.csv",
                            "text/csv")
                                    
                        st.markdown("# Important genes in pathway: ")
                        selected_pathway = st.selectbox(
                            "🔍 Select a pathway to see the top genes",
                            top_pathways["Term"].tolist(),
                            key="pathway"
                        )

                        # Get full row of selected pathway
                        selected_pathway_row = top_pathways[top_pathways["Term"] == selected_pathway].iloc[0]

                        # Get overlapping genes for that pathway
                        #default file path for csv files of top 10 pathways
                        folder_path = Path("../demoData/all_pathway_genes_csvs") 

                        # Display pathway name
                        st.subheader(f"{selected_pathway}")

                        file_name = f"{selected_pathway.replace(' ', '_')}.csv"
                        matching_file = folder_path / file_name
                        pathway_genes_newNames = pd.read_csv(matching_file, sep=",")
                        
                        # Render HTML table with scroll
                        html_genespathway_table = pathway_genes_newNames.to_html(escape=False, index=False, table_id="pathwayGeneTable")


                        # HTML and DataTables JS
                        html_code_pathway = f"""
                        <link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
                        <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
                        <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>

                        <style>
                            .dataTables_wrapper, .dataTables_wrapper * {{
                                font-family: "Segoe UI", "Helvetica", "Arial", sans-serif !important;
                                color: #31333f !important;
                                font-size: 15px !important;
                            }}

                            #pathwayGeneTable {{
                                width: 100% !important;
                            }}

                            #pathwayGeneTable, #pathwayGeneTable thead, #pathwayGeneTable tbody, #pathwayGeneTable tr, #pathwayGeneTable td, #pathwayGeneTable th {{
                                background-color: white !important;
                                color: #31333f !important;
                                border-color: #ccc !important;
                            }}

                            #pathwayGeneTable thead th {{
                                font-weight: 600 !important;
                            }}

                            #pathwayGeneTable input {{
                                width: 100%;
                                box-sizing: border-box;
                                background-color: white !important;
                                color: #31333f !important;
                                font-size: 14px !important;
                                border: 1px solid #ccc !important;
                            }}

                            .dataTables_length select {{
                                background-color: white !important;
                                color: #31333f !important;
                                border: 1px solid #ccc !important;
                                font-size: 14px !important;
                            }}

                            .dataTables_info {{
                                color: #31333f !important;
                            }}

                            .dataTables_paginate a {{
                                color: #31333f !important;
                                font-size: 14px !important;
                                font-weight: normal !important;
                                padding: 6px 12px;
                                border-radius: 6px;
                                text-decoration: none;
                                margin: 0 2px;
                            }}

                            .dataTables_paginate a.current {{
                                background-color: #f0f0f0 !important;
                                font-weight: bold !important;
                                border: 1px solid #aaa !important;
                            }}

                            .dataTables_paginate a:hover {{
                                background-color: #e3e3e3 !important;
                            }}
                        </style>

                        <script>
                            $(document).ready(function() {{
                                var table = $('#pathwayGeneTable').DataTable({{
                                    scrollCollapse: true,
                                    paging: true,
                                    orderCellsTop: true,
                                    fixedHeader: false,
                                    autoWidth: false
                                }});

                                $('#pathwayGeneTable thead tr').clone(true).appendTo('#pathwayGeneTable thead');
                                $('#pathwayGeneTable thead tr:eq(1) th').each(function(i) {{
                                    var title = $(this).text();
                                    $(this).html('<input type="text" placeholder="Search ' + title + '" />');
                                }});

                                $('#pathwayGeneTable thead').on('keyup change', 'input', function () {{
                                    let i = $(this).parent().index();
                                    table.column(i).search(this.value).draw();
                                }});
                            }});
                        </script>

                        <div style="overflow-x:auto">
                            {html_genespathway_table}
                        </div>
                        """

                        height = estimate_table_height(pathway_genes_newNames)
                        components.html(html_code_pathway, height=height, scrolling=True)

                        st.download_button(
                            "📥 Download Important Genes CSV",
                            pathway_genes_newNames.to_csv(index=False),
                            "genes_pathway.csv",
                            "text/csv"
                        )

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

                    #default path for csv files of top 10 pathways
                    folder_path = Path("../demoData/all_drug_csvs") 
                    file_name = f"{selected_pathway.replace(' ', '_')}_drugs.csv"
                    matching_file = folder_path / file_name
                    drug_df = pd.read_csv(matching_file, sep=",")

                if not drug_df.empty:                    
                    # HTML table
                    html_drug_table = drug_df.to_html(escape=False, index=False, table_id="drugTable")

                    html_code_drug = f"""
                    <link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
                    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
                    <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>

                    <style>
                        .dataTables_wrapper, .dataTables_wrapper * {{
                            font-family: "Segoe UI", "Helvetica", "Arial", sans-serif !important;
                            color: #31333f !important;
                            font-size: 15px !important;
                        }}

                        #drugTable {{
                            width: 100% !important;
                        }}

                        #drugTable, #drugTable thead, #drugTable tbody, #drugTable tr, #drugTable td, #drugTable th {{
                            background-color: white !important;
                            color: #31333f !important;
                            border-color: #ccc !important;
                        }}

                        #drugTable thead th {{
                            font-weight: 600 !important;
                        }}

                        #drugTable input {{
                            width: 100%;
                            box-sizing: border-box;
                            background-color: white !important;
                            color: #31333f !important;
                            font-size: 14px !important;
                            border: 1px solid #ccc !important;
                        }}

                        .dataTables_length select {{
                            background-color: white !important;
                            color: #31333f !important;
                            border: 1px solid #ccc !important;
                            font-size: 14px !important;
                        }}

                        .dataTables_info {{
                            color: #31333f !important;
                        }}

                        .dataTables_paginate a {{
                            color: #31333f !important;
                            font-size: 14px !important;
                            font-weight: normal !important;
                            padding: 6px 12px;
                            border-radius: 6px;
                            text-decoration: none;
                            margin: 0 2px;
                        }}

                        .dataTables_paginate a.current {{
                            background-color: #f0f0f0 !important;
                            font-weight: bold !important;
                            border: 1px solid #aaa !important;
                        }}

                        .dataTables_paginate a:hover {{
                            background-color: #e3e3e3 !important;
                        }}
                    </style>

                    <script>
                        $(document).ready(function() {{
                            var table = $('#drugTable').DataTable({{
                                scrollCollapse: true,
                                paging: true,
                                orderCellsTop: true,
                                fixedHeader: false,
                                autoWidth: false
                            }});

                            $('#drugTable thead tr').clone(true).appendTo('#drugTable thead');
                            $('#drugTable thead tr:eq(1) th').each(function(i) {{
                                var title = $(this).text();
                                $(this).html('<input type="text" placeholder="Search ' + title + '" />');
                            }});

                            $('#drugTable thead').on('keyup change', 'input', function () {{
                                let i = $(this).parent().index();
                                table.column(i).search(this.value).draw();
                            }});
                        }});
                    </script>

                    <div style="overflow-x:auto">
                        {html_drug_table}
                    </div>
                    """
                    height = estimate_table_height(drug_df)
                    components.html(html_code_drug, height=height, scrolling=True)                    
        
                    st.download_button("📥 Download Drug Interactions CSV", drug_df.to_csv(index=False), "drug_interactions.csv", "text/csv")
                    
                    # Extract plain gene names into a new column
                    drug_df['Gene_name'] = drug_df['Gene'].apply(lambda x: re.search(r'>(.*?)<', x).group(1))

                    # Gene selection
                    unique_genes = drug_df['Gene_name'].unique()
                    selected_gene = st.selectbox("Select a gene to view its drug interaction scores", unique_genes)

                    # Filter using the plain gene name column
                    gene_df = drug_df[drug_df['Gene_name'] == selected_gene]
                    gene_df = gene_df.sort_values(by="Interaction Score", ascending=False)

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

                        drug_df['Interaction Type'] = drug_df['Interaction Type'].fillna("N/A")
                        drug_df['Drug_name_only'] = drug_df['Drug'].apply(
                            lambda x: re.search(r'>(.*?)<', x).group(1) if pd.notnull(x) else ""
                        )

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

                        # Group by type of interaction
                        interaction_grouped = (
                            drug_df.groupby('Interaction Type', group_keys=False)
                            .apply(lambda df: [f"{g} → {d}" for g, d in zip(df['Gene_name'], df['Drug_name_only'])], include_groups=False)
                        )

                        # Convert to DataFrame with columns according to interaction type
                        max_len = interaction_grouped.map(len).max()

                        interaction_wide = pd.DataFrame({
                            interaction_type: list(values) + [""] * (max_len - len(values))
                            for interaction_type, values in interaction_grouped.items()
                        })

                        # Layout: 2 columns
                        col1, col2 = st.columns([1.5, 1.8])

                        with col1:
                            st.plotly_chart(pie_fig)

                        with col2:
                            st.markdown("**Gene–Drug Interactions by Interaction Type**")
                            st.dataframe(interaction_wide, use_container_width=True)              
                        
                    else:
                        st.info("No 'Interaction Type' column found.")

                else:
                    st.warning("No drug-gene interactions found.")
                    
            else:
                st.warning("No enriched pathways found.")
        # -- Merge all data into a full report table --
        with st.spinner("Creating summary table..."):
            if True:
                
                # Apply to your merged table
                #default file
                merged_with_links = pd.read_csv("../demoData/full_results_table.csv", sep=",")
                merged_with_links = merged_with_links.fillna("NaN")

                st.markdown("## 📦 Download Full Results Table")
                
                html_final_table = merged_with_links.to_html(escape=False, index=False, table_id="finalTable")

                html_code_final = f"""
                    <link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
                    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
                    <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>

                    <style>
                        .dataTables_wrapper, .dataTables_wrapper * {{
                            font-family: "Segoe UI", "Helvetica", "Arial", sans-serif !important;
                            color: #31333f !important;
                            font-size: 12px !important;
                        }}

                        #finalTable {{
                            width: 100% !important;
                        }}

                        #finalTable, #finalTable thead, #finalTable tbody, #finalTable tr, #finalTable td, #finalTable th {{
                            background-color: white !important;
                            color: #31333f !important;
                            border-color: #ccc !important;
                        }}

                        #finalTable thead th {{
                            font-weight: 600 !important;
                        }}

                        #finalTable input {{
                            width: 100%;
                            box-sizing: border-box;
                            background-color: white !important;
                            color: #31333f !important;
                            font-size: 14px !important;
                            border: 1px solid #ccc !important;
                        }}

                        .dataTables_length select {{
                            background-color: white !important;
                            color: #31333f !important;
                            border: 1px solid #ccc !important;
                            font-size: 14px !important;
                        }}

                        .dataTables_info {{
                            color: #31333f !important;
                        }}

                        .dataTables_paginate a {{
                            color: #31333f !important;
                            font-size: 14px !important;
                            font-weight: normal !important;
                            padding: 6px 12px;
                            border-radius: 6px;
                            text-decoration: none;
                            margin: 0 2px;
                        }}

                        .dataTables_paginate a.current {{
                            background-color: #f0f0f0 !important;
                            font-weight: bold !important;
                            border: 1px solid #aaa !important;
                        }}

                        .dataTables_paginate a:hover {{
                            background-color: #e3e3e3 !important;
                        }}
                    </style>

                    <script>
                        $(document).ready(function() {{
                            var table = $('#finalTable').DataTable({{
                                scrollCollapse: true,
                                paging: true,
                                orderCellsTop: true,
                                fixedHeader: false,
                                autoWidth: false
                            }});

                            $('#finalTable thead tr').clone(true).appendTo('#finalTable thead');
                            $('#finalTable thead tr:eq(1) th').each(function(i) {{
                                var title = $(this).text();
                                $(this).html('<input type="text" placeholder="Search ' + title + '" />');
                            }});

                            $('#finalTable thead').on('keyup change', 'input', function () {{
                                let i = $(this).parent().index();
                                table.column(i).search(this.value).draw();
                            }});
                        }});
                    </script>

                    <div style="overflow-x:auto">
                        {html_final_table}
                    </div>
                    """
                height = estimate_table_height(merged_with_links)
                components.html(html_code_final, height=height, scrolling=True)

                st.download_button(
                        "📥 Download Full Combined CSV",
                        merged_with_links.to_csv(index=False),
                        "full_results_table.csv",
                        "text/csv"
                )
                # --- ZIP with chosen tables (as they are) ---
                # Build the list of available dataframes (desc, filename, df)
                candidates = []
                try: candidates.append(("Genes with links",     "genes_with_links.csv", df_genes))
                except NameError: pass
                try: candidates.append(("Open Targets",         "openTargets_df.csv",   openTargets_df))
                except NameError: pass
                try: candidates.append(("Top pathways",         "top_pathways.csv",     top_pathways))
                except NameError: pass
                try: candidates.append(("Pathway genes",        "pathway_genes.csv",    pathway_genes_newNames))
                except NameError: pass
                try: candidates.append(("Drug interactions",    "drug_df.csv",          drug_df))
                except NameError: pass

                if candidates:
                    st.markdown("### Select tables to include in ZIP")

                    # Table headers
                    h1, h2 = st.columns([0.12, 0.88])
                    h1.markdown("** **")
                    h2.markdown("**Dataframe**")

                    # Individual checkboxes
                    for desc, fname, df in candidates:
                        c1, c2 = st.columns([0.12, 0.88])
                        c1.checkbox(
                            label=f"Include {desc}",
                            key=f"_dl_{fname}",
                            value=True,  # selected by default
                            label_visibility="collapsed"
                        )
                        c2.write(f"{desc} ({len(df)} rows)")

                    # Gather selected
                    selected = [(fname, df) for _, fname, df in candidates if st.session_state.get(f"_dl_{fname}", False)]

                # Create columns for the two download buttons
                col1, col2 = st.columns(2)

                # --- Column 1: Download selected tables ---
                with col1:
                    if not selected:
                        st.warning("Select at least one table to enable the download.")
                    else:
                        zip_buf = io.BytesIO()
                        with zipfile.ZipFile(zip_buf, "w", zipfile.ZIP_DEFLATED) as zf:
                            for fname, df in selected:
                                zf.writestr(fname, df.to_csv(index=False))
                        zip_buf.seek(0)

                        st.download_button(
                            label="📥 Download selected (.zip)",
                            data=zip_buf.getvalue(),
                            file_name=f"selected_tables_{datetime.date.today().isoformat()}.zip",
                            mime="application/zip",
                            key="zip_download_btn"
                        )

                # --- Column 2: Download all tables ---
                with col2:
                    # Define folders containing csvs
                    genes_folder = "all_pathway_genes_csvs"
                    drugs_folder = "all_drug_csvs"

                    # Optional extra CSVs
                    extra_files = [
                        "genes.csv",
                        "openTargets_genes.csv",
                        "topPathways.csv",
                        "genes_pathway.csv",
                        "drug_interactions.csv",
                        "full_results_table.csv"
                    ]

                    def count_csvs(folder):
                        return len([f for f in os.listdir(folder) if f.endswith(".csv")])

                    st.markdown(
                        "⚠️ **Note:** Clicking the button below will start the process of generating and packaging all tables "
                        "(it will generate all the tables for each pathway). This may take a few minutes depending on your system and file sizes."
                    )

                    if st.button("Generate and Download Tables"):
                        with st.spinner("Generating all CSVs..."):
                            drug_csvs = save_drug_csvs(df_genes, top_pathways)
                            pathway_csvs = save_pathway_csvs(df_genes, top_pathways)

                            # Build ZIP in memory
                            zip_buf = io.BytesIO()
                            with zipfile.ZipFile(zip_buf, "w", zipfile.ZIP_DEFLATED) as zf:
                                # Add pathway CSVs
                                for fname, data in pathway_csvs.items():
                                    zf.writestr(fname, data)
                                # Add drug CSVs
                                for fname, data in drug_csvs.items():
                                    zf.writestr(fname, data)
                                # Add any extra files
                                for f in extra_files:
                                    if os.path.exists(f):
                                        zf.write(f, arcname=os.path.basename(f))

                            zip_buf.seek(0)

                        st.download_button(
                            label="📥 Download All Tables",
                            data=zip_buf.getvalue(),
                            file_name="all_tables.zip",
                            mime="application/zip"
                        )