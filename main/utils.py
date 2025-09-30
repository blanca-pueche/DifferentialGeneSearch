import os, json, pandas as pd

def load_precomputed(demo_root="demoData"):
    base = os.path.join(demo_root)
    results = {}

    results["genes"] = pd.read_csv(os.path.join(base, "genes.csv"))
    results["top_pathways"] = pd.read_csv(os.path.join(base, "topPathways.csv"))
    results["drug_interactions"] = pd.read_csv(os.path.join(base, "drug_interactions.csv"))
    results["openTargets"] = pd.read_csv(os.path.join(base, "openTargets.csv"))
    results["full_table"] = pd.read_csv(os.path.join(base, "full_results_table.csv"))

    if os.path.exists(os.path.join(base, "metadata.json")):
        with open(os.path.join(base, "metadata.json")) as f:
            results["metadata"] = json.load(f)

    if os.path.exists(os.path.join(base, "pathway_csvs")):
        results["pathway_csvs_folder"] = os.path.join(base, "pathway_csvs")

    if os.path.exists(os.path.join(base, "drug_csvs")):
        results["drug_csvs_folder"] = os.path.join(base, "drug_csvs")

    return results
