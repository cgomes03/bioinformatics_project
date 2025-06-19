import pandas as pd
import numpy as np
from cobra.io import read_sbml_model
from EA_rMTA.utils.aux import get_reactions_from_geneID
from Preprocessing.preprocessing import DifferentialExpression
from MTAs.rmta import rMTA

# --- File paths (adjust if needed) ---
model_path = "/home/catarinagomes/catarina_stuff/Recon3D_consistent.xml"
vref_path = "/home/catarinagomes/catarina_stuff/Sampler_ACHR_Recon3D_EOL.csv"
dge_path = "/home/catarinagomes/catarina_stuff/DGE_results_limma_voom_EOL.csv"
conv_dict_path = "/home/catarinagomes/catarina_stuff/dict_conv_final.txt"
trace_path = "/home/catarinagomes/catarina_stuff/EAs/result_merged/simplification_trace_all.csv"

# --- Load model and samples ---
model = read_sbml_model(model_path)
samples = pd.read_csv(vref_path, header=0).T
idx = np.zeros((samples.shape[0],), dtype=int)
for i in range(len(idx)):
    idx[i] = model.reactions.index(samples.index[i])
samples_temp = np.zeros((len(model.reactions), samples.shape[1]))
samples_temp[idx, :] = samples
Vref = pd.DataFrame(samples_temp.mean(axis=1))
Vref.index = [x.id for x in model.reactions]

# --- Load differential expression and map gene IDs ---
dge = pd.read_csv(dge_path)
conv_dict = pd.read_csv(conv_dict_path, sep='\t', header=None).set_index(0).T.to_dict('records')[0]
dge['gene_id'] = dge['gene_id'].replace(conv_dict)
dge = dge[dge['gene_id'].isin([x.id for x in model.genes])]

# --- Create DifferentialExpression and rMTA instance ---
diff_expr = DifferentialExpression(model, dge, pd.DataFrame(Vref))
diff_expr.get_differentially_expressed_genes(1, 0.05)
diff_expr.get_differentially_expressed_reactions_from_genes()
diff_expr.samples_filter(samples)
rMTA_instance = rMTA(model, diff_expr.rxnFBS, Vref, 0.66, samples_temp)

# --- Load simplification trace and extract core genes ---
df = pd.read_csv(trace_path)
best_id = df[df["Step"] == 0].sort_values(by="Fitness", ascending=False).iloc[0]["Individual"]
core_genes_line = df[(df["Individual"] == best_id) & (df["Genes Count"] == 5)].iloc[0]
core_genes = core_genes_line["Remaining Genes"].split(";")

# --- Evaluate each gene independently ---
results = []
for gene in core_genes:
    try:
        reactions = get_reactions_from_geneID([gene], rMTA_instance.model)
        if not reactions:
            raise ValueError("No reactions found for this gene.")
        value, _, _, _ = rMTA_instance.evaluate_KO(reactions)
        score = float(value.rTS.values[0])
    except Exception as e:
        score = f"error: {e}"
    results.append({"Gene": gene, "Single-gene rTS": score})

# --- Build and save Table  ---
table = pd.DataFrame(results)
print("\nTable  — Individual rTS for Core Genes:")
print(table)
table.to_csv("result_merged\table_individual_rts.csv", index=False)
