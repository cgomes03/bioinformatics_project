from EA_rMTA.simplification.simplifier import SolutionSimplifier
from cobra.io import read_sbml_model
import pandas as pd
import numpy as np
from Preprocessing.preprocessing import DifferentialExpression
from MTAs.rmta import rMTA

run_folder = "/home/catarinagomes/catarina_stuff/EAs/result_merged"
model_path = "/home/catarinagomes/catarina_stuff/Recon3D_consistent.xml"
vref_path = "/home/catarinagomes/catarina_stuff/Sampler_ACHR_Recon3D_EOL.csv"
dge_path = "/home/catarinagomes/catarina_stuff/DGE_results_limma_voom_EOL.csv"
conv_dict_path = "/home/catarinagomes/catarina_stuff/dict_conv_final.txt"

model = read_sbml_model(model_path)
samples = pd.read_csv(vref_path, header=0).T
idx = np.zeros((samples.shape[0],), dtype=int)
for i in range(len(idx)):
    idx[i] = model.reactions.index(samples.index[i])
samples_temp = np.zeros((len(model.reactions), samples.shape[1]))
samples_temp[idx, :] = samples
Vref = pd.DataFrame(samples_temp.mean(axis=1))
Vref.index = [x.id for x in model.reactions]
diff_express_genes = pd.read_csv(dge_path, header=0)
conv_dict = pd.read_csv(conv_dict_path, sep='\t', header=None)
conv_dict = conv_dict.set_index(0).T.to_dict('records')[0]
diff_express_genes['gene_id'] = diff_express_genes['gene_id'].replace(conv_dict)
diff_express_genes = diff_express_genes[diff_express_genes['gene_id'].isin([x.id for x in model.genes])]
diff_expr = DifferentialExpression(model, diff_express_genes, pd.DataFrame(Vref))
diff_expr.get_differentially_expressed_genes(1, 0.05)
diff_expr.get_differentially_expressed_reactions_from_genes()
diff_expr.samples_filter(samples)
rMTA_instance = rMTA(model, diff_expr.rxnFBS, Vref, 0.66, samples_temp)

# simplificação só do top 10 já presente no final.csv
simplifier = SolutionSimplifier(rMTA_instance)
simplifier.simplify(folder=run_folder, top_n=10, max_workers=20, save=True)

print("Simplification plot and CSV generated.")
