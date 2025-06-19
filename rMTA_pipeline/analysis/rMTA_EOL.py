import cobra
import pandas as pd
import scipy
import os
from Preprocessing.preprocessing import DifferentialExpression
model = cobra.io.read_sbml_model('Recon3D_consistent.xml')
model
samples = pd.read_csv('Sampler_ACHR_Recon3D_EOL.csv', header=0)
samples = samples.T
samples
import numpy as np
idx = np.zeros((samples.shape[0],), dtype=int)
len(idx)

for i in range(len(idx)):
    idx[i] = model.reactions.index(samples.index[i])
import numpy as np



rxnInactive = set(range(len(model.reactions))) - set(idx)  # inactive reactions
rxnInactive.__len__()




samples_temp = np.zeros((len(model.reactions), samples.shape[1]))
samples_temp
samples_temp[idx, :] = samples
#samples = samples.T
len(samples_temp)
Vref = samples_temp.mean(axis=1)
Vref
Vref = pd.DataFrame(Vref)
Vref.index = [x.id for x in model.reactions]
Vref
diff_express_genes = pd.read_csv('DGE_results_limma_voom_EOL.csv', header=0)
diff_express_genes
#read the dict tat is in txt format
conv_dict = pd.read_csv('dict_conv_final.txt', sep='\t', header=None)
# convert to dict
conv_dict = conv_dict.set_index(0).T.to_dict('records')[0]
conv_dict
# replace in the diff_express_genes the gene names with the corresponding ids
diff_express_genes['gene_id'] = diff_express_genes['gene_id'].replace(conv_dict)
diff_express_genes
# if the gene is not in the diff_express_genes remove it
diff_express_genes = diff_express_genes[diff_express_genes['gene_id'].isin([x.id for x in model.genes])]
diff_express_genes

# Test
Vref
len(model.reactions)
"""
    'balanced_stringency': {'logFC_required': 1.0,  'pval_required': 0.05},  # 2-fold, FDR≤5%   P 120718
    'increased_sensitivity': {'logFC_required': 0.58, 'pval_required': 0.10},  # 1.5-fold, FDR≤10%   P 121592
    'high_specificity': {'logFC_required': 2.0,  'pval_required': 0.05},    P 122007
"""
logfc_thresh = 2
pval_thresh = 0.05


# Test
diff_expr = DifferentialExpression(model, diff_express_genes, pd.DataFrame(Vref))
diff_expr.get_differentially_expressed_genes(logfc_thresh, pval_thresh)
diff_expr.get_differentially_expressed_reactions_from_genes()
diff_expr.samples_filter(samples)
from MTAs.rmta import rMTA
rMTA = rMTA(model, diff_expr.rxnFBS, Vref, 0.66, samples_temp)
rMTA.run()

rMTA.TScores.to_csv('Recon3D_EOL_delete_2.csv')

outfile = f"Recon3D_EOL_logFC{logfc_thresh}_pval{pval_thresh}.csv"
rMTA.TScores.to_csv(outfile)