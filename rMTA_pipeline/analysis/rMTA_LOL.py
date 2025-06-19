import cobra
import pandas as pd
import scipy
import os
from Preprocessing.preprocessing import DifferentialExpression
model = cobra.io.read_sbml_model('Recon3D_consistent.xml')
model
samples = pd.read_csv('Sampler_ACHR_Recon3D_LOL.csv', header=0)
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


len(samples_temp)
Vref = samples_temp.mean(axis=1)
Vref
Vref = pd.DataFrame(Vref)
Vref.index = [x.id for x in model.reactions]
Vref


diff_express_genes = pd.read_csv('DGE_results_limma_voom_LOL.csv', header=0)
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
# Test
diff_expr = DifferentialExpression(model, diff_express_genes, pd.DataFrame(Vref))
diff_expr.get_differentially_expressed_genes(1, 0.05)
diff_expr.get_differentially_expressed_reactions_from_genes()
diff_expr.samples_filter(samples)
from MTAs.rmta import rMTA
rMTA = rMTA(model, diff_expr.rxnFBS, Vref, 0.66, samples_temp)
print(rMTA.evaluate_KO('6542_AT1'))
#rMTA.run()
#rMTA.TScores.to_csv('Recon3D_LOL.csv')