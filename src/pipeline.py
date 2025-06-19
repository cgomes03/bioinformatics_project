import os
from EA_rMTA.ea.engine import EAEngine
from EA_rMTA.plotting.plotter import ResultPlotter
from EA_rMTA.simplification.simplifier import SolutionSimplifier


class Pipeline:
    """
    Orchestrates the EA run, plotting, and solution simplification.
    """

    def __init__(self, rMTA_obj, output_folder: str):
        self.rMTA_obj = rMTA_obj
        self.output_folder = output_folder

    def run_ea(self, **ea_params):
        engine = EAEngine(self.rMTA_obj, **ea_params)
        run_folder = engine.run(self.output_folder)
        return run_folder

    def plot_results(self, run_folder: str, top_n: int = None, save: bool = False):
        plotter = ResultPlotter(run_folder)
        plotter.plot_all(top_n=top_n, save=save)

    def simplify_solutions(self, run_folder: str, top_n: int = 1, max_workers: int = 4, save: bool = False):
        simplifier = SolutionSimplifier(self.rMTA_obj)
        simplifier.simplify(run_folder, top_n=top_n, max_workers=max_workers, save=save)


if __name__ == "__main__":

    from cobra.io import read_sbml_model
    import pandas as pd
    import numpy as np
    from Preprocessing.preprocessing import DifferentialExpression
    from MTAs.rmta import rMTA


    model = read_sbml_model('/home/catarinagomes/catarina_stuff/Recon3D_consistent.xml')
    model

    samples = pd.read_csv('/home/catarinagomes/catarina_stuff/Sampler_ACHR_Recon3D_EOL.csv', header=0)
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
    # samples = samples.T
    len(samples_temp)
    Vref = samples_temp.mean(axis=1)
    Vref
    Vref = pd.DataFrame(Vref)
    Vref.index = [x.id for x in model.reactions]
    Vref
    diff_express_genes = pd.read_csv('/home/catarinagomes/catarina_stuff/DGE_results_limma_voom_EOL.csv', header=0)
    diff_express_genes
    # read the dict tat is in txt format
    conv_dict = pd.read_csv('/home/catarinagomes/catarina_stuff/dict_conv_final.txt', sep='\t', header=None)
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
    logfc_thresh = 1
    pval_thresh = 0.05

    # Test
    diff_expr = DifferentialExpression(model, diff_express_genes, pd.DataFrame(Vref))
    diff_expr.get_differentially_expressed_genes(logfc_thresh, pval_thresh)
    diff_expr.get_differentially_expressed_reactions_from_genes()
    diff_expr.samples_filter(samples)
    from MTAs.rmta import rMTA

    rMTA_instance = rMTA(model, diff_expr.rxnFBS, Vref, 0.66, samples_temp)

    output_folder = "/home/catarinagomes/catarina_stuff/EAs/outputs"

    pipeline = Pipeline(rMTA_instance, output_folder)

    # adjust the EA params
    ea_params = {
        "pop_size": 400,
        "ngen": 100,
        "refresh_fraction": 0.75,
        "cxpb": 0.7,
        "mutpb": 0.2,
        "n_cores": 50,
        "elitism_size": 0,
        "repair": True,
        "remove_duplicates_bool": False,
        "max_knockouts": 10,
        "indpb": 0.01,
        "crossover_indpb": 0.5
    }

    run_folder = pipeline.run_ea(**ea_params)
    #run_folder = './outputs/RUN_YYYY-MM-DD_HH-MM'
    pipeline.plot_results(run_folder, top_n=25, save=True)
    pipeline.simplify_solutions(run_folder, top_n=2, max_workers=32, save=True)
