import os
import argparse
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing

from cobra.io import read_sbml_model

from MTA.MTAs.rmta import rMTA
from MTA.Preprocessing.preprocessing import DifferentialExpression

def get_reactions_from_geneID(gene_list, model):
    """
    Given a list of gene IDs, return a list of associated reaction IDs.
    """
    reactions = []
    for gene in gene_list:
        try:
            gene_obj = model.genes.get_by_id(gene)
            reactions.extend([rxn.id for rxn in gene_obj._reaction])
        except Exception as e:
            logging.getLogger(__name__).exception(f"Error processing gene {gene}")
    return reactions


GLOBAL_MODEL = None
GLOBAL_RMTA = None

def init_worker(model, rMTA_instance):
    """
    Initializer for worker processes.
    """
    global GLOBAL_MODEL, GLOBAL_RMTA
    GLOBAL_MODEL = model
    GLOBAL_RMTA = rMTA_instance


def process_individual(item):
    """
    Process one individual solution.
    'item' is a tuple (ind_idx, knocked_genes_str), where knocked_genes_str
    is a comma-separated string of gene IDs.
    """
    ind_idx, knocked_genes_str = item
    logger = logging.getLogger(__name__)
    logger.debug(f"Processing individual {ind_idx} with knocked_genes_str: {knocked_genes_str!r}")

    if knocked_genes_str is None or knocked_genes_str.strip() == "":
        logger.error(f"Individual {ind_idx} has an invalid 'Knocked Genes' value: {knocked_genes_str!r}")
        return (ind_idx, None)

    # Convert the string to a list of gene IDs.
    initial_genes = [g.strip() for g in knocked_genes_str.split(",") if g.strip()]
    if not initial_genes:
        logger.error(f"Individual {ind_idx} resulted in an empty gene list after splitting.")
        return (ind_idx, None)

    def evaluate_knockout(gene_list):
        try:
            model_copy = GLOBAL_MODEL.copy()
        except Exception as e:
            logger.exception("Error copying model in worker.")
            return -1e6
        associated_reactions = get_reactions_from_geneID(gene_list, model_copy)
        try:
            # Evaluate knockout using the rMTA instance.
            value, _, _, _ = GLOBAL_RMTA.evaluate_KO(associated_reactions)
            return float(value.rTS.values[0])
        except Exception as e:
            logger.exception("Error evaluating knockout in worker.")
            return -1e6


    trace = []
    step = 0
    current_genes = initial_genes.copy()
    current_fitness = evaluate_knockout(current_genes)
    trace.append({
        "Step": step,
        "Remaining Genes": ";".join(current_genes),
        "Genes Count": len(current_genes),
        "Removed Gene": None,
        "Fitness": current_fitness,
        "Delta": 0.0
    })

    while len(current_genes) > 1:
        best_candidate = None
        best_new_fitness = -np.inf
        best_delta = None
        for gene in current_genes:
            candidate_set = [g for g in current_genes if g != gene]
            new_fitness = evaluate_knockout(candidate_set)
            delta = current_fitness - new_fitness
            if new_fitness > best_new_fitness:
                best_new_fitness = new_fitness
                best_candidate = gene
                best_delta = delta
        current_genes.remove(best_candidate)
        step += 1
        trace.append({
            "Step": step,
            "Remaining Genes": ";".join(current_genes),
            "Genes Count": len(current_genes),
            "Removed Gene": best_candidate,
            "Fitness": best_new_fitness,
            "Delta": best_delta
        })
        current_fitness = best_new_fitness

    return (ind_idx, trace)


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] : %(message)s")
    logger = logging.getLogger(__name__)
    parser = argparse.ArgumentParser(
        description="Parallel iterative simplification of knockout solutions using final.csv. "
                    "The input CSV must have columns: 'Knocked Genes', 'Knocked Count', 'Fitness'."
    )
    parser.add_argument("folder", help="Folder containing final.csv and where outputs will be saved.")
    parser.add_argument("--topn", type=int, default=1,
                        help="Number of top individuals (by Fitness) to simplify (default: 1).")
    parser.add_argument("--max_workers", type=int, default=4, help="Maximum number of parallel workers.")
    parser.add_argument("--save", action="store_true", help="Save plot to file instead of showing interactively.")
    args = parser.parse_args()

    folder = os.path.abspath(args.folder)
    final_csv_path = os.path.join(folder, "final.csv")
    if not os.path.isfile(final_csv_path):
        logger.error(f"'final.csv' not found in {folder}")
        return

    model_path = '../examples/Case Study/iCEL1314_consistent.xml'
    diff_expr_path = '../examples/Case Study/limma_unc_results.csv'
    samples_mat_path = '../examples/Case Study/Sampler_ACHR_unc_imat.csv'

    model = read_sbml_model(model_path)
    samples_mat = pd.read_csv(samples_mat_path).T
    idx = np.zeros((samples_mat.shape[0],), dtype=int)
    for i in range(len(idx)):
        idx[i] = model.reactions.index(samples_mat.index[i])
    samples = np.zeros((len(model.reactions), samples_mat.shape[1]))
    samples[idx, :] = samples_mat.values
    vref = samples.mean(axis=1)
    vref = pd.DataFrame(vref)
    vref.index = [reaction.id for reaction in model.reactions]
    diff_express_genes = pd.read_csv(diff_expr_path)
    gene_dict = {}
    for gene in model.genes:
        gene_dict[gene.name] = gene.id
        dot = gene.name.split('.')[0]
        gene_dict[dot] = gene.id
    diff_express_genes['GENE'] = diff_express_genes['GENE'].apply(lambda x: gene_dict.get(x))
    diff_express_genes = diff_express_genes.dropna().iloc[:, [0, 1, 4]]
    de = DifferentialExpression(model, diff_express_genes, vref)
    de.get_differentially_expressed_genes()
    de.get_differentially_expressed_reactions_from_genes()
    de.samples_filter(samples_mat)
    rMTA_instance = rMTA(model, de.rxnFBS, vref, 0.66, samples)

    df_final = pd.read_csv(final_csv_path)
    if df_final.empty:
        logger.error("final.csv is empty.")
        return
    df_final = df_final.sort_values(by="Fitness", ascending=False).head(args.topn)
    selected_items = []
    for i, row in df_final.iterrows():
        knocked_genes_str = row["Knocked Genes"]
        selected_items.append((i, knocked_genes_str))
    total_individuals = len(selected_items)
    logger.info(f"Evaluating {total_individuals} individuals for simplification.")

    results = []
    with multiprocessing.Pool(processes=args.max_workers,
                              initializer=init_worker,
                              initargs=(model, rMTA_instance)) as pool:
        for i, res in enumerate(pool.imap(process_individual, selected_items), start=1):
            results.append(res)
            logger.info(f"Processed {i}/{total_individuals} individuals.")

    all_traces = []
    for ind_idx, trace in results:
        if trace is None:
            logger.error(f"Individual {ind_idx} processing failed.")
            continue
        trace_df = pd.DataFrame(trace)
        trace_df["Individual"] = ind_idx
        all_traces.append(trace_df)
    if not all_traces:
        logger.error("No valid individual traces were produced.")
        return

    combined_trace = pd.concat(all_traces, ignore_index=True)
    out_csv = os.path.join(folder, "simplification_trace_all.csv")
    combined_trace.to_csv(out_csv, index=False)
    logger.info(f"Combined simplification trace saved to: {out_csv}")

    plt.figure(figsize=(10, 6))
    for ind, df_ind in combined_trace.groupby("Individual"):
        plt.plot(df_ind["Genes Count"], df_ind["Fitness"], marker='o', label=f"Ind {ind}")
    plt.xlabel("Solution Length")
    plt.ylabel("Score")
    plt.title("Simplification Trace for All Individuals")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    if args.save:
        plot_path = os.path.join(folder, "simplification_plot_all.png")
        plt.savefig(plot_path)
        logger.info(f"Combined plot saved to: {plot_path}")
    else:
        plt.show()

if __name__ == "__main__":
    main()
