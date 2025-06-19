import os
import random
import time
import multiprocessing
import logging
from datetime import datetime
import pandas as pd
from deap import base, creator, tools
from aux import init_worker, evaluate_individual, init_individual, repair_individual, remove_duplicates
from typing import List, Dict, Any

# Setup logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


class MTA_EA:
    """
    Encapsulates the evolutionary algorithm (EA) for gene knockout analysis.

    This class requires a pre-built rMTA instance (which holds the model and related data).
    It applies elitism and, after each generation, updates a single CSV file that aggregates
    the top individuals from all generations. This allows you to plot the evolution at any time.

    All output files are saved into a dynamically generated subfolder (named with the run ID)
    inside a given base output folder.
    """

    def __init__(self, rMTA_obj: Any, pop_size: int, ngen: int, refresh_fraction: float,
                 cxpb: float, mutpb: float, n_cores: int, elitism_size: int,
                 base_output_folder: str,
                 repair: bool = False, remove_duplicates_bool: bool = False,
                 max_knockouts: int = 5, indpb: float = 0.01, crossover_indpb: float = 0.5) -> None:
        self.rMTA_obj = rMTA_obj
        self.pop_size = pop_size
        self.ngen = ngen
        self.refresh_fraction = refresh_fraction
        self.cxpb = cxpb
        self.mutpb = mutpb
        self.n_cores = n_cores
        self.repair = repair
        self.remove_duplicates_bool = remove_duplicates_bool
        self.max_knockouts = max_knockouts
        self.base_output_folder = base_output_folder
        self.indpb = indpb
        self.crossover_indpb = crossover_indpb
        self.elitism_size = elitism_size if elitism_size is not None else 0
        self.stop_early = False

        self.genes = list(self.rMTA_obj.model.genes)
        logger.info(f"Using rMTA instance with {len(self.genes)} genes.")
        self.top_individuals_all: List[Dict[str, Any]] = []

    def run_ea(self) -> None:
        run_id = "RUN_" + datetime.now().strftime("%Y-%m-%d_%H-%M")
        output_folder = os.path.join(self.base_output_folder, run_id)
        os.makedirs(output_folder, exist_ok=True)

        timing_log_path = os.path.join(output_folder, "timing.txt")
        final_results_path = os.path.join(output_folder, "final.csv")
        top_individuals_path = os.path.join(output_folder, "top.csv")
        script_name = os.path.basename(__file__) if "__file__" in globals() else "unknown_script"

        logger.info(f"Starting EA, run_id={run_id}, script={script_name}")
        logger.info(
            f"Repair={self.repair}, Remove Duplicates={self.remove_duplicates_bool}, MaxKnock={self.max_knockouts}")

        import signal
        def signal_handler(sig, frame):
            logger.info("Ctrl+C detected, will finish current generation then stop.")
            self.stop_early = True

        signal.signal(signal.SIGINT, signal_handler)

        total_start = time.time()
        pool = multiprocessing.Pool(
            processes=self.n_cores,
            initializer=init_worker,
            initargs=(self.rMTA_obj,)
        )
        toolbox = self.setup_toolbox(pool)
        pop = toolbox.population(n=self.pop_size)
        self.evaluate_population(pop, toolbox)

        gen_times: List[float] = []
        for gen in range(self.ngen):
            if self.stop_early:
                logger.info(f"Stopping early before generation {gen}")
                break
            gen_start = time.time()
            pop = self.ea_generation(pop, toolbox)
            best_ind = tools.selBest(pop, 1)[0]
            logger.info(f"Gen {gen}, Best fitness: {best_ind.fitness.values[0]}")
            gen_times.append(time.time() - gen_start)
            self.append_generation_top_individuals(gen, tools.selBest(pop, self.pop_size))
            self.save_all_generation_top_individuals(top_individuals_path)
            if self.stop_early:
                logger.info(f"Stopping early after generation {gen}")
                break

        pool.close()
        pool.join()

        best_ind = tools.selBest(pop, 1)[0]
        logger.info(f"Final Best individual: {best_ind}, Fitness: {best_ind.fitness.values[0]}")
        self.save_results(pop, final_results_path)
        self.save_timing_log(timing_log_path, script_name, run_id, gen_times, total_start, best_ind)
        logger.info(f"Run complete. All files saved in folder: {output_folder}")

    def setup_toolbox(self, pool: multiprocessing.Pool) -> Any:
        creator.create("FitnessMax", base.Fitness, weights=(1.0,))
        creator.create("Individual", list, fitness=creator.FitnessMax)
        toolbox = base.Toolbox()
        toolbox.register("individual", init_individual, creator.Individual,
                         size=len(self.genes), max_knockouts=self.max_knockouts)
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)
        toolbox.register("evaluate", evaluate_individual)
        toolbox.register("mate", tools.cxUniform, indpb=self.crossover_indpb)
        toolbox.register("mutate", tools.mutFlipBit, indpb=self.indpb)
        toolbox.register("select", tools.selTournament, tournsize=3)
        toolbox.register("map", pool.map)
        return toolbox

    def evaluate_population(self, population: List[Any], toolbox: Any) -> None:
        invalid_ind = [ind for ind in population if not ind.fitness.valid]
        fitnesses = list(toolbox.map(toolbox.evaluate, invalid_ind))
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

    def ea_generation(self, pop: List[Any], toolbox: Any) -> List[Any]:
        offspring = toolbox.select(pop, len(pop))
        offspring = list(map(toolbox.clone, offspring))
        # Crossover.
        for child1, child2 in zip(offspring[::2], offspring[1::2]):
            if random.random() < self.cxpb:
                toolbox.mate(child1, child2)
                del child1.fitness.values
                del child2.fitness.values
        # Optional repair.
        if self.repair:
            for ind in offspring:
                repair_individual(ind, max_knockouts=self.max_knockouts)
        # Mutation.
        for mutant in offspring:
            if random.random() < self.mutpb:
                toolbox.mutate(mutant)
                del mutant.fitness.values
        self.evaluate_population(offspring, toolbox)
        offspring.sort(key=lambda ind: ind.fitness.values[0], reverse=True)
        keep_count = int((1 - self.refresh_fraction) * self.pop_size)
        kept = offspring[:keep_count]
        new_inds = [toolbox.individual() for _ in range(self.pop_size - keep_count)]
        self.evaluate_population(new_inds, toolbox)
        next_gen = kept + new_inds
        if self.remove_duplicates_bool:
            next_gen = remove_duplicates(next_gen, toolbox, max_attempts=10)
            self.evaluate_population(next_gen, toolbox)
        # Apply elitism.
        elites = tools.selBest(pop, self.elitism_size)
        next_gen.extend(elites)
        next_gen.sort(key=lambda ind: ind.fitness.values[0], reverse=True)
        next_gen = next_gen[:self.pop_size]
        return next_gen

    def append_generation_top_individuals(self, gen: int, top_individuals: List[Any]) -> None:
        for ind in top_individuals:
            knocked_genes = [self.genes[i].id for i, val in enumerate(ind) if val == 0]
            self.top_individuals_all.append({
                "Generation": gen,
                "Individual": ''.join(map(str, ind)),
                "Knocked Genes": ','.join(knocked_genes),
                "Knocked Count": len(knocked_genes),
                "Fitness": ind.fitness.values[0]
            })

    def save_all_generation_top_individuals(self, top_individuals_path: str) -> None:
        df = pd.DataFrame(self.top_individuals_all)
        df.to_csv(top_individuals_path, index=False)
        logger.info(f"Aggregated top individuals updated in {top_individuals_path}")

    def save_results(self, pop: List[Any], final_results_path: str) -> None:
        records = []
        for ind in pop:
            knocked_genes = [self.genes[i].id for i, val in enumerate(ind) if val == 0]
            records.append({
                "Knocked Genes": ','.join(knocked_genes),
                "Knocked Count": len(knocked_genes),
                "Fitness": ind.fitness.values[0]
            })
        df = pd.DataFrame(records)
        df.to_csv(final_results_path, index=False)
        logger.info(f"Final results saved to {final_results_path}")

    def save_timing_log(self, timing_log_path: str, script_name: str, run_id: str,
                        gen_times: List[float], total_start: float, best_ind: Any) -> None:
        total_run_time = time.time() - total_start
        mean_gen_time = sum(gen_times) / len(gen_times) if gen_times else 0.0
        with open(timing_log_path, 'w') as f:
            f.write("=== EA Timing & Parameters Log ===\n\n")
            f.write(f"Script Name: {script_name}\n")
            f.write(f"Run ID: {run_id}\n\n")
            f.write("=== Parameters Used ===\n")
            f.write(f"Number of Genes: {len(self.genes)}\n")
            f.write(f"Population Size: {self.pop_size}\n")
            f.write(f"Generations: {self.ngen}\n")
            f.write(f"Refresh Fraction: {self.refresh_fraction}\n")
            f.write(f"Crossover Probability: {self.cxpb}\n")
            f.write(f"Crossover indpb: {self.crossover_indpb}\n")
            f.write(f"Mutation Probability: {self.mutpb}\n")
            f.write(f"n_cores: {self.n_cores}\n")
            f.write(f"Repair: {self.repair}\n")
            f.write(f"Remove Duplicates: {self.remove_duplicates_bool}\n")
            f.write(f"Max Knockouts: {self.max_knockouts}\n")
            f.write(f"indpb: {self.indpb}\n")
            f.write(f"Elitism Size: {self.elitism_size}\n")
            f.write(f"\nFinal Fitness: {best_ind.fitness.values[0]}\n")
            f.write("\n=== Timing Info ===\n")
            f.write(f"Total Run Time: {total_run_time:.4f} seconds\n")
            f.write(f"Mean Generation Time: {mean_gen_time:.4f} seconds\n")
            f.write(f"Number of Generations Run: {len(gen_times)}\n")


def main() -> None:
    from MTA.MTAs.rmta import rMTA
    from MTA.Preprocessing.preprocessing import DifferentialExpression
    from cobra.io import read_sbml_model
    import numpy as np

    model_path = '../examples/Case Study/iCEL1314_consistent.xml'
    diff_expr_path = '../examples/Case Study/limma_unc_results.csv'
    samples_mat_path = '../examples/Case Study/Sampler_ACHR_unc_imat.csv'

    model = read_sbml_model(model_path)

    samples_mat = pd.read_csv(samples_mat_path)
    samples_mat = samples_mat.T

    idx = np.zeros((samples_mat.shape[0],), dtype=int)
    for i in range(len(idx)):
        idx[i] = model.reactions.index(samples_mat.index[i])
    rxnInactive = set(range(len(model.reactions))) - set(idx)

    samples = np.zeros((len(model.reactions), samples_mat.shape[1]))
    samples[idx, :] = samples_mat

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
    diff_express_genes = diff_express_genes.dropna()
    diff_express_genes = diff_express_genes.iloc[:, [0, 1, 4]]

    de = DifferentialExpression(model, diff_express_genes, vref)
    de.get_differentially_expressed_genes()
    de.get_differentially_expressed_reactions_from_genes()
    de.samples_filter(samples_mat)

    rMTA_instance = rMTA(model, de.rxnFBS, vref, 0.66, samples)

    base_output_folder = "trash_outputs"

    ea = MTA_EA(
        rMTA_obj=rMTA_instance,
        pop_size=256,
        ngen=250,
        refresh_fraction=0.50,
        cxpb=0.7,
        mutpb=0.2,
        n_cores=40,
        elitism_size=0,
        repair=True,
        remove_duplicates_bool=False,
        max_knockouts=20,
        base_output_folder=base_output_folder,
        indpb=0.01,
        crossover_indpb=0.5
    )

    ea.run_ea()


if __name__ == "__main__":
    main()
