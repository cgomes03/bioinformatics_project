import re
import cobra
import pandas as pd

# Import omics datasets and required Troppo modules
from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ModelBasedWrapper
from troppo.omics.integration import ContinuousScoreIntegrationStrategy
from troppo.methods.reconstruction.imat import IMAT, IMATProperties

# ------------------------------------------------------------------------------
# 1) GPR Cleanup: Remove '__COBAMPGPRDOT__X' artifacts from GPR strings
# ------------------------------------------------------------------------------
patt = re.compile(r'__COBAMPGPRDOT__[0-9]{1}')

def replace_alt_transcripts(gpr_string: str) -> str:
    return patt.sub('', gpr_string)

# ------------------------------------------------------------------------------
# 2) Load Model and Omics Data
# ------------------------------------------------------------------------------
config = cobra.Configuration()
config.solver = "glpk"

model = cobra.io.read_sbml_model('Recon3D_consistent.xml')
omics_data = pd.read_csv('early_onset_expression.tsv', sep="\t",  index_col=0)

# Create a TabularReader container from the omics DataFrame
omics_container = TabularReader(
    path_or_df=pd.DataFrame(omics_data).T,
    nomenclature='ensemble_gene_id',  # Adjust if necessary
    omics_type='transcriptomics'
).to_containers()[0]

# ------------------------------------------------------------------------------
# 3) Map Omics Data to the Model
# ------------------------------------------------------------------------------
model_wrapper = ModelBasedWrapper(
    model=model,
    ttg_ratio=9999,
    gpr_gene_parse_function=replace_alt_transcripts
)

data_map = omics_container.get_integrated_data_map(
    model_reader=model_wrapper.model_reader,
    and_func=min,  # For AND rules in GPR
    or_func=sum    # For OR rules in GPR
)

# ------------------------------------------------------------------------------
# 4) Integrate Reaction Scores
# ------------------------------------------------------------------------------
def score_apply(reaction_map_scores):
    """Replace None with 0 and return a list of reaction-level scores."""
    return [0 if score is None else score for _, score in reaction_map_scores.items()]

continuous_integration = ContinuousScoreIntegrationStrategy(score_apply=score_apply)
scores = continuous_integration.integrate(data_map=data_map)
print(scores)

from collections import Counter
print(Counter(scores))

data = {'rxnExpr': scores}


# Save the reaction expression data to a .mat file
import scipy.io
scipy.io.savemat('expressionData_imat_early.mat', data)
print("Expression data saved to expressionData_imat_early.mat")

