import cobra
from cobra.io import load_matlab_model
import pandas as pd


# Load the original metabolic model 
model = cobra.io.read_sbml_model("Recon3D_consistent.xml")

# Load the context specific metabolic model from MATLAB format
model_context = load_matlab_model("tissueModel_imat_early.mat")



# Count elements from the original model
nreactions_og = len(model.reactions)
ngenes_og = len(model.genes)
nmetabolites_og = len(model.metabolites)

#Count elements from context specific model
nreactions_cont = len(model_context.reactions)
ngenes_cont = len(model_context.genes)
nmetabolites_cont = len(model_context.metabolites)


# Create dictionary with data 
data = {
    "Element": ["Reactions", "Genes", "Metabolites"],
    "Before integration": [nreactions_og, ngenes_og, nmetabolites_og],
    "After integration": [nreactions_cont, ngenes_cont, nmetabolites_cont]
}

# Criar DataFrame
table_comp = pd.DataFrame(data)

# Exportar como CSV
table_comp.to_csv("model_comparison_table.csv", index=False)

# Exibir a tabela
print(table_comp)
