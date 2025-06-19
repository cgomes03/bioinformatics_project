addpath(genpath('cobratoolbox'))

% if initcibra is not initialized, run the following command
initCobraToolbox(false)

%initCobraToolbox(false)
changeCobraSolver ('gurobi', 'all');


% change the current folder to  C:\Users\bsa\Desktop\troppo_validation
cd 'C:\Users\bsa\Desktop\DATASET_CATARINA\context_models'

% Load the model
model = readCbModel('Recon3D_consistent.xml');