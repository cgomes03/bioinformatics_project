%% iMAT tissue?specific model extraction with objective protection

% 1. Load COBRA model and expression data
load('expressionData_imat_late.mat');    % Contains rxnExpr
% Assumes `model` is already in your workspace, otherwise:
% model = readCbModel('yourModel.mat');

% 2. Map expression values to reactions
%    Create a map from reaction IDs to expression scores
expressionMap = containers.Map(model.rxns, rxnExpr);

%    Initialize array of length nRxns
nRxns = numel(model.rxns);
expressionRxns = zeros(nRxns,1);

for i = 1:nRxns
    rxnID = model.rxns{i};
    if isKey(expressionMap, rxnID)
        expressionRxns(i) = expressionMap(rxnID);
    else
        % Mark reactions without data as “low” (or unknown)
        expressionRxns(i) = -1;
    end
end

% 3. Set your expression thresholds
threshold_lb = 0.2;   % below or equal = “lowly expressed”
threshold_ub = 0.8;   % above or equal = “highly expressed”

% 4. Protect the objective (e.g. biomass) reaction
%    Find the index of the reaction with a nonzero objective coefficient
biomassIdx = find(model.c ~= 0, 1);
if isempty(biomassIdx)
    error('No reaction in model.c has a nonzero objective coefficient.');
end

%    Force its expression score just above the high?expr threshold
expressionRxns(biomassIdx) = threshold_ub + eps;

runtimeLimit = 300;

% 5. Run iMAT
tissueModel = iMAT(model, expressionRxns, threshold_lb, threshold_ub,[], [], [] , runtimeLimit);

% 6. Remove any orphan metabolites (optional)
if exist('removeUnusedMets','file')
    tissueModel = removeUnusedMets(tissueModel);
    fprintf('Orphan metabolites removed from the tissue?specific model.\n');
else
    warning('removeUnusedMets not found. Skipping orphan metabolite removal.');
end

% 7. Save the resulting tissue?specific model
outputFile = 'tissueModel_imat_late.mat';
save(outputFile, 'tissueModel');
fprintf('iMAT extraction complete. Tissue?specific model saved to %s\n', outputFile);

% 8. Display basic stats
fprintf('Extracted model contains %d reactions and %d metabolites.\n', ...
    numel(tissueModel.rxns), numel(tissueModel.mets));
