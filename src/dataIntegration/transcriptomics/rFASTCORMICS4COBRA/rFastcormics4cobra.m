function [model,Afinal] = rFastcormics4cobra(model, discretized, rownames, colnames,dico, consensusProportion, epsilon, optionalSettings, biomassReactionName)
% The rFASTCORMICS is a context-specific building algorithm for
% reconstructing a tissue, a cell-specific, or any context-specific model from RNAseq data
% and a generic reconstruction (Pacheco et al. 2019)

% USAGE:
%
%    A = rFastcormics4cobra(model, discretized, rownames, colnames,dico, consensusProportion, epsilon, optionalSettings, biomassReactionName, printFigures, path)
% REQUIREMENTS:    Matlab
%                         * Statistics and Machine Learning Toolbox
%                         * Curve fitting toolbox
% INPUTS:
% model:                 (the following fields are required - others can be supplied)
%                         * S  - `m x 1` Stoichiometric matrix
%                         * lb - `n x 1` Lower bounds
%                         * ub - `n x 1` Upper bounds
%                         * rxns   - `n x 1` cell array of reaction abbreviations
%                         * metFormulas m*1 metabolite Formulas
% discretized:                  discretized values for the samples, size(discretized,1) =
%                        length(geneIDs), size(discretized,2)= length(colnames)
% rownames:              cell array with the gene IDs
% colnames:              cell array with the sample names
% dico:                  table which contains corresponding gene identifier information. Needed
%                        to map the rownames to the genes in the model.

% OPTIONAL INPUTS:
% consensusProportion :  the rate of samples that have to express or not to express a gene for the
%                        gene to be considered expressed or not in the context of interest
% epsilon:       smallest flux that is considered nonzero
% optionalSettings:       object
%                        * Func - cell array of reaction abbreviations that should carry a flux
%                        * medium - cell array of metabolites abbreviation that defines metabolites
%                          in the growth medium of cells to constrain the model, see example medium_example.mat
%                        * not_medium_constrained- react
% biomassReactionName:      cell array or all other functions that might require some exchange reactions to
% printFigures:           1 print all plots, 0 do not print plots
% path:                    Folder in which the plots should be saved




% OUTPUT:
% A            indices of the active reactions

% .. Authors:
%       - Maria Pires Pacheco, Thomas Sauter, 2016, University of Luxembourg
%       - Maria Pires Pacheco, Thomas Sauter, 2023, adaptation of the code to the Cobra toolbox

if nargin <11
    p = genpath(pwd);
    p = strtok(p,';');
    p = strcat(p,'/');    
end

if nargin < 10
    printFigures = [];
end


if nargin < 9
    biomassReactionName = [];
end
if nargin < 8
    optionalSettings = [];
end
if nargin < 7 || isempty(epsilon)
    epsilon = getCobraSolverParams('LP', 'feasTol')*100;
end
if nargin<6
    consensusProportion = 0.9;
end


% Data discretization
% The discretized values will be discretized into expressed, not expressed, and unknown
% expression status.


if ischar(model.subSystems{1})
else
    model.subSystems = vertcat(model.subSystems{:});
end
% fix reversibilities
modelOrg=model;

model = fixIrrRFASTCORMICS(model);
[A, ~, ~] = fastccCobra(model, epsilon, 0, 0, 'original');
model=generateRules(model);
% create aconsistent models
consistentModel = removeRxns(model, model.rxns(setdiff(1:numel(model.rxns),A)));
clear model

mapping = mapExpressionToModel(consistentModel, discretized, dico, rownames);
mapping = sparse(mapping);


if sum(mapping) == 0
    disp('no genes were mapped, check again')
    return
end

if numel(rownames) ~= size(discretized,1)
    disp('data and iDs do not correspond')
    return
end
if size(mapping,1)~= numel(consistentModel.rxns)
    disp('when the option already_mapped is used the size of the data has to correpond to the number of reaction of the model')
    return
end

requiredToolboxes = {'curve_fitting_toolbox', 'optimization_toolbox'};
if sum(isstruct(prepareTest('requiredToolboxes', requiredToolboxes)))
%% Optionally the model can be constrained in function of the medium
%composition
if ~isempty(optionalSettings) && isfield(optionalSettings, 'medium')
    if isfield(optionalSettings, 'notMediumConstrained')
        notMediumConstrained = optionalSettings.notMediumConstrained;
    else
        notMediumConstrained = [];
    end
    mediumMets = optionalSettings.medium;
    if ~isfield(optionalSettings,'func')
        optionalSettings.func = '';
    end
    
    mediumConstrainedmodel =  constrainModelRFASTCORMICS(consistentModel, mediumMets, notMediumConstrained, biomassReactionName, optionalSettings.func);
    
elseif isempty(optionalSettings)
    warning('No optional settings detected')
else
    warning('No medium set')
end
%% Identify reactions that are under the control of expressed genes
numberOfArrayPerModel = size(discretized,2);
C =  find(sum(mapping,2) >= (consensusProportion * numberOfArrayPerModel));

%% Additions of the reactions needed for a given function to carry a flux
% to the core set

if isfield(optionalSettings, 'func')
    functionKeep = optionalSettings.func;
else
    functionKeep = '';
end

if ~isempty(functionKeep)
    
    
    B = find(ismember(mediumConstrainedmodel.rxns,functionKeep));
    if isempty(B)
        warning('no functions set to be kept')
    elseif numel(B) ~= numel(functionKeep)
        warning('Not all functions set to be kept were found in the model')
    end
    BiomassRelatedRxns = fastcoreCobraNonPen(mediumConstrainedmodel, B, epsilon, 0, C);

    C = union(C,BiomassRelatedRxns); % add the ATP and biomass reactions to the core set
    
else
    BiomassRelatedRxns = [];
end

%% Identification of the inactive reactions set (medium composition and
% reactions under control of unexpressed genes the or  and removing of
% inactive branches
notexpressed = find(sum(mapping,2) <= (- consensusProportion * numberOfArrayPerModel));
notexpressed = setdiff(notexpressed,BiomassRelatedRxns);
mediumConstrainedmodel.lb(notexpressed) = 0;
mediumConstrainedmodel.ub(notexpressed) = 0;

% Building of a consistent model
[A, ~, ~] = fastcc(mediumConstrainedmodel, epsilon, 0, 0, 'original');
ConsistentmediumConstrainedmodel = removeRxns(mediumConstrainedmodel,mediumConstrainedmodel.rxns(setdiff(1:numel(mediumConstrainedmodel.rxns),A)));

%% Establishment of the reactions in the core set

[~, ~, IB] = intersect(consistentModel.rxns(C), ConsistentmediumConstrainedmodel.rxns);
C = IB;
%Removal of transporters from the core set

NonPenalised = [];
if ~isempty(optionalSettings) && isfield(optionalSettings, 'unpenalized')
    NonPenalised = find(ismember(ConsistentmediumConstrainedmodel.rxns(C),optionalSettings.unpenalized));
    NonPenalised = C(unique(NonPenalised));
end
C = setdiff(C,NonPenalised);
%% building of the context-specific model
A2 = fastcoreCobraNonPen(ConsistentmediumConstrainedmodel, C, 1e-4, 0, NonPenalised);
else
    disp( 'Statistics and Machine Learning and Curve fitting toolbox must be installed')                         
end
%% building of the context-specific model
model   = removeRxns(ConsistentmediumConstrainedmodel,ConsistentmediumConstrainedmodel.rxns(setdiff(1:numel(ConsistentmediumConstrainedmodel.rxns),A2)));
Afinal = find(ismember(modelOrg.rxns,ConsistentmediumConstrainedmodel.rxns(A2)));
end
