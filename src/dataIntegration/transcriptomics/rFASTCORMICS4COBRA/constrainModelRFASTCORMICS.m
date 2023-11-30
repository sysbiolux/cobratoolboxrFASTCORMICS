function [model] = constrainModelRFASTCORMICS(model, mediumMets, notMediumConstrained, biomassReactionID, ObjectiveFunction)
% The constrainModelRFASTCORMICS closes the inputs of metabolites that are not present in the medium 

% USAGE:
%
%[model] = constrainModelRFASTCORMICS(model, mediumMets, notMediumConstrained, biomassReactionID, ObjectiveFunction)
% 

%INPUTS:
% model:                 (the following fields are required - others can be supplied)
%                         * S  - `m x 1` Stoichiometric matrix
%                         * lb - `n x 1` Lower bounds
%                         * ub - `n x 1` Upper bounds
%                         * rxns   - `n x 1` cell array of reaction abbreviations
%

% mediumMets:              cell array with the identifiers of the metabolites

%OPTIONAL INPUTS:

% notMediumConstrained:   cell array with the metabolites that should not be constrained
% biomassReactionID:              name of the biomass reaction
% ObjectiveFunction:            cell array with the identifiers of objective functions

% OUTPUTS:
% model                     constrained cobra model

% .. Authors:
%       - Maria Pires Pacheco, Thomas Sauter, 2016, University of
%       Luxembourg modified by Tamara Bintener
%       - Maria Pires Pacheco, Thomas Sauter, 2022, adaptation of the code to the Cobra toolbox
             

if nargin <3
      notMediumConstrained=[];
end

if nargin < 4
    biomassReactionID=[];
end
if nargin < 5
    ObjectiveFunction=[];
end

if ~isempty(notMediumConstrained)
    notMediumConstrained = find(ismember(model.rxns,notMediumConstrained));
    lb = model.lb(notMediumConstrained);
    ub = model.ub(notMediumConstrained);
end
[exRxns, Ex_orgaInd] = findEXRxnsRFASTCORMICS(model, biomassReactionID, ObjectiveFunction);
notConstrained = intersect(findRxnsFromMets(model,mediumMets),exRxns);
notConstrained = find(ismember(model.rxns,notConstrained));
lb2 = model.lb(notConstrained);
ub2 = model.ub(notConstrained);
[~,match]= find(model.S(:,Ex_orgaInd)<0);
if ~isempty(match)
model.lb(Ex_orgaInd(match)) = 0; % close all the carbon sources
end
[~,match]= find(model.S(:,Ex_orgaInd)>0);
if ~isempty(match)
model.ub(Ex_orgaInd(match)) = 0; % close all the carbon sources
end
if ~isempty(notMediumConstrained)
    model.lb(notMediumConstrained) = lb;
    model.ub(notMediumConstrained) = ub;
end
if  ~isempty(notConstrained)
    model.lb(notConstrained) = lb2;
    model.ub(notConstrained) = ub2;
end
end
