function [exRxns, ExOrgaInd, ExOrga] = findEXRxnsRFASTCORMICS(model, biomassReactionID, ObjectiveFunction)
%The function identifies exchange reactions

% USAGE:
%
%[exRxns, ExOrgaInd, ExOrga] = findEXRxnsRFASTCORMICS(model, biomassReactionID, ObjectiveFunction)
%INPUTS:
%    model:             (the following fields are required - others can be supplied)
%                         * S  - `m x 1` Stoichiometric matrix
%                         * lb - `n x 1` Lower bounds
%                         * ub - `n x 1` Upper bounds
%                         * rxns   - `n x 1` cell array of reaction
%                         abbreviations
%OPTIONAL INPUTS:

% notMediumConstrained:   cell array with the metabolites that should not be constrained
% biomassReactionID:              name of the biomass reaction
% ObjectiveFunction:            cell array with the identifiers of objective functions

% OUTPUTS:
% exRxns                    cell array with the exchange reactions
% ExOrgaInd                indices of the exchnage reactions for organic
%                          compounds
% ExOrga                   cell array of the exchange reactions for organic
%                          compounds
% .. Authors:
%       - Maria Pires Pacheco, Thomas Sauter, 2016, University of Luxembourg
%       - Maria Pires Pacheco, Thomas Sauter, 2022, adaptation of the code to the Cobra toolbox

exchange_mets   = [];
exRxnsInd       = find(sum(abs(model.S),1)==1);
biomass_id      = find(ismember(model.rxns, biomassReactionID));


if isempty(biomass_id) && isempty(ObjectiveFunction)
    warning('No biomass set for the model, please verify if the medium constraints do not affect the biomass production')
elseif isempty(biomass_id) && ~isempty(ObjectiveFunction)
    warning('No biomass set for the model, please verify if the medium constraints do not affect the biomass production')
    function_id = find(ismember(model.rxns, ObjectiveFunction));
    exRxnsInd = setdiff(exRxnsInd, function_id);
else
    function_id = find(ismember(model.rxns, ObjectiveFunction));
    function_id = unique([biomass_id; function_id]);
    exRxnsInd = setdiff(exRxnsInd, function_id);
end

for i=1:numel(exRxnsInd)
    exmets  = find(model.S(:,exRxnsInd(i)));
    exchange_mets(end+1) = exmets;
    if model.S(exmets,exRxnsInd(i))==1
        model.S(exmets,exRxnsInd(i))=-1;
    end
end
exRxns  = model.rxns(exRxnsInd);
model.mets(exchange_mets);
ex_mets_carbon  = (regexp(model.metFormulas(exchange_mets),'C'));
is_organic      = ~cellfun('isempty', ex_mets_carbon);
mets_ex_orga    = exchange_mets(is_organic);
mets_ex_inorga  = setdiff(exchange_mets,mets_ex_orga);
model.mets(mets_ex_orga);
model.mets(mets_ex_inorga);
ExOrgaInd  = exRxnsInd(is_organic);
ExOrga     = model.rxns(ExOrgaInd);
end