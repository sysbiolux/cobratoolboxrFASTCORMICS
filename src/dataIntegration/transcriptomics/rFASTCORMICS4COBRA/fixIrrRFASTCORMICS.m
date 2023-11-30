function[model] = fixIrrRFASTCORMICS(model)
% The function converts irreversible backwards reactions into irreversible
% forward reactions

% USAGE:
%
%   [model] = fixIrrRFASTCORMICS(model)

%INPUTS:
%    model:             (the following fields are required - others can be supplied)
%                         * S  - `m x 1` Stoichiometric matrix
%                         * lb - `n x 1` Lower bounds
%                         * ub - `n x 1` Upper bounds
%                         * rxns   - `n x 1` cell array of reaction
%                         abbreviations

% OUTPUT:                 
% model:                 model with corrected reversibilties

% .. Authors:
%       - Maria Pires Pacheco, Thomas Sauter, 2016, University of Luxembourg
%       - Maria Pires Pacheco, Thomas Sauter, 2022, adaptation of the code to the Cobra toolbox


model.rev = zeros(numel(model.rxns),1);
model.rev(model.lb <0 & model.ub> 0) = 1;
Irr =(model.lb >=0 & model.ub >0| model.ub <=0 & model.lb <0);
model.rev(Irr) = 0;
FakeIrr= model.ub <=0 & model.lb<0;
model.S(:, FakeIrr) = -model.S(:,FakeIrr);
model.ub(FakeIrr) = -model.lb(FakeIrr);
model.lb(FakeIrr) = zeros(sum(FakeIrr),1);
end