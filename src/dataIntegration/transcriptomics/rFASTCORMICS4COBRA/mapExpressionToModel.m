function [mapping] = mapExpressionToModel(model, data, dico, rownames)
% The mapExpressionToModel function map the expression to a model following
% the GPR rules (Pacheco et al,2019)

% USAGE:
%
% [mapping] = mapExpressionToModel(model, data, dico, rownames)


% INPUTS:
% model:                 (the following fields are required - others can be supplied)
%                         * S  - `m x 1` Stoichiometric matrix
%                         * lb - `n x 1` Lower bounds
%                         * ub - `n x 1` Upper bounds
%                         * rxns   - `n x 1` cell array of reaction abbreviations
%

% data:                  expression or dicretized values for the samples, size(expression,1) =
%                        lenghth(gene IDs size(expression,2) = number of
%                        samples
% rownames:              cell aray with the gene IDs
% dico:                  table which contains corresponding gene identifier information. Needed
%                        to map the rownames to the genes in the model.

% OUTPUTS:
% mapping                 matrix containing the expression values that were mapped to the reactions acoording to the GPR rules
                         %size(mapping,1) is equal to the number of reactions and size (mapping,2) is equal to size(data,2)


% .. Authors:
%       - Maria Pires Pacheco, Thomas Sauter, 2016, University of
%       Luxembourg modified by Tamara Bintener
%       - Maria Pires Pacheco, Thomas Sauter, 2022, adaptation of the code to the Cobra toolbox
                         
%% search the dataIds in the dictionnary
if istable(dico)
    dico = table2array(dico);
end


if any(cellfun(@isempty, rownames))
    disp('rownames contains empty entries, please check')
    return
%     rownames(cellfun(@isempty, rownames)) = cellstr('not found');
end


col = (sum(ismember(dico,rownames)) == max(sum(ismember(dico,rownames)))); % find matching column, usually the one with the highest number of matches
[~,idico,irownames] = intersect(dico(:,col),rownames); % get indices

if isempty(irownames)
    'the dico does match the dataIds';
    return
end

mapped = data(irownames,:); %get data for matching rownames

mapped2(:,1) = rownames(irownames,1); %create dico with same order than input
for i = 1:size(dico,2)
    try        mapped2(:,i+1) = dico(idico,i);
end
    try         mapped2(:,i+1) = cellstr(dico(idico,i));
    end %might need conversion to cell
end

% deal with transcripts (if any)
if any(contains(model.genes,'.')) && ~contains(model.description, 'recon','IgnoreCase',true) % if possible transcripts but not recon
    warning('Does your input model contain transcripts? If not, please remove any dots . in the model.genes field')
    disp('Temporarily removing transcripts...')
    model.genes = regexprep(model.genes, '\.[0-9]*','');
elseif contains(model.description, 'recon','IgnoreCase',true) %if Recon model
    model.genes = regexprep(model.genes, '\.[0-9]*','');
end


mappedToGenes = zeros(numel(model.genes), size(mapped,2)); %initialise variable

genesMatched = 0;

% maps the discretized expression data to the genes

for i=1:numel(model.genes)
    [match,~] = find(ismember(mapped2,model.genes(i))); %find row of match
    if numel(match)==1
        mappedToGenes(i,:) = mapped(match,:);
    elseif isempty(match)
    else
        mappedToGenes(i,:) = max(mapped(match,:),[],1); % take the highest value if
        %more probeIDs correspond to one modelID
    end
    if ~isempty(match)
        genesMatched =  genesMatched + 1;
    end
end

fprintf('%i of %i genes matched\n', genesMatched, numel(model.genes))

mapping = zeros(numel(model.rxns), size(mapped,2)); % maps the expression data to the reactions
for j= 1:size(mappedToGenes,2)
    x = mappedToGenes(:,j);
    for k = 1:numel(model.rxns)       
        mapping(k,j) = GPRrulesMapperRFASTCORMICS(cell2mat(model.rules(k)),x);
    end
end

