function pred = VSL2B(seq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function pred = VSL2B(seq)
%
% VSL2 (Baseline) Protein Disorder Predictor
% 
% Inputs:
%        seq:    The amino acid equence - only the 20 standard amino acids 
%                are allowed
%
% Output:
%        pred:   A length(seq)-by-1 vector of per-residue disorder probabilities
% 
% Author:
%        Kang Peng
%        Temple University
%        2003-2006
%
% References:
%        1. Peng K., Radivojac P., Vucetic S., Dunker A.K., and Obradovic Z., Length-Dependent 
%           Prediction of Protein Intrinsic Disorder, BMC Bioinformatics 7:208, 2006. 
%        2. Obradovic Z., Peng K., Vucetic S., Radivojac P., and Dunker A.K., Exploiting 
%           Heterogeneous Sequence Properties Improves Prediction of Protein Disorder, 
%           Proteins 61(S7):176-182, 2005. 
%
% Note: 
%        The model used is the VSL2B predictor, as described in Peng et al. 2006. That is, 
%        attributes are constructed using the amino acid sequence only.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global CURRDIR;

pred = [];

if length(seq) < 25
    error('VSL2: input sequence should have at least 25 residues!!');
end

% check input sequence
seq = upper(seq);
if isempty(setdiff(seq, 'WCFIYVLHMATRGQSNPDEK')) == 0
    error('VSL2: input sequence contains non-standard amino acid codes !!');
end

% fprintf('\tMaking VSL2 prediction ...\n');

data15 = makeattr(seq, 15);
data41 = makeattr(seq, 41);
data61 = makeattr(seq, 61);

% "All-in-one" linear SVM models - the weights combine (a) attribute normalization, 
% (b) linear SVM forward calculation, and (c) logisitic regression to convert   
% the SVM outputs into probability scores

load(strcat(CURRDIR, filesep, 'all_models', filesep, 'VSL2B.mat'), 'weightL', 'biasL', 'weightS', 'biasS', 'weightM', 'biasM');

% Long disorder prediction (VSL2-L)
predL = [ones(size(data41, 1), 1) data41] * [biasL weightL]';
predL = 1 ./ (1 + exp(-predL));
predL = smoothing(predL, 31);

% Short disorder prediction (VSL2-S)
predS = [ones(size(data15, 1), 1) data15] * [biasS weightS]';
predS = 1 ./ (1 + exp(-predS));
predS = smoothing(predS, 5);

% Meta prediction (VSL2-M)
predM = [ones(size(data61, 1), 1) data61] * [biasM weightM]';
predM = 1 ./ (1 + exp(-predM));

% Final VSL2 prediction
pred = predL .* predM + predS .* (1 - predM);

return
