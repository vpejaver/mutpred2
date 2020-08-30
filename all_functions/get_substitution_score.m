% Vikas Pejaver
% October 2012

% Function that takes in a wild type and mutant amino acid and
% returns the scores for corresponding substitutions from matrices
% such as BLOSUM and PAM
% NOTE: Returns a vector of values from all possible matrices in MATLAB

function [scores, feature_names] = get_substitution_score(W, M)

% Constants and defaults
scores = [];
feature_names = {};
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY';
blosumgen = [30:5:90, 62, 100];
pamgen = [10:10:500];

% Check error
if ~ismember(W, AMINO_ACIDS) || ~ismember(M, AMINO_ACIDS)
    error('Illegal substitution! Allowed characters are ACDEFGHIKLMNPQRSTVWY!');
end

% Get all BLOSUM scores
for i = 1:length(blosumgen)
    tmp_mat = blosum(blosumgen(i), 'order', AMINO_ACIDS);
    scores = [scores tmp_mat(AMINO_ACIDS == W, AMINO_ACIDS == M)];
    feature_names = [feature_names strcat('BLOSUM', num2str(blosumgen(i)))];
end

% Get all PAM scores
for i = 1:length(pamgen)
    tmp_mat = pam(pamgen(i), 'order', AMINO_ACIDS);
    scores = [scores tmp_mat(AMINO_ACIDS == W, AMINO_ACIDS == M)];
    feature_names = [feature_names strcat('PAM', num2str(pamgen(i)))];
end

return
