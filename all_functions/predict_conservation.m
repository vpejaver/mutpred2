% Vikas Pejaver
% December 2014

% Function to predict conservation features (primates, mammals and vertebrates), given a sequence and a
% PSSM (both indices and frequencies)

function [U, N, F] = predict_conservation(sequence, PSSM)

% Constants and defaults
U = [];
N = [];
F = [];
amino_acids = 'ARNDCQEGHILKMFPSTWYV'; % ordering as per PSI-BLAST matrix
breakpoints = [1 : 20 : 120];
offset = length(amino_acids) - 1;

% Extract sequence-based features
bin_mat = zeros(length(sequence), length(amino_acids));
for i = 1:length(sequence)
    bin_mat(i, amino_acids == sequence(i)) = 1;
end

% Get positions
abs_pos = [1:length(sequence)]';
rel_pos = abs_pos / length(sequence);

% Make feature matrix
X = [bin_mat PSSM abs_pos rel_pos];

% Predict
U = predict_cindex_features(X);
[meanv, stdv, N] = normalize(U, [], []);
raw_F = predict_freqs_features(X);

% Addional step - normalize so that frequencies sum to one
raw_F(raw_F < 0) = 0;
for i = 1:length(breakpoints)
    start = breakpoints(i);
    finish = breakpoints(i) + offset;
    M = raw_F(:, start : finish);
    F = [F, M ./ repmat(sum(M, 2), 1, size(M, 2))];
end

return
