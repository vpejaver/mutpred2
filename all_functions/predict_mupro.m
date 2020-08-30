% Vikas Pejaver
% October 2013

% Function that predicts change in stability through the MuPro
% method, given a sequence and a cell array of mutations

function [ri, ddg] = predict_mupro(sequence, mutations)

% Constants and defaults
%class_model = 'model_class_final.mat';
%reg_model = 'model_regr_final.mat';
ri = zeros(length(mutations), 1);
ddg = zeros(length(mutations), 1);
window = 7;
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY';

global G_RI ALPHA_RI BIAS_RI SV_RI G_DDG ALPHA_DDG BIAS_DDG SV_DDG;

% Make feature matrix
for i = 1:length(mutations)
    wild = mutations{i}(1);
    mut = mutations{i}(end);
    position = str2num(mutations{i}(2:end-1));

    % Make vector for mutated position
    mutated = zeros(1, length(AMINO_ACIDS));
    mutated(AMINO_ACIDS == wild) = -1;
    mutated(AMINO_ACIDS == mut) = 1;

    % Get left and right fragments
    bound = (window - 1)/2;
    left_margin = [max(1, position-bound) : position-1];
    right_margin = [position+1 : min(position+bound, length(sequence))];

    left_flank = sequence(left_margin);
    right_flank = sequence(right_margin);

    left_bin = [];
    for j = 1:length(left_flank)
        left_bin = [left_bin ismember(AMINO_ACIDS, left_flank(j))];
    end

    right_bin = [];
    for j = 1:length(right_flank)
        right_bin = [right_bin ismember(AMINO_ACIDS, right_flank(j))];
    end

    %%% BUG FIX for positions at the termini %%%
    if length(left_flank) < bound
        padding = zeros(1, (bound - length(left_flank)) * length(AMINO_ACIDS));
	left_bin = [padding left_bin];
    end
    if length(right_flank) < bound
        padding = zeros(1, (bound - length(right_flank)) * length(AMINO_ACIDS));
	right_bin = [right_bin padding];
    end
    %%%%%

    % Combine
    X(i, :) = [mutated left_bin right_bin];
end

% Predict classification model
%load(class_model);
ri = rbf_svm_predict(X, ALPHA_RI, SV_RI, G_RI, BIAS_RI);

% Predict regression model
%load(reg_model);
ddg = rbf_svm_predict(X, ALPHA_DDG, SV_DDG, G_DDG, BIAS_DDG);

return;
