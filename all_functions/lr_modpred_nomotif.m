% Vikas Rao Pejaver
% February 2012

% Function that loads the appropriate random forest model and runs the data matrix through it
% Input is the data matrix, PTM type, motif and pssm flag
% If motif is '.', assumes the PTM doesn't have motif and uses regular model, otherwise uses a model containing 'nomotif'

function [posteriors, confidences] = lr_modpred_nomotif(X, ptm, m, pflag)

global CURRDIR;

% Define which model file to use
if m == '.'
    model_name = ptm;
    if pflag == 1
        load(strcat(CURRDIR, filesep, 'all_models', filesep, 'thresholds_PSSM_lr.mat'));
        model_mat = strcat(CURRDIR, filesep, 'all_models', filesep, ptm, '_', 'PSSM.mat');
    else
        load(strcat(CURRDIR, filesep, 'all_models', filesep, 'thresholds_noPSSM_lr.mat'));
        model_mat = strcat(CURRDIR, filesep, 'all_models', filesep, ptm, '_', 'noPSSM.mat');
    end
else
    model_name = strcat(ptm, '_nomotif');
    if pflag == 1
        load(strcat(CURRDIR, filesep, 'all_models', filesep, 'thresholds_PSSM_lr.mat'));
        model_mat = strcat(CURRDIR, filesep, 'all_models', filesep, ptm, '_', 'nomotif_PSSM.mat');
    else
        load(strcat(CURRDIR, filesep, 'all_models', filesep, 'thresholds_noPSSM_lr.mat'));
        model_mat = strcat(CURRDIR, filesep, 'all_models', filesep, ptm, '_', 'nomotif_noPSSM.mat');
    end
end

% Load the model
load(model_mat);

% Get posterior probabilities
posteriors = zeros(1, size(X, 1));
for b = 1 : length(lr_model)
    test = X(:, to_keep{b});
    norm_X = (test - repmat(mu{b}, size(test, 1), 1)) ./ repmat(sigma{b}, size(test, 1), 1);
    test = norm_X * proj_mats{b};
    posteriors = posteriors + glmval(lr_model{b}, test, 'logit')';
end
posteriors = posteriors/length(lr_model);

% Set confidence values based on probability cutoffs (1 - Low
% confidence, 2 - Medium confidence and 3 - High confidence)
confidences = zeros(1, size(X, 1));
thr = cutoffs(find(strcmp(model_list, model_name)), :);
rounded = round(posteriors * 100) / 100; % Round posteriors to nearest 2nd decimal place (fixes 'boundary' bug)
confidences(find(rounded >= thr(1) & rounded < thr(2))) = 1;
confidences(find(rounded >= thr(2) & rounded < thr(3))) = 2;
confidences(find(rounded >= thr(3))) = 3;

return
