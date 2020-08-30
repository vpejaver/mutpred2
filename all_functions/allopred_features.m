% Vikas Rao Pejaver
% April 2014

% Function that predicts allosteric sites from features

function posteriors = allopred_features(X)

% Constants and defaults
%model_mat = 'allorf_model_042814.mat';
global ALLO_MODEL;

% Load the model
%load(model_mat, 'allo_model');

% Get posterior probabilities
posteriors = zeros(1, size(X, 1));
for b = 1:length(ALLO_MODEL)
    posteriors = posteriors + eval(ALLO_MODEL{b}, X)';
end
posteriors = posteriors/length(ALLO_MODEL);

return
