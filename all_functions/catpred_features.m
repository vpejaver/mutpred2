% Vikas Pejaver
% March 2015

% Function to predict catalytic residues from feature matrix

function fx = catpred_features(X)

% Global variables
global CAT_MODELS CAT_PROJS CAT_MUS CAT_SIGMAS;

% Constants and defaults
fx = [];
%model_mat = 'cat_model_040115.mat';

% Load model
%load(model_mat);

% Loop through each model
test_prediction = [];
for i = 1:length(CAT_MODELS)
    X_ts = X; %X(:, cat_feats{i});
    [~, ~, X_ts] = normalize(X_ts, CAT_MUS{i}, CAT_SIGMAS{i});
    X_ts = X_ts * CAT_PROJS{i};
    X_ts = [ones(size(X_ts, 1), 1) X_ts];
    preds(:, i) = logsig(X_ts * CAT_MODELS{i});
end

% Predict
fx = mean(preds, 2);

return
