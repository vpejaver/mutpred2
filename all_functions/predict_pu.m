% Vikas Pejaver
% June 2016

% Function to predict using PU models

function fx = predict_pu(X, q)

% Global variables
global PU_MODELS

% Constants and defaults
fx = [];
%model_mat = sprintf('~/updated_mutpred2/src/retrain_pu/model_%s_061616.mat', prop_name);

% Load model
%load(model_mat);

% Set up model
features = PU_MODELS(q).features;
mus = PU_MODELS(q).mus;
sigmas = PU_MODELS(q).sigmas;
projections = PU_MODELS(q).projections;
nn_model = PU_MODELS(q).nn_model;
%PU_MODELS(q).this_prop

% Loop through each model
test_prediction = [];
for i = 1:length(nn_model)
    X_ts = X(:, features{i});
    [meanv, stdv, X_ts] = normalize(X_ts, mus{i}, sigmas{i});
    X_ts = X_ts * projections{i};
    test_prediction = [test_prediction predict_nn(nn_model{i}, X_ts)];
end
fx = mean(test_prediction, 2);

return
