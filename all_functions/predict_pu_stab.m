% Vikas Pejaver
% June 2016

% Function to predict using PU models

function fx = predict_pu_stab(X)

% Global variables
global PU_MODELS;

% Constants and defaults
fx = [];
%model_mat = sprintf('~/updated_mutpred2/src/retrain_pu/model_%s_080116.mat', prop_name);
%model_mat = sprintf('~/updated_mutpred2/src/retrain_pu/model_Stability_081416.mat');


% Load model
%load(model_mat);
%PU_MODELS(end).this_prop

% Loop through each model
test_prediction = [];
for i = 1:length(PU_MODELS(end).nn_model)
    X_ts = X(:, PU_MODELS(end).features{i});
    [meanv, stdv, X_ts] = normalize(X_ts, PU_MODELS(end).mus{i}, PU_MODELS(end).sigmas{i});
    %test_prediction = [test_prediction; sim(rsa_model{i}, X_ts')];
    test_prediction = [test_prediction predict_nn(PU_MODELS(end).nn_model{i}, X_ts)];
end
fx = mean(test_prediction, 2);

return
