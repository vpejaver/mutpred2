% Vikas Pejaver
% February 2014

% Function to predict SIGTM from features

function fx = predict_sigtm_feat(X)

% Constants and defaults
fx = [];
global SIGTM_MODELS SIGTM_PROJS SIGTM_MEANS SIGTM_STDS;

%global sigtm_models sigtm_projs sigtm_means sigtm_stds;
%model_mat = 'sig_tm_model_030914.mat';

% Load model
%load(model_mat, 'models', 'mus', 'sigmas', 'projections');

% Loop through each model
test_prediction = [];
for i = 1:length(SIGTM_MODELS)
    [meanv, stdv, X_ts] = normalize(X, SIGTM_MEANS{i}, SIGTM_STDS{i});
    X_ts = X_ts * SIGTM_PROJS{i};
    test_prediction(:, :, i) = predict_nn(SIGTM_MODELS{i}, X_ts);
end

% Average predictions
fx = mean(test_prediction, 3);

return
