% Vikas Pejaver
% February 2014

% Function to predict RSA from feature matrix

function fx = predict_ppi_features(X)

% Constants and defaults
fx = [];
%model_mat = 'rsa_nnmodel_021014.mat';
%model_mat = 'ppbind_model_042814.mat';

global PPI_MODELS PPI_PROJS PPI_MEANS PPI_STDS PPI_F;

% Load model
%load(model_mat);

% Loop through each model
test_prediction = [];
%for i = 1:length(ppi_models)
%    X_ts = X(:, features{i});
%    [meanv, stdv, X_ts] = normalize(X_ts, mus{i}, sigmas{i});
%    X_ts = X_ts * projections{i};
%    test_prediction = [test_prediction predict_nn(ppi_models{i}, X_ts)];
%end
for i = 1:length(PPI_MODELS)
    X_ts = X(:, PPI_F{i});
    [meanv, stdv, X_ts] = normalize(X_ts, PPI_MEANS{i}, PPI_STDS{i});
    X_ts = X_ts * PPI_PROJS{i};
    test_prediction = [test_prediction predict_nn(PPI_MODELS{i}, X_ts)];
end

fx = mean(test_prediction, 2);

return
