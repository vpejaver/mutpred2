% Vikas Pejaver
% January 2015

% Function to predict frequencies from feature matrix

function fx = predict_freqs_features(X)

global CURRDIR;

% Constants and defaults
fx = [];
%model_mat = 'cindex_freqs_model_122714.mat';
model_mat = strcat(CURRDIR, filesep, 'all_models', filesep, 'cindex_freqs_model_012015.mat');

%global ppi_models ppi_projs ppi_means ppi_stds ppi_f;

% Load model
load(model_mat);

% Loop through each model
test_prediction = [];
for i = 1:length(freq_models)
    [~, ~, X_ts] = normalize(X, mus{i}, sigmas{i});
    test_prediction(:, :, i) = predict_nn(freq_models{i}, X_ts) .* repmat(sig_target{i}, size(X_ts, 1), 1) + repmat(mu_target{i}, size(X_ts, 1), 1);
end

fx = mean(test_prediction, 3);

return
