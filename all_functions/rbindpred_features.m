% Vikas Pejaver
% May 2014

% Function to predict RNA-binding sites from feature matrix

function fx = rbindpred_features(X)

% Constants and defaults
fx = [];
%model_mat = 'rbind_model_052714.mat';
global RBIND_MODELS RBIND_PROJS RBIND_MUS RBIND_SIGMAS RBIND_FEATS;

% Load model
%load(model_mat);

% Loop through each model
test_prediction = [];
for i = 1:length(RBIND_MODELS)
    X_ts = X(:, RBIND_FEATS{i});
    [meanv, stdv, X_ts] = normalize(X_ts, RBIND_MUS{i}, RBIND_SIGMAS{i});
    X_ts = X_ts * RBIND_PROJS{i};
    test_prediction = [test_prediction predict_nn(RBIND_MODELS{i}, X_ts)];
end
fx = mean(test_prediction, 2);

return
