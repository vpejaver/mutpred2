% Vikas Pejaver
% May 2014

% Function to predict hotspot residues from feature matrix

function fx = predict_hotspot_features(X)

% Constants and defaults
fx = [];
%model_mat = 'hotspot_logreg_052814.mat';
global HOTSPOT_MODELS HOTSPOT_PROJS HOTSPOT_MUS HOTSPOT_SIGMAS HOTSPOT_FEATS;

% Load model
%load(model_mat);

% Transform
for i = 1:length(HOTSPOT_MODELS)
    X_ts = X(:, HOTSPOT_FEATS{i});
    [~, ~, X_ts] = normalize(X_ts, HOTSPOT_MUS{i}, HOTSPOT_SIGMAS{i});
    X_ts = X_ts * HOTSPOT_PROJS{i};
    X_ts = [ones(size(X_ts, 1), 1) X_ts];
    preds(:, i) = logsig(X_ts * HOTSPOT_MODELS{i});
end

% Predict
fx = mean(preds, 2);

return
