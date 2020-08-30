% Vikas Pejaver
% May 2014

% Function to predict DNA-binding from feature matrix

function fx = predict_dbind_features(X)

% Constants and defaults
fx = [];
%model_mat = 'dbind_model_043014.mat';

global DBIND_MODEL DBIND_PROJ DBIND_MEAN DBIND_STD DBIND_F;

% Load model
%load(model_mat);

% Transform
%X_ts = X(:, features);
%[~, ~, X_ts] = normalize(X_ts, mu, sig);
%X_ts = X_ts * proj;
%X_ts = [ones(size(X_ts, 1), 1) X_ts];
X_ts = X(:, DBIND_F);
[~, ~, X_ts] = normalize(X_ts, DBIND_MEAN, DBIND_STD);
X_ts = X_ts * DBIND_PROJ;
X_ts = [ones(size(X_ts, 1), 1) X_ts];


% Predict
fx = logsig(X_ts * DBIND_MODEL);

return
