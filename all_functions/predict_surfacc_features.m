% Vikas Pejaver
% February 2014

% Function to predict RSA from feature matrix

function fx = predict_surfacc_features(X)

% Constants and defaults
fx = [];
global RSA_MODELS RSA_PROJS RSA_MEANS RSA_STDS RSA_F;

%global rsa_models rsa_projs rsa_means rsa_stds rsa_f;
%model_mat = 'rsa_nnemodel_030914.mat';

% Load model
%load(model_mat);

% Loop through each model
test_prediction = [];
for i = 1:length(RSA_MODELS)
    X_ts = X(:, RSA_F{i});
    [meanv, stdv, X_ts] = normalize(X_ts, RSA_MEANS{i}, RSA_STDS{i});
    X_ts = X_ts * RSA_PROJS{i};
    %test_prediction = [test_prediction; sim(rsa_model{i}, X_ts')];
    test_prediction = [test_prediction predict_nn(RSA_MODELS{i}, X_ts)];
end
fx = mean(test_prediction, 2);

return
