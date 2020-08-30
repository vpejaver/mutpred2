function [prediction] = predictBfactors (sequence, Win, Wout)

global CURRDIR;

load(strcat(CURRDIR, filesep, 'all_models', filesep, 'linpred_ns.mat'));
load(strcat(CURRDIR, filesep, 'all_models', filesep, 'FHC.mat'));

y{1} = sequence; % have to do this, because of make_attribute()
x = make_attribute(y, Win, FHC);
x = x(:, 1 : 24);
x = x(:, [1 : 19, 21 : 24]);

data = x(:, [1 : 18 20 23]);

n_features = size(data, 2);

[meanv, stdv, data(:, 1 : n_features)] = normalize(data(:, 1 : n_features), meanv, stdv);

% add additional column of ones so that LR can be done
data = [data(:, 1 : n_features) ones(size(data, 1), 1)];

% get prediction on test data
raw_prediction = data(:, 1 : n_features + 1) * results.beta;

% transform results to 0-1 interval (logsig function)
raw_prediction = 1 ./ (1 + exp(-raw_prediction));

prediction = moving_average(raw_prediction, Wout);

return