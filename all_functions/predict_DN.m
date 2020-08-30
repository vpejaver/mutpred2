% Vikas Pejaver
% July 2015

% Function to predict disease/neutral based on feature set - gives
% the main MutPred2 scores

function fx = predict_DN(X, markers)

%%%%%% Global variables %%%%%%
global NN_MODELS NN_MUS NN_PROJS NN_SELECTED NN_STDS;
global PNN_MODELS PNN_MUS PNN_PROJS PNN_SELECTED PNN_STDS;


%%%%%% Constants and defaults %%%%%%
fx = zeros(size(X, 1), 1);


%%%%%% Separate out rows %%%%%%
inds1 = find(markers == 1);
X1 = X(inds1, :);
inds2 = find(markers == 2);
X2 = X(inds2, :);


%%%%%% Predict %%%%%%
if ~isempty(inds1)
    test_prediction = [];
    for i = 1:length(NN_MODELS)
        X_ts = X1(:, NN_SELECTED{i});
	[meanv, stdv, X_ts] = normalize(X_ts, NN_MUS{i}, NN_STDS{i});
	X_ts = X_ts * NN_PROJS{i};
	test_prediction = [test_prediction predict_nn(NN_MODELS{i}, X_ts)];
    end
    fx(inds1) = mean(test_prediction, 2);
end
if ~isempty(inds2)
    test_prediction = [];
    for i = 1:length(PNN_MODELS)
        X_ts = X2(:, PNN_SELECTED{i});
	[meanv, stdv, X_ts] = normalize(X_ts, PNN_MUS{i}, PNN_STDS{i});
	X_ts = X_ts * PNN_PROJS{i};
	test_prediction = [test_prediction predict_nn(PNN_MODELS{i}, X_ts)];
    end
    fx(inds2) = mean(test_prediction, 2);
end

return
