% Vikas Pejaver
% October 2013

% Function to make predictions given models from SVMLight - written
% mainly to use MuPro models
% Takes in the prediction set matrix and the model parameters -
% alphaY, support vectors, gamma and bias

function scores = rbf_svm_predict(X, ay, sv, g, b);

% Constants and defaults
scores = [];

% Predict
scores = (exp(-g * (pdist2(sv, X) .^ 2))' * ay) - b;

return
