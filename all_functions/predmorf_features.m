% Vikas Rao Pejaver
% March 2014

% Function that predicts MoRFs from features

function posteriors = predmorf_features(X)

% Constants and defaults
%model_mat = 'morf_model_032614.mat';
global MORF_MODEL;

% Load the model
%load(model_mat, 'MORF_MODEL');

% Get posterior probabilities
posteriors = zeros(1, size(X, 1));
for b = 1:length(MORF_MODEL)
    posteriors = posteriors + eval(MORF_MODEL{b}, X)';
end
posteriors = posteriors/length(MORF_MODEL);

return
