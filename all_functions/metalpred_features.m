% Vikas Pejaver
% May 2014

% Function to predict metal-binding sites from a feature matrix

function scores = metalpred_features(X)

global CURRDIR;

% Constants and defaults
scores = zeros(size(X, 1), 11);
model_names = strcat(CURRDIR, filesep, 'all_models', filesep, {'cadmium_model_052614.mat', 'calcium_model_050514.mat', 'cobalt_model_050514.mat', 'copper_model_052614.mat', 'iron_model_052614.mat', 'magnesium_model_050514.mat', 'manganese_model_052614.mat', 'nickel_model_050514.mat', 'potassium_model_050514.mat', 'sodium_model_050514.mat', 'zinc_model_050514.mat'});
rf_flag = [0 1 1 0 0 1 0 1 1 1 1];

% Make all predictions
for i = 1:length(model_names)

    % load MAT file
    load(model_names{i});
    
    % predict based on whether random forest or neural network
    if rf_flag(i)
        posteriors = zeros(size(X, 1), 1);
        for b = 1:length(metal_model)
	    posteriors = posteriors + eval(metal_model{b}, X);
	end
	posteriors = posteriors/length(metal_model);
    else
        posteriors = [];
        for b = 1:length(metal_model)
	    X_ts = X(:, features{b});
	    [meanv, stdv, X_ts] = normalize(X_ts, mus{b}, sigmas{b});
	    X_ts = X_ts * projections{b};
	    posteriors = [posteriors predict_nn(metal_model{b}, X_ts)];
	end
	posteriors = mean(posteriors, 2);
    end

    % Add to final matrix
    scores(:, i) = posteriors;%/length(metal_model);

end

return
