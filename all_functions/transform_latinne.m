% Vikas Pejaver
% June 2016

% Function to transform prediction scores using Latinne's formula -
% based on code that Jose gave me

function adjusted_posterior = transform_latinne(predictions, training_prior, true_prior)

% Constants and defaults
training_priors = [(1 - training_prior) training_prior];
true_priors = [(1 - true_prior) true_prior];

% Do adjustment for class priors
for i = 1 : length(predictions)
    adjusted_posterior(:, i) = (true_priors(2)/training_priors(2) * predictions(i)) / sum((true_priors ./ training_priors) .* [1 - predictions(i) predictions(i)]);
end

%prior = (1.0/length(predictions)) * sum(adjusted_posterior);

return
