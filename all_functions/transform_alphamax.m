% Vikas Pejaver
% June 2016

% Function to transform raw scores to posteriors using Shantanu's
% transformation formulae

function posteriors = transform_alphamax(tau, C, ab, method, padding)

% Constants and defaults
posteriors = [];

% Transform
alpha = ab(1);
int_trans = C * (tau ./ (1 - tau));

% Correct for noise
if strcmp(method, 'noisy')
    beta = ab(2);
    posteriors = ((alpha * (1 - alpha)) / (beta - alpha)) * (int_trans - ((1 - beta) / (1 - alpha)));
else
    posteriors = alpha * int_trans;
end

% Pad values between zero and one
if padding
    posteriors(posteriors < 0) = 0;
    posteriors(posteriors > 1) = 1;
end

return
