% Vikas Pejaver
% October 2012

% Function to calculate loss and gain for a feature, given two
% vectors (wild and mutant in that order) - Basically calculates
% P(loss) = P(p|sw).(1 - P(p|sm)) and P(gain) = P(p|sm).(1 -
% P(p|sw))

function [loss, gain] = calculate_gain_loss(wild, mutant)

% Constants and defaults
loss = zeros(1, length(wild));
gain = zeros(1, length(mutant));

% Check error
if length(wild) ~= length(mutant)
    error('Number of scores from wild and mutant are not the same!')
end

% Get loss
loss = wild .* (1 - mutant);
gain = (1 - wild) .* mutant;

return
