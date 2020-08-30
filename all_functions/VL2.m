function [p1, p2, p3, p4] = VL2 (sequence)

% function [p1, p2, p3, p4] = VL2 (sequence)
% returns VL2 predictions
% p1 is general predictor
% p2-p4 are flavors

global CURRDIR;

load(strcat(CURRDIR, filesep, 'all_models', filesep, 'models.mat'));

pfilt = predict_linear(sequence, models);

p1 = pfilt(:, 4)';
p2 = pfilt(:, 1)';
p3 = pfilt(:, 2)';
p4 = pfilt(:, 3)';

return