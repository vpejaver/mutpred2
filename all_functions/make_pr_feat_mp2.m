% Vikas Pejaver
% March 2014

% Function to run all relevant predictors on a given the wild and
% mutant predictions, positions and length of sequence. The
% function also calculates gain and loss of properties over a
% window and creates a feature matrix based of these

function [feat_mat, locations] = make_pr_feat_mp2(s_wild, s_mut, positions, len)

% Constants and defaults
window = 11;
feat_mat = [];
locations = [];

% Define window over which predictions will be used
bound = (window - 1)/2;

% Loop
for i = 1:length(positions)
    position = positions(i);
    left_margin = [max(1, position-bound) : position-1];
    right_margin = [position+1 : min(position+bound, len)];
    left_flank = [];
    right_flank = [];
    if length(left_margin) < bound
        left_flank = zeros(1, bound - length(left_margin));
    end
    if length(right_margin) < bound
        right_flank = zeros(1, bound - length(right_margin));
    end
    inds = [left_margin position right_margin];

    % Calculate loss and gain of properties and add four
    % features (wild, mutant, loss and gain)
    this_features = [];
    this_loc = [];
    for j = 1:length(s_wild)
        [L, G] = calculate_gain_loss(s_wild{j}(inds), s_mut{j}(inds));
        [loss_max, I1] = max(L);
        [gain_max, I2] = max(G);
        [~, m] = max([loss_max, gain_max]);
	if m == 1
	    this_features = [this_features s_wild{j}(inds(I1)) s_mut{j}(inds(I1)) L(I1) G(I1)];
	    this_loc = [this_loc inds(I1)];
	else
	    this_features = [this_features s_wild{j}(inds(I2)) s_mut{j}(inds(I2)) L(I2) G(I2)];
	    this_loc = [this_loc inds(I2)];
	end
    end

    % Add to main matrix
    feat_mat = [feat_mat; this_features];
    locations = [locations; this_loc];
end

return
