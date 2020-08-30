% Vikas Pejaver
% October 2014

% Function to assign p-values to each hypothesis
% NOTE: Reorders the features to have stability at the end

function [comb_scores, comb_types, comb_pvals] = generate_hypotheses(P, option)

%%%%%% Global variable %%%%%%
global NEUTRAL_PROPERTIES PROPERTY_NAMES;
global NEUTRAL_PROPERTIES_PU PROPERTY_NAMES_PU;
global NAMES_REORDERED;


%%%%%% Constants and defaults %%%%%%
comb_scores = [];
comb_types = [];
comb_pvals = [];


%%%%%% Set mode %%%%%%
if strcmp(option, 'pu')
    this_props = NEUTRAL_PROPERTIES_PU;
    this_names = PROPERTY_NAMES_PU;
else
    this_props = NEUTRAL_PROPERTIES;
    this_names = PROPERTY_NAMES;
end
N = size(this_props);


%%%%%% Sanity check %%%%%%
if size(P, 2) ~= N
    error('The number of predicted properties does not match!');
end


%%%%%% Separate out features %%%%%% 
loss_inds = find(~cellfun(@isempty, regexp(this_names, '_(loss)', 'match')));
gain_inds = find(~cellfun(@isempty, regexp(this_names, '_(gain)', 'match')));

predicted_losses = P(:, loss_inds);
loss_neutrals = this_props(:, loss_inds);

predicted_gains = P(:, gain_inds);
gain_neutrals = this_props(:, gain_inds);
NAMES_REORDERED = strrep(this_names(loss_inds), '_loss', '');


%%%%%% Take care of stability %%%%%%
if strcmp(option, 'pu')
    stab_idx = find(~cellfun(@isempty, regexp(this_names, 'Stability', 'match')));

    predicted_stabs = P(:, stab_idx);
    stab_losses = this_props(:, stab_idx);
    stab_gains = stab_losses;
else
    stab_idx = find(~cellfun(@isempty, regexp(this_names, 'Stability_reliability', 'match')));

    predicted_stabs = P(:, stab_idx);
    stab_losses = this_props(this_props(:, stab_idx) < 0, stab_idx);
    stab_gains = this_props(this_props(:, stab_idx) >= 0, stab_idx);
end
predicted_losses = [predicted_losses predicted_stabs];
predicted_gains = [predicted_gains predicted_stabs];

NAMES_REORDERED = [NAMES_REORDERED; 'Stability'];


%%%%%% Assign P-values %%%%%%
Nsl = length(stab_losses);
Nsg = length(stab_gains);
for i = 1:size(P, 1)
    for j = 1:size(loss_neutrals, 2)
        p_vals_loss(i, j) = length(find(loss_neutrals(:, j) >= predicted_losses(i, j))) / N(1);
    end
    for j = 1:size(gain_neutrals, 2)
        p_vals_gain(i, j) = length(find(gain_neutrals(:, j) >= predicted_gains(i, j))) / N(1);
    end
    if strcmp(option, 'pu')
        p_vals_loss(i, j+1) = length(find(stab_losses >= predicted_stabs(i))) / Nsl;
    else
        p_vals_loss(i, j+1) = length(find(stab_losses < predicted_stabs(i))) / Nsl;
    end
    p_vals_gain(i, j+1) = length(find(stab_gains >= predicted_stabs(i))) / Nsg;

    % Combine scores (assign loss or gain)
    this_pvals = [p_vals_loss(i, :); p_vals_gain(i, :)];
    this_scores = [predicted_losses(i, :); predicted_gains(i, :)];
    [comb_scores(i, :), type_inds] = max(this_scores);
    
    comb_types(i, :) = type_inds;
    comb_pvals(i, type_inds == 1) = this_pvals(1, type_inds == 1);
    comb_pvals(i, type_inds == 2) = this_pvals(2, type_inds == 2);
end
return
