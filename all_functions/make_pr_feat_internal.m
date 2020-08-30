% Vikas Pejaver
% November 2012

% Function to make features from predictors internally for MutPred2

function [vect_pr, features, feat_wins] = make_pr_feat_internal(seq, residues, predictions, windows_pr, F)

%windows_pr = [1 7 11 21];

vect_pr = [];
features = {};
feat_wins = [];
for residue = 1:length(residues)
    tmp = [];
    for w = 1 : length(windows_pr)
        % determine margins for the windows (important for the regions near ends
	left_margin = max(1, residues(residue) - (windows_pr(w) - 1) / 2);
	right_margin = min(residues(residue) + (windows_pr(w) - 1) / 2, length(seq));
    
	% average predictions for different models
	for j = 1 : length(predictions)
	    tmp = [tmp mean(predictions{j}(left_margin : right_margin)) std(predictions{j}(left_margin : right_margin)) max(predictions{j}(left_margin : right_margin))];
	    if residue == 1
	        features = [features strcat(F{j}, '_mean') strcat(F{j}, '_std') strcat(F{j}, '_max')];
		feat_wins = [feat_wins windows_pr(w) windows_pr(w) windows_pr(w)];
	    end
	end  
    end
    vect_pr(residue, :) = tmp;
end		

return
