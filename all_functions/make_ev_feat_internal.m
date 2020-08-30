% Vikas Pejaver
% November 2012

% Function to make evolutionary features for the internal
% predictors used in MutPred2

function vect_ev = make_ev_feat_internal(seq, residues, pssm, windows_ev)

% Constant
%windows_ev = [1 3 11 21];

vect_ev = [];
for residue = 1:length(residues)
    tmp = [];
    for w = 1 : length(windows_ev)
        % determine margins for the windows (important
	% for the regions near ends
	left_margin = max(1, residues(residue) - (windows_ev(w) - 1) / 2);
	right_margin = min(residues(residue) + (windows_ev(w) - 1) / 2, length(seq));
    
	for j = 1 : 42 % Last column excluded
	    tmp = [tmp mean(pssm(j, left_margin : right_margin))];
	end
    end
    vect_ev(residue, :) = tmp;
end		
return
