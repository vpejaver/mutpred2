% Vikas Pejaver
% November 2012

% Function to make evolutionary features given sequence, positions
% and a data matrix (either PSSMs or conservation indices)

function [vect_ev, feat_wins] = make_ev_feat_mp2(seq, residues, evmat)

% Constants
windows_ev = [1 5 11 21];

vect_ev = [];
feat_wins = [];
for residue = 1:length(residues)
    tmp = [];
    for w = 1:length(windows_ev)
        % determine margins for the windows (important
	% for the regions near ends
	left_margin = max(1, residues(residue) - (windows_ev(w) - 1) / 2);
	right_margin = min(residues(residue) + (windows_ev(w) - 1) / 2, length(seq));
    
	for j = 1 : size(evmat, 1)
	    tmp = [tmp mean(evmat(j, left_margin : right_margin))];
	    if residue == 1
	        feat_wins = [feat_wins windows_ev(w)];
	    end
	end
    end
    vect_ev(residue, :) = tmp;
end
return
