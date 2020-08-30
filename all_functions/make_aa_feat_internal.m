% Vikas Pejaver
% November 2012

% Function to make amino acid features internally for MutPred2

function [vect_aa, features, feat_wins] = make_aa_feat_internal(seq, residues, windows_aa)

% Constants
AminoAcids = 'ACDEFGHIKLMNPQRSTVWY';
%windows_aa = [3 7 11 21];
vect_aa = [];
features = {};
feat_wins = [];

for residue = 1:length(residues)
    for w = 1 : length(windows_aa)
        % determine margins for the windows (important
	% for the regions near ends
	left_margin = max(1, residues(residue) - (windows_aa(w) - 1) / 2);
	right_margin = min(residues(residue) + (windows_aa(w) - 1) / 2,length(seq));
	fragment = seq(left_margin : right_margin);
    
	% 20 aa relative frequencies
	for j = 1 : 20
	    vect_aa(residue, 28 * (w - 1) + j) = length(find(fragment == AminoAcids(j))) / length(fragment);
	    features{28 * (w - 1) + j} = 'AA_frequency';
	    feat_wins(28 * (w - 1) + j) = windows_aa(w);
	end
	% entropy and generalized entropy
	vect_aa(residue, 28 * (w - 1) + 21) = entropy(vect_aa(residue, 28 * (w - 1) + 1 : 28 * (w - 1) + 20));
	features{28 * (w - 1) + 21} = 'Entropy';
	feat_wins(28 * (w - 1) + 21) = windows_aa(w);

	vect_aa(residue, 28 * (w - 1) + 22) = beta_entropy(vect_aa(residue, 28 * (w - 1) + 1 : 28 * (w - 1) + 20), 1.25);
	features{28 * (w - 1) + 22} = 'Beta_entropy_1.25';
	feat_wins(28 * (w - 1) + 22) = windows_aa(w);
	
	vect_aa(residue, 28 * (w - 1) + 23) = beta_entropy(vect_aa(residue, 28 * (w - 1) + 1 : 28 * (w - 1) + 20), 1.50);
	features{28 * (w - 1) + 23} = 'Beta_entropy_1.50';
	feat_wins(28 * (w - 1) + 23) = windows_aa(w);

	vect_aa(residue, 28 * (w - 1) + 24) = beta_entropy(vect_aa(residue, 28 * (w - 1) + 1 : 28 * (w - 1) + 20), 1.75);
	features{28 * (w - 1) + 24} = 'Beta_entropy_1.75';
	feat_wins(28 * (w - 1) + 24) = windows_aa(w);
	
	% net and total charge
	pos_charge = length(find(fragment == 'K' | fragment == 'R'));
	neg_charge = length(find(fragment == 'D' | fragment == 'E'));
	vect_aa(residue, 28 * (w - 1) + 25) = (pos_charge - neg_charge)/ length(fragment);
	features{28 * (w - 1) + 25} = 'Net_charge';
        feat_wins(28 * (w - 1) + 25) = windows_aa(w);

	vect_aa(residue, 28 * (w - 1) + 26) = (pos_charge + neg_charge)/ length(fragment);
	features{28 * (w - 1) + 26} = 'Total_charge';
	feat_wins(28 * (w - 1) + 26) = windows_aa(w);

	% aromatics
	vect_aa(residue, 28 * (w - 1) + 27) = length(find(fragment == 'F' | fragment == 'Y' | fragment == 'W')) / length(fragment);
	features{28 * (w - 1) + 27} = 'Aromatics';
	feat_wins(28 * (w - 1) + 27) = windows_aa(w);
	
	% charge-hydrophobicity ratio
	vect_aa(residue, 28 * (w - 1) + 28) = charge_hydrophobicity(fragment);
	features{28 * (w - 1) + 28} = 'Charge_hydrophobicity';
	feat_wins(28 * (w - 1) + 28) = windows_aa(w);
    end
end
return
