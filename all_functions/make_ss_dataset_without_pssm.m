function [D] = make_ss_dataset_without_pssm(seq)
% Args:
%   sequence: a string of length n
% 
% Last modified: Wed Feb 23 17:05:31 2011

AminoAcids = 'ACDEFGHIKLMNPQRSTVWY';
windows_aa = [1 3 7 11 15 21];  % amino acid compositions
windows_ev = [1 3 7 11 15 21];  % evolutionary information
D = [];

vect_aa = [];
vect_ev = [];
for residue = 1 : length(seq)
    % features made from amino acid sequence
    for w = 1 : length(windows_aa)
        % determine margins for the windows (important for the regions near ends)
        e = (windows_aa(w) - 1) / 2;
        left_margin = max(1, residue - e);
        right_margin = min(residue + e, length(seq));
        fragment = seq(left_margin : right_margin);
        
        % 20 aa relative frequencies
        for j = 1 : 20
            vect_aa(residue, 28 * (w - 1) + j) = length(find(fragment == AminoAcids(j))) / ...
                length(fragment);
        end
        % entropy and generalized entropy
        vect_aa(residue, 28 * (w - 1) + 21) = entropy(vect_aa(residue, ...
                                                          28 * (w - 1) + 1 : 28 * (w - 1) + 20));
        vect_aa(residue, 28 * (w - 1) + 22) = beta_entropy(vect_aa(residue, ...
                                                          28 * (w - 1) + 1 : 28 * (w - 1) + 20), 1.25);
        vect_aa(residue, 28 * (w - 1) + 23) = beta_entropy(vect_aa(residue, ...
                                                          28 * (w - 1) + 1 : 28 * (w - 1) + 20), 1.50);
        vect_aa(residue, 28 * (w - 1) + 24) = beta_entropy(vect_aa(residue, ...
                                                          28 * (w - 1) + 1 : 28 * (w - 1) + 20), 1.75);
        
        % net and total charge
        pos_charge = length(find(fragment == 'K' | fragment == 'R'));
        neg_charge = length(find(fragment == 'D' | fragment == 'E'));
        vect_aa(residue, 28 * (w - 1) + 25) = (pos_charge - neg_charge) / length(fragment);
        vect_aa(residue, 28 * (w - 1) + 26) = (pos_charge + neg_charge) / length(fragment);
        
        % aromatics
        vect_aa(residue, 28 * (w - 1) + 27) = length(find(fragment == 'F' | fragment == 'Y' | ...
                                                          fragment == 'W')) / length(fragment);
        
        % charge-hydrophobicity ratio
        vect_aa(residue, 28 * (w - 1) + 28) = charge_hydrophobicity(fragment);
    end
    
    % tmp = [];
    % for w = 1 : length(windows_ev)
    %     % determine margins for the windows (important for the regions near ends)
    %     e = (windows_ev(w) - 1) / 2;
    %     left_margin = max(1, residue - e);
    %     right_margin = min(residue + e, length(seq));
        
    %     for j = [1 : 20] % 'only' 20 evolutionary attributes are used
    %         tmp = [tmp mean(PSSM(j, left_margin : right_margin))];
    %     end
    % end
    % vect_ev(residue, :) = tmp;          
end

D = [vect_aa vect_ev];


return