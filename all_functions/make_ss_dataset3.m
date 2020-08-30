function [D] = make_ss_dataset3(seq, PSSM)
% Create PSSM and AA frequency features for input proteins. This feature set combines
% Pedja's and those from psipred and works best among all tried three feature sets.
% 
% Args:
%   sequence: a string of length n
%   PSSM: a 42-by-n matrix
% 
% Last modified: Wed Feb 23 17:05:49 2011

AminoAcids = 'ACDEFGHIKLMNPQRSTVWY';
windows_aa = 15;                      % amino acid compositions
windows_ev = 15;                      % evolutionary information
D = [];

vect_aa = [];
vect_ev = [];
n = length(seq);
for residue = 1 : n
    % features made from amino acid sequence
    % determine margins for the windows (important for the regions near ends)
    e = (windows_aa - 1) / 2;
    left_margin = max(1, residue - e);
    right_margin = min(residue + e, n);
    fragment = seq(left_margin : right_margin);
        
    % 20 aa relative frequencies
    for j = 1 : 20
        vect_aa(residue, j) = length(find(fragment == AminoAcids(j))) / length(fragment);
    end

    vect_aa(residue, 21) = entropy(vect_aa(residue, 1 : 20));
    vect_aa(residue, 22) = beta_entropy(vect_aa(residue, 1 : 20), 1.25);
    vect_aa(residue, 23) = beta_entropy(vect_aa(residue, 1 : 20), 1.50);
    vect_aa(residue, 24) = beta_entropy(vect_aa(residue, 1 : 20), 1.75);
        
    % net and total charge
    pos_charge = length(find(fragment == 'K' | fragment == 'R'));
    neg_charge = length(find(fragment == 'D' | fragment == 'E'));
    vect_aa(residue, 25) = (pos_charge - neg_charge) / length(fragment);
    vect_aa(residue, 26) = (pos_charge + neg_charge) / length(fragment);
        
    % aromatics
    vect_aa(residue, 27) = length(find(fragment == 'F' | fragment == 'Y' | ...
                                       fragment == 'W')) / length(fragment);
        
    % charge-hydrophobicity ratio
    vect_aa(residue, 28) = charge_hydrophobicity(fragment);
    
    % 'only' 20 evolutionary attributes are used
    tmp = zeros(1, 20 * windows_ev);
    % determine margins for the windows (important for the regions near ends)
    e = (windows_ev + 1) / 2;
    for j = 1:windows_ev
        rel_ix = j + residue - e;
        if rel_ix < 1 || rel_ix > n
            continue
        end
        tmp(1, (1+(j-1)*20):(j*20)) = PSSM(1:20, rel_ix);
    end
    vect_ev(residue, :) = tmp;
end

D = [vect_aa vect_ev];

return