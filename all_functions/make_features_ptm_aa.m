function vect_aa = make_features_ptm_aa(seq,residues)

AminoAcids = 'ACDEFGHIKLMNPQRSTVWY';
windows_aa = [3 7 11 21];

vect_aa = [];
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
    end
    % entropy and generalized entropy
    vect_aa(residue, 28 * (w - 1) + 21) = entropy(vect_aa(residue, 28 * (w - 1) + 1 : 28 * (w - 1) + 20));
    vect_aa(residue, 28 * (w - 1) + 22) = beta_entropy(vect_aa(residue, 28 * (w - 1) + 1 : 28 * (w - 1) + 20), 1.25);
    vect_aa(residue, 28 * (w - 1) + 23) = beta_entropy(vect_aa(residue, 28 * (w - 1) + 1 : 28 * (w - 1) + 20), 1.50);
    vect_aa(residue, 28 * (w - 1) + 24) = beta_entropy(vect_aa(residue, 28 * (w - 1) + 1 : 28 * (w - 1) + 20), 1.75);
    
    % net and total charge
    pos_charge = length(find(fragment == 'K' | fragment == 'R'));
    neg_charge = length(find(fragment == 'D' | fragment == 'E'));
    vect_aa(residue, 28 * (w - 1) + 25) = (pos_charge - neg_charge)/ length(fragment);
    vect_aa(residue, 28 * (w - 1) + 26) = (pos_charge + neg_charge)/ length(fragment);
    
    % aromatics
    vect_aa(residue, 28 * (w - 1) + 27) = length(find(fragment == 'F' | fragment == 'Y' | fragment == 'W')) / length(fragment);
    
    % charge-hydrophobicity ratio
    vect_aa(residue, 28 * (w - 1) + 28) = charge_hydrophobicity(fragment);
  end
end		
return