function vect_ev = make_features_ptm_ev(seq,residues,pssm)

windows_ev = [1 3 11 21];

vect_ev = [];
for residue = 1:length(residues)
  tmp = [];
  for w = 1 : length(windows_ev)
    % determine margins for the windows (important
    % for the regions near ends
    left_margin = max(1, residues(residue) - (windows_ev(w) - 1) / 2);
    right_margin = min(residues(residue) + (windows_ev(w) - 1) / 2, length(seq));
    
    %for j = [1 : 20 41] % 'only' 21 evolutionary attributes are used
    for j = 1 : 41 % Last column excluded
      tmp = [tmp mean(pssm(j, left_margin : right_margin))];
    end
  end
  vect_ev(residue, :) = tmp;
end		
return