function vect_pr = make_features_ptm_pr(seq,residues,predictions)

windows_pr = [1 7 11 21];

vect_pr = [];
for residue = 1:length(residues)
  tmp = [];
  for w = 1 : length(windows_pr)
    % determine margins for the windows (important for the regions near ends
    left_margin = max(1, residues(residue) - (windows_pr(w) - 1) / 2);
    right_margin = min(residues(residue) + (windows_pr(w) - 1) / 2, length(seq));
    
    % average predictions for 18 different models
    for j = 1 : length(predictions)
      tmp = [tmp mean(predictions{j}(left_margin : right_margin)) std(predictions{j}(left_margin : right_margin)) max(predictions{j}(left_margin : right_margin))];
    end  
  end
  vect_pr(residue, :) = tmp;
end		
return