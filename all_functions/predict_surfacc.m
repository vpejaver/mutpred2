% Vikas Pejaver
% February 2014

% Function to predict surface accessibility (2-class) from sequence

function scores = predict_surfacc(sequence, PSSM)

% Constants and defaults
scores = [];
bin_win = 7;
win_aa = [3 7 11 21];
win_pr = [1 7 11 21];
win_ev = [1 3 11 21];

% Make amino acid features
%vect_aa = make_features_ptm_aa(sequence, [1:length(sequence)]);
vect_aa = make_aa_feat_internal(sequence, [1:length(sequence)], win_aa);

% Make binary features
bin_vec = [];
padding = repmat('*', 1, floor(bin_win/2));
tmp_seq = [padding sequence padding];
for j = 1:length(sequence)  
    start = j;
    finish = j + bin_win - 1;
    peptide = tmp_seq(start : finish);
    bin_vec = [bin_vec; my_make_spacer(peptide)];
end
vect_aa = [vect_aa bin_vec];	

% Make evolutionary features
%vect_ev = make_features_rsa_ev(sequence, [1:length(sequence)], PSSM');
vect_ev = make_ev_feat_internal(sequence, [1:length(sequence)], PSSM', win_ev);

% Clean protein
protein = clean_protein(sequence);

% Make predicted features
preds = rsa_run_predictors(protein);
%vect_pr = make_features_ptm_pr(protein, [1:length(protein)], preds);
vect_pr = make_pr_feat_internal(sequence, [1:length(sequence)], preds, win_pr, repmat({'Dummy'}, 1, length(preds)));

% Combine and predict
scores = predict_surfacc_features([vect_aa vect_ev vect_pr]);

return;
