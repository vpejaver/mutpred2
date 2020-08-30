% Vikas Pejaver
% March 2014

% Function to predict surface accessibility and signal-transmembrane
% Categories = {'RSA', 'N-terminal signal', 'Signal helix', 'C-terminal signal', 'Cleavage site', 'Inside', 'Membrane', 'Outside', 'Negative'};

function [rsa, sigtm, X_rsa, X_sigtm] = run_int_struct_predictors(sequence, PSSM)

% Constants and defaults
scores = [];
bin_win = 7;
win_aa = [3 7 11 21];
win_pr = [1 7 11 21];
win_ev_rsa = [1 3 11 21];
win_ev_sigtm = [1 7 11 21];

residues = [1:length(sequence)];

% Make amino acid features
vect_aa = make_aa_feat_internal(sequence, residues, win_aa);

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

% Make evolutionary features
vect_ev_rsa = make_ev_feat_internal(sequence, residues, PSSM', win_ev_rsa);
vect_ev_sigtm = make_ev_feat_internal(sequence, residues, PSSM', win_ev_sigtm);

% Clean protein
protein = clean_protein(sequence);

% Get predicted features
preds = use_predictors(protein, [], 0);
vect_pr_rsa = make_pr_feat_internal(protein, residues, preds(1:12), win_pr, repmat({'Dummy'}, 1, length(preds(1:12))));
vect_pr_sigtm = make_pr_feat_internal(protein, residues, preds, win_pr, repmat({'Dummy'}, 1, length(preds)));

% Combine and predict
X_rsa = [vect_aa bin_vec vect_ev_rsa vect_pr_rsa];
X_sigtm = [vect_aa vect_pr_sigtm vect_ev_sigtm];
rsa = predict_surfacc_features(X_rsa);
sigtm = predict_sigtm_feat(X_sigtm);

return;
