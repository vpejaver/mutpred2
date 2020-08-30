% Vikas Pejaver
% February 2014

% Function to predict topology and status of signal and TM helices,
% given sequence and PSSM
% Categories = {'N-terminal signal', 'Signal helix', 'C-terminal signal', 'Cleavage site', 'Inside', 'Membrane', 'Outside', 'Negative'};

function scores = sigtm_pred(sequence, PSSM)

% Constants and defaults
scores = [];
win_aa = [3 7 11 21];
win_pr = [1 7 11 21];
win_ev = [1 7 11 21];

% Clean protein
seq = clean_protein(sequence);
residues = [1:length(seq)];

% Get amino acid features
vect_aa = make_aa_feat_internal(seq, residues, win_aa);

% Get predicted features
preds = use_predictors(seq, [], 0);
vect_pr = make_pr_feat_internal(seq, residues, preds, win_pr, repmat({'Dummy'}, 1, length(preds)));

% Get evolutionary features
vect_ev = make_ev_feat_internal(seq, residues, PSSM, win_ev);

% Combine features and add to feature matrix
X = [vect_aa vect_pr vect_ev];

% Predict from features
scores = predict_sigtm_feat(X);

return
