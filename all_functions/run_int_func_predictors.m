% Vikas Pejaver
% April 2015

% Function to predict catalytic residue, DNA-binding, RNA-binding, PPI, hotspot,
% metal-binding, ligand-binding, hotspot and PTM sites

function [catres, dbind, rbind, ppi, hs, allo, morf, metal, X_common, X_hs, X_morf, X_metal] = run_int_func_predictors(sequence, PSSM)

% Constants and defaults
scores = [];
bin_win = 7;
win_aa = [3 7 11 21]; % dna, rna, ppi, hs, allo, morf, metal
win_pr = [1 7 11 21]; % dna, rna, ppi, hs, allo, morf, metal
win_ev = [1 3 11 21]; % dna, rna, ppi, hs, allo, morf, metal
%win_ev_sigtm = [1 7 11 21];
hotspot_feats = [1:28, 113:116, 117:120, 121:124, 125, 126, 127, 128, 129, 130, 131:133, 134:136, 137:139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 341:382];

% Set residues
residues = [1:length(sequence)];

% Make amino acid features
vect_aa = make_aa_feat_internal(sequence, residues, win_aa);

% Make evolutionary features
vect_ev = make_ev_feat_internal(sequence, residues, PSSM', win_ev);

% Clean protein
%protein = clean_protein(sequence);

% Get predicted features
preds = use_predictors_func(sequence, PSSM', 0);
vect_pr = make_pr_feat_internal(sequence, residues, preds, win_pr, repmat({'Dummy'}, 1, length(preds)));

% Make binary features
bin_vec_mms = [];
bin_vec_mes = [];
padding = repmat('*', 1, floor(bin_win/2));
tmp_seq = [padding sequence padding];
for j = 1:length(sequence)  
    start = j;
    finish = j + bin_win - 1;
    peptide = tmp_seq(start : finish);
    bin_vec_mms = [bin_vec_mms; my_make_spacer(peptide)];
    bin_vec_mes = [bin_vec_mes; make_example_spacer(peptide([1:3 5:7]))];
end
metal_bin_pos = [1:63 85:size(bin_vec_mms, 2)];
bin_vec_metal = bin_vec_mms(:, metal_bin_pos);

% Combine and predict
X_common = [vect_aa vect_pr vect_ev bin_vec_mes];
catres = catpred_features(X_common);
dbind = predict_dbind_features(X_common);
rbind = rbindpred_features(X_common);
ppi = predict_ppi_features(X_common);
allo = allopred_features(X_common)';
X_metal = [vect_aa vect_pr vect_ev bin_vec_metal];
metal = metalpred_features(X_metal);
X_morf = [vect_aa vect_ev vect_pr bin_vec_mms];
morf = predmorf_features(X_morf)';

% Finish the hotspot prediction
preds{end + 1} = ppi;
mutations = strcat(cellstr(sequence'), arrayfun(@num2str, residues, 'uniformoutput', false)', 'A');
[ri, ddg] = predict_mupro(sequence, mutations);
preds{end + 1} = ri';
preds{end + 1} = ddg';
vect_pr_hs = make_pr_feat_internal(sequence, residues, preds, win_pr, repmat({'Dummy'}, 1, length(preds)));
X_hs = [vect_aa vect_pr_hs vect_ev];
X_hs = X_hs(:, hotspot_feats);
hs = predict_hotspot_features(X_hs);

return;
