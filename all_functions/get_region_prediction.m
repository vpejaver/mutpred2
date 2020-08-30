function reg_pred = get_region_prediction (res_pred, sequence)

% construct dataset for the second stage
[D region_prediction rg_pred] = get_region_dataset(res_pred, sequence);

% if a residue based predictor did not predict any CaM binding residues
if isempty(D)
    reg_pred = region_prediction;
    return
end

% otherwise, load the predictor
load CaMPredReg.mat

% features after feature selection are hard coded
features = [5 6 7 8 13 14 15 16 21 22 23 24 29 30 31 32 37 38 39 ...
        40 61 62 63 64 69 70 71 72 86 87 90];

% do feature selection
D = D(:, features);

% predict using 15 LR classifiers
N = length(CaMPredReg.mn1);
for n = 1 : N
    [tmp1, tmp2, Dn] = normalize(D, CaMPredReg.mn1{n}, CaMPredReg.st1{n});
    Dn = Dn * CaMPredReg.T{n};
    [tmp1, tmp2, Dn] = normalize(Dn, CaMPredReg.mn2{n}, CaMPredReg.st2{n});
    Dn = [Dn ones(size(Dn, 1), 1)];
    tmp_prediction(:, n) = Dn * CaMPredReg.beta{n};
    tmp_prediction(:, n) = 1 ./ (1 + exp(-tmp_prediction(:, n)));
end

% average prediction
prediction = mean(tmp_prediction(:, 1 : N), 2);

% assign regions
reg_pred = zeros(1, length(region_prediction));

% this is order (0) vs. Swiss-Prot prediction (1)
protein_prediction_score = get_sequence_prediction(sequence);

% set the threshold depending on the type of proteins
if protein_prediction_score == 1
    threshold = 0.1;
else
    threshold = 0.4;
end

% set the region prediction (a region is  predicted as CaMBT if its score is above threshold)
ln = length(rg_pred) / 2;
for i = 1 : ln
    if prediction(i) > threshold
        reg_pred(rg_pred(2 * i - 1) : rg_pred(2 * i)) = 1;
    end
end

return
