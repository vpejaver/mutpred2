function prediction = get_residue_prediction (sequence)

global CURRDIR;

% selected features according to the t-test with p-value threshold 0.01
features = [1 2 3 4 6 7 8 9 10 12 13 14 15 16 17 19 20 21 22 23 24 25 26 27 28 29 30 ...
        31 33 34 36 37 38 39 40 41 43 44 45 46 47 49 50 51 52 53 54 55 56 57 58 59 60 ...
        61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 ...
        87 88 89 90 91 92];

% construct sample for the first predictor
D = make_sample_stage_1(sequence);
% select features
D = D(:, features);

load(strcat(CURRDIR, filesep, 'all_models', filesep, 'CaMPredRes.mat'));

N = length(CaMPredRes.mn1);
for n = 1 : N
    [tmp1, tmp2, Dn] = normalize(D, CaMPredRes.mn1{n}, CaMPredRes.st1{n});
    Dn = Dn * CaMPredRes.T{n};
    [tmp1, tmp2, Dn] = normalize(Dn, CaMPredRes.mn2{n}, CaMPredRes.st2{n});
    Dn = [Dn ones(size(Dn, 1), 1)];
    tmp_prediction(:, n) = Dn * CaMPredRes.beta{n};
    tmp_prediction(:, n) = 1 ./ (1 + exp(-tmp_prediction(:, n)));
end

prediction = mean(tmp_prediction(:, 1 : N), 2)';
prediction = moving_average(prediction, 13);

return



    

    
