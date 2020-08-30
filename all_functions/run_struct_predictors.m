% Vikas Pejaver
% February 2014

% Function to run the different structural predictors for MutPred on a
% sequence

function [stp, stpu, features, status] = run_struct_predictors(seq, pssm)

% Global variables
global PRIOR_PROPS ALPHAS_BETAS RATIOS;
global PU_MODELS;

% Constants and defaults
stp = {};
stpu = {};
features = {};
status = ones(1, 15);

% Get VSL2B disorder
stp{1} = VSL2B(seq)';
features{1} = 'VSL2B_disorder';
idx = find(strcmp(PRIOR_PROPS, features{1}));
stpu{1} = transform_alphamax(stp{1}, RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);

% Get B factors
stp{2} = predictBfactors(seq, 3, 5)';
features{2} = 'B_factor';
idx = find(strcmp(PRIOR_PROPS, features{2}));
stpu{2} = transform_latinne(stp{2}, RATIOS(idx), ALPHAS_BETAS(idx, 1));

% Get combined structural predictions
if ~isempty(pssm)
    [surfacc, sigtm, fsurfacc, fsigtm] = run_int_struct_predictors(seq, pssm');
    idx = find(strcmp(PRIOR_PROPS, 'Surface_accessibility'));
    surfaccu = predict_pu(fsurfacc, idx);
    surfaccu = transform_alphamax(surfaccu, RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);

    idx = find(strcmp(PRIOR_PROPS, 'N-terminal_signal'));
    sigtmu(:, 1) = predict_pu(fsigtm, idx);
    sigtmu(:, 1) = transform_alphamax(sigtmu(:, 1), RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);

    idx = find(strcmp(PRIOR_PROPS, 'Signal_helix'));
    sigtmu(:, 2) = predict_pu(fsigtm, idx);
    sigtmu(:, 2) = transform_alphamax(sigtmu(:, 2), RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);

    idx = find(strcmp(PRIOR_PROPS, 'C-terminal_signal'));
    sigtmu(:, 3) = predict_pu(fsigtm, idx);
    sigtmu(:, 3) = transform_alphamax(sigtmu(:, 3), RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);

    idx = find(strcmp(PRIOR_PROPS, 'Signal_cleavage'));
    sigtmu(:, 4) = predict_pu(fsigtm, idx);
    sigtmu(:, 4) = transform_alphamax(sigtmu(:, 4), RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);

    idx = find(strcmp(PRIOR_PROPS, 'Intracellular_loop'));
    sigtmu(:, 5) = predict_pu(fsigtm, idx);
    sigtmu(:, 5) = transform_alphamax(sigtmu(:, 5), RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);

    idx = find(strcmp(PRIOR_PROPS, 'Transmembrane_region'));
    sigtmu(:, 6) = predict_pu(fsigtm, idx);
    sigtmu(:, 6) = transform_alphamax(sigtmu(:, 6), RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);

    idx = find(strcmp(PRIOR_PROPS, 'Extracellular_loop'));
    sigtmu(:, 7) = predict_pu(fsigtm, idx);
    sigtmu(:, 7) = transform_alphamax(sigtmu(:, 7), RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);

    idx = find(strcmp(PRIOR_PROPS, 'Non_transmembrane'));
    sigtmu(:, 8) = predict_pu(fsigtm, idx);
    sigtmu(:, 8) = transform_alphamax(sigtmu(:, 8), RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);
else
    surfaccu = zeros(length(seq), 1);
    sigtmu = zeros(length(seq), 8);
    surfacc = zeros(length(seq), 1);
    sigtm = zeros(length(seq), 8);
    status(3) = 0;
    status(7:14) = 0;
end

% Get surface accessibility
%stp{3} = predict_surfacc(seq, pssm')';
stp{3} = surfacc';
stpu{3} = surfaccu';
features{3} = 'Surface_accessibility';

% Get secondary structure (ignore disorder)
if ~isempty(pssm)
    ss = ssp4(seq, pssm);
else
    ss = zeros(length(seq), 4);
    status(4:6) = 0;
end
stp{4} = ss(:, 1)';
stp{5} = ss(:, 2)';
stp{6} = ss(:, 3)';
features{4} = 'Helix';
features{5} = 'Strand';
features{6} = 'Loop';
idx = find(strcmp(PRIOR_PROPS, features{4}));
stpu{4} = transform_latinne(stp{4}, RATIOS(idx), ALPHAS_BETAS(idx, 1));
idx = find(strcmp(PRIOR_PROPS, features{5}));
stpu{5} = transform_latinne(stp{5}, RATIOS(idx), ALPHAS_BETAS(idx, 1));
idx = find(strcmp(PRIOR_PROPS, features{6}));
stpu{6} = transform_latinne(stp{6}, RATIOS(idx), ALPHAS_BETAS(idx, 1));

% Get signal and transmembrane regions
%sigtm = sigtm_pred(seq, pssm);
stp{7} = sigtm(:, 1)';
stp{8} = sigtm(:, 2)';
stp{9} = sigtm(:, 3)';
stp{10} = sigtm(:, 4)';
stp{11} = sigtm(:, 5)';
stp{12} = sigtm(:, 6)';
stp{13} = sigtm(:, 7)';
stp{14} = sigtm(:, 8)';
stpu{7} = sigtmu(:, 1)';
stpu{8} = sigtmu(:, 2)';
stpu{9} = sigtmu(:, 3)';
stpu{10} = sigtmu(:, 4)';
stpu{11} = sigtmu(:, 5)';
stpu{12} = sigtmu(:, 6)';
stpu{13} = sigtmu(:, 7)';
stpu{14} = sigtmu(:, 8)';
features{7} = 'N-terminal_signal';
features{8} = 'Signal_helix';
features{9} = 'C-terminal_signal';
features{10} = 'Signal_cleavage';
features{11} = 'Intracellular_loop';
features{12} = 'Transmembrane_region';
features{13} = 'Extracellular_loop';
features{14} = 'Non_transmembrane';

% Get coiled-coil regions
stp{15} = ccpred(seq, 0)';
features{15} = 'Coiled_coil';
idx = find(strcmp(PRIOR_PROPS, features{15}));
stpu{15} = transform_alphamax(stp{15}, RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);

return
