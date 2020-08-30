% Vikas Pejaver
% March 2014

% Function to run the different functional predictors for MutPred on a
% sequence

function [fup, fupu, features, motinfo, motpos, status] = run_func_predictors(seq, pssm)

% Global variables
global PRIOR_PROPS ALPHAS_BETAS RATIOS;
global PU_MODELS;

% Constants and defaults
fup = {};
fupu = {};
features = {};
motinfo = {};
motpos = [];
status = ones(1, 42);
modnames = {'Acetylation', 'ADP-ribosylation', 'Amidation_motif', 'Amidation_nomotif', 'C-linked_glycosylation', 'Carboxylation', 'Disulfide_linkage', 'Farnesylation', 'Geranylgeranylation', 'GPI_anchor_amidation', 'Hydroxylation', 'Methylation_K', 'Methylation_R', 'Myristoylation', 'N-linked_glycosylation_motif', 'N-linked_glycosylation_nomotif', 'N-terminal_acetylation_A', 'N-terminal_acetylation_G', 'N-terminal_acetylation_M', 'N-terminal_acetylation_S', 'N-terminal_acetylation_T', 'O-linked_glycosylation_S', 'O-linked_glycosylation_T', 'Palmitoylation', 'Phosphorylation_S', 'Phosphorylation_T', 'Phosphorylation_Y', 'Proteolytic_cleavage', 'Pyrrolidone_carboxylic_acid', 'Sulfation', 'SUMOylation_motif', 'SUMOylation_nomotif', 'Ubiquitylation'};
motifs = {'.', '.', '(.G[RK][RK])', '(.G[RK][RK])', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '(N[^P][ST][^P])', '(N[^P][ST][^P])', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '([AVFPILM]K.(E|D))', '([AVFPILM]K.(E|D))', '.'};
modres = {'K', 'ER', 'ACDEFGHIKLMNPQRSTVWY', 'ACDEFGHIKLMNPQRSTVWY', 'W', 'E', 'C', 'C', 'C', 'N', 'KPY', 'K', 'R', 'G', 'N', 'N', 'A', 'G', 'M', 'S', 'T', 'S', 'T', 'C', 'S', 'T', 'Y', 'ACDEFGHIKLMNPQRSTVWY', 'Q', 'Y', 'K', 'K', 'K'};
offsets = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0];
umods = unique(regexprep(modnames, '(_(motif|nomotif|.))$', ''), 'stable');

if ~isempty(pssm)
    [catres, dnabind, rnabind, pint, hotspot, allosite, morfres, mbind, fcommon, fhotspot, fmorf, fmetal] = run_int_func_predictors(seq, pssm);
    idx = find(strcmp(PRIOR_PROPS, 'Catalytic_site'));
    catresu = predict_pu(fcommon, idx);
    catresu = transform_alphamax(catresu, RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);
    
    idx = find(strcmp(PRIOR_PROPS, 'DNA_binding'));
    dnabindu = predict_pu(fcommon, idx);
    dnabindu = transform_alphamax(dnabindu, RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);

    idx = find(strcmp(PRIOR_PROPS, 'RNA_binding'));
    rnabindu = predict_pu(fcommon, idx);
    rnabindu = transform_alphamax(rnabindu, RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);

    idx = find(strcmp(PRIOR_PROPS, 'PPI_residue'));
    pintu = predict_pu(fcommon, idx);
    pintu = transform_alphamax(pintu, RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);

    idx = find(strcmp(PRIOR_PROPS, 'PPI_hotspot'));
    hotspotu = predict_pu(fhotspot, idx);
    hotspotu = transform_alphamax(hotspotu, RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);

    idx = find(strcmp(PRIOR_PROPS, 'Allosteric_site'));
    allositeu = predict_pu(fcommon, idx);
    allositeu = transform_alphamax(allositeu, RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);

    idx = find(strcmp(PRIOR_PROPS, 'MoRF'));
    morfresu = predict_pu(fmorf, idx);
    morfresu = transform_alphamax(morfresu, RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);

    idx = find(strcmp(PRIOR_PROPS, 'Cadmium_binding'));
    mbindu(:, 1) = predict_pu(fmetal, idx);
    mbindu(:, 1) = transform_alphamax(mbindu(:, 1), RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);

    idx = find(strcmp(PRIOR_PROPS, 'Calcium_binding'));
    mbindu(:, 2) = predict_pu(fmetal, idx);
    mbindu(:, 2) = transform_alphamax(mbindu(:, 2), RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);

    idx = find(strcmp(PRIOR_PROPS, 'Cobalt_binding'));
    mbindu(:, 3) = predict_pu(fmetal, idx);
    mbindu(:, 3) = transform_alphamax(mbindu(:, 3), RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);

    idx = find(strcmp(PRIOR_PROPS, 'Copper_binding'));
    mbindu(:, 4) = predict_pu(fmetal, idx);
    mbindu(:, 4) = transform_alphamax(mbindu(:, 4), RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);

    idx = find(strcmp(PRIOR_PROPS, 'Iron_binding'));
    mbindu(:, 5) = predict_pu(fmetal, idx);
    mbindu(:, 5) = transform_alphamax(mbindu(:, 5), RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);

    idx = find(strcmp(PRIOR_PROPS, 'Magnesium_binding'));
    mbindu(:, 6) = predict_pu(fmetal, idx);
    mbindu(:, 6) = transform_alphamax(mbindu(:, 6), RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);

    idx = find(strcmp(PRIOR_PROPS, 'Manganese_binding'));
    mbindu(:, 7) = predict_pu(fmetal, idx);
    mbindu(:, 7) = transform_alphamax(mbindu(:, 7), RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);

    idx = find(strcmp(PRIOR_PROPS, 'Nickel_binding'));
    mbindu(:, 8) = predict_pu(fmetal, idx);
    mbindu(:, 8) = transform_alphamax(mbindu(:, 8), RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);

    idx = find(strcmp(PRIOR_PROPS, 'Potassium_binding'));
    mbindu(:, 9) = predict_pu(fmetal, idx);
    mbindu(:, 9) = transform_alphamax(mbindu(:, 9), RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);

    idx = find(strcmp(PRIOR_PROPS, 'Sodium_binding'));
    mbindu(:, 10) = predict_pu(fmetal, idx);
    mbindu(:, 10) = transform_alphamax(mbindu(:, 10), RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);

    idx = find(strcmp(PRIOR_PROPS, 'Zinc_binding'));
    mbindu(:, 11) = predict_pu(fmetal, idx);
    mbindu(:, 11) = transform_alphamax(mbindu(:, 11), RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);
else
    catres = zeros(length(seq), 1);
    dnabind = zeros(length(seq), 1);
    rnabind = zeros(length(seq), 1);
    pint = zeros(length(seq), 1);
    hotspot = zeros(length(seq), 1);
    allosite = zeros(length(seq), 1);
    morfres = zeros(length(seq), 1);
    mbind = zeros(length(seq), 11);
    catresu = zeros(length(seq), 1);
    dnabindu = zeros(length(seq), 1);
    rnabindu = zeros(length(seq), 1);
    pintu = zeros(length(seq), 1);
    hotspotu = zeros(length(seq), 1);
    allositeu = zeros(length(seq), 1);
    morfresu = zeros(length(seq), 1);
    mbindu = zeros(length(seq), 11);
    status([1, 3:19]) = 0;
end

% Get catalytic residue
fup{1} = catres; %predictCR(seq);
fupu{1} = catresu;
features{1} = 'Catalytic_site';

% Get calmodulin-binding
fup{2} = CaMBTPredictor(seq)';
features{2} = 'Calmodulin_binding';
idx = find(strcmp(PRIOR_PROPS, features{2}));
fupu{2} = transform_alphamax(fup{2}, RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);

% Get DNA-binding
fup{3} = dnabind;
fupu{3} = dnabindu;
features{3} = 'DNA_binding';

% Get RNA-binding
fup{4} = rnabind;
fupu{4} = rnabindu;
features{4} = 'RNA_binding';

% Get protein-protein interaction
fup{5} = pint;
fupu{5} = pintu;
features{5} = 'PPI_residue';

% Get PPI hotspot
fup{6} = hotspot;
fupu{6} = hotspotu;
features{6} = 'PPI_hotspot';

% Get MoRF prediction
fup{7} = morfres;
fupu{7} = morfresu;
features{7} = 'MoRF';

% Get allosteric residues
fup{8} = allosite;
fupu{8} = allositeu;
features{8} = 'Allosteric_site';

% Get metal-binding
fup{9} = mbind(:, 1);
fup{10} = mbind(:, 2);
fup{11} = mbind(:, 3);
fup{12} = mbind(:, 4);
fup{13} = mbind(:, 5);
fup{14} = mbind(:, 6);
fup{15} = mbind(:, 7);
fup{16} = mbind(:, 8);
fup{17} = mbind(:, 9);
fup{18} = mbind(:, 10);
fup{19} = mbind(:, 11);
fupu{9} = mbindu(:, 1);
fupu{10} = mbindu(:, 2);
fupu{11} = mbindu(:, 3);
fupu{12} = mbindu(:, 4);
fupu{13} = mbindu(:, 5);
fupu{14} = mbindu(:, 6);
fupu{15} = mbindu(:, 7);
fupu{16} = mbindu(:, 8);
fupu{17} = mbindu(:, 9);
fupu{18} = mbindu(:, 10);
fupu{19} = mbindu(:, 11);
features = [features 'Cadmium_binding', 'Calcium_binding', 'Cobalt_binding', 'Copper_binding', 'Iron_binding', 'Magnesium_binding', 'Manganese_binding', 'Nickel_binding', 'Potassium_binding', 'Sodium_binding', 'Zinc_binding'];

% Get PTM sites
bin_vec = ones(1, 23);
bin_vec(19) = 0;
if ~isempty(pssm)
    [ptms, ~, ~, ~, ~, ~, fptms] = predict_ptms(seq, pssm', bin_vec);
    for i = 1:length(modnames)
        pos = find(ismember(seq, modres{i}));
        if ~strcmp(motifs{i}, '.')
	    r = regexp(seq, motifs{i}, 'start') + offsets(i);
	    if ~isempty(strfind(modnames{i}, '_motif'))
	        pos = intersect(pos, r);
	    elseif ~isempty(strfind(modnames{i}, '_nomotif'))
	        pos = setdiff(pos, r);
	    end
	end
	ptmsu{i} = zeros(length(seq), 1);
	idx = find(strcmp(PRIOR_PROPS, modnames{i}));
	ptmsu{i}(pos) = predict_pu(fptms(pos, :), idx);
	ptmsu{i} = transform_alphamax(ptmsu{i}, RATIOS(idx), ALPHAS_BETAS(idx, [1, 3]), 'noisy', 1);
    end
    % combine
    for i = 1:length(umods)
        inds = find(~cellfun(@isempty, regexp(modnames, umods{i}, 'match')));
	tmp_ptm{i} = sum([ptmsu{inds}], 2);
    end
    ptmsu = tmp_ptm;
else
    ptms = predict_ptms(seq, bin_vec);
    ptmsu = repmat({zeros(length(seq), 1)}, 1, length(umods))
    status(20:41) = 0;
end
tmp = cellfun(@transpose, num2cell(ptms, 2), 'uniformoutput', false)';
fup = [fup tmp];
fupu = [fupu ptmsu];
features = [features umods];

% Get motifs
[motscores, motids, motnames, motpatts, motpos] = check_motifs(seq);
%motinfo = strcat(motids, ':', motpatts, ':', motnames);
motinfo = strcat(motids, ':', regexprep(motnames, ':', '-'));
%motinfo = strcat(motids, ':', motpatts);
fup{end + 1} = motscores';
fupu{end + 1} = fup{end};
features{end + 1} = 'Motifs';
return
