% Vikas Pejaver
% October 2014

% Function to make the feature matrix for a set of mutations from a single sequence
% Version - alpha

function [data_matrix, position_matrix, property_matrix, position_matrix_pu, motif_list, model_markers, remarks] = make_feature_matrix(sequence, mutations)

%%%%%% Global variables %%%%%%
global GRANTHAM;
global AMINO_ACIDS;
global SKIP_PSIBLAST;
global PREDICT_CONS;
global DO_HOMOLOGY;
global CURRDIR;
global PRIOR_PROPS ALPHAS_BETAS RATIOS;


%%%%%% Constants and defaults %%%%%%
if DO_HOMOLOGY
    Num_features = 1345;
else
    Num_features = 1325;
end
Num_props = 57;
Num_remarks = 9;
Win_aa = [3, 7, 11, 21];
Win_pr = [1, 7, 11, 21];
Breakpoints = [1 : 20 : 120];
Pssm_aas = 'ARNDCQEGHILKMFPSTWYV';
Offset = length(Pssm_aas) - 1;
PU_stab_features = [4, 5, 6, 1];

data_matrix = zeros(length(mutations), Num_features);
property_matrix = zeros(length(mutations), (Num_props * 4) + 1);
position_matrix = zeros(length(mutations), Num_props);
position_matrix_pu = zeros(length(mutations), Num_props);
motif_list = repmat({'None'}, 1, length(mutations));
remarks = zeros(length(mutations), Num_remarks);
model_markers = zeros(length(mutations), 1);


%%%%%% Validate input %%%%%%
[validated, idx, logs] = validate_input(sequence, mutations);
if isempty(validated)
    remarks = [logs zeros(length(mutations), 2)]
    return;
end


%%%%%% Parse mutations %%%%%%
positions = cellfun(@str2num, regexprep(validated, '[ACDEFGHIKLMNPQRSTVWY]', ''));


%%%%%% Sequence-based features %%%%%%
% Get amino acid features
vect_aa = make_aa_feat_internal(sequence, positions, Win_aa);

% Clean sequence for predictions
cleaned = clean_protein(sequence);

% Get hydrophobicity
p{1} = hydrophobicity(cleaned, 1);
pf{1} = 'Hydrophobicity';

% Get hydrophobic moments
p{2} = hydrophobic_moment(cleaned, 11, 100);
pf{2} = 'Hydrophobic_moments_100';
p{3} = hydrophobic_moment(cleaned, 11, 160);
pf{3} = 'Hydrophobic_moments_160';
p{4} = hydrophobic_moment(cleaned, 11, 120);
pf{4} = 'Hydrophobic_moments_120';

% Get Vihinen flexibility
p{5} = vihinen(cleaned);
pf{5} = 'Vihinen_flexibility';

% Get amino acid volumes
p{6} = get_aa_volumes(cleaned);
pf{6} = 'Amino_acid_volume';

% Extract features at positions    
vect_pr = make_pr_feat_internal(cleaned, positions, p, Win_pr, pf);

% Add to main feature matrix
this_matrix = [vect_aa vect_pr];


%%%%%% Mutation-based features %%%%%%
vect_submat = [];
for j = 1:length(validated)
    wild = validated{j}(1);
    mut = validated{j}(end);
    sub_scores = get_substitution_score(wild,  mut);
    vect_submat = [vect_submat; [sub_scores GRANTHAM(AMINO_ACIDS == wild, AMINO_ACIDS == mut)]];
end

% Add binary encoding
pairings = regexprep(validated, '[0-9]+', '');
vect_bin = encode_mutation(pairings);

% Add to main feature matrix
this_matrix = [this_matrix vect_submat vect_bin];


%%%%%% PSSM-based features %%%%%%
% Fetch or generate PSSM
PSSM = get_pssm(sequence)';

pssm_flag = zeros(length(mutations), 1);
if isempty(PSSM) & SKIP_PSIBLAST
    vect_pssm = zeros(length(validated), 168);
    vect_bin = zeros(length(validated), 1);
    pssm_flag = ones(length(mutations), 1);
elseif isempty(PSSM) & ~SKIP_PSIBLAST
    tmp_file = strcat(tempname(pwd), '.faa');
    [~, tmp_head] = fileparts(tmp_file);
    pssm_file = strcat(tmp_head, '.pssm');
    
    writeFASTA({sequence}, {tmp_head}, tmp_file);
    system(sprintf('%s/blast-2.2.18/bin/blastpgp -i %s -d %s/nr062813/nr -Q %s -a 3 -h 0.0001 -j 3 > /dev/null', CURRDIR, tmp_file, CURRDIR, pssm_file));
    PSSM = PSSMread(pssm_file)';
    %system(sprintf('rm %s*', tmp_head));
    delete(tmp_file, pssm_file);

    vect_pssm = make_ev_feat_mp2(sequence, positions, PSSM);
    vect_bin = ones(length(validated), 1);
else
    % Make PSSM features
    vect_pssm = make_ev_feat_mp2(sequence, positions, PSSM);
    vect_bin = ones(length(validated), 1);
end

% Add to main feature matrix
this_matrix = [this_matrix vect_pssm vect_bin];


%%%%%% Alignment-based features %%%%%%
% Intialize
vect_trans = zeros(length(validated), 6);
vect_freq = zeros(length(validated), 12);
vect_cindex = zeros(length(validated), 216);
vect_bin = zeros(length(validated), 3);
aln_flag = ones(length(mutations), 1);
pred_aln_flag = zeros(length(mutations), 1);
model_markers(:, 1) = 1;

% Fetch alignment
[headers, alignment] = get_alignment(sequence);

% Fetch conservation scores
cons_mats = get_cindex(sequence);    

% Check if any of the matrices are empty
if PREDICT_CONS | ~strcmp(alignment, '') 
    % Make transition frequency features
    vect_trans = extrans(sequence, validated);

    if PREDICT_CONS & strcmp(alignment, '') 
        if ~isempty(PSSM)

	    % Predict conservation scores
	    cons_mats = {};
	    [cons_mats{1}, cons_mats{2}, fmat] = predict_conservation(sequence, PSSM');
	
	    % Get column frequency
	    vect_freq = [];
	    for j = 1:length(validated)
	        W = validated{j}(1);
		M = validated{j}(end);
		tmp_freq = [];
		for k = 1:length(Breakpoints)
		    start = Breakpoints(k);
		    finish = Breakpoints(k) + Offset;
		    this_mat = fmat(:, start : finish);
		    tmp_freq = [tmp_freq this_mat(positions(j), Pssm_aas == W) this_mat(positions(j), Pssm_aas == M)];
		end
		vect_freq = [vect_freq; tmp_freq];
	    end
	
	    % Make mutation features
	    vect_cindex = [];
	    vect_bin = [];
	    for j = 1:length(cons_mats)
	        if isempty(cons_mats{j})
		    thisc = zeros(length(validated), 72);
		    if j == 1
		        vect_bin = [vect_bin zeros(length(validated), 3)];
		    end
		else
		    thisc = make_ev_feat_mp2(sequence, positions, cons_mats{j}');
		    if j == 1
		        vect_bin = [vect_bin ones(length(validated), 3)];
		    end
		end
		vect_cindex = [vect_cindex thisc];
	    end
	    aln_flag = zeros(length(mutations), 1);
	    pred_aln_flag = ones(length(mutations), 1);
	end

        % Assign model
      	model_markers(:, 1) = 2;

    elseif ~strcmp(alignment, '')

        % Get column frequency
	vect_freq = get_aafreq_aln(alignment, headers, validated);
	
	% Update positions
	new_positions = [];
	for j = 1:length(positions)
	    new_positions(j, :) = update_position(alignment, headers, positions(j));
	end

	% Make mutation features
	vect_cindex = [];
	vect_bin = [];
	for j = 1:length(cons_mats)
	    if isempty(cons_mats{j})
	        thisc = zeros(length(validated), 72);
		vect_bin = [vect_bin zeros(length(validated), 1)];
	    else
	        thisc = make_ev_feat_mp2(sequence, new_positions(:, j), cons_mats{j}');
		vect_bin = [vect_bin ones(length(validated), 1)];
	    end
	    vect_cindex = [vect_cindex thisc];
	end
	aln_flag = zeros(length(mutations), 1);

    end
end

% Add to main feature matrix
this_matrix = [this_matrix vect_trans vect_freq vect_cindex vect_bin];

% Finalize remarks matrix
remarks = [logs pssm_flag aln_flag pred_aln_flag];


%%%%%% Homology count features if set to on %%%%%%
if DO_HOMOLOGY
    % Get paralogy profiles
    vect_hprofile = get_para_prof(sequence);

    % Add to main feature matrix
    this_matrix = [this_matrix repmat(vect_hprofile, length(validated), 1)];
end


%%%%%% Predicted property-based features %%%%%%
% Run predictors on wild-type sequence
[str_wild_preds, str_wild_pu] = run_struct_predictors(cleaned, PSSM);
[fnc_wild_preds, fnc_wild_pu, ~, motif_wild_names, motif_wild_pos] = run_func_predictors(cleaned, PSSM');

% Loop through each mutation and predict on mutated sequence
vect_str = [];
vect_fnc = [];
vect_str_pu = [];
vect_fnc_pu = [];
vect_stab_pu = [];
vect_bin_str = [];
vect_bin_fnc = [];
vect_loc_str = [];
vect_loc_fnc = [];
vect_loc_str_pu = [];
vect_loc_fnc_pu = [];
for j = 1:length(validated)
    % Mutate
    mut = validated{j}(end);
    mutated = cleaned;
    mutated(positions(j)) = mut;

    % Run predictors on mutated sequence
    [str_mut_preds, str_mut_pu, ~, vb] = run_struct_predictors(mutated, PSSM);
    vect_bin_str = [vect_bin_str; vb];
    [fnc_mut_preds, fnc_mut_pu, ~, motif_mut_names, motif_mut_pos, vb] = run_func_predictors(mutated, PSSM');
    vect_bin_fnc = [vect_bin_fnc; vb];

    % Get motif information
    this_mot = unique([motif_wild_names(motif_wild_pos(:, positions(j)) == 1), motif_mut_names(motif_mut_pos(:, positions(j)) == 1)]);
    if ~isempty(this_mot)
        this_mot = regexprep(this_mot, ',|;', ' -');
        motif_list{j} = regexprep(sprintf('%s; ', this_mot{:}), '; $', '');
    end

    % Make gain-loss features
    [this_str, locs] = make_pr_feat_mp2(str_wild_preds, str_mut_preds, positions(j), length(cleaned));
    vect_str = [vect_str; this_str];
    vect_loc_str = [vect_loc_str; locs];
    
    [this_str_pu, locs] = make_pr_feat_mp2(str_wild_pu, str_mut_pu, positions(j), length(cleaned));
    vect_str_pu = [vect_str_pu; this_str_pu];
    vect_loc_str_pu = [vect_loc_str_pu; locs];    

    [this_fnc, locs] = make_pr_feat_mp2(fnc_wild_preds, fnc_mut_preds, positions(j), length(cleaned));
    vect_fnc = [vect_fnc; this_fnc];
    vect_loc_fnc = [vect_loc_fnc; locs];
    
    [this_fnc_pu, locs] = make_pr_feat_mp2(fnc_wild_pu, fnc_mut_pu, positions(j), length(cleaned));
    vect_fnc_pu = [vect_fnc_pu; this_fnc_pu];
    vect_loc_fnc_pu = [vect_loc_fnc_pu; locs];    

    % Pull features for stability (PU)
    tmp_stab = [];
    for k = 1:length(PU_stab_features)
        K = PU_stab_features(k);
	tmp_stab = [tmp_stab, str_wild_preds{K}(positions(j)) str_mut_preds{K}(positions(j))];
    end
    vect_stab_pu = [vect_stab_pu; tmp_stab];
end

% Add stability
[vect_ri, vect_ddg] = predict_mupro(cleaned, validated);
vect_stab_pu = [vect_ri, vect_ddg, vect_stab_pu];
vect_ri_pu = predict_pu_stab(vect_stab_pu);
pidx = find(strcmp(PRIOR_PROPS, 'Stability'));
vect_ri_pu = transform_alphamax(vect_ri_pu, RATIOS(pidx), ALPHAS_BETAS(pidx, [1, 3]), 'noisy', 1);
vect_bin_str = [vect_bin_str ones(length(validated), 2)];

% Add to main feature matrix
this_matrix = [this_matrix vect_str vect_ri vect_ddg vect_bin_str vect_fnc vect_bin_fnc];

% Update matrices
data_matrix(idx, :) = this_matrix;
position_matrix(idx, :) = [vect_loc_str vect_loc_fnc];
property_matrix(idx, :) = [vect_str_pu vect_fnc_pu vect_ri_pu]; 
position_matrix_pu(idx, :) = [vect_loc_str_pu vect_loc_fnc_pu];

return
