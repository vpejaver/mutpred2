% Vikas Pejaver
% April 2017

% The main function for MutPred2
% Run this to generate predictions on a FASTA file
% Version - final (added posterior-based ranking,
% ontology output format and ANNOVAR intake)

function mutpred2(varargin)

%%%%%% Global variables %%%%%%
% Re-initialize
clearvars -global

% Command-line arguments
global IN_FILE OID DO_HOMOLOGY PREDICT_CONS SKIP_PSIBLAST CSV_FORM ONTO_FORM PTHR; % SCORES_ONLY TOPK;

% Data from files
global GRANTHAM;
global RSA_MODELS RSA_PROJS RSA_MEANS RSA_STDS RSA_F;
global SIGTM_MODELS SIGTM_PROJS SIGTM_MEANS SIGTM_STDS;
global EST_EMISSION EST_TRANSITION HIDDEN_STATES OBSERVED_STATES;
global G_RI ALPHA_RI BIAS_RI SV_RI G_DDG ALPHA_DDG BIAS_DDG SV_DDG;
global CAT_MODELS CAT_PROJS CAT_MUS CAT_SIGMAS;
global DBIND_MODEL DBIND_PROJ DBIND_MEAN DBIND_STD DBIND_F;
global RBIND_MODELS RBIND_PROJS RBIND_MUS RBIND_SIGMAS RBIND_FEATS;
global PPI_MODELS PPI_PROJS PPI_MEANS PPI_STDS PPI_F;
global HOTSPOT_MODELS HOTSPOT_PROJS HOTSPOT_MUS HOTSPOT_SIGMAS HOTSPOT_FEATS;
global MORF_MODEL;
global ALLO_MODEL;
global MOTIFS MOTIF_DESCS MOTIF_NAMES MOTIF_PATTS;
global PRIOR_PROPS ALPHAS_BETAS RATIOS;
global G NODE_LABELS;
global PU_MODELS;
global NEUTRAL_PROPERTIES NN_MODELS NN_MUS NN_PROJS NN_SELECTED NN_STDS PROPERTY_NAMES PROP_INDS;
global PNN_MODELS PNN_MUS PNN_PROJS PNN_SELECTED PNN_STDS;

% Other constants
global AMINO_ACIDS;
global CURRDIR;
global MINDS;
global NAMES_REORDERED;


%%%%%% Constants and defaults %%%%%%
OID = 1;
DO_HOMOLOGY = 0;
PREDICT_CONS = 1;
SKIP_PSIBLAST = 0;
CSV_FORM = 1;
ONTO_FORM = 1;
%SCORES_ONLY = 0;
%TOPK = 5;
PTHR = 0.05;
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY';
DISP_ITER = 100;
tokens = regexp(fileparts(mfilename('fullpath')), filesep, 'split');
%CURRDIR = strjoin(tokens(1:end-2), filesep);
CURRDIR = strjoin(tokens, filesep);
Batch_size = 4;


%%%%%% Parse command-line %%%%%%
parse_status = parse_args(varargin);
if ~parse_status
    return;
end


%%%%%% Set warnings off if std. output %%%%%%
if OID == 1
    warning('off', 'all');
end


%%%%%% Load dependency files %%%%%%
if OID ~= 1
    fprintf(1, 'Initializing ...\n');
end
get_data();


%%%%%% Read in file and pre-process %%%%%%
[ids, substitutions, sequences] = preprocess();


%%%%%% Process %%%%%%
if OID ~= 1
    fprintf(1, 'Predicting ...\n');
end
% Set display iterations
if length(sequences) <= 10
    DISP_ITER = 1;
elseif length(sequences) <= 100
    DISP_ITER = 10;
elseif length(sequences) <= 1000
    DISP_ITER = 100;
else
    DISP_ITER = 1000;
end

% Print header in output file
fprintf(OID, 'ID,Substitution,MutPred2 score,Molecular mechanisms with Pr >= 0.01 and P < %.2f,Motif information,Remarks\n', PTHR);

% Loop through and predict
for i = 1:length(sequences)

    % Make features
    [feats, positions, propX_pu, positions_pu, motif_info, models, notes] = make_feature_matrix(sequences{i}, substitutions{i});

    % Predict pathogenicity
    S = predict_DN(feats, models);

    % Get molecular mechanisms
    if CSV_FORM == 1 | ONTO_FORM == 1
	[prop_scores_pu, prop_types_pu, prop_pvals_pu] = generate_hypotheses(propX_pu, 'pu');
    else
        prop_scores_pu = [];
	prop_types_pu = [];
	prop_pvals_pu = [];
	loss_inds = find(~cellfun(@isempty, regexp(PROPERTY_NAMES, '_(loss)', 'match')));
	NAMES_REORDERED = strrep(PROPERTY_NAMES(loss_inds), '_loss', '');
	NAMES_REORDERED = [NAMES_REORDERED; 'Stability'];
    end

    % Print output
    print_output(ids{i}, sequences{i}, substitutions{i}, notes, positions_pu, S, prop_scores_pu, prop_types_pu, prop_pvals_pu, motif_info)

    % Update user
    if OID ~= 1
        if ~mod(length(sequences), DISP_ITER)
	    fprintf(1, 'Finished %d of %d sequences\n', i, length(sequences));
	end
    end
end
if OID ~= 1
    fclose(OID);
end

return
