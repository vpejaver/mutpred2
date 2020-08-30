% Vikas Pejaver
% October 2014

% Function to load up all dependency files

function [] = get_data()

%%%%%% Global variables %%%%%%
global CURRDIR;
global DO_HOMOLOGY;
global PREDICT_CONS;
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
global PU_MODELS;
global NEUTRAL_PROPERTIES NN_MODELS NN_MUS NN_PROJS NN_SELECTED NN_STDS PROPERTY_NAMES PROP_INDS;
global NEUTRAL_PROPERTIES_PU PROPERTY_NAMES_PU;
global NODE_LABELS G;
global PNN_MODELS PNN_MUS PNN_PROJS PNN_SELECTED PNN_STDS;
global MINDS


%%%%%% Constants and defaults %%%%%%
Grantham_file = strcat(CURRDIR, filesep, 'data', filesep, 'grantham.mat');
Rsa_file = strcat(CURRDIR, filesep, 'all_models', filesep, 'rsa_nnemodel_030914.mat');
Sigtm_file = strcat(CURRDIR, filesep, 'all_models', filesep, 'sig_tm_model_030914.mat');
Cc_file = strcat(CURRDIR, filesep, 'all_models', filesep, 'marcoil_like_072313.mat');
Stab_cfile = strcat(CURRDIR, filesep, 'all_models', filesep, 'model_class_final.mat');
Stab_rfile = strcat(CURRDIR, filesep, 'all_models', filesep, 'model_regr_final.mat');
Cat_file = strcat(CURRDIR, filesep, 'all_models', filesep, 'cat_model_040115.mat');
Dbind_file = strcat(CURRDIR, filesep, 'all_models', filesep, 'dbind_model_043014.mat');
Rbind_file = strcat(CURRDIR, filesep, 'all_models', filesep, 'rbind_model_052714.mat');
Ppi_file = strcat(CURRDIR, filesep, 'all_models', filesep, 'ppbind_model_042814.mat');
Hotspot_file = strcat(CURRDIR, filesep, 'all_models', filesep, 'hotspot_logreg_052814.mat');
%Morf_file = strcat(CURRDIR, filesep, 'predmorf', filesep, 'morf_model_032614.mat');
Morf_file = strcat(CURRDIR, filesep, 'all_models', filesep, 'morf_model_030815.mat');
Allo_file = strcat(CURRDIR, filesep, 'all_models', filesep, 'allorf_model_042814.mat');
Motif_file = strcat(CURRDIR, filesep, 'data', filesep, 'elm_0214_prosite_0713.mat');
Priors_file = strcat(CURRDIR, filesep, 'data', filesep, 'property_priors_082016.mat');%'property_priors_082016.txt';
Pu_model_files = strcat(CURRDIR, filesep, 'all_models', filesep, 'model_#_061616.mat');
Pu_file = strcat(CURRDIR, filesep, 'data', filesep, 'pu_features_null_distributions.mat');
Graph_file = strcat(CURRDIR, filesep, 'data', filesep, 'updated_asymm_mech_graph.mat');
%Model_file = strcat(CURRDIR, filesep, 'model', filesep, 'alpha_nn_model_100914.mat');
Model_files = strcat(CURRDIR, filesep, 'all_models', filesep, ...
		     {'beta_nn_model_all_050415.mat', ...
		    'beta_nn_model_noparalog_050415.mat', ...
		    'beta_nn_model_predaln_paralog_050415.mat', ...
		    'beta_nn_model_predaln_noparalog_050415.mat'});
To_skip = {'VSL2B_disorder', 'B_factor', 'Helix', 'Strand', 'Loop', 'Coiled_coil', 'Calmodulin_binding', 'Motifs'};
PU_MODELS = [];


%%%%%% Set model indices %%%%%%
if DO_HOMOLOGY
    if PREDICT_CONS
        MINDS = [1, 3];
    else
        MINDS = [1];
    end
else
    if PREDICT_CONS
        MINDS = [2, 4];
    else
        MINDS = [2];
    end
end


%%%%%% Load files %%%%%%
% Grantham matrix
load(Grantham_file);
GRANTHAM = grantham_mat;

% Solvent accessibility
load(Rsa_file);
RSA_MODELS = rsa_model;
RSA_PROJS = projection;
RSA_MEANS = mu(:)';
RSA_STDS = sig;
RSA_F = features;

% Signal/TM
load(Sigtm_file);
SIGTM_MODELS = models;
SIGTM_PROJS = projections;
SIGTM_MEANS = mus;
SIGTM_STDS = sigmas;

% Coiled-coil
load(Cc_file);
EST_EMISSION = est_emission;
EST_TRANSITION = est_transition;
HIDDEN_STATES = hidden_states;
OBSERVED_STATES = observed_states;

% Stability
load(Stab_cfile);
G_RI = G;
ALPHA_RI = alpha_Y;
BIAS_RI = bias;
SV_RI = support_vectors;

load(Stab_rfile);
G_DDG = G;
ALPHA_DDG = alpha_Y;
BIAS_DDG = bias;
SV_DDG = support_vectors;

% Catalytuc residue
load(Cat_file);
CAT_MODELS = cat_models;
CAT_PROJS = cat_projs;
CAT_MUS = cat_mus;
CAT_SIGMAS = cat_sigmas;

% DNA
load(Dbind_file);
DBIND_MODEL = dbind_model;
DBIND_PROJ = proj;
DBIND_MEAN = mu(:)';
DBIND_STD = sig;
DBIND_F = features;

% RNA
load(Rbind_file);
RBIND_MODELS = rbind_models;
RBIND_PROJS = rbind_projs;
RBIND_MUS = rbind_mus;
RBIND_SIGMAS = rbind_sigmas;
RBIND_FEATS = rbind_feats;

% PPI
load(Ppi_file);
PPI_MODELS = ppi_models;
PPI_PROJS = projections;
PPI_MEANS = mus;
PPI_STDS = sigmas;
PPI_F = features;

% Hotspot
load(Hotspot_file);
HOTSPOT_MODELS = hotspot_models;
HOTSPOT_PROJS = hotspot_projs;
HOTSPOT_MUS = hotspot_mus;
HOTSPOT_SIGMAS = hotspot_sigmas;
HOTSPOT_FEATS = hotspot_feats;

% MoRF
load(Morf_file);
MORF_MODEL = morf_model;

% Allosteric
load(Allo_file);
ALLO_MODEL = allo_model;

% Motif
load(Motif_file);
MOTIFS = motifs;
MOTIF_DESCS = motif_descs;
MOTIF_NAMES = motif_names;
MOTIF_PATTS = motif_patts;

% Priors
%fid = fopen(Priors_file, 'r');
%if fid == -1
%    error('ERROR: Cannot open priors file!');
%end
%rows = textscan(fid, '%s %.4f %.6f %.6f %.6f %.6f %.2f', 'delimiter', '\t');
%PRIOR_PROPS = rows{1};
%ALPHAS_BETAS = [rows{3} rows{4} rows{6}];
%RATIOS = rows{7};
%fclose(fid);
load(Priors_file);
PRIOR_PROPS = prior_props;
ALPHAS_BETAS = ab;
RATIOS = c;

% Load in PU model files
for i = 1:length(PRIOR_PROPS)
    if ~ismember(PRIOR_PROPS{i}, To_skip)
        in_mat = strrep(Pu_model_files, '#', PRIOR_PROPS{i});
        if strcmp(PRIOR_PROPS{i}, 'Stability')
  	    PU_MODELS = [PU_MODELS; load(strrep(in_mat, '061616', '081416'))];
	else
  	    PU_MODELS = [PU_MODELS; load(in_mat)];
	end
    else
        PU_MODELS = [PU_MODELS; struct('this_prop', PRIOR_PROPS{i}, 'mus', NaN, 'sigmas', NaN, 'projections', NaN, 'features', NaN, 'nn_model', NaN, 'auc', NaN)];
    end
end

% PU null distributions
load(Pu_file);
NEUTRAL_PROPERTIES_PU = neut_ntrans_feats;
PROPERTY_NAMES_PU = neut_properties;

% Load graph file
load(Graph_file);
NODE_LABELS = node_labels;
G = onto_graph;

% MP2 Model/s
load(Model_files{MINDS(1)});
NEUTRAL_PROPERTIES = neutral_properties;
PROPERTY_NAMES = property_names;
PROP_INDS = other_inds;
NN_MODELS = nn_models;
NN_MUS = nn_mus;
NN_PROJS = nn_projs;
NN_SELECTED = nn_selected;
NN_STDS = nn_stds;
if length(MINDS) > 1
    load(Model_files{MINDS(2)});
    PNN_MODELS = nn_models;
    PNN_MUS = nn_mus;
    PNN_PROJS = nn_projs;
    PNN_SELECTED = nn_selected;
    PNN_STDS = nn_stds;
end

return
