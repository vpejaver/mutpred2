% Vikas Rao Pejaver
% February 2012

% Function that predicts PTM sites, given a sequence
% Optional arguments - binary vector of which PTMs to predict (default: all PTMs)
%                    - PSSM matrix from PSI-BLAST (default: no PSSMs used)
% Order of PTM list: Acetylation, ADP-ribosylation, Amidation, 
% C-linked glycosylation, Carboxylation,
% Disulfide linkage, 
% Farnesylation,
% Geranylgeranylation, GPI anchor amidation,
% Hydroxylation,
% Methylation, Myristoylation, 
% N-linked glycosylation, 
% O-linked glycosylation, 
% Palmitoylation, Phosphorylation, Proteolytic cleavage, PUPylation, Pyrrolidone Carboxylic Acid, 
% Sulfation, SUMOylation, 
% Ubiquitination

function [post_probs, mod_states, ptms, positions, fragments, motif_states, ptm_feats] = predict_ptms(sequence, varargin)

% Constants and defaults
ptm_list = {'Acetylation', 'ADP-ribosylation', 'Amidation', 'C-linked_glycosylation', 'Carboxylation', 'Disulfide_linkage', 'Farnesylation', 'Geranylgeranylation', 'GPI_anchor_amidation', 'Hydroxylation', 'Methylation', 'Myristoylation', 'N-linked_glycosylation', 'N-terminal_acetylation', 'O-linked_glycosylation', 'Palmitoylation', 'Phosphorylation', 'Proteolytic_cleavage', 'PUPylation', 'Pyrrolidone_carboxylic_acid', 'Sulfation', 'SUMOylation', 'Ubiquitination'};
res_list = {'K', 'ER', 'ACDEFGHIKLMNPQRSTVWY', 'W', 'E', 'C', 'C', 'C', 'N', 'KPY', 'KR', 'G', 'N', 'AGMST', 'ST', 'C', 'STY', 'ACDEFGHIKLMNPQRSTVWY', 'K', 'Q', 'Y', 'K', 'K'};
motifs = {'.', '.', '(.G[RK][RK])', '.', '.', '.', '.', '.', '.', '.', '.', '.', '(N[^P][ST][^P])', '.', '.', '.', '.', '.', '.', '.', '.', '([AVFPILM]K.(E|D))', '.'};
offsets = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0];
to_predict = ones(1, 23);
pssm = [];
pssm_flag = 0;
min_len = 80;

% Check argument list and change defaults
if length(varargin) == 1
    if size(varargin{1}, 1) == 1
        to_predict = varargin{1};	
    else
        pssm = varargin{1};
        pssm_flag = 1;
    end
elseif length(varargin) == 2    
    if size(varargin{1}, 1) == 1
        to_predict = varargin{1};
        pssm = varargin{2};
    else
        to_predict = varargin{2};
        pssm = varargin{1};
    end
    pssm_flag = 1;
end     

% Make features (Only one pass through all relevant residues)
all_res = sprintf('%s', res_list{to_predict == 1});
uniq_res = unique(all_res);
[dataset, peptides, residues] = make_features_ptm(sequence, uniq_res, pssm);
ptm_feats = dataset;

% Run through random forest models
post_probs = [];
mod_states = [];
positions = {};
fragments = {};
motif_states = {};

ptms = ptm_list(to_predict == 1);
res_set = res_list(to_predict == 1);
motif_set = motifs(to_predict == 1);
offs = offsets(to_predict == 1);
for i = 1:length(ptms) % For each PTM
    probs = zeros(1, size(sequence, 2));
    mstate = zeros(1, size(sequence, 2));
    conf = zeros(1, size(sequence, 2));
    
    if ~strcmp(ptms{i}, 'ADP-ribosylation') && ~strcmp(ptms{i}, 'Amidation') && ~strcmp(ptms{i}, 'Hydroxylation') && ~strcmp(ptms{i}, 'Proteolytic_cleavage')
        % Above 'if' statement checks if the PTMs are adp-ribosylation, amidation, hydroxylation,
        % proteolytic cleavage - these PTMs have models trained with amino acids together
	
	peps = repmat({''}, 1, size(sequence, 2));
	if strcmp(ptms{i}, 'Myristoylation') || strcmp(ptms{i}, 'N-terminal_acetylation')
	    % This 'if' statement exclusively works for N-terminal
            % PTMs - although it uses make features again, it is
            % not expensive as only one site is considered for the
            % given sequence (since these PTMs have separately
            % trained models for each amino acid - only this part
            % of the preceding 'if' block has this part of the code)
	    X = [];
	    frags = {};
	    aas = [];
	    pos = [];
	    if ismember(sequence(1), res_set{i})
                [X, frags, aas] = make_terminal_features_ptm(sequence, res_set{i}, pssm, 'N');
		pos = [pos 1];
	    end
	    if ismember(sequence(2), res_set{i})
	        tmp_seq = sequence(2 : end);
		if ~isempty(pssm)
		    tmp_pssm = pssm(:, 2:end);
		else
		    tmp_pssm = [];
		end
		[tmpX, tmp_frags, tmp_aas] = make_terminal_features_ptm(tmp_seq, res_set{i}, tmp_pssm, 'N');
		tmp_frags{1}(12) = sequence(1); % Correct the fragment for output purposes
		X = [X; tmpX];
		frags = [frags tmp_frags];
		aas = [aas; tmp_aas];
		pos = [pos 2];
	    end

	    % Predict (no motif mode - not applicable to these PTMs)
	    for j = 1:size(X, 1)
	        %[probs(pos(j)), conf(pos(j))] = rf_modpred_nomotif(X(j, :), strcat(ptms{i}, '_', aas(j)), motif_set{i}, pssm_flag);
		[probs(pos(j)), conf(pos(j))] = lr_modpred_nomotif(X(j, :), strcat(ptms{i}, '_', aas(j)), motif_set{i}, pssm_flag);
		peps(pos(j)) = frags(j);
	    end
	elseif strcmp(ptms{i}, 'Farnesylation') || strcmp(ptms{i}, 'Geranylgeranylation')
	    if ~isempty(regexp(sequence, 'C...$', 'match'))
	        pos = length(sequence) - 3;
	        [X, frags, aas] = make_terminal_features_ptm(sequence, res_set{i}, pssm, 'C');
		% Predict (no motif mode - not applicable to these PTMs)
		%[probs(pos), conf(pos)] = rf_modpred_nomotif(X(1, :), strcat(ptms{i}, '_', aas(1)), motif_set{i}, pssm_flag);
		[probs(pos), conf(pos)] = lr_modpred_nomotif(X(1, :), strcat(ptms{i}, '_', aas(1)), motif_set{i}, pssm_flag);
		peps(pos) = frags(1);
	    end
	else % Standard prediction block
	    for j = 1:length(res_set{i}) % For each residue
		peps(sequence == res_set{i}(j)) = peptides(residues == res_set{i}(j));
		X = dataset(find(residues == res_set{i}(j)), :);
		if strcmp(motif_set{i}, '.')
		    if ~isempty(X)
		        %[probs(sequence == res_set{i}(j)), conf(sequence == res_set{i}(j))] = rf_modpred_nomotif(X, strcat(ptms{i}, '_', res_set{i}(j)), motif_set{i}, pssm_flag);
			[probs(sequence == res_set{i}(j)), conf(sequence == res_set{i}(j))] = lr_modpred_nomotif(X, strcat(ptms{i}, '_', res_set{i}(j)), motif_set{i}, pssm_flag);
		    end
		else % For motif PTMs further split datasets
		    [ind_m, dm] = intersect(find(sequence == res_set{i}(j)), regexp(sequence, motif_set{i}, 'start') + offs(i));
		    [ind_n, dn] = setdiff(find(sequence == res_set{i}(j)), regexp(sequence, motif_set{i}, 'start') + offs(i));
		    Xm = X(dm, :);
		    Xn = X(dn, :);
                              
		    if ~isempty(Xm)
		        %[probs(ind_m), conf(ind_m)] = rf_modpred_motif(Xm, strcat(ptms{i}, '_', res_set{i}(j)), pssm_flag);
			[probs(ind_m), conf(ind_m)] = lr_modpred_motif(Xm, strcat(ptms{i}, '_', res_set{i}(j)), pssm_flag);
			mstate(ind_m) = 1;
		    end
                
		    if ~isempty(Xn)
		        %[probs(ind_n), conf(ind_n)] = rf_modpred_nomotif(Xn, strcat(ptms{i}, '_', res_set{i}(j)), motif_set{i}, pssm_flag);
			[probs(ind_n), conf(ind_n)] = lr_modpred_nomotif(Xn, strcat(ptms{i}, '_', res_set{i}(j)), motif_set{i}, pssm_flag);
		    end
		end
	    end
	end	
	peps = peps(~strcmp(peps, ''));
    else
        % Joint amino acid models
        X = [];
        Xm = [];
        inds = [];
        indr = [];
        for j = 1:length(res_set{i}) % For each residue
            inds = [inds find(sequence == res_set{i}(j))];
            indr = [indr find(residues == res_set{i}(j))'];
        end
        
        peps = repmat({''}, 1, size(sequence, 2));
        if strcmp(motif_set{i}, '.')
            X = dataset(indr, :);
	    if ~isempty(X)
	        %[probs(inds), conf(inds)] = rf_modpred_nomotif(X, strcat(ptms{i}, '_', res_set{i}), motif_set{i}, pssm_flag);
		[probs(inds), conf(inds)] = lr_modpred_nomotif(X, strcat(ptms{i}, '_', res_set{i}), motif_set{i}, pssm_flag);
            end
        else % For motif PTMs further split datasets
            [ind_m, dm] = intersect(inds, regexp(sequence, motif_set{i}, 'start') + offs(i));
            [ind_n, dn] = setdiff(inds, regexp(sequence, motif_set{i}, 'start') + offs(i));
            Xm = dataset(dm, :);
            X = dataset(dn, :);
            if ~isempty(Xm)
	        %[probs(ind_m), conf(ind_m)] = rf_modpred_motif(Xm, strcat(ptms{i}, '_', res_set{i}), pssm_flag);
                [probs(ind_m), conf(ind_m)] = lr_modpred_motif(Xm, strcat(ptms{i}, '_', res_set{i}), pssm_flag);
		mstate(ind_m) = 1;
            end
            
            if ~isempty(X)
	        %[probs(ind_n), conf(ind_n)] = rf_modpred_nomotif(X, strcat(ptms{i}, '_', res_set{i}), motif_set{i}, pssm_flag);
		[probs(ind_n), conf(ind_n)] = lr_modpred_nomotif(X, strcat(ptms{i}, '_', res_set{i}), motif_set{i}, pssm_flag);
	    end
        end
        peps(inds) = peptides(indr);
        peps = peps(~strcmp(peps, ''));
    end
 
    post_probs = [post_probs; probs];
    mod_states = [mod_states; conf];
    positions = [positions find(probs)];
    fragments = [fragments {peps}];
    motif_states = [motif_states mstate(find(probs))];
end

return
