% Vikas Rao Pejaver
% February 2012

% Function that takes in a sequence, residues to consider and a
% PSSM (optional) and returns a feature matrix, peptides and amino acids
function [siteFeat, fragments, aminos] = make_features_ptm(protein, res_considered, PSSM)

fragment_half_size = 12;

siteFeat = []; % the big matrix
aminos = [];
fragments = {};
k = 1;
res = [];
    
% DEFINE RESIDUES %
q = [];
for i = 1:length(res_considered)
    q = [q find(protein == res_considered(i))];
end

if ~isempty(q)
    res = sort(q);

    % REMOVE AMBIGUOUS RESIDUES %
    newSeq = clean_protein(protein);

    % GET PREDICTION VALUES FROM OTHER PREDICTORS %
    predValues = use_predictors_ptm(newSeq, PSSM, 0);

    % MAKE FEATURES %
    vect = make_features_ptm_aa(protein, res);
    vect = [vect make_features_ptm_pr(protein, res, predValues)];

    % GET PSSM only if PSSM is not empty %
    if ~isempty(PSSM)
        vect = [vect make_features_ptm_ev(protein, res, PSSM)];
    end
    siteFeat = [siteFeat; vect];

    % BUILD FRAGMENTS %
    for j = 1 : length(res)
        %residues(k) = res(j);
	aminos = [aminos; protein(res(j))];
	left_margin = max(1, res(j) - fragment_half_size);
	right_margin = min(res(j) + fragment_half_size, length(protein));
	if res(j) - fragment_half_size < 1
	    a = blanks(fragment_half_size - res(j) + 1);
	    a = strrep(a,' ', '*');
	    fragments{k} = [a protein(1 : right_margin)];
	elseif res(j) + fragment_half_size > length(protein)
	    a = blanks(res(j) + fragment_half_size - length(protein));
	    a = strrep(a,' ', '*');
	    fragments{k} = [protein(left_margin : length(protein)) a];
	else
	    fragments{k} = protein(left_margin : right_margin);
	end
	k = k + 1;
    end

    % ADD BINARY FEATURES
    for i = 1 : length(fragments)
        feat_bin(i, :) = make_example_spacer(fragments{i}([10:12 14:16]));
    end
    siteFeat = [siteFeat feat_bin];
end
return
