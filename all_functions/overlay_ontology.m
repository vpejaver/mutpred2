% Vikas Pejaver
% January 2017

% Function to overlay predictions on ontology and print

function onto_str = overlay_ontology(path_probs, path_pvals, property_paths, mech_types, posinfo, motinfo, prefix)

% Global variables
global PTHR;

% Constants and defaults
onto_str = '';
num_rows = 79;
num_cols = 7;
link_len = 6;
%str_len = 200;
padding = 10;

nodes = {'Effect', 'Structure_and_dynamics', 'Functional_residue', ...
	 'Special_structural_signatures', 'Secondary_structure', ...
	 'Stability_and_conformational_flexibility', 'Metal_binding', ...
	 'Macromolecular_binding', 'PTM_site', 'Transmembrane_protein', ...
	 'Signal_peptide', 'Conformational_flexibility', 'Protein_binding', ...
	 'Ordered_interface', 'Disordered_interface', 'Intrinsic_disorder', ...
	 'B-factor', 'Relative_solvent_accessibility', 'Helix', 'Strand', ...
	 'Loop', 'N-terminal_signal', 'Signal_helix', 'C-terminal_signal', ...
	 'Signal_cleavage', 'Cytoplasmic_loop', 'Transmembrane_region', ...
	 'Non-cytoplasmic_loop', 'Coiled_coil', ...
	 'Catalytic_site', 'Calmodulin_binding', 'DNA_binding', ...
	 'RNA_binding', 'PPI_residue', 'PPI_hotspot', 'MoRF', ...
	 'Allosteric_site', 'Cadmium_binding', 'Calcium_binding', ...
	 'Cobalt_binding', 'Copper_binding', 'Iron_binding', ...
	 'Magnesium_binding', 'Manganese_binding', 'Nickel_binding', ...
	 'Potassium_binding', 'Sodium_binding', 'Zinc_binding', ...
	 'Acetylation', 'ADP-ribosylation', 'Amidation', 'C-linked_glycosylation','Carboxylation', ...
	 'Disulfide_linkage', 'Farnesylation', ...
	 'Geranylgeranylation', 'GPI-anchor_amidation', ...
	 'Hydroxylation', 'Methylation', 'Myristoylation', ...
	 'N-linked_glycosylation', 'N-terminal_acetylation', ...
	 'O-linked_glycosylation', 'Palmitoylation', ...
	 'Phosphorylation', 'Proteolytic_cleavage', ...
	 'Pyrrolidone_carboxylic_acid', 'Sulfation', ...
	 'SUMOylation', 'Ubiquitylation', 'Motifs', 'Stability'};
%inds = [1, 3, 30, 17, 5, 10, 45, 33, 58, 24, 19, 12, 36, 41, 37, 13, ...
%      14, 4, 6, 7, 8, 20, 21, 22, 23, 26, 27, 28, 18, 31, 38, 34, ...
%	35, 42, 43, 39, 32, 46:56, 59:63, 15, 64:79, 2, 11];
inds = [1, 3, 30, 17, 5, 10, 45, 33, 58, 25, 19, 13, 36, 41, 37, 14, ...
      15, 4, 6, 7, 8, 20, 21, 22, 23, 26, 27, 28, 18, 31, 38, 34, ...
	35, 42, 43, 39, 32, 46:56, 59:63, 12, 64:79, 2, 11];

startlink = ['|' repmat('-', 1, link_len-1)];
midlink = repmat('-', 1, link_len+3);
endlink = [repmat('-', 1, link_len-1) '>'];
singlink = ['|' repmat('-', 1, link_len-2) '>'];
vertlink = '|';

Altered = {'Effect', 'Structure_and_dynamics', 'Functional_residue', 'Special_structural_signatures', 'Secondary_structure', 'Stability_and_conformational_flexibility', 'Metal_binding', 'Macromolecular_binding', 'PTM_site', 'Transmembrane_protein', 'Signal_peptide', 'Conformational_flexibility', 'Protein_binding', 'Ordered_interface', 'Disordered_interface', 'N-terminal_signal', 'Signal_helix', 'C-terminal_signal', 'Signal_cleavage', 'Cytoplasmic_loop', 'Transmembrane_region', 'Non-cytoplasmic_loop', 'Non_transmembrane', 'Coiled_coil', 'Calmodulin_binding', 'DNA_binding', 'RNA_binding', 'PPI_residue', 'PPI_hotspot', 'MoRF', 'Stability'};
Region = [Altered, 'Helix', 'Strand', 'Loop', 'Intrinsic_disorder', 'B-factor', 'Relative_solvent_accessibility'];
Fixed = {'Effect', 'Structure_and_dynamics', 'Functional_residue', 'Special_structural_signatures', 'Secondary_structure', 'Stability_and_conformational_flexibility', 'Metal_binding', 'Macromolecular_binding', 'PTM_site', 'Transmembrane_protein', 'Signal_peptide', 'Conformational_flexibility', 'Protein_binding', 'Ordered_interface', 'Disordered_interface'};
Types = {'Loss', 'Gain'};

% Loop through and print
tmp_vec = repmat({''}, num_rows, num_cols);
to_show = ones(num_rows, 1);
strfunc_pr = 0;
for i = 1:length(path_probs)
    for j = 1:length(path_probs{i})
        if strcmp(property_paths{i}{j}, 'Non_transmembrane')
	    continue;
	end
	if strcmp(property_paths{i}{j}, 'Structure_dynamics_and_functional_residue')
	    strfunc_pr = path_probs{i}(j);
	    continue;
	end

	% make adjustments for properties that are both structural
        % and functional in nature
	if strcmp(property_paths{i}{j}, 'Structure_and_dynamics') | strcmp(property_paths{i}{j}, 'Functional_residue')
	    path_probs{i}(j) = max(strfunc_pr, path_probs{i}(j));
	end

	c = j;
	idx = find(strcmp(nodes, property_paths{i}{j}));
	r = inds(idx);
	if path_pvals{i}(j) < 0.01
	    formatted_pval = sprintf('P = %.1e', path_pvals{i}(j));
	else
	    formatted_pval = sprintf('P = %.2f', path_pvals{i}(j));
	end
	if ~ismember(property_paths{i}{j}, Fixed) & ~strcmp(property_paths{i}{j}, 'Motifs') & (path_pvals{i}(j) >= PTHR | path_probs{i}(j) < 0.01)
	    to_show(r) = 0;
	end

	if ismember(property_paths{i}{j}, Fixed)
	    tmp_vec{r, c} = sprintf('Altered %s (Pr = %.2f)', property_paths{i}{j}, path_probs{i}(j));
	elseif strcmp(property_paths{i}{j}, 'Motifs')
	    this_mot = regexp(motinfo, ':|(; )', 'split');
	    mot_str = regexprep(sprintf('%s|', this_mot{1 : 2 : end}), '\|$', '');
	    tmp_vec{r, c} = sprintf('Altered %s (%s)', property_paths{i}{j}, mot_str);
	elseif ismember(property_paths{i}{j}, Altered) & ismember(property_paths{i}{j}, Region)
	    tmp_vec{r, c} = sprintf('Altered %s (Pr = %.2f | %s)', property_paths{i}{j}, path_probs{i}(j), formatted_pval);
	elseif ismember(property_paths{i}{j}, Altered)
	    tmp_vec{r, c} = sprintf('Altered %s at %s (Pr = %.2f | %s)', property_paths{i}{j}, strrep(sorted_pos{j}, ' ', ''), path_probs{i}(j), formatted_pval);
	elseif ismember(property_paths{i}{j}, Region)
	    tmp_vec{r, c} = sprintf('%s of %s (Pr = %.2f | %s)', Types{mech_types(i)}, property_paths{i}{j}, path_probs{i}(j), formatted_pval);
	else
	    tmp_vec{r, c} = sprintf('%s of %s at %s (Pr = %.2f | %s)', Types{mech_types(i)}, property_paths{i}{j}, strrep(posinfo{i}, ' ', ''), path_probs{i}(j), formatted_pval);
	end

	% add links (thresholds hardcoded)
	if c ~= 1
	    if (r > 3 & r < 30) | (r > 33 & r ~= 45 & r < 58)
	        tmp_vec{r, c-1} = regexprep(singlink, '^\|', '     \|');
	    else
	        tmp_vec{r, c-1} = singlink;
	    end
	end
    end
end

% Add other links in a hardcoded fashion
tmp_vec(4:29, 1) = repmat({vertlink}, length(4:29), 1); 
tmp_vec(31:end, 1) = repmat({'      '}, length(31:size(tmp_vec, 1)), 1); 
tmp_vec([34:44, 46:57], 2) = repmat({vertlink}, length([34:44, 46:57]), 1); 

vertlink = regexprep(vertlink, '\|', '     \|');
tmp_vec([6:9, 11:16], 2) = repmat({vertlink}, length([6:9, 11:16]), 1);
tmp_vec([18:29, 59:end], 2) = repmat({'      '}, length([18:29, 59:size(tmp_vec, 1)]), 1);
tmp_vec(20:24, 3) = repmat({vertlink}, length(20:24), 1);
tmp_vec([14:15, 26:28, 37:44], 3) = repmat({'      '}, length([14:15, 26:28, 37:44]), 1);
tmp_vec(38:40, 4) = repmat({vertlink}, length(38:40), 1);
tmp_vec(42:43, 4) = repmat({'      '}, length(42:43), 1); 

% Make string
for j = 1:size(tmp_vec, 2)
     q = find(strcmp(tmp_vec(:, j), ''));
     tmp_vec(q, j) = regexprep(repmat({'#'}, length(q), 1), '#', '     ');
end
for i = 1:size(tmp_vec, 1)
    tmp_vec2{i} = sprintf('%s,#  %s%s', prefix, sprintf('%s', tmp_vec{i, :}));
end

k = 1;
str_len = max(cellfun(@length, tmp_vec2(3:45))) + padding;
for i = 1:size(tmp_vec2, 2)
    if ~to_show(i)
        continue;
    end
    
    if i == 3 
        tmp_vec3{k} = regexprep(tmp_vec2{i}, ' +$', '-');
	tmp_vec3{k} = pad(tmp_vec3{k}, str_len+1, 'right', '-');
	tmp_vec3{k} = regexprep(tmp_vec3{k}, '-$', '\|');
    elseif i == 32 | i == 45
        tmp_vec3{k} = regexprep(tmp_vec2{i}, ' +$', '<');
	tmp_vec3{k} = pad(tmp_vec3{k}, str_len+1, 'right', '-');
	tmp_vec3{k} = regexprep(tmp_vec3{k}, '-$', '\|');
    elseif i > 2 & i < 46
        tmp_vec3{k} = regexprep(tmp_vec2{i}, ' +$', '');
	tmp_vec3{k} = pad(tmp_vec3{k}, str_len+1, 'right', ' ');
	tmp_vec3{k} = regexprep(tmp_vec3{k}, ' $', '\|');
    else
        tmp_vec3{k} = regexprep(tmp_vec2{i}, ' +$', '');
    end
    k = k + 1;
end

% Save
onto_str = sprintf('%s\n', tmp_vec3{:});

return
