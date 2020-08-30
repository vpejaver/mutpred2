% Vikas Pejaver
% October 2015

% Function to reformat P-values based on graph structure

function [path_probs, path_pvals, property_paths] = propagate_posteriors(probs, pvals)

global NAMES_REORDERED NODE_LABELS G;

%probs = [0.4397, 0.9508, 0.3091, 0.0117, 0.4170, 0.4775, 0.1584, 0.4725, 0.8385, 0.5314, 0.6683, 0.6047, 0.8697, 0.8778, 0.0854, 0.0592, 0.0318, 0.1193, 0.3151, 0.7620, 0.2572, 0.1681, 0.1596, 0.7414, 0.1515, 0.3070, 0.2303, 0.6302, 0.8566, 0.7897, 0.4482, 0.6252, 0.8912, 0.5191, 0.8840, 0.5012, 0.5505, 1.0000, 0.0838, 0.4706, 1.0000, 1.0000, 1.0000, 0.2949, 0.1165, 1.0000, 1.0000, 1.0000, 1.0000, 0.6407, 0.1873, 0.1213, 1.0000, 0.7762, 0.2802, 0.7526, 0.5795, 0.0558];
%pvals = [0.4397, 0.9508, 0.3091, 0.0117, 0.4170, 0.4775, 0.1584, 0.4725, 0.8385, 0.5314, 0.6683, 0.6047, 0.8697, 0.8778, 0.0854, 0.0592, 0.0318, 0.1193, 0.3151, 0.7620, 0.2572, 0.1681, 0.1596, 0.7414, 0.1515, 0.3070, 0.2303, 0.6302, 0.8566, 0.7897, 0.4482, 0.6252, 0.8912, 0.5191, 0.8840, 0.5012, 0.5505, 1.0000, 0.0838, 0.4706, 1.0000, 1.0000, 1.0000, 0.2949, 0.1165, 1.0000, 1.0000, 1.0000, 1.0000, 0.6407, 0.1873, 0.1213, 1.0000, 0.7762, 0.2802, 0.7526, 0.5795, 0.0558];
%NAMES_REORDERED = {'Intrinsic_disorder', 'B-factor', 'Relative_solvent_accessibility', 'Helix', 'Strand', 'Loop', 'N-terminal_signal', 'Signal_helix', 'C-terminal_signal', 'Signal_cleavage', 'Cytoplasmic_loop', 'Transmembrane_region', 'Non-cytoplasmic_loop', 'Non_transmembrane', 'Coiled_coil', 'Catalytic_site', 'Calmodulin_binding', 'DNA_binding', 'RNA_binding', 'PPI_residue', 'PPI_hotspot', 'MoRF', 'Allosteric_site', 'Cadmium_binding', 'Calcium_binding', 'Cobalt_binding', 'Copper_binding', 'Iron_binding', 'Magnesium_binding', 'Manganese_binding', 'Nickel_binding', 'Potassium_binding', 'Sodium_binding', 'Zinc_binding', 'Acetylation', 'ADP-ribosylation', 'Amidation', 'C-linked_glycosylation', 'Carboxylation', 'Disulfide_linkage', 'Farnesylation', 'Geranylgeranylation', 'GPI_anchor_amidation', 'Hydroxylation', 'Methylation', 'Myristoylation', 'N-linked_glycosylation', 'N-terminal_acetylation', 'O-linked_glycosylation', 'Palmitoylation', 'Phosphorylation', 'Proteolytic_cleavage', 'Pyrrolidone_carboxylic_acid', 'Sulfation', 'SUMOylation', 'Ubiquitylation', 'Motifs', 'Stability'};

% Constants and defaults
property_paths = {};
path_probs = {};
path_pvals = {};
%graph_file = 'updated_asymm_mech_graph.txt';

% Read in graph file
%[g, N] = tblread(graph_file, '\t');
%G = sparse(g);
%NODE_LABELS = cellstr(N);
[~, parents] = graphtraverse(G, 1, 'method', 'bfs');
to_assign = setdiff(NODE_LABELS, NAMES_REORDERED, 'stable');

% Loop through each property, set P-values for shallow nodes and
% generate tree trace
j = 1;
tmp_probs = [zeros(1, length(to_assign)) probs];
tmp_pvals = [ones(1, length(to_assign)) pvals];
for i = 1:length(NODE_LABELS)
    if ismember(NODE_LABELS{i}, to_assign) %isempty(idx)
        traversed = graphtraverse(G, i, 'method', 'bfs');
	[tp, m] = max(tmp_probs(traversed));
	[tmp_probs(i), m] = max(tmp_probs(traversed));
	tmp_pvals(i) = tmp_pvals(traversed(m));
    else
        k = i;
	path_trace = [];
	while k ~= 0
	    path_trace = [k path_trace];
	    k = parents(k);
	end

	property_paths{j} = NODE_LABELS(path_trace);
	path_probs{j} = tmp_probs(path_trace);
	path_pvals{j} = tmp_pvals(path_trace);
	j = j + 1;
    end
end

return
