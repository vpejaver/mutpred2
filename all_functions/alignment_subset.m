% Vikas Pejaver
% October 2012

% Function to extract a subset of a multiple sequence alignment,
% given the alignment, an NCBI taxa list corresponding to the order
% of the sequences and distance cutoff between 0 and 1
% Basically, it removes any sequences OLDER than the given taxon cutoff
% and returns a cleaned multiple sequence alignment with all-gap
% columns removed
% NOTE: The distances are measured with respect to the human sequence

function [new_aln, new_head] = alignment_subset(aln, headers, cutoff)

global CURRDIR

% Constants and defaults
dist_file = strcat(CURRDIR, filesep, 'data', filesep, 'ucsc_alignments', filesep, 'ucsc_taxa_dist_real.mat');
new_aln = '';
new_head = '';

% Load distance file
load(dist_file);

% Identify retained taxa
to_keep = taxa_list(distances(1, :) < cutoff);
to_keep = cellfun(@str2num, to_keep);

% Find these remaining taxa in the given taxa list
inds = find(ismember(headers, to_keep));
new_head = headers(inds);
new_aln = aln(inds, :);

% Remove all-gap columns
new_aln = clean_allgap_cols(new_aln);

return
