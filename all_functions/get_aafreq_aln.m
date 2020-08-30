% Vikas Pejaver
% July 2013

% Function to get amino acid frequency in alignment column for a
% given cell array of mutations (both wild and mutant residues)
% Returns a 12-element vector (wild-mutant paired, with and without
% gaps over 3 thresholds)

function [freq_vec] = get_aafreq_aln(aln, head, mutations)

% Constants and defaults
thresholds = [0.28, 0.84, 2];
species = 9606;
freq_vec = zeros(length(mutations), 6);

% Loop through all mutations
for i = 1:length(mutations)
    wild = mutations{i}(1);
    mut = mutations{i}(end);
    pos = str2num(mutations{i}(2:end-1));

    % Loop through different threshold cutoffs
    k = 1;
    for j = 1:length(thresholds)
    
        % Get sub-alignment
	[this_aln, this_head] = alignment_subset(aln, head, thresholds(j));

	% Get the human sequence
	this_seq = this_aln(this_head == species, :);
    
	% Get positions without gaps and update
	gapless = find(this_seq ~= '-');
	new_pos = gapless(pos);
	
	% Get alignment column and calculate frequencies
	aln_col = this_aln(:, new_pos);
	freq_vec(i, k) = length(find(aln_col == wild)) / length(aln_col);
	freq_vec(i, k+1) = length(find(aln_col == mut)) / length(aln_col);
	freq_vec(i, k+2) = length(find(aln_col == wild)) / length(find(aln_col ~= '-'));
	freq_vec(i, k+3) = length(find(aln_col == mut)) / length(find(aln_col ~= '-'));

	k = k + 4;
    end
end

return
