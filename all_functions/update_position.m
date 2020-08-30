% Vikas Pejaver
% November 2012

% Function to get conservation index matrices for different
% thresholds and to revise positions based on alignment

function new_pos = update_position(aln, head, pos)

% Constants and defaults
thresholds = [0.28, 0.84, 2];
species = 9606;

% Loop through different threshold cutoffs
for i = 1:length(thresholds)
    
    % Get sub-alignment
    [this_aln, this_head] = alignment_subset(aln, head, thresholds(i));

    % Get the human sequence
    this_seq = this_aln(this_head == species, :);
    
    % Get positions without gaps
    gapless = find(this_seq ~= '-');

    new_pos(i) = gapless(pos);
end

return
