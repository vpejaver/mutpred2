% Vikas Pejaver
% October 2012

% Function to remove columns from a multiple sequence alignment
% that contain only gaps

function aln = clean_allgap_cols(msa)

% Constants and defaults
aln = '';

% Loop through all columns (No easier way to do this)
inds = [];
for j = 1:size(msa, 2)
    if length(find(msa(:, j) == '-')) ~= size(msa, 1)
        inds = [inds j];
    end
end
aln = msa(:, inds);

return