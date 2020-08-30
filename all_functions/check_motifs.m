% Vikas Pejaver
% June 2013

% Function to assign scores to each position of a sequence
% depending on how many PROSITE motifs occur at that position
% Returns a vector of numbers normalized between 0 and 1 (where 0
% means no motifs are seen in a given position and 1 means all
% motifs are seen)

function [scores, ids, desc, patt, positions] = check_motifs(sequence)

% Global variables
global MOTIFS MOTIF_DESCS MOTIF_NAMES MOTIF_PATTS;

% Constants and defaults
scores = [];
ids = {};
desc = {};
patt = {};
positions = [];

% Load PROSITE
%load(motif_mat);

% Loop through all patterns 
bin_mat = zeros(length(MOTIFS), length(sequence));
for i = 1:length(MOTIFS)
    [m, n] = regexp(sequence, MOTIFS{i}, 'start', 'end');
    if ~isempty(m)
        for j = 1:length(m)
	    bin_mat(i, m(j) : n(j)) = 1;
	end
    end
end

% Sum up counts
total = sum(bin_mat);

% Calculate score
%scores = total / length(MOTIFS)
%scores = (2 * total) ./ (length(MOTIFS) + total);
scores = logsig(total);
%scores(scores == 0.5) = 0; % Set non-motif positions to zero

% Get motif information
idx = find(any(bin_mat, 2));
ids = MOTIF_NAMES(idx);
desc = MOTIF_DESCS(idx);
patt = MOTIFS(idx); %MOTIF_PATTS(idx);
positions = bin_mat(idx, :);

return
