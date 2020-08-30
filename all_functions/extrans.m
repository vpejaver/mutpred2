% Vikas Pejaver
% November 2012

% Function to calculate the change in transition frequencies
% because of mutations given a sequence and a list of mutations
% Returns an N x 6 matrix where N is the number of mutations. Every
% odd column contains a zero or a one (depending on whether the
% difference is positive or not) and every even column contains the
% absolute difference

function trans_diffs = extrans(sequence, mutations)

% Global variable
global CURRDIR;

% Constants and defaults
mat_file = strcat(CURRDIR, filesep, 'all_models', filesep, 'TransitionsBIG.mat');
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWYU';
window = 5;
trans_diffs = zeros(length(mutations), 6);

% Load MAT file
load(mat_file)

% Loop through mutations
for i = 1:length(mutations)
    wild = mutations{i}(1);
    mut = mutations{i}(end);
    index = str2num(mutations{i}(2:end-1));

    % Get left and right fragments
    bound = (window - 1)/2;
    left_margin = [max(1, index-bound) : index-1];
    right_margin = [index+1 : min(index+bound, length(sequence))];
    left_flank = sequence(left_margin);
    right_flank = sequence(right_margin);

    % Add 'U' for terminii
    left_diff = abs(length(left_flank) - bound);
    right_diff = abs(length(right_flank) - bound);
    if left_diff ~= 0
        left_flank = strcat(repmat('U', 1, left_diff), left_flank);
    end
    if right_diff ~= 0
        right_flank = strcat(right_flank, repmat('U', 1, right_diff));
    end
    
    % Make fragments
    wild_frag = [left_flank wild right_flank];
    mut_frag = [left_flank mut right_flank];

    % Loop through positions and get triplet frequencies
    for j = 1:3
        diff_scores(j) = frequencies(strcmp(triplets, wild_frag(j:j+2))) - frequencies(strcmp(triplets, mut_frag(j:j+2)));
    end
    binaries = zeros(1, 3);
    binaries(diff_scores >= 0) = 1;
    
    % Add to output matrix
    for j = 1:6
        if mod(j, 2) == 0
	    trans_diffs(i, j) = abs(diff_scores(round(j/2)));
	else
	    trans_diffs(i, j) = binaries(round(j/2));
	end
    end
end
return
