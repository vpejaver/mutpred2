% Vikas Pejaver
% November 2012

% Function that represents a mutation as a 21-element vector where
% the wild type residue is represented by -1, the mutant residue by
% 1 and all other elements as zero. This can specified over a
% window to indicate a presence or absence

function bin_vec = make_mutation_binary(sequence, mutation, win)

% Constants and defaults
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY*';
mutated = zeros(1, length(AMINO_ACIDS));
bin_vec = [];

% Parse mutation
wild = mutation(1);
mut = mutation(end);
index = str2num(mutation(2:end-1));

% Get left and right fragments
bound = (win - 1)/2;
left_margin = [max(1, index-bound) : index-1];
right_margin = [index+1 : min(index+bound, length(sequence))];
left_flank = sequence(left_margin);
right_flank = sequence(right_margin);

% Add '*' for terminii
left_diff = abs(length(left_flank) - bound);
right_diff = abs(length(right_flank) - bound);
if left_diff ~= 0
    left_flank = strcat(repmat('*', 1, left_diff), left_flank);
end
if right_diff ~= 0
    right_flank = strcat(right_flank, repmat('*', 1, right_diff));
end

% Make vector for left and right fragments
left_bin = my_make_spacer(left_flank);
right_bin = my_make_spacer(right_flank);

% Make vector for mutated position
mutated(AMINO_ACIDS == wild) = -1;
mutated(AMINO_ACIDS == mut) = 1;

% Combine
bin_vec = [left_bin mutated right_bin];

return
