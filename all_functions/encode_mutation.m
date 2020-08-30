% Vikas Pejaver
% April 2013

% Function to encode mutations as 380-element binary vectors for (20
% x 19) possible mutations
% Input a vector containing substitutions as follows: {'AB', 'CD',
% ... etc.}

function [binmat, possibilities] = encode_mutation(mut_list)

% Constants and defaults
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY';

% Get all possible combinations
possibilities = combnk(AMINO_ACIDS, 2);
possibilities = cellstr([possibilities; fliplr(possibilities)]);

% Loop through each mutation in the list and encode
for i = 1:length(mut_list)
    binvec = zeros(1, length(possibilities));  
    binvec(strcmp(possibilities, mut_list{i})) = 1;
    if ~any(binvec)
        binvec(strcmp(possibilities, [mut_list{i}(2) mut_list{i}(1)])) = 1;
    end
    binmat(i, :) = binvec;
end

return
