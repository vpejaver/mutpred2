% Vikas Pejaver
% December 2012

% Function to mutate sequences in silico, given the wild type
% sequence and an array of mutation positions in the following
% format (K43C, T54R - wild residue, position and new residue) 
% Returns a list of sequences that are mutated at the given position

function new_sequence = mutate_sequence2(sequence, mutations)

% Constants and defaults
new_sequence = '';

% Loop through sequences and return mutated ones
for i = 1:length(mutations)
    res = mutations{i}(end);
    pos = str2num(mutations{i}(2:end-1));
    sequence(pos) = res;
    new_sequence = sequence;
end

return
