% Vikas Pejaver
% October 2012

% Function to mutate sequences in silico, given the wild type
% sequence and an array of mutation positions in the following
% format (K43C, T54R - wild residue, position and new residue) 
% Returns a list of sequences that are mutated at the given position

function sequences = mutate_sequence(sequence, mutations)

% Constants and defaults
sequences = {};

% Loop through sequences and return mutated ones
for i = 1:length(mutations)
    wild = sequence;
    res = mutations{i}(end);
    pos = str2num(mutations{i}(2:end-1));
    wild(pos) = res;
    sequences = [sequences wild];
end

return
