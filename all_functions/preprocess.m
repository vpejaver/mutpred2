% Vikas Pejaver
% April 2017

% Function to read in the FASTA file and pre-process (check for
% ANNOVAR format and convert if required

function [prots, subs, seqs] = preprocess()

% Global variables
global IN_FILE OID;

% Constants and defaults
prots = {};
subs = {};
seqs = {};

% Read in file
[sequences, headers] = readFASTA(IN_FILE);

% If ANNOVAR, then convert
if ~isempty(regexp(headers{1}, 'WILDTYPE', 'match'))
    [sequences, headers] = av2mp2(sequences, headers);
    if OID ~= 1
        fprintf(1, 'ANNOVAR input format detected. Converting ...\n');
    end
end

% Merge duplicate sequences (to be more efficient)
[U, ~, iseq] = unique(sequences, 'stable');
for i = 1:length(U)
    idx1 = find(iseq == i);

    mids = {};
    mmuts = {};
    for j = 1:length(idx1)
        tokens = regexp(headers{idx1(j)}, ' ', 'split');
        mmuts = [mmuts; tokens(2:end)'];
        mids = [mids; repmat(strrep(tokens(1), '>', ''), length(tokens)-1, 1)];
    end

    [subs{i}, ~, im] = unique(mmuts, 'stable');
    prots{i} = {};
    for j = 1:length(subs{i})
        idx2 = find(im == j);
        prots{i} = [prots{i}; regexprep(sprintf('%s ', mids{idx2}), ' $', '')];
    end
    seqs{i} = U{i};
end
if OID ~= 1
    fprintf(1, 'After merging identical sequences, %d unique substitutions in %d sequences were found.\n', sum(cellfun(@length, subs)), length(seqs));
end

return
