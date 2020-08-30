% Vikas Pejaver
% March 2016

% Function to convert ANNOVAR coding change file to MutPred2 input format

function [mp2_seqs, mp2_headers] = av2mp2(seqs, heads)

% Constants and defaults
mp2_headers = {};
mp2_seqs = {};

% Loop through sequences and extract only substitutions
for i = 2 : 2 : length(seqs)
    if isempty(regexp(heads{i}, 'protein-altering', 'match'))
        continue;
    end
    if length(seqs{i}) ~= length(seqs{i-1}) % make sure WT and MUT sequences are of the same lengths
        continue;
    end

    if length(find(seqs{i} ~= seqs{i-1})) > 1
        continue;
    end
    
    tokens = regexp(heads{i}, ' +', 'split');
    if tokens{end-2} == '*' | strrep(tokens{end}, ')', '') == '*'
        continue;
    end
    mp2_headers = [mp2_headers; sprintf('%s %s', regexprep(sprintf('%s|', tokens{1:3}), '\|$', '')), [tokens{end-2}, tokens{end-5}, strrep(tokens{end}, ')', '')]];
    mp2_seqs = [mp2_seqs; strrep(seqs{i-1}, '*', '')];
end

return
