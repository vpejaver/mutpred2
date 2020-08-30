function [flag] = writeFASTA (sequence, header, filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [flag] = writeFASTA (sequences, headers, filename)
% 
% Writes preteins into a file in FastA format. 
%
% headers = cell array of sequence headers (strings)
%         = symbol '>' will be inserted at the beginning of each header
%           if it is not already there
% sequences = cell array of sequences (strings)
%           = first sequence corresponds to the first header sequence etc.
% filename = string, name of the output file
% flag = integer, 0 if writing was successfully completed
%
% The function will NOT check if the sequences are valid protein sequences.
%
% Example:
%         headers{1} = '>1BDO:_ ACETYL-COA CARBOXYLASE'
%         headers{2} = '>1BOV:A VEROTOXIN-1'
%         sequences{1} = 'EISGHIVRSPMVGTFYRTPSP'
%         sequences{2} = 'PDCVTGKVEYTKYND'
%
%         writeFASTA(sequences, headers, 'outfile.txt');
%
% Predrag Radivojac, 2001-2002
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(header) ~= length(sequence)
    error('Header and sequence cell arrays differ in length!');
end

fid = fopen(filename, 'wt');
if (fid == -1)
    error(['File ' filename ' cannot be opened!']);
end

len = length(header);

for i = 1 : len
    % add '>' if it is not at position 1 in header
    if header{i}(1) ~= '>'
        header{i} = ['>' header{i}];
    end

    fprintf(fid, '%s', header{i});
    fprintf(fid, '\n');

    if length(sequence{i}) <= 80
        fprintf(fid, '%s', sequence{i});
        fprintf(fid, '\n\n');
    else
        j = 1;
        while j <= length(sequence{i})
            fprintf(fid, '%s', sequence{i}(j : min(j + 79, length(sequence{i}))));
            fprintf(fid, '\n');
            j = j + 80;
        end
        fprintf(fid, '\n');
    end
end

flag = fclose(fid);

return