function [sequence, header] = readFASTA (filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [sequence, header] = readFASTA (filename)
% 
% Reads file containing multiple proteins in FastA format. 
%
% header = cell array of sequence headers
% sequence = cell array of sequences
%          = first sequence corresponds to the first header sequence etc.
% filename = text file in FastA format
%
% Note: if header extends to more than the first line this function will
%       not work properly. Make sure the input file is in appropriate format.
% 
% Predrag Radivojac, 2001-2002
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(filename, 'rt');
if (fid == -1)
   error(['The file ' filename ' cannot be opened!']);
end

sequence = [];
header = [];

j = 0;
while feof(fid) == 0
    line = fgetl(fid);
    if ~isempty(line) & ischar(line)
        if line(1) == '>'
            j = j + 1;
            header{j} = line;
            sequence{j} = '';
        else
            sequence{j} = strcat(sequence{j}, line);
        end
    end
end

fclose(fid);

return