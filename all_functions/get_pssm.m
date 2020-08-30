% Vikas Pejaver
% November 2012

% Function to retrieve PSSMs from the local NCBI Refseq PSSM database
% database, given a sequence
% Returns a N x 42 matrix where N is length of the protein

function [pssm] = get_pssm(sequence)

% Global variable
global CURRDIR;

% Constants and defaults
pssm = [];
file_tag = strcat(CURRDIR, filesep, 'data', filesep, 'pssms', filesep, 'refseq_sp_082713_pssms_');
thresholds = [30 40 47 54 60 66 72 78 84 90 95 101 106 111 116 121 126 131 136 140 146 151 156 161 166 172 178 184 190 196 204 211 219 227 236 245 254 263 272 283 293 304 314 325 336 348 360 372 385 399 413 427 442 457 471 487 503 520 537 556 576 598 620 645 671 699 729 762 801 842 891 942 1002 1072 1154 1248 1373 1581 1880 2437 35991];

% Get length of the sequence
len = length(sequence);

% Get correct MAT file and load it
I = find(thresholds > len, 1, 'first');
if len == thresholds(end) % Last sequence is a special case
    mat_file = strcat(file_tag, num2str(thresholds(end-1)), '_', num2str(thresholds(end)), '.mat');
elseif I ~= 1 & ~isempty(I)
    mat_file = strcat(file_tag, num2str(thresholds(I-1)), '_', num2str(thresholds(I)), '.mat');
else
    warning('No match found in database. Returning empty matrix.');
    return
end
load(mat_file);

% Get indices in terms of lengths and the actual matches
l_inds = find(lengths == len);
inds = find(strcmp(sequences(l_inds), sequence));

% Get the PSSM
if length(inds) >= 1
    index = inds(1); % If more than one match, return the first one
                     % as the pssms would be the same anyways
    pssm = pssms{l_inds(index)};
else
    warning('No match found in database. Returning empty matrix.');
end

return
