% Vikas Pejaver
% October 2014

% Function to retrieve alignments from the local UCSC alignment
% database, given a sequence
% Returns an array of integers that represent the NCBI taxa numbers
% for the 46 species in these alignments and the alignment
% (character array)

function [headers, aln] = get_alignment(sequence)

% Global variable
global CURRDIR;

% Constants and defaults
headers = '';
aln = '';
file_tag = strcat(CURRDIR, filesep, 'data', filesep, 'ucsc_alignments', filesep, 'hg19_refGene_alignments_');
thresholds = [17 121 161 198 233 268 305 332 364 399 437 474 516 563 625 699 792 914 1109 1512 33423];

% Get length of the sequence
len = length(sequence);

% Get correct MAT file and load it
I = find(thresholds > len, 1, 'first');
if len == thresholds(end) % Last sequence is a special case
    mat_file = strcat(file_tag, num2str(thresholds(end-1)), '_', num2str(thresholds(end)), '.mat');
elseif I ~= 1 & ~isempty(I)
    mat_file = strcat(file_tag, num2str(thresholds(I-1)), '_', num2str(thresholds(I)), '.mat');
else
    warning('No match found in database. Returning empty string.');
    return
end
load(mat_file);

% Get indices in terms of lengths and the actual matches
l_inds = find(lengths == len);
inds = find(strcmp(sequences(l_inds), sequence));

% Get the header info and alignment
if length(inds) > 1
    index = inds(1); % If more than one match, return the first one
                     % as the alignments would be the same anyways
    headers = aln_headers{l_inds(index)};
    aln = alignments{l_inds(index)};
elseif isempty(inds)
    warning('No match found in database. Returning empty string.');
else
    headers = aln_headers{l_inds(inds)};
    aln = alignments{l_inds(inds)};
end

return
