% Vikas Pejaver
% February 2014

% Function to split up a large sequence into smaller sequences so
% that processing can be sped up

function [intervals, newpos] = split_sequence2(positions, maxlen, padding)

% Constants and defaults
intervals = [];
newpos = [];

% Loop through each position and get intervals
for i = 1:length(positions)
    intervals(i, 1) = max([1, positions(i) - padding]);
    intervals(i, 2) = min([maxlen, positions(i) + padding]);
    tmp = [intervals(i, 1) : intervals(i, 2)];
    newpos(i) = find(tmp == positions(i));
end

return
