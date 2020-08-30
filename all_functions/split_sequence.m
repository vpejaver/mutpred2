% Vikas Pejaver
% February 2014

% Function to split up a large sequence into smaller sequences so
% that processing can be sped up

function [subseqs, newpos, idx, clusters] = split_sequence(sequence, positions, padding)

% Constants and defaults
subseqs = {};
newpos = {};
idx = {};
clusters = {};

% Sort
[sorted, I] = sort(positions);

% Loop through
j = 1;
clusters{1} = [sorted(1)];
idx{1} = [I(1)];
for i = 2:length(sorted)
    if sorted(i) - sorted(i-1) < padding
        clusters{j} = [clusters{j} sorted(i)];
	idx{j} = [idx{j} I(i)];
    else
        this_start = max([1, min(clusters{j}) - padding]);
        this_end = min([length(sequence), max(clusters{j}) + padding]);
	subseqs{j} = sequence(this_start : this_end);
	% convert to new positions
	orig_vec = this_start : this_end;
	new_vec = 1 : length(subseqs{j});
	newpos{j} = new_vec(ismember(orig_vec, clusters{j}));

        j = j + 1;

	clusters{j} = [sorted(i)];
	idx{j} = [I(i)];
    end
end

% Add last subsequence and updated positions
this_start = max([1, min(clusters{j}) - padding]);
this_end = min([length(sequence), max(clusters{j}) + padding]);
subseqs{j} = sequence(this_start : this_end);
orig_vec = this_start : this_end;
new_vec = 1 : length(subseqs{j});
newpos{j} = new_vec(ismember(orig_vec, clusters{j}));

return
