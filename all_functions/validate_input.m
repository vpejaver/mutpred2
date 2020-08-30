% Vikas Pejaver
% October 2014

% Function to check if the input sequences and mutations are in the
% correct formats and flag them

function [valid_mutations, indices, markers] = validate_input(sequence, mutations)

%%%%%% Constants and defaults %%%%%%
valid_mutations = {};
indices = [];
markers = zeros(length(mutations), 6);

Mut_format = '(^[ACDEFGHIKLMNPQRSTVWY][0-9]+[ACDEFGHIKLMNPQRSTVWY]$)';
Nonstd = '[BJOUXZ]';
Window = 21; % TO DO: Double-check to see if this is the largest window size


%%%%%% Check sequence length %%%%%%
if length(sequence) < 30
    markers(:, 1) = repmat(1, length(mutations), 1);
    return;
end


%%%%%% Loop through each mutation %%%%%%
j = 1;
L = length(sequence);
for i = 1:length(mutations)
    if isempty(regexp(mutations{i}, Mut_format, 'match'))
        markers(i, 2) = 1;
    else
        wild = mutations{i}(1);
	mut = mutations{i}(end);
	position = str2num(mutations{i}(2:end-1));
	start = max(1, position - floor(Window / 2));
	finish = min(L, position + floor(Window / 2));

	% Check is position is within the length, wild-type residue
        % and presence of non-std amino acids in window
	if position > L
	    markers(i, 3) = 1;
        elseif sequence(position) ~= wild
	    markers(i, 4) = 1;
	elseif wild == mut
	    markers(i, 5) = 1;
	elseif ~isempty(regexp(sequence(start:finish), Nonstd, 'match'))
	    markers(i, 6) = 1;
	else
	    valid_mutations{j} = mutations{i};
	    indices(j) = i;
	    j = j + 1;
	end
    end
end

return
