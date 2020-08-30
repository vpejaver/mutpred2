function [run] = find_runs(pos, min_dist)

  % pos should be an array of numbers (row vector)
  % min_dist should be a positive constant (scalar)

%pos = [10 40 40 100 150 500]
  %min_dist = 10;

%% Check inputs

if ~isnumeric(pos) %|| ~isrow(pos)
error('Illegal format for the array of positions')
end

if ~isnumeric(min_dist) || ~isscalar(min_dist) || min_dist <= 0
  error('Illegal format for the minimum distance');
end

%% Here's where it starts

% elements taken in previous runs (as we go)
taken = zeros(1, length(pos));

% define unique positions in the array (sort them at the same time)
uqpos = unique(pos);

% get counts for each position
for i = 1 : length(uqpos)
    cnt(i) = length(find(pos == uqpos(i)));
end

run = {};
i = 1;

while ~isempty(uqpos)
    % grab first available element
    
    % x contains the location of the first element (we actually don't need
% it, but that's okay)
    x = uqpos(1);
    
    % y contains the index of the first free element in the original array
    qy = find(pos == uqpos(1) & taken == 0);
    y = qy(1);
    taken(qy(1)) = 1;
    
    % if count for this location is 1, remove location
    % otherwise, reduce count by 1 and leave the element in
    if cnt(1) == 1
        uqpos = uqpos(2 : length(uqpos));
        cnt = cnt(2 : length(cnt));
    else
        cnt(1) = cnt(1) - 1;
    end

    % start adding other elements based on the minimum distance constraint
    while true
        q = find(uqpos >= x(length(x)) + min_dist);
        
        if isempty(q)
            break
        else
            x = [x uqpos(q(1))]; % not needed but so what

            qy = find(pos == uqpos(q(1)) & taken == 0);
            y = [y qy(1)];
            taken(qy(1)) = 1;
            
            if cnt(q(1)) == 1
                uqpos = [uqpos(1 : q(1) - 1) uqpos(q(1) + 1 : length(uqpos))];
                cnt = [cnt(1 : q(1) - 1) cnt(q(1) + 1 : length(cnt))];
            else
                cnt(q(1)) = cnt(q(1)) - 1;
            end
        end
    end
    run{i} = y;
    i = i + 1;
end    
       
return
