function seq = stdaa(seq)
% Replace non-standard amino acids with standard ones
% 
% Last modified: Wed Feb 23 17:04:14 2011

seq = upper(seq);

q = find(seq == 'X');
if ~isempty(q)
    seq(q) = 'A';
end
q = find(seq == 'B');
if ~isempty(q)
    seq(q) = 'D';
end
q = find(seq == 'Z');
if ~isempty(q)
    seq(q) = 'E';
end
q = find(seq == 'J');
if ~isempty(q)
    seq(q) = 'A';
end
q = find(seq == 'O');
if ~isempty(q)
    seq(q) = 'A';
end
q = find(seq == 'U');
if ~isempty(q)
    seq(q) = 'C';
end

return
