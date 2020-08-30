function [p] = composition_profile (s)

AA = 'ACDEFGHIKLMNPQRSTVWY';

p = zeros(1, length(AA));

for i = 1 : length(AA)
    p(i) = length(find(s == AA(i)));
end

return
