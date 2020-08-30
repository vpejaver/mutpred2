% Vikas Pejaver
% May 2013

% FDR procedure with original BH procedure for the case when
% the p values are based on statistically independent tests.
% input: P is a vector of p-values, q is a scalar which controls 
% the false discovery rate. output: F the indices of the element 
% that pass the procedure
% NOTE: Code modified from Broad Institute's FDR package (cannot
% have NaNs)

function F = benjamini_hochberg(P, q)

% Constants and defaults
N = length(P);
F = [];

% Prepare probability vector
[m, n] = size(P);
if n > m
   P = P'; 
end
P_sor = sort(P);

% Calculate differences between p-values and corrected q-value
P_diff = P_sor - ([1:N]' * q / N);

% Find the largest 'k' such that P_k <= corrected q-value
k = find(P_diff <= 0, 1, 'last');
if k > 0
    p_d = P_sor(k);
    F = find(P <= p_d);
end

return
