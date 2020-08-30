function y = moving_average (x, window)

% function y = moving_average (x, window)
%
% x      = a vector of real numbers
% window = a two-sided window of an odd length
% y      = output, a filtered sequence
%
% Predrag Radivojac
% Temple University, 2002

if mod(window, 2) == 0 % window has to be an odd number
    y = [];
    return
% x has to be a vector, not a matrix
elseif size(x, 1) > 1 & size(x, 2) > 1
    y = [];
    return
% if a sequence length is comparable to the window size (or shorter) we do quadratic loop
% the assumption is that both window and the sequence x are short 
elseif window >= length(x)
    % set output to zero vector and adjust y to be of the same size as x
    y = zeros(size(x, 1), size(x, 2)); 
    len = length(x);
    for i = 1 : length(x)
        y(i) = mean(x(max(i - floor(window / 2), 1) : min(i + floor(window / 2), len)));
    end
    return
end

% otherwise, we do linear loop

y = zeros(size(x, 1), size(x, 2)); % set output to zero vector
    
one_sided_window = floor(window / 2) + 1;

t = sum(x(1 : one_sided_window)); % t is a running sum

% first, we expand the window at the left end
for i = 1 : one_sided_window - 1
    y(i) = t / (one_sided_window + i - 1);
    t = t + x(one_sided_window + i);
end

% second, we work with the full-sized window ('middle' of the sequence x)
len = length(x);
for i = one_sided_window : len - one_sided_window
    y(i) = t / window;
    t = t - x(i - one_sided_window + 1) + x(one_sided_window + i);
end

% third, we collapse the window at the right end
j = 0;
for i = length(x) - one_sided_window + 1 : len
    y(i) = t / (window - j);
    t = t - x(i - one_sided_window + 1);
    j = j + 1;
end
    
return