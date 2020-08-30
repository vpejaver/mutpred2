function output = sim2(X, w1, w2, xmin, xmax)
% Implement 'sim' function for NN model
% 
% Args:
%   X      data set for test; each row represents one data point
%   w1, w2        weight matrices for the first and second layers
%   xmin, xmax    minimum and maximum values for each column
% 
% Dimension:
%   m     batch size
%   n     number of input variables, including constant 1
%   hn    number of hidden nodes for network 1
%   X     [m   x        n]
%   w1    [hn  x        n]
%   w2    [4   x (hn + 1)]
%   xmin  [1   x  (n - 1)]
%   xmax  [1   x  (n - 1)]
% 
% Returns:
%   prediction with each row representing one data point.
% 
% Last modified: Tue Feb 22 17:16:56 2011

% preprocess the data
n = size(X, 1);
keep = ones(n, 1);
rangex = xmax - xmin;
rangex(rangex == 0) = 1;                % avoid divisions by zero
% normalize test to [-1, 1] interval (minmax normalization)
X = (X - xmin(keep, :)) ./ rangex(keep, :) * 2 - 1;

% produce outputs after hidden neurons
h = tanh([X keep] * (w1'));

% produce outputs
o = tanh([h keep] * (w2'));

% normalize output to [0, 1] interval (minmax normalization)
output = (o + 1) ./ 2;

return