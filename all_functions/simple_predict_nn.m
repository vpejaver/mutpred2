% Vikas Pejaver
% March 2014

% Function to make predictions using a reduced neural network data
% structure as obtained from extract_NF_weights.m
% NOTE: Assumes single-layer network

function fx = simple_predict_nn(net, P)

% Constants and defaults
fx = [];
N = size(P, 1);

% Set min and maxes for normalization
%ymax = net.ymax;
%ymin = net.ymin;
%xmax = repmat(net.xmax, N, 1);
%xmin = repmat(net.xmin, N, 1);
%xomax = net.xomax;
%xomin = net.xomin;
%yomin = net.yomin;
%yomax = net.yomax;

% Min-max normalization
%Pc = (((P - xmin) ./ (xmax - xmin)) * (ymax - ymin)) + ymin;

% Add column of ones
Pc = P;
Pc = [ones(size(Pc, 1), 1) Pc];

% Hidden layer calculation    
if strcmp(net.activation{1}, 'tansig')
    Z = tansig(Pc * net.IW');
elseif strcmp(net.activation{1}, 'logsig')
    Z = logsig(Pc * net.IW');
elseif strcmp(net.activation{1}, 'purelin') 
    Z = Pc * net.IW';
else
    error('ERROR: Incorrect activation function for hidden layer!');
end

% Output layer calculation
Zc = [ones(size(Z, 1), 1) Z];
if strcmp(net.activation{2}, 'tansig')
    output = tansig(Zc * net.LW');
elseif strcmp(net.activation{2}, 'logsig')
    output = logsig(Zc * net.LW');
elseif strcmp(net.activation{2}, 'purelin') 
    output = Zc * net.LW';
else
    error('ERROR: Incorrect activation function for output layer!');
end

% Normalize between 0 and 1
%multi = diag(xomax' - xomin');
%addition = repmat(xomin', N, 1);
%fx = (((output - yomin) / (yomax - yomin)) * multi) + addition;
fx = output;
    
return
