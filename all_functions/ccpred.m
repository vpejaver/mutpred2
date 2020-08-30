% Vikas Pejaver
% July 2013

% Function to predict coiled-coil regions in a given protein
% Default mode = 1: outputs a n x 64 matrix (the MARCOIL model has 64
% states) and also provides that state names (Refer to Delorenzi
% and Speed, 2002 for state definitions)
% When mode = 0: outputs an n x 1 vector (probability of
% coiled-coil occurring at position)

function [posteriors, states] = ccpred(sequence, varargin)

% Global variables
global EST_EMISSION EST_TRANSITION HIDDEN_STATES OBSERVED_STATES;

% Constants and defaults
posteriors = [];

%model_mat = 'marcoil_like_072313.mat';
if isempty(varargin)
    out = 1;
elseif varargin{1} ~= 0 && varargin{1} ~= 1
    error('ERROR: Output mode can only be 0 or 1!');
else
    out = varargin{1};
end

% Load model
%load(CC_MODEL);

% Clean non-standard amino acids
seq = clean_protein(sequence);

% Make prediction
preds = hmmdecode(seq, EST_TRANSITION, EST_EMISSION, 'symbols', char(OBSERVED_STATES)');

% Depending on output style, output prediction 
if out == 1
    posteriors = preds';
    states = HIDDEN_STATES';
else
    posteriors = 1 - preds(1, :)';
    states = {};
end

return
