function [soft_prediction] = VL3 (sequence, boundary_type, window_out, clip_type);

% sequence = string of amino acids (capital letters)
% boundary_type = 'unidirectional' - forward prediction (for boundaries)
% boundary_type = 'bidirectional' - forward-backward prediction (for boundaries)
% boundary_type = 'no boundary' - no boundary prediction is used
% window_out = an odd numbered two-sided window for postfiltering
% clip_type = 'clip' - prediction is bounded to [0, 1]
% clip_type = 'no clip' - prediction go below 0 and above 1

global CURRDIR;

examples = make_sample_vl3(sequence, 41);

%if strmatch(version('-release'), '2006b', 'exact')
%    load(strcat(CURRDIR, filesep, 'all_models', filesep, 'model41vl3_R2006b.mat'));
%else
%    %load model41vl3v6.mat
%    load(strcat(CURRDIR, filesep, 'all_models', filesep, 'model41vl3.mat'));
%end
load(strcat(CURRDIR, filesep, 'all_models', filesep, 'model41vl3_compilable.mat'));

T = length(net);

for t = 1 : T
    % normalize and calculate prediction
    %tmp_prediction(:, t, :) = sim(net{t}, normalize_dataset(examples, meanv{t}, stdv{t})')';
    tmp_prediction(:, t, :) = simple_predict_nn(net{t}, normalize_dataset(examples, meanv{t}, stdv{t}));
end

prediction = sum(tmp_prediction(:, 1 : T, :), 2) / T;

X(:, 1) = moving_average(prediction(:, 1), window_out);
X(:, 2) = moving_average(prediction(:, 2), window_out);

% construct soft_prediction
soft_prediction = 0.5 + 0.5 * (X(:, 2) - X(:, 1))';

% rescale it into 0-1 interval (from 0.1-0.9)
soft_prediction = (soft_prediction - 0.1) / 0.8;

if ~isempty(strmatch(boundary_type ,'unidirectional','exact')) | ~isempty(strmatch(boundary_type ,'bidirectional','exact'))
    soft_prediction = boundary_augment(soft_prediction, sequence, boundary_type);
elseif isempty(strmatch(boundary_type ,'no boundary','exact'))
    error('Wrong format of variable boundary_type');
end

if ~isempty(strmatch(clip_type ,'clip','exact'))
    % clip if > 1 to 1 or if < 0 to 0
    q = find(soft_prediction > 1);
    if ~isempty(q)
        soft_prediction(q) = 1;
    end
    
    q = find(soft_prediction < 0);
    if ~isempty(q)
        soft_prediction(q) = 0;
    end
elseif isempty(strmatch(clip_type ,'no clip','exact'))
    error('Invalid type of clipping');
end

return