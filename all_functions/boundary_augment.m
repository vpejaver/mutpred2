function [soft_prediction] = boundary_augment (soft_prediction, sequence, boundary_type)

% boundary_type = 'unidirectional' - forward prediction (for boundaries)
% boundary_type = 'bidirectional' - forward-backward prediction (for boundaries)
% boundary_type = 'no boundary' - no boundary prediction is used

K = 0.06; N = 12;
% K = 0.16; N = 5;

% form quantized prediction
qpred = zeros(1, length(soft_prediction));
q = find(soft_prediction >= 0.5);
qpred(q) = 1;

% determine disorder boundaries
pos = 1;
dis_pos = [];
while pos < length(qpred)
    while pos < length(qpred) & qpred(pos) == 0
        pos = pos + 1;
    end
    if pos >= length(qpred) 
        break
    end
    dis_pos(length(dis_pos) + 1) = pos;
    while pos < length(qpred) & qpred(pos) == 1
        pos = pos + 1;
    end
    dis_pos(length(dis_pos) + 1) = pos;
end

% calculate where to reverse inputs for FB boundary predictor
reversing = zeros(1, length(soft_prediction));
for j = 2 : 2 : length(dis_pos)
    if j == length(dis_pos)
        reversing(floor(dis_pos(j - 1) + .67 * (dis_pos(j) - dis_pos(j - 1))) : length(soft_prediction)) = 1;
    else
        reversing(floor(dis_pos(j - 1) + .67 * (dis_pos(j) - dis_pos(j - 1))) : floor(dis_pos(j) + .5 * (dis_pos(j + 1) - dis_pos(j)))) = 1;
    end
end

% predict boundaries
[boundF, boundR, boundFR] = predict_boundaries(sequence, reversing);

% fit soft_prediction so that it crosses 0.5 at the predicted boundary
[ba_soft_predictionF ba_soft_predictionFR] = adjust_lin_bagg(soft_prediction, boundF, boundFR, K, N);

% output what user wants
if ~isempty(strmatch(boundary_type ,'unidirectional','exact'))
    soft_prediction = ba_soft_predictionF;
elseif ~isempty(strmatch(boundary_type ,'bidirectional','exact'))
    soft_prediction = ba_soft_predictionFR;
elseif isempty(strmatch(boundary_type ,'no boundary','exact'))
    error('Wronf format of input parameters...');
end
    
return
