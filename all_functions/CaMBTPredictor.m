function [res_pred, reg_pred] = CaMBTPredictor (sequence)

% function [res_pred, reg_pred] = CaMBTPredictor (sequence)
%
% The function predicts calmodulin (CaM) binding regions (targets) from amino acid sequence
%
% Input:
%        sequence = amino acid sequence (string)
%
% Output:
%        res_pred = approximated probability that a residue is calmodulin binding (vector)
%        reg_pred = consecutive 1's for predicted CaM binding regions, 0's otherwise (vector)
%
% Note: if only a residue based prediction is needed, the function should be called as
%
%       res_pred = CaMBTPredictor(sequence);
%
%       to save on detecting the CaM binding regions, which approximately doubles prediction time.
%
% Predrag Radivojac
% Indiana university
% June 2004

sequence = upper(sequence);

% allow only 20 amino acid codes together with Z, B, and X (ZBX will later be converted into T's)
if ~isempty(find(sequence ~= 'A' & sequence ~= 'C' & sequence ~= 'D' & sequence ~= 'E' & sequence ~= 'F' & ...
                 sequence ~= 'G' & sequence ~= 'H' & sequence ~= 'I' & sequence ~= 'K' & sequence ~= 'L' & ...
                 sequence ~= 'M' & sequence ~= 'N' & sequence ~= 'P' & sequence ~= 'Q' & sequence ~= 'R' & ...
                 sequence ~= 'S' & sequence ~= 'T' & sequence ~= 'V' & sequence ~= 'W' & sequence ~= 'Y' & ... 
                 sequence ~= 'X' & sequence ~= 'B' & sequence ~= 'Z'))
    res_pred = [];
    reg_pred = [];
    warning('Unwanted characters in input sequence: no prediction is made');
    return
end

% change all B's into D's
q = find(sequence == 'B');
sequence(q) = 'D';
% change all Z' into E's
q = find(sequence == 'Z');
sequence(q) = 'E';
% randomly assign residues to all positions where X occurs
q = find(sequence == 'X');
if ~isempty(q)
    array = [65 67 68 69 70 71 72 73 75 76 77 78 80 81 82 83 84 86 87 89];
    qq = rand(1, length(q));
    qq = 1 + round(qq * 19);
    for i = 1 : length(q)
        sequence(q(i)) = char(array(qq(i)));
    end
end

% predict the sites on a residue basis
res_pred = get_residue_prediction(sequence);

% if the number of output arguments is 2
if nargout == 2
    % get final region prediction if a user wants it
    reg_pred = get_region_prediction(res_pred, sequence);
end

return
