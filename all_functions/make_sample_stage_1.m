function [D] = make_sample_stage_1 (sequence)

windows_aa = [11 21];
windows_pr = [1 7 11 21];
AminoAcids = 'ACDEFGHIKLMNPQRSTVWY';

q = find(sequence == 'X');
if ~isempty(q)
    sequence(q) = 'T';
end

[p{1}, p{2}, p{3}, p{4}] = VL2(sequence);
q = find(p{1} > 1); p{1}(q) = 1; q = find(p{1} < 0); p{1}(q) = 0;
q = find(p{2} > 1); p{2}(q) = 1; q = find(p{2} < 0); p{2}(q) = 0;
q = find(p{3} > 1); p{3}(q) = 1; q = find(p{3} < 0); p{3}(q) = 0;
q = find(p{4} > 1); p{4}(q) = 1; q = find(p{4} < 0); p{4}(q) = 0;
p{5} = vihinen(sequence);
p{6} = hydrophobic_moment(sequence, 11, 100);
p{7} = hydrophobic_moment(sequence, 11, 160);
p{8} = hydrophobic_moment(sequence, 11, 120);
p{9} = predictBfactors(sequence, 5, 5);
if length(sequence) >= 30 % vlxt doesn't work for short sequences
    p{10} = VLXT(sequence);
else
    p{10} = p{1};
end
p{11} = VL3(sequence, 'no boundary', 1, 'clip');

vect_aa = [];
vect_pr = [];

for residue = 1 : length(sequence)
    % features made from amino acid sequence
    for w = 1 : length(windows_aa)
        % determine margins for the windows (important for the regions near ends
        left_margin = max(1, residue - (windows_aa(w) - 1) / 2);
        right_margin = min(residue + (windows_aa(w) - 1) / 2, length(sequence));
        
        fragment = sequence(left_margin : right_margin);
        
        % 20 aa relative frequencies
        for j = 1 : 20
            vect_aa(residue, 24 * (w - 1) + j) = length(find(fragment == AminoAcids(j))) / length(fragment);
        end
        % entropy
        vect_aa(residue, 24 * (w - 1) + 21) = entropy(vect_aa(residue, 1:20));
        
        % net and total charge
        pos_charge = length(find(fragment == 'K' | fragment == 'R'));
        neg_charge = length(find(fragment == 'D' | fragment == 'E'));
        vect_aa(residue, 24 * (w - 1) + 22) = (pos_charge - neg_charge) / length(fragment);
        vect_aa(residue, 24 * (w - 1) + 23) = (pos_charge + neg_charge) / length(fragment);
        
        % aromatics
        vect_aa(residue, 24 * (w - 1) + 24) = length(find(fragment == 'F' | fragment == 'Y' | fragment == 'W')) / length(fragment);
    end
    
    % features made from different predictors
    tmp = [];
    for w = 1 : length(windows_pr)
        % determine margins for the windows (important for the regions near ends
        left_margin = max(1, residue - (windows_pr(w) - 1) / 2);
        right_margin = min(residue + (windows_pr(w) - 1) / 2, length(sequence));
        
        for j = 1 : length(p)
            tmp = [tmp mean(p{j}(left_margin : right_margin))];
        end
    end
    vect_pr(residue, :) = tmp;
end

D = [vect_aa vect_pr];

return
