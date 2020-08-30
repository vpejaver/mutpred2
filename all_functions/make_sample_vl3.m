function sample = make_sample_vl3 (sequence, window)

AminoAcids = 'ACDEFGHIKLMNPQRSTVWY';

sample = [];

for i = 1 : length(sequence)
    x = [];
    for k = 1 : 20
        q = find(sequence(max(1, i - floor(window / 2)) : min(length(sequence), i + floor(window / 2))) == AminoAcids(k));
        x(k) = length(q) / length(sequence(max(1, i - floor(window / 2)) : min(length(sequence), i + floor(window / 2))));
    end
    x(21) = entropy(x);
    sample(i, :) = x;
end

% don't use 20th amino acid frequency
sample = sample(:, [1 : 19 21]);

return