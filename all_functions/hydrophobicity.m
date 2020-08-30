function h = hydrophobicity (sequence, window)

AminoAcids = 'ACDEFGHIKLMNPQRSTVWYX';

% Kyte-Doolittle scale
KDscale = [1.8 2.5 -3.5 -3.5 2.8 -0.4 -3.2 4.5 -3.9 3.8 1.9 -3.5 -1.6 -3.5 -4.5 -0.8 -0.7 4.2 -0.9 -1.3 -0.49];

h = zeros(1, length(sequence));
for i = 1 : 21
    h(find(sequence == AminoAcids(i))) = KDscale(i);
end

if window > 1
    h = moving_average(h, window);
end

return