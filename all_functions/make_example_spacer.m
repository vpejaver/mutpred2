function example = make_example_spacer (sequence)

AminoAcids = 'ACDEFGHIKLMNPQRSTVWY ';

for j = 1 : length(sequence)
    q = find(AminoAcids == sequence(j));
    w = zeros(1, 21);
    w(q) = 1;
    example(1, (j - 1) * 21 + 1 : j * 21) = w;
end

return

