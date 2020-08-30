function data = makeattr(seq, win)

AA = 'WCFIYVLHMATRGQSNPDEK';   

% Hydropathicity., Author:  Kyte J., Doolittle R.F., Reference: J. Mol. Biol. 157:105-132(1982).
AAS = [-0.900, 2.500, 2.800, 4.500, -1.300, 4.200, 3.800, -3.200, 1.900, 1.800, -0.700, -4.500, -0.400, -3.500, -0.800, -3.500, -1.600, -3.500, -3.500, -3.9];
% mi = min(AAS);
% ma = max(AAS);
% AAS = (AAS - mi) / (ma - mi);

N = length(seq);
W = round((win - 1) / 2);

freqseq = zeros(N, 20);
aaseq = zeros(N, 1);
for a = 1 : 20
    idx = find(seq == AA(a)); 
    freqseq(idx, a) = 1;
    aaseq(idx) = repmat(AAS(a), length(idx), 1);
end

spacer = [ones(W, 1); zeros(N, 1); ones(W, 1)];

data = [];

for n = 1 : N
    
    s = max(1, n - W);
    e = min(N, n + W);
    L = e - s + 1;

    % 1-20 : AA frequencies
    xf = mean(freqseq(s : e, :));
    
    % 21 : K2-entropy
    i = find(xf ~= 0);
    xent = - sum(xf(i) .* log2(xf(i)));

    % 22 : spacer frequency
    xsp = mean(spacer(n : n + win - 1));

    % 24 : mean hydropathy
    xs = mean(aaseq(s : e, :));

    % 25 : net charge
    nc = xf(20) + xf(12) - xf(18) - xf(19);
    
    % 25 : charge/hydropathy ratio
    if abs(xs) > 0.01
        nh = nc / xs;
    else
        nh = 100;
    end
    
    data = [data; [xf xent xsp nc xs nh]];
end

% 26 : Vihinen's flexibility (B-factor) prediction
bp = vhBfactor(seq);
data = [data bp];

return

