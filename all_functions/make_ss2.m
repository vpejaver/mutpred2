function D1 = make_ss2(D, w1, w2, xmin, xmax)
% create features for the second network

windows_aa = 15;
e = (windows_aa + 1) / 2;
p = sim2(D, w1, w2, xmin, xmax);
[nr nc] = size(p);                      % nc == 4
D1 = zeros(nr, windows_aa * nc);
for residue = 1:nr
    for j = 1:windows_aa
        rel_ix = j + residue - e;
        if rel_ix < 1 || rel_ix > nr
            continue
        end
        D1(residue, (1+(j-1)*nc):(j*nc)) = p(rel_ix, :);
    end
end
