function pred = vhBfactor(seq)

AA = 'ACDEFGHIKLMNPQRSTVWY';

% David Smith's flexibility scale
% FS = [-0.605 -0.693 -0.279 -0.160 -0.719 -0.537 -0.662 -0.682 -0.043 -0.631 -0.626 -0.381 -0.271 -0.369 -0.448 -0.423 -0.525 -0.669 -0.727 -0.721];
% Vihinen's flexibility scale
FS = [0.984 0.906 1.068 1.094 0.915 1.031 0.950 0.927 1.102 0.935 0.952 1.048 1.049 1.037 1.008 1.046 0.997 0.931 0.904 0.929];

% David Smith's 
% coef = [0.2 0.4 0.6 0.8 1 0.8 0.6 0.4 0.2];
% Vihinen's
coef = [0.25 0.4375 0.625 0.8125 1 0.8125 0.625 0.4275 0.25];

N = length(seq);
W = length(coef); 

flexseq = [];
for i = 1 : 20
    idx = find(seq == AA(i)); 
    flexseq(idx) = FS(i);
end

pred = zeros(N, 1);

for j = 1 : length(seq)
    
    s = max(1, j - 4);
    e = min(N, j + 4);
    
    if j <= 4
        c = coef(5 - j + 1 : W);
    elseif j > N - 4  
        c = coef(1 : W - 4 + N - j);
    else
        c = coef;
    end
    
    c = c * sum(coef) / sum(c);        
    
    pred(j) = c * flexseq(s : e)';
end

return
