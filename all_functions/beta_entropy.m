function H = beta_entropy (p, beta)

% function H = beta_entropy (p, beta)
%
% Daroczy's beta-entropy (converges to Shannon's entropy when beta goes to 1)
%
% Input:
%        p = vector of probabilities (non-negative numbers that sum to one; not checked)
%        beta = parameter, within (1, 2]
%
% Output:
%         H = beta-entropy
%
% Reference:
%            Z. Daroczy. Generalized information functions. Information and Control, 16: 36-51, 1970.
%
% Predrag Radivojac
% Indiana university School of Medicine
% July 2004

H = 1 / (1 - 2 ^ (1 - beta)) * (1 - sum(p .^ beta));

return
