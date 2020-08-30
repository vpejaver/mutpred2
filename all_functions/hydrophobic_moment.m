function [hm] = hydrophobic_moment (sequence, window, angle)

% function [hm] = hydrophobic_moment (sequence, window, angle)
% 
% Calculates the hydrophobic (amphiphilic) moment for a protein sequence
%
% Input: sequence: protein sequnce (string from the alphabet of 20 symbols)
%        window  : window size (a number)
%        angle   : rotation angle in degrees (100 degrees for helix)
%
% Output: hm     : vector of hydrophobic moments
%
% Comment: the code does not handle illegal characters in the input sequence
%
% Predrag Radivojac
% Indiana University
% May 2004

AminoAcids = 'ACDEFGHIKLMNPQRSTVWY';

% normalized concensus scale values for 20 amino acids (from Eisenberg's Nature paper 1982)
% A =  0.25; L =  0.53; R = -1.76; K = -1.10;
% N = -0.64; M =  0.26; D = -0.72; F =  0.61;
% C =  0.04; P = -0.07; Q = -0.69; S = -0.26;
% E = -0.62; T = -0.18; G =  0.16; W =  0.37;
% H = -0.40; Y =  0.02; I =  0.73; V =  0.54;

% initialize stuff
hydrophobic_indices = [.25 .04 -.72 -.62 .61 .16 -.40 .73 -1.10 .53 .26 -.64 -0.07 -.69 -1.76 -.26 -.18 .54 .37 .02];
half_window = floor(window / 2);

sequence = upper(sequence); % convert a sequence into uppercase sequence
angle = angle * pi / 180;   % convert an angle into radians

% q is an array of indices corresponding to particular amino acids in the whole sequence
q = zeros(1, length(sequence));
for i = 1 : 20
    qq = find(sequence == AminoAcids(i)); %find all positions with amino acid i
    q(qq) = i;
end

% calculate hydrophobic indices over the input sequence
for i = 1 : length(sequence)
    hm_sin = 0;
    hm_cos = 0;
    for j = max(1, i - half_window) : min(i + half_window, length(sequence))
        hm_sin = hm_sin + sin(angle * j) * hydrophobic_indices(q(j));
        hm_cos = hm_cos + cos(angle * j) * hydrophobic_indices(q(j));
    end
    
    % hydrophobic moment is normalized with the window length
    if min(i + half_window, length(sequence)) - max(1, i - half_window) + 1 < window
        hm(i) = sqrt(hm_sin ^ 2 + hm_cos ^ 2) / (min(i + half_window, length(sequence)) - max(1, i - half_window) + 1); 
    else
        hm(i) = sqrt(hm_sin ^ 2 + hm_cos ^ 2) / window;
    end
end

return