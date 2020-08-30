% function [d seq1 seq2] = protglobal (p1, p2, S, g, fg)
%
% Funtion calculates optimal gobal alignment of two proteins, p1 and p2.
%
% Input arguments:
%       p1 = string of amino acids
%       p2 = string of amino acids
%       S = 20 x 20 scoring matrix
%       g = gap penalty (should be a negative number)
%       fg = penalty for the first gap (should be a negative number)
%
% Output arguments:
%        d = score of the optimal alignment
%        seq1 = p1, aligned to p2 (possibly containing gaps)
%        seq2 = p2, aligned to p1 (possibly containing gaps)
%
% Amino acids corresponding to rows and columns of S must be in the following order:
% "CSTPAGNDEQHRKMILVFYW". Matrix S should be symmetric. The function does NOT test symmetry.
% Case of the characters in p1 and p2 does not influence the final alignment.
% Letters 'B' and 'Z' are treated as 'D' and 'E', respectively. Every amino acid
% scores -1 if aligned with X. If there are more than one alignments with the same
% score, the function returnes only one of them.
%
% Example: p1 = 'ACCALI'
%          p2 = 'ACQVI'
%          S = BLOSUM62
%          g = -4
%          fg = -12
%
%          The function returns:
%
%          d = 5
%          seq1 = ACCALI
%          seq2 = A-CQVI
%
% For information about BLOSUM62 matrix see S. Henikoff, J. Henikoff, "Amino acid
% substitution matrices from protein blocks," Proc. Natl. Acad. Sci. USA, Vol 89.,
% pp. 10915-10919, November 1992, Biochemistry.
%
% Predrag Radivojac (predrag@ist.temple.edu)
% Temple University, Philadelphia, Pennsylvania
