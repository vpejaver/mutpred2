function v = get_aa_volumes (sequence)

% Funtion calculates the volume of each residue in the protein sequence
%  
% Input arguments:
%       sequence = string of amino acids
% 
% Output arguments:
%        v = vector of numbers, each corresponding to the volume of a
%        residue
%  
% Residues 'B' and 'Z' are treated as 'D' and 'E', respectively. For residue
% 'X' an average over 20 residue volumes is used.
%  
% 
% Predrag Radivojac (predrag@indiana.edu)
% Indiana University, Bloomington, Indiana
% November 2007

sequence == upper(sequence);

if ~isempty(find(sequence ~= 'A' & sequence ~= 'C' & sequence ~= 'D' & sequence ~= 'E' & sequence ~= 'F' & ...
        sequence ~= 'G' & sequence ~= 'H' & sequence ~= 'I' & sequence ~= 'K' & sequence ~= 'L' & ...
        sequence ~= 'M' & sequence ~= 'N' & sequence ~= 'P' & sequence ~= 'Q' & sequence ~= 'R' & ...
        sequence ~= 'S' & sequence ~= 'T' & sequence ~= 'V' & sequence ~= 'W' & sequence ~= 'Y' & ...
        sequence ~= 'B' & sequence ~= 'Z' & sequence ~= 'X'))
    error('Illegal character in amino acid sequence');
end

AminoAcids = 'ARNDCQEGHILKMFPSTWYVBZX';

%  Residue volume (Goldsack-Chalifoux, 1973)
Scales = [88.3 181.2 125.1 110.8 112.4 148.7 140.5 60.0 152.6 168.5 168.5 175.6 ...
    162.2 189.0 122.2 88.7 118.2 227.0 193.0 141.4 110.8 140.5 143.7];

v = zeros(1, length(sequence));

len = length(AminoAcids);
for i = 1 : len
    q = find(sequence == AminoAcids(i));
    v(q) = Scales(i);
end

return
