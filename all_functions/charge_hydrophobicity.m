function ratio = charge_hydrophobicity (sequence)

% function ratio = charge_hydrophobicity (sequence)
%
% The function determines ratio between net charge and hydrophobicity
% for an input amino acid sequence
%
% Predrag Radivojac
% Indiana University
% May 2004


AminoAcids = 'ACDEFGHIKLMNPQRSTVWY';

% Kyte-Doolittle scale
KDscale = [1.8 2.5 -3.5 -3.5 2.8 -0.4 -3.2 4.5 -3.9 3.8 1.9 -3.5 -1.6 -3.5 -4.5 -0.8 -0.7 4.2 -0.9 -1.3];

% mean net charge
net_charge = (length(find(sequence == 'K' | sequence == 'R')) - length(find(sequence == 'D' | sequence == 'E'))) / length(sequence);

% mean hydrophobicity using KD scale
hydrophobicity = 0;
for i = 1 : 20
    q = find(sequence == AminoAcids(i));
    hydrophobicity = hydrophobicity + length(q) * KDscale(i);
end
hydrophobicity = hydrophobicity / length(sequence);

% determine unfoldedness
%if hydrophobicity * 2.785 < net_charge + 1.151
%    unfolded = 1;
%else
%    unfolded = 0;
%end

if abs(hydrophobicity) > 0.01
    ratio = net_charge / hydrophobicity;
else
    ratio = 100;
end

return
