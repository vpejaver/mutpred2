% Vikas Pejaver
% January 2015

% Function to return remarks given a binary vector

function comments = get_remarks(vec)

%%%%%% Constants and defaults %%%%%%
comments = '';
Remark_list = {'Sequence length less than 30 residues', 'Substitution in incorrect format', 'Position lies outside the sequence', 'Reference residue does not match the one at the position in the sequence', 'Reference residue and substituted residue are the same', 'Non-standard amino acids found within 10 residues', 'No PSSM', 'No conservation scores', 'Predicted conservation scores'};


%%%%%% Make remarks for output %%%%%%
comments = regexprep(sprintf('%s | ', Remark_list{find(vec)}), ' \| $', '');

return
