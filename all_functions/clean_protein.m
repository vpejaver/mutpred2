% Function to get rid of ambiguous amino acids
% Date: 09/21/2010

function seq = clean_protein(seq)
    q = find(seq == 'X');
    if ~isempty(q)
        seq(q) = 'A';
    end
    q = find(seq == 'B');
    if ~isempty(q)
        seq(q) = 'D';
    end
    q = find(seq == 'Z');
    if ~isempty(q)
        seq(q) = 'E';
    end
    q = find(seq == 'U');
    if ~isempty(q)
        seq(q) = 'C';
    end
    q = find(seq == 'O');
    if ~isempty(q)
        seq(q) = 'K';
    end
    return
