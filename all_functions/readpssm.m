function pssm = readpssm(pssmFileName)
% readpssm - read PSSM file and extract the 42 x length matrix
% Input:
%     pssmFileName  Name of PSSM file; PSSM block for each sequence begins with
%         a blank line and a comment line followed by actual result and ends with a blank line
%         followed by 5 lines of comments.
% Output:
%     pssm  A 42 x length matrix
%
% Last modified: Wed Apr 14 10:10:50 2010

pssm = [];
if ~exist(pssmFileName, 'file')
    return
end

ix = 1;
fid = fopen(pssmFileName, 'r');
while ~feof(fid)
    tline = fgetl(fid);
    % skip the first 3 lines
    if ix < 4
        ix = ix + 1;
        continue
    end

    stripLine = regexp(tline, '\s+', 'split');
    n = length(stripLine);
    if n == 1 % blank line
        fclose(fid);
        return
    end

    stripLine = str2num(char(stripLine{n-41:n}));
    pssm = [pssm stripLine];
end
