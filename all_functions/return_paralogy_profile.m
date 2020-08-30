function [profile] = return_paralogy_profile (s, db)

% profile = return_paralogy_profile (s, db)
% Input: sequence s, string
%        database db, cell array of strings
%        Both s and db must come from a 20-amino acid alphabet
%
% Output: profile, an array of integers containint numbers of sequences
%         over a series of thresholds
%
% Predrag Radivojac
% Indiana University
% October 2014

lng = length(s);
ldb = length(db);
thrs = [.5 .55 .6 .65 .7 .75 .8 .85 .9 .95];

% load BLOSUM62 matrix for sequence alignments
b62 = load('b62.txt');

% minimum threshols below which one need not align proteins
mnthr = min(thrs);


qq = find(s == 'U');
s(qq) = 'C';
qq = find(s == 'X');
s(qq) = 'A';
qq = find(s == 'B');
s(qq) = 'D';
qq = find(s == 'Z');
s(qq) = 'E';


% vector of sequence identites above minimum threshold
sid = zeros(1, ldb);

%c=0;

% perform sequence alignments that are necessary
for i = 1 : ldb
    if min([lng length(db{i})]) >= mnthr * max([lng length(db{i})]) & sum(min(composition_profile(s), composition_profile(db{i}))) / max([lng length(db{i})]) >= mnthr
    %if min([lng length(db{i})]) >= mnthr * max([lng length(db{i})])
        %c = c + 1;
        ss = db{i};
        
        qq = find(ss == 'U');
        ss(qq) = 'C';
        qq = find(ss == 'X');
        ss(qq) = 'A';
        qq = find(ss == 'B');
        ss(qq) = 'D';
        qq = find(ss == 'Z');
        ss(qq) = 'E';

        [~, sq1, sq2] = protglobal(s, ss, b62, -1, -11);
        sid(i) = length(find(sq1 == sq2)) / max([lng length(db{i})]);
    end
end

%fprintf(1, '%d\t', c);

% vector to be returned
profile = zeros(1, length(thrs));
for i = 1 : length(thrs)
    profile(i) = length(find(sid >= thrs(i)));
end

return
