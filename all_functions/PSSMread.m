function [pssm, seq] = readPSSM(fname)

if exist(fname) == 0
    error('file %s doesn''t exist!!!', fname);
end

f = fopen(fname, 'rb');
ln = fgetl(f);
ln = fgetl(f);
ln = fgetl(f);
buf = fread(f);
fclose(f);
buf = char(buf');

i = findstr(buf, '                      K         Lambda');
buf = buf(1 : i - 2);

p = find(buf == 10);
p = [0 p];

seq = [];
pssm = [];
for m = 1 : length(p) - 1
    seq = [seq buf(p(m) + 7)];
    ln = buf(p(m) + 9 : p(m + 1) - 1);
    a = str2num(ln);
    pssm = [pssm; a];
end

return
