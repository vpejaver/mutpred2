function [ss] = ssp4(seq, pssm)
% Make secondary structure prediction with given protein sequence with or without PSSM.
% 
% Args:
%   seq    input protein sequence with n (>= 30) amino acids
%   pssm   matrix of profile from running blastpgp; dimension 42 x n
% 
% Returns:
%   A n x 4 matrix representing probabilities for 4 possible structures; From left to
%   right the columns correspond to Helix (H), Strand (E), Loop (L) and Disorder (D), 
%   respectively.
% 
% Example (PDB ID: 1JKV_A):
% test_seq = ['MFKHTRKLQYNAKPDRSDPIMARRLQESLGGQWGETTGMMSYLSQGWASTGAEKYKDLLL'...
%             'DTGTEEMAHVEMISTMIGYLLEDAPFGPEDLKRDPSLATTMAGMDPEHSLVHGLNASLNN'...
%             'PNGAAWNAGYVTSSGNLVADMRFNVVRESEARLQVSRLYSMTEDEGVRDMLKFLLARETQ'...
%             'HQLQFMKAQEELEEKYGIIVPGDMKEIEHSEFSHVLMNFSDGDGSKAFEGQVAKDGEKFT'...
%             'YQENPEAMGGIPHIKPGDPRLHNHQG'];
% ss = ssp4(test_seq);                       % without PSSM
% The first three rows look like:
%     0.0693    0.1204    0.1119    0.8046
%     0.4342    0.1270    0.1077    0.2859
%     0.4769    0.1034    0.1269    0.2834
% 
% ss = spp4(test_seq, [])                     % call blastpgp to get PSSM
% 
% ss = ssp4(test_seq, test_pssm);             % with PSSM
% Again the first three rows look like:
%     0.0650    0.1031    0.0769    0.8098
%     0.0693    0.1920    0.0802    0.7001
%     0.1935    0.4140    0.1283    0.1577
% 
% And the secondary structure from PDB for the first three positions is DEE.
% 
% When reading sequences and PSSMs from disk the following functions can help
% 
%   [seq, h] = readFASTA('seq_file_name');
%   pssm = readpssm('pssm_file_name');
% 
% Last modified: Wed Feb 23 17:14:57 2011

global CURRDIR;

% input sequence must be more than 30 AAs.
if length(seq) < 30
    fprintf(1, ['Input sequence must contain at least 30 amino acids. Exit without ' ...
                'prediction.\n'])
    return
end

% standardize protein sequence
seq = stdaa(seq);

% --**************************************************************--
% --* change the following line when moving ssp4 to other places *--
% --**************************************************************--
localpath = strcat(CURRDIR, filesep, 'all_models'); %fileparts(which('ssp4'));

% read network parameters
sprintf('%s%s%s', localpath, filesep, 'ssp4_para.bin');
fid = fopen(sprintf('%s%s%s', localpath, filesep, 'ssp4_para.bin'));
par = fread(fid, [194 329], 'double');
fclose(fid);

% parameters for individual networks
nn1_nopssm_xmin = par(1, 1:164);
nn1_nopssm_xmax = par(2, 1:164);
nn1_nopssm_wb1  = par(3:37, 1:165);
nn1_nopssm_wb2  = par(38:41, 1:36);

nn2_nopssm_xmin = par(42, 1:60);
nn2_nopssm_xmax = par(43, 1:60);
nn2_nopssm_wb1  = par(44:53, 1:61);
nn2_nopssm_wb2  = par(54:57, 1:11);

nn3_pssm_xmin = par(58, 1:328);
nn3_pssm_xmax = par(59, 1:328);
nn3_pssm_wb1  = par(60:124, 1:329);
nn3_pssm_wb2  = par(125:128, 1:66);

nn4_pssm_xmin = par(129, 1:60);
nn4_pssm_xmax = par(130, 1:60);
nn4_pssm_wb1  = par(131:190, 1:61);
nn4_pssm_wb2  = par(191:194, 1:61);

% create features for the first network
if nargin == 1
    % make prediction without pssm
    f = make_ss_dataset_without_pssm(seq); % dimension: n x 168

    % remove rows with constant values. From
    %   net1.inputs{1}.processSettings{2}.remove
    keep = setdiff(1:168, [21 22 23 24]);
    f = f(:, keep);
    
    % output for the first network
    f = make_ss2(f, nn1_nopssm_wb1, nn1_nopssm_wb2, nn1_nopssm_xmin, nn1_nopssm_xmax);

    % prediction from the second network
    ss = sim2(f, nn2_nopssm_wb1, nn2_nopssm_wb2, nn2_nopssm_xmin, nn2_nopssm_xmax);
    
    return
elseif isempty(pssm)
    % call blastpgp to generate PSSM
    sfn = sprintf('%s/pssm/seq.txt', localpath);
    fid = fopen(sfn, 'w');
    fprintf(fid, '>seq\n%s\n', seq);
    fclose(fid);
    cmd1 = ['/tools/BLAST/blast-2.2.18/bin/blastpgp -d /data_magic/filternr/filtnr ' ...
            '-h 0.0001 -j 3 -o /dev/null'];
    pfn = sprintf('%s/pssm/pssm.txt', localpath);
    s = system([sprintf('%s -i %s -Q %s', cmd1, sfn, pfn)]);
    if s ~= 0
        error('Can not run blastpgp to generate PSSM.')
    end
    pssm = readpssm(pfn);
    delete(sfn);
    delete(pfn);
else
    % nothing to do
end

% make prediction with pssm
f = make_ss_dataset3(seq, pssm); % dimension: n x 328

% output for the first network
f = make_ss2(f, nn3_pssm_wb1, nn3_pssm_wb2, nn3_pssm_xmin, nn3_pssm_xmax);

% prediction from the second network
ss = sim2(f, nn4_pssm_wb1, nn4_pssm_wb2, nn4_pssm_xmin, nn4_pssm_xmax);

return
