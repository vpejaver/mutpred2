% Vikas Pejaver
% October 2014

% Function 

function pprof = get_para_prof(sequence)

% Global variable
global CURRDIR;

% Constants and defaults
pprof = [];
run_flag = 0;
mat_file =  strcat(CURRDIR, filesep, 'data', filesep, 'fulldata_human_mouse_sp_r2014_09.mat');

% Load MAT file
load(mat_file);

% Get the profile
inds = find(strcmp(sequences, sequence));
if length(inds) >= 1
    index = inds(1); % If more than one match, return the first one
                     % as the profiles would be the same anyways
    %pprof = para_profs(l_inds(index), :);
    pprof = para_profs(index, :);
else
    run_flag = 1;
end

% If not found perform sequence similarity search
if run_flag == 1
    %'RUNNIN'
    pprof = return_paralogy_profile(sequence, human_seqs);
    pprof = [pprof return_paralogy_profile(sequence, mouse_seqs)];
end

return
