% Vikas Pejaver
% October 2014

% Function to parse command-line arguments

function ps = parse_args(args)

% Global variables
global IN_FILE OID DO_HOMOLOGY PREDICT_CONS SKIP_PSIBLAST CSV_FORM ONTO_FORM PTHR; %SCORES_ONLY TOPK;

% Constants and defaults
ps = 1;
Valid_p = [0, 1];
Valid_c = [0, 1];
Valid_b = [0, 1];
%Valid_s = [0, 1];
%Valid_k = [1, 57];
Valid_t = [0, 1];
Valid_f = [1:4];
%Flag_list = 'iopcbsk';
Flag_list = 'iopcbstf';
Sing_tab = sprintf('\t');
Padding1 = sprintf('\t\t\t\t\t\t\t\t');
Padding2 = sprintf('\t\t\t');

% Convert everything to strings
args = cellfun(@num2str, args, 'uniformoutput', false); % TODO : not necessary for the final executable

% Count number of arguments
total_count = length(args);
flag_count = length(find(~cellfun(@isempty, regexp(args, '^-', 'match'))));
help_flag = any(~cellfun(@isempty, regexp(args, '^-h$', 'match')));

% Help message
if flag_count < 1 | flag_count > length(Flag_list) | (total_count / flag_count) ~= 2 | help_flag
    fprintf(1, 'USAGE: mutpred2 arguments (see below)\n');
    fprintf(1, '%s-i%s<FASTA file name (String)>%sSubstitutions must be in headers, delimited by space\n', Sing_tab, Sing_tab, Padding2);
    fprintf(1, '%s-o%s[Output file name (String)]%sDefaults to the standard output\n', Sing_tab, Sing_tab, Padding2);
    fprintf(1, '%s-p%s[Use model with homology profiles (0 or 1)]%sIf 0 (default), no human and mouse proteome homolog counts are computed\n', Sing_tab, Sing_tab, Sing_tab);
    fprintf(1, '%sIf 1, these counts are computed (much slower but more accurate)\n', Padding1);
    fprintf(1, '%s-c%s[Predict conservation features (0 or 1)]%sIf 0, for substitutions from proteins where conservation scores are not available, these features are marked as zeros\n', Sing_tab, Sing_tab, Sing_tab);
    fprintf(1, '%sIf 1 (default), predicted conservation scores are used (more accurate than not using conservation features but less accurate than using known conservation scores)\n', Padding1);
    fprintf(1, '%s-b%s[Skip PSI-BLAST (0 or 1)]%sIf 0 (default and also when "-c" is 1), in cases where precomputed PSSMs are not available, PSI-BLAST is run\n', Sing_tab, Sing_tab, Padding2);
    fprintf(1, '%sIf 1, they are treated as missing features (much faster but scores not as reliable for such proteins)\n', Padding1);
    fprintf(1, '%s-f%s[Output file format]%s%sIf 1 (default), loss/gain of structural and functional properties are output in both row and ontological format\n', Sing_tab, Sing_tab, Padding2, Sing_tab);
    fprintf(1, '%sIf 2, loss/gain of structural and functional properties are output in row format only\n', Padding1);
    fprintf(1, '%sIf 3, loss/gain of structural and functional properties are output in ontological format only\n', Padding1);
    fprintf(1, '%sIf 4, only MutPred2 general pathogenicity scores are output (smaller output file)\n', Padding1);
    %fprintf(1, '%s-s%s[Scores only (0 or 1]%s%sIf 0 (default), both MutPred2 scores and loss/gain of structural and functional properties are output\n', Sing_tab, Sing_tab, Padding2, Sing_tab);
    %fprintf(1, '%sIf 1, only MutPred2 scores are output (smaller output file)\n', Padding1);
    %fprintf(1, '%s-k%s[Top "k" properties (>= 1 and <= 57)]%s%sNumber of best-ranking property losses/gains to be output\n', Sing_tab, Sing_tab, Sing_tab, Sing_tab);
    fprintf(1, '%s-t%s[P-value threshold (>= 0 and <= 1)]%s%sShow only those mechanisms with P-value < this value (default: 0.05) (for Bonferroni correction, set "t" to 0.0009)\n', Sing_tab, Sing_tab, Sing_tab, Sing_tab);
    ps = 0;
    return
end

% Loop through arguments and store in structure (this will help
% accept arguments in any order)
args_struct = struct();
for j = 1:2:length(args)-1
    flag_name = strrep(args{j}, '-', '');
    if ~any(flag_name == Flag_list)
        error('Invalid argument flag "-%s"! Type "mutpred2 -h" for a list of valid arguments.', flag_name);
    end
    args_struct.(flag_name) = args{j+1};
end

% Set global variables and validate
if isfield(args_struct, 'i')
    IN_FILE = args_struct.i;
else
    error('Input file (-i) is a mandatory argument!');
end
if isfield(args_struct, 'o')
    OID = fopen(args_struct.o, 'w');
    if OID == -1
        error('Cannot open the output file!');
    end
end
if isfield(args_struct, 'p')
    DO_HOMOLOGY = str2num(args_struct.p);
    if ~ismember(DO_HOMOLOGY, Valid_p)
        error('Argument -p must be 0 or 1!');
    end
end
if isfield(args_struct, 'c')
    PREDICT_CONS = str2num(args_struct.c);
    if ~ismember(PREDICT_CONS, Valid_c)
        error('Argument -c must be 0 or 1!');
    end
end
if isfield(args_struct, 'b')
    SKIP_PSIBLAST = str2num(args_struct.b);
    if ~ismember(SKIP_PSIBLAST, Valid_b)
        error('Argument -b must be 0 or 1!');
    end
end
%if isfield(args_struct, 's')
%    SCORES_ONLY = str2num(args_struct.s);
%    if ~ismember(SCORES_ONLY, Valid_s)
%        error('Argument -s must be 0 or 1!');
%    end
%end
%if isfield(args_struct, 'k')
%    TOPK = str2num(args_struct.k);
%    if TOPK < Valid_k(1) | TOPK > Valid_k(2)
%        error('Argument -k must be between 1 and 57');
%    end
%end
if isfield(args_struct, 't')
    PTHR = str2num(args_struct.t);
    if PTHR < Valid_t(1) | PTHR > Valid_t(2)
        error('Argument -t must be between 0 and 1');
    end
end
if isfield(args_struct, 'f')
    out_form = str2num(args_struct.f);
    if ~ismember(out_form, Valid_f)
        error('Argument -f must be 1, 2, 3 or 4!');
    end
    if out_form == 1
        CSV_FORM = 1;
	ONTO_FORM = 1;
    elseif out_form == 2
        ONTO_FORM = 0;
    elseif out_form == 3
        CSV_FORM = 0;
    elseif out_form == 4
        CSV_FORM = 0;
	ONTO_FORM = 0;
    end
end

return;
