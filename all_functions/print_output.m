% Vikas Pejaver
% October 2015

% Function to print the output in the desired format
% Updated February 2017 to change what the output looks like

function print_output(id, sequence, mutations, annotations, pos, fx, prs, pt, pp, minfo)


%%%%%% Global variables %%%%%%
global OID NAMES_REORDERED CSV_FORM ONTO_FORM PTHR; % TOPK


%%%%%% Constants and defaults %%%%%%
To_skip = '(Non_transmembrane)';
Altered = {'N-terminal_signal', 'Signal_helix', 'C-terminal_signal', 'Signal_cleavage', 'Cytoplasmic_loop', 'Transmembrane_region', 'Non-cytoplasmic_loop', 'Non_transmembrane', 'Coiled_coil', 'Calmodulin_binding', 'DNA_binding', 'RNA_binding', 'PPI_residue', 'PPI_hotspot', 'MoRF', 'Stability'};
Region = [Altered, 'Helix', 'Strand', 'Loop', 'Intrinsic_disorder', 'B-factor', 'Relative_solvent_accessibility'];
Types = {'Loss', 'Gain'};
%Act_gthr = 0.68;
%Conf_gthr = 0.79;
%Vconf_gthr = 0.79;
%Act_pthr = 0.05;
%Conf_pthr = 0.05 / 57;
%Vconf_pthr = 0.05 / 57;


%%%%%% Update property list %%%%%%
NAMES_REORDERED = regexprep(NAMES_REORDERED, 'Intracellular', 'Cytoplasmic');
NAMES_REORDERED = regexprep(NAMES_REORDERED, 'Extracellular', 'Non-cytoplasmic');
NAMES_REORDERED = regexprep(NAMES_REORDERED, 'VSL2B_disorder', 'Intrinsic_disorder');
NAMES_REORDERED = regexprep(NAMES_REORDERED, 'Surface_accessibility', 'Relative_solvent_accessibility');
NAMES_REORDERED = regexprep(NAMES_REORDERED, 'GPI_anchor_amidation', 'GPI-anchor_amidation');
NAMES_REORDERED = regexprep(NAMES_REORDERED, 'B_factor', 'B-factor');


%%%%%% Loop through each substitution %%%%%%
for i = 1:length(fx)

    % Initialize
    out_str = '';
    action = '';
    conf = '';
    vconf = '';
    remark = '-';
    %new_pvals = zeros(1, length(lp(i, :)));

    %% Print basic info
    %idmut = sprintf('%s,%s', id, mutations{i});
    %fprintf(OID, '%s,', idmut);

    % Get special remarks
    if any(annotations(i, :))
        remark = get_remarks(annotations(i, :));
    end

    % Get residues at position
    if CSV_FORM == 1 | ONTO_FORM == 1
        %if ~any(pos(i, :))
	%    %fprintf(OID, '-,%s,%s,%s,%s,%s\n', action, conf, vconf, out_str, remark);
	%    fprintf(OID, '-,%s,%s\n', out_str, remark);
	%    continue;
	%end
	res = sequence(pos(i, :));
	pos_info = cellstr(strcat(res', num2str(pos(i, :)')));
	pos_info = [pos_info; mutations{i}(1:end-1)]; % add dummy position for stability

	% Hard-set motif scores and P-values to be zero
	prs(i, end-1) = 0;
	pp(i, end-1) = 1;

	% Run through ontology graph and map
	[new_scores, new_pvals, paths] = propagate_posteriors(prs(i, :), pp(i, :));    

	% Sort everything based on posteriors
	%[sorted_pvals, S] = sort(comb_pvals);
	[sorted_scores, S] = sort(prs(i, :), 'descend');
	sorted_pos = pos_info(S);
	sorted_types = Types(pt(i, S));
	sorted_names = NAMES_REORDERED(S);
	%sorted_pvals = pp(i, S);
	sorted_gpvals = new_pvals(S);
	%sorted_gscores = new_scores(S);
	%sorted_paths = paths(S);
	
	% benjamini-hochberg correction
	%bh_flags5 = zeros(length(sorted_pvals), 1);
	%bh_flags1 = zeros(length(sorted_pvals), 1);
	%bh_flags5(benjamini_hochberg(sorted_pvals', Alpha5)) = 1;
	%bh_flags1(benjamini_hochberg(sorted_pvals', Alpha1)) = 1;

	% Print rest of the results
	%k = 1;
	for j = 1:length(sorted_names)
	    if isempty(regexp(sorted_names{j}, To_skip, 'match'))
	        this_mech = '';
	        %path_mech = '';
		%for m = 1:length(sorted_paths{j})-1
		%	if sorted_gpvals{j}(m) < 0.01
		%	    formatted_pval = sprintf('P = %.1e', sorted_gpvals{j}(m));
		%	else
		%	    formatted_pval = sprintf('P = %.2f', sorted_gpvals{j}(m));
		%	end
		%    path_mech = sprintf('%s -> %s (Pr = %.2f | %s)', path_mech, sorted_paths{j}{m}, sorted_gscores{j}(m), formatted_pval);
		%end
		%path_mech = regexprep(path_mech, '^( -> )', '');

		if (sorted_scores(j) < 0.01 | sorted_gpvals{j}(end) >= PTHR) & ~strcmp(sorted_names{j}, 'Motifs')
		    continue;
		end

		if sorted_gpvals{j}(end) < 0.01
		    formatted_pval = sprintf('P = %.1e', sorted_gpvals{j}(end));
		else
		    formatted_pval = sprintf('P = %.2f', sorted_gpvals{j}(end));
		end
		if strcmp(sorted_names{j}, 'Motifs')
		    motif_mech = sprintf('Altered Motifs (%s)', minfo{i});
		    %this_mech = sprintf('Altered Motifs (%s)', minfo{i});
	        %elseif strcmp(sorted_names{j}, 'Stability')
		%    this_mech = sprintf('%s of Stability (Pr = %.2f | %s)', sorted_types{j}, sorted_scores(j), formatted_pval);
		elseif ismember(sorted_names{j}, Altered) & ismember(sorted_names{j}, Region)
		    this_mech = sprintf('Altered %s (Pr = %.2f | %s)', sorted_names{j}, sorted_scores(j), formatted_pval);
		elseif ismember(sorted_names{j}, Altered)
		    this_mech = sprintf('Altered %s at %s (Pr = %.2f | %s)', sorted_names{j}, strrep(sorted_pos{j}, ' ', ''), sorted_scores(j), formatted_pval);
		elseif ismember(sorted_names{j}, Region)
		    this_mech = sprintf('%s of %s (Pr = %.2f | %s)', sorted_types{j}, sorted_names{j}, sorted_scores(j), formatted_pval);
		else
		    this_mech = sprintf('%s of %s at %s (Pr = %.2f | %s)', sorted_types{j}, sorted_names{j}, strrep(sorted_pos{j}, ' ', ''), sorted_scores(j), formatted_pval);
		end
		%path_mech = sprintf('%s -> %s', path_mech, this_mech);
		
		%if sorted_W(j) >= 0.5 | sorted_M(j) >= 0.5
	            %if fx(i) >= 0.5
		         %if fx(i) >= Vconf_gthr & sorted_pvals(j) < Vconf_pthr %bh_flags1(j) == 1 
			 %    vconf = sprintf('%s; %s', vconf, path_mech);
			 %elseif fx(i) >= Conf_gthr & sorted_pvals(j) < Conf_pthr  %bh_flags5(j) == 1
			 %    conf = sprintf('%s; %s', conf, path_mech);
			 %elseif fx(i) >= Act_gthr & sorted_pvals(j) < Act_pthr %bh_flags5(j) == 1
			 %    action = sprintf('%s; %s', action, path_mech);
			 %end
			 %end
		
		    %if k <= TOPK
		    out_str = sprintf('%s; %s', out_str, this_mech);
		    %out_str = sprintf('%s; %s', out_str, path_mech);
		%    k = k + 1;
	 	     %end
		%end
	    end
	end
    end

    %action = regexprep(action, '^(; )', '');
    %conf = regexprep(conf, '^(; )', '');
    %vconf  = regexprep(vconf, '^(; )', '');
    out_str = regexprep(out_str, '(^(; ))|((; )$)', '');

    %fprintf(OID, '%.3f,%s,%s\n', fx(i), out_str, remark);
    tokens = regexp(id{i}, ' ', 'split');
    for j = 1:length(tokens)
        if CSV_FORM
	    fprintf(OID, '%s,%s,%.3f,%s,%s,%s\n', tokens{j}, mutations{i}, fx(i), out_str, motif_mech, remark);
	else
	    fprintf(OID, '%s,%s,%.3f,%s\n', tokens{j}, mutations{i}, fx(i), remark);
	end
	if ONTO_FORM
	    pfix = sprintf('%s,%s', tokens{j}, mutations{i});
	    formatted_onto = overlay_ontology(new_scores, new_pvals, paths, pt(i, :), pos_info, minfo{i}, pfix);
	    fprintf(OID, '%s', formatted_onto);
	end
    end
end

return
