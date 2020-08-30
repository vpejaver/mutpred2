function p = use_predictors_func(seq, pssm, flag)
    % get VL2 prediction
    [p{1}, p{2}, p{3}, p{4}] = VL2(seq);

    p{5} = VSL2B(seq);

    % get Vihinen's flexibility
    p{6} = vihinen(seq);

    % get hydrophobic moments
    p{7} = hydrophobic_moment(seq, 11, 100);
    p{8} = hydrophobic_moment(seq, 11, 160);
    p{9} = hydrophobic_moment(seq, 11, 120);

    p{10} = predictBfactors(seq, 3, 5);

    p{11} = get_aa_volumes(seq);

    p{12} = hydrophobicity(seq, 1);

    % get secondary structure predictions (ignore disorder)
    if isempty(pssm)
        tmp = ssp4(seq);
    elseif ~isempty(pssm) && flag == 0
        tmp = ssp4(seq, pssm);    
    elseif ~isempty(pssm) && flag == 1
        tmp = [ssp4(seq) ssp4(seq, pssm)];
	p{13} = tmp(:, 1)';
	p{14} = tmp(:, 2)';
	p{15} = tmp(:, 3)';
	p{16} = tmp(:, 5)'; 
	p{17} = tmp(:, 6)';
	p{18} = tmp(:, 7)';    
    end
    if flag == 0
        p{13} = tmp(:, 1)';
        p{14} = tmp(:, 2)';
        p{15} = tmp(:, 3)';
    end

    % get surface accessibility
    p{end+1} = predict_surfacc(seq, pssm');

    return
