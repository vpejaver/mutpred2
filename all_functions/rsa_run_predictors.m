function p = rsa_run_predictors(seq)
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

return
