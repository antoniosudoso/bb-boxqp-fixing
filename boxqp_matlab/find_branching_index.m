function index = find_branching_index(N_set, x)

    bin_var = x(N_set);
    pairs = sortrows([N_set, min(bin_var, 1-bin_var)], 2, 'descend'); % pairs (id_x, viol_bin)
    %disp(pairs)
    value = pairs(1, 2);
    if value <= 1e-4
        index = -1;
    else
        index = pairs(1, 1);
    end

end