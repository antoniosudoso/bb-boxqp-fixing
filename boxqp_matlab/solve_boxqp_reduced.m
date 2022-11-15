function result = solve_boxqp_reduced(original_n, new_Q, new_c, const, xfix, xin, best_ub, best_x_gub, parent_Bcell, parent_l, params)

    shr_idx = setdiff(1:original_n, xfix(:, 1)', 'stable');
        
    % update triangle inequalities for the reduced problem
    [new_B_TRI, new_l_TRI] = reduce_triangle(xfix, original_n, parent_Bcell, parent_l);
    % update xin
    [~, new_xin] = ismember(xin, shr_idx);
    
    result = solve_boxqp_sdpnal(new_Q, new_c, best_ub-const, best_x_gub(shr_idx), [], new_xin, new_B_TRI, new_l_TRI, params);

    result.new_xin = new_xin;
    result.best_lb = result.best_lb + const;
    result.best_ub = result.best_ub + const;    
    result.lb_list = result.lb_list + const;
    result.ub_list = result.ub_list + const;
    
    result.best_Bcell_shr = result.best_Bcell;
    result.best_l_shr = result.best_l;
    % save the reduced solution
    result.best_x_lb_shr = result.best_x_lb;
    
    % recover triangle inequalities with original indices
    [or_B_TRI, or_l_TRI] = recover_triangle(xfix, original_n, result.best_Bcell_shr, result.best_l_shr);
    result.best_Bcell = or_B_TRI;
    result.best_l = or_l_TRI;

    % recover original solution ub
    original_x_ub = zeros(original_n, 1);
    original_x_ub(xfix(:, 1)) = xfix(:, 2);
    original_x_ub(shr_idx) = result.best_x_ub;
    result.best_x_ub = original_x_ub;
    
    % recover original solution lb
    original_x_lb = zeros(original_n, 1);
    original_x_lb(xfix(:, 1)) = xfix(:, 2);
    original_x_lb(shr_idx) = result.best_x_lb;
    result.best_x_lb = original_x_lb;
    
    % keep track of the mapping
    result.shr_idx = shr_idx;

end
