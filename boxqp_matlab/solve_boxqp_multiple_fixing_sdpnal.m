function result = solve_boxqp_multiple_fixing_sdpnal(Q, c, fix0, fix1, xfix, xin, init_B_TRI, init_l_TRI, params)
    
    % fix0: variables we are trying to fix to 0
    % fix1: variables we are trying to fix to 1
    % xfix: variables already fixed
        
    maxNumCompThreads(params.n_threads);
    
    result = struct();
    
    n = size(Q, 1);
    
    C_full = zeros(n+1, n+1);
    C_full(2:n+1, 1) = c;
    C_full(1, 2:n+1) = c';
    C_full(2:n+1, 2:n+1) = Q; 
    C_full = 0.5 * C_full;
    C = cell(1);
    C{1} = sparse(C_full);
        
    blk = cell(1);
    blk{1,1} = 's';
    blk{1,2} = n+1;
    
    B_FIX = cell(1);
    n0 = size(fix0, 1); % -
    n1 = size(fix1, 1); % +
    e0 = ones(n0, 1);
    e1 = ones(n1, 1);
    B_FIX{1} = -sparse([e0; e1; fix0+1; fix1+1], [fix0+1; fix1+1; e0; e1], ...
        [-0.5.*e0; 0.5.*e1; -0.5.*e0; 0.5.*e1], n+1, n+1);
    
    L = cell(1);
    L{1} = 0; % we want Y >= 0

    A_one = cell(1, 1);
    A_one{1} = sparse(1, 1, 1, n+1, n+1);

    if isempty(xfix) && isempty(xin)
        At = svec(blk, A_one, 1);
        At{1} = At{1}';
        b = 1;
    else
        % add xfix
        n_fix = size(xfix, 1);
        A_fix = cell(1, n_fix);
        b_fix = zeros(n_fix, 1);
        for v=1:n_fix
            id = xfix(v, 1);
            val = xfix(v, 2);
            A_fix{v} = sparse([1, id+1], [id+1, 1], [0.5, 0.5], n+1, n+1);
            b_fix(v) = val;
        end
        % add xin
        n_xin = length(xin);
        A_xin = cell(1, n_xin);
        b_xin = zeros(n_xin, 1);
        for i=1:n_xin
            Qi = Q(xin(i), :);
            A_xin{i} = sparse([ones(1, n), 2:(n+1)], [2:(n+1), ones(1, n)], [Qi/2, Qi/2], n+1, n+1);
            b_xin(i) = -c(xin(i));
        end
        At = svec(blk, [A_one, A_fix, A_xin], 1);
        b = [1; b_fix, b_xin];
    end

    [B_RLT, l_RLT] = add_RLT(n);
    Bt = svec(blk, [B_RLT, init_B_TRI, B_FIX], 1);
    l = [l_RLT; init_l_TRI; -(n1-1+0.01)];
        
    options.printlevel = params.sdp_verbose;
    options.tol = params.sdp_tol;
    options.stopoption = 0;
    options.maxiter = 100000;
    options.maxtime = 10800;

    % solve SDP
    [~, ~, ~, ~, Z1, ~, ~, ~, info, ~] = sdpnalplus(blk, At, C, b, L, [], Bt, l, [], options);
    [~, lb, timeLP] = safe_bound_lp(blk, Z1, At, Bt, C, b, l, params.n_threads);
    
    result.lb = lb;
    result.timeSDP = info.totaltime;
    result.timeLP = timeLP;
    
    
end