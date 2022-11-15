function result = solve_boxqp_sdpnal(Q, c, gub, x_gub, xfix, xin, init_B_TRI, init_l_TRI, params)

% This script solves the SDP relaxation of the following nonconvex qudratic programming problem:
%
%    min      1/2*x'*Q*x + c'*x
%    s.t.     0 <=  x <= 1
%
% * Q, c: data of the QP formulation above
% * gub, x_gub: global upper bound
% * xfix: indices of variables fixed to 0 or 1
% * xin: inceces of variables inside the unit box
% * init_B_TRI, init_l_TRI: initial set of triangle inequalities
% * params: a structure with the following possible fields
%   - n_threads: number of threads
%   - opt_tol: optimality tolerance of the B&B algorithm (relative gap)
%   - sdp_tol: accuracy of SDPNAL+ (relative KKT residual)
%   - sdp_verbose: 0 do not display SDPNAL+ log, 1 display SDPNAL+ log
%   - cp_maxiter: maximum number of cutting-plane iterations
%   - cp_tol: relative tolerance of the cutting-plane algorithm
%   - cp_maxineq: maximum number of triangle inequalities to separate
%   - cp_epsineq: tolerance for checking the violation of triangle inequalities
%   - cp_percineq: fraction of violated triangle inequalities to add
%   - cp_activeineq: tolerance for finding active triangle inequalities
    
    maxNumCompThreads(params.n_threads);
    
    result = struct();
    result.message = "Maximum number of iterations";
    
    result.ineq_list = [];
    result.timeSDP_list = [];
    result.timeLP_list = [];
    result.lb_list = [];
    result.ub_list = [];
    result.gap_list = [];
    
    result.best_lb = inf;
    result.best_ub = gub;
    result.best_gap = inf;
    result.best_x_ub = x_gub;
    result.best_x_lb = [];
    result.best_X_lb = [];
    result.cp_flag = -1;

    n = size(Q, 1);
    result.n = n;
    
    C_full = zeros(n+1, n+1);
    C_full(2:n+1, 1) = c;
    C_full(1, 2:n+1) = c';
    C_full(2:n+1, 2:n+1) = Q; 
    C_full = 0.5 * C_full;
        
    blk = cell(1);
    blk{1,1} = 's';
    blk{1,2} = n+1;
    
    C = cell(1);
    C{1} = sparse(C_full);

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
        b = [1; b_fix; b_xin];
    end
    
    n_rlt = (nchoosek(n, 2)+n)*3;
    [B_RLT, l_RLT] = add_RLT(n);
    Bt = svec(blk, [B_RLT, init_B_TRI], 1);
    l = [l_RLT; init_l_TRI];
        
    options.printlevel = params.sdp_verbose;
    options.tol = params.sdp_tol;
    options.stopoption = 0;
    options.maxiter = 100000;
    options.maxtime = 10800;

    % solve SDP
    [~, Yopt, ~, ~, Z1, ~, ~, ~, info, ~] = sdpnalplus(blk, At, C, b, L, [], Bt, l, [], options);
    [~, lb, timeLP] = safe_bound_lp(blk, Z1, At, Bt, C, b, l, params.n_threads);
    
    result.lb_list = [result.lb_list; lb];
    result.best_Bcell = init_B_TRI;
    result.best_l = init_l_TRI;
    
    % local search
    options_local.screen = 'off';
    Y = Yopt{1};
    x_SDP = Y(2:n+1, 1);
    [x_ub, ub] = snsolve(@(x) qp_obj(x, Q, c), x_SDP, [], [], [], [], zeros(n, 1), ones(n, 1), options_local);
    if ub < result.best_ub
        result.best_ub = ub;
        result.best_x_ub = x_ub;
    end
    result.best_x_lb = x_SDP;
    result.best_X_lb = Y(2:n+1, 2:n+1);
    result.ub_list = [result.ub_list; ub];
    result.ineq_list = [result.ineq_list; size(init_B_TRI, 2)+n_rlt];
    result.timeSDP_list = [result.timeSDP_list; info.totaltime];
    result.timeLP_list = [result.timeLP_list; timeLP];
    result.best_lb = lb;
    result.best_gap = (result.best_ub - result.best_lb)/abs(result.best_ub);
    result.gap_list = [result.gap_list; result.best_gap];
    result.cp_iter = 0;
    if result.best_gap <= params.opt_tol
        result.message = "Pruning";
        result.cp_flag = 0;
        return;
    end
    
    fprintf('\nRelative gap = %d \n', result.best_gap);
    
    B_TRI = init_B_TRI;
    l_TRI = init_l_TRI;
    
    for i=1:params.cp_maxiter
        
        if ~isempty(B_TRI) && ~isempty(l_TRI)
            Bt = svec(blk, B_TRI, 1);
            Yvec = svec(blk, Yopt{1}, 1);
            active_id = abs(Bt{1}' * Yvec - l_TRI) <= params.cp_activeineq;
            B_TRI = B_TRI(active_id);
            l_TRI = l_TRI(active_id);
        end
        
        [B_triangle, l_triangle, ~] = separate_triangle(Yopt{1}, params.cp_epsineq, params.cp_maxineq, params.cp_percineq);
        n_triangle = size(B_triangle, 2);
        if n_triangle < n*10
            result.message = "Less than n*10 triangle";
            result.cp_flag = 1;
            break;
        end
        
        B_TRI = [B_TRI, B_triangle];
        l_TRI = [l_TRI; l_triangle];
        Bt = svec(blk, [B_RLT, B_TRI], 1);
        l = [l_RLT; l_TRI];
        
        % solve SDP
        [~, Yopt, ~, ~, Z1, ~, ~, ~, info, ~] = sdpnalplus(blk, At, C, b, L, [], Bt, l, [], options);
        [~, lb, timeLP] = safe_bound_lp(blk, Z1, At, Bt, C, b, l, params.n_threads);       
        
        result.lb_list = [result.lb_list; lb];
        result.cp_iter = result.cp_iter + 1;
        
        % local search
        Y = Yopt{1};
        x_SDP = Y(2:n+1, 1);
        [x_ub, ub] = snsolve(@(x) qp_obj(x, Q, c), x_SDP, [], [], [], [], zeros(n, 1), ones(n, 1), options_local);
        if ub < result.best_ub
            result.best_ub = ub;
            result.best_x_ub = x_ub;
        end
        
        result.ub_list = [result.ub_list; ub];
        result.ineq_list = [result.ineq_list; size(B_TRI, 2)+n_rlt];
        result.timeSDP_list = [result.timeSDP_list; info.totaltime];
        result.timeLP_list = [result.timeLP_list; timeLP];
        
        if abs(lb - result.best_lb)/abs(lb) <= params.opt_tol
            result.best_lb = lb;
            result.best_x_lb = x_SDP;
            result.best_X_lb = Y(2:n+1, 2:n+1);
            result.best_gap = (result.best_ub-result.best_lb)/abs(result.best_ub);
            result.gap_list = [result.gap_list; result.best_gap];
            result.best_Bcell = B_TRI;
            result.best_l = l_TRI;
            result.message = "Bound not significantly improved";
            result.cp_flag = 2;
            if result.best_gap <= params.opt_tol
                result.message = "Pruning";
                result.cp_flag = 0;
            end
            break;
        end
        
        if (lb - result.best_lb)/abs(lb) <= params.cp_tol
            result.message = "Bound not improved";
            result.cp_flag = 3;
            if result.best_gap <= params.opt_tol
                result.message = "Pruning";
                result.cp_flag = 0;
            end
            break;
        end
        
        result.best_lb = lb;
        result.best_x_lb = x_SDP;
        result.best_X_lb = Y(2:n+1, 2:n+1);
        if ub < result.best_ub
            result.best_ub = ub;
            result.best_x_ub = x_ub;
        end
        
        result.best_Bcell = B_TRI;
        result.best_l = l_TRI;
        
        result.best_gap = (result.best_ub-result.best_lb)/abs(result.best_ub);
        result.gap_list = [result.gap_list; result.best_gap];
        if result.best_gap <= params.opt_tol
            result.message = "Pruning";
            result.cp_flag = 0;
            break;
        end
        
        fprintf('\nRelative gap = %d \n', result.best_gap);
        
    end
    

end
