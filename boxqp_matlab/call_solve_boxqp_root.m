function result = call_solve_boxqp_root(Q, c, params)

% The B&B algorithm solves the nonconvex qudratic programming problem:
%
%    min      1/2*x'*Q*x + c'*x
%    s.t.     0 <=  x <= 1
%
% This script is called by the C++ program at the root node
%
% * Q, c: data of the QP formulation above
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

    disp(params)
    
    n = size(Q, 1);
    % local search with multistart
    options.screen = 'off';
    best_ub = inf;
    best_x_gub = [];
    for i=1:200
        [x, fval] = snsolve(@(x) qp_obj(x, Q, c), rand(n, 1), [], [], [], [], zeros(n, 1), ones(n, 1), options);
        if fval < best_ub
            best_ub = fval;
            best_x_gub = x;
            disp(best_ub)
        end
    end
    
    result = solve_boxqp_sdpnal(Q, c, best_ub, best_x_gub, [], [], cell(0), [], params);
    
    if ~params.fixing
        
        % do not fix variables
        result.time_fix = 0;
        result.sdp_fix = 0;
        result.init_xfix = [];
        result.n_fixed = 0;
        
    else
        
        % trying to fix variables
        if result.best_gap <= params.opt_tol
            result.time_fix = 0;
            result.sdp_fix = 0;
            result.init_xfix = [];
            result.n_fixed = 0;
        elseif result.best_gap < 0.01
            % fixing
            result_fix = multiple_fixing(Q, c, 0, result.best_x_lb, result.best_ub, [], [], result.best_Bcell, result.best_l, params);
            result.time_fix = result_fix.time_fixing;
            result.sdp_fix = result_fix.sdp_fixing;
            result.init_xfix = result_fix.xfix;
            result.n_fixed = length(result_fix.xfix);
        else
            result.time_fix = 0;
            result.sdp_fix = 0;
            result.init_xfix = [];
            result.n_fixed = 0;   
        end
    
    end
    
    % branch on variables in N set
    N_set = find(diag(Q) <= 0);
    result.N_bin = length(N_set);
    result.branching_type = 2;

    if result.cp_flag == 0 || result.best_gap <= params.opt_tol
        result.idx_i = -1;
    else
        if ~isempty(result.init_xfix)
            N_set = setdiff(N_set, result.init_xfix(:, 1), 'stable');
        end
        if isempty(N_set)
            % branch on variables in P set
            result.branching_type = 3;
            P_set = find(diag(Q) > 0);
            i = find_branching_index_P(P_set, [], Q, result.best_x_lb, result.best_X_lb);
            result.idx_i = i;
        else
            i = find_branching_index(N_set, result.best_x_lb);
            result.idx_i = i;
        end
    end
    
    %disp(result)
    

end
