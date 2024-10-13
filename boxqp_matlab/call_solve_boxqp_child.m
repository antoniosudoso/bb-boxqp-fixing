function result = call_solve_boxqp_child(Q, c, xfix, xin, best_ub, best_x_gub, parent_Bcell, parent_l, params)

% The B&B algorithm solves the following nonconvex qudratic programming problem:
%
%    min      1/2*x'*Q*x + c'*x
%    s.t.     0 <=  x <= 1
%
% This script is called by the C++ program at the child node
%
% * Q, c: data of the QP formulation above
% * xfix: indices of variables fixed to 0 or 1
% * xin: inceces of variables inside the unit box
% * best_ub, best_x_ub: global upper bound
% * parent_Bcell, parent_l: initial set of triangle inequalities
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

    %disp(xfix)
    %disp(xin)
    
    if ~params.fixing
        
        % solve BoxQP
        result = solve_boxqp_sdpnal(Q, c, best_ub, best_x_gub, xfix, xin, parent_Bcell, parent_l, params);
        
        % do not fix variables
        result.time_fix = 0;
        result.sdp_fix = 0;
        result.init_xfix = [];
        result.n_fixed = 0;
        
    else
        
        % solve reduced BoxQP
        or_n = size(Q, 1);
        [new_Q, new_c, const] = reduce_problem_data(Q, c, xfix);
        result = solve_boxqp_reduced(or_n, new_Q, new_c, const, xfix, xin, best_ub, best_x_gub, parent_Bcell, parent_l, params);
    
        if result.best_gap <= params.opt_tol
            result.time_fix = 0;
            result.sdp_fix = 0;
            result.init_xfix = [];
            result.n_fixed = 0;
        % trying to fix variables
        elseif result.best_gap < 0.01
            % fixing on the reduced problem
            result_fix = multiple_fixing(new_Q, new_c, const, result.best_x_lb_shr, result.best_ub, [], result.new_xin, ...
                result.best_Bcell_shr, result.best_l_shr, params);
            result.time_fix = result_fix.time_fixing;
            result.sdp_fix = result_fix.sdp_fixing;
            result.init_xfix = result_fix.xfix;
            % recover original indices
            if ~isempty(result.init_xfix)
                shr_idx_fix = result.init_xfix(:, 1);
                result.init_xfix(:, 1) = result.shr_idx(shr_idx_fix);
            end
            result.n_fixed = length(result_fix.xfix);
        else
            result.time_fix = 0;
            result.sdp_fix = 0;
            result.init_xfix = [];
            result.n_fixed = 0;
        end
    
    end
    
    % update N and branch
    N_set = find(diag(Q) <= 0);
    N_set = setdiff(N_set, xfix(:, 1), 'stable');
    result.N_bin = length(N_set);
    result.branching_type = 2;
    
    result.best_ub =  min(result.ub_list);
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
            P_set = setdiff(P_set, xin, 'stable'); % update P
            P_set = setdiff(P_set, xfix(:, 1), 'stable'); % also remove from xfix
            i = find_branching_index_P(P_set, xfix, Q, result.best_x_lb, result.best_X_lb);
            result.idx_i = i;
        else
            i = find_branching_index(N_set, result.best_x_lb);
            result.idx_i = i;
        end
    end
    
    %disp(result)
    

end
