function result = multiple_fixing(Q, c, const, x_lb, gub, xfix, xin, B_TRI, l_TRI, params)

% This script implements the multiple fixing strategy
% * Q, c, const: data of the QP formulation
% * x_lb: optimal solution of the SDP relaxation
% * xfix: indices of variables fixed to 0 or 1
% * xin: inceces of variables inside the unit box
% * gub: global upper bound
% * B_TRI, l_TRI: set of triangle inequalities
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
%   - fixing_tol: tolerance for fixing variables

    N_set = find(diag(Q) <= 0);
    if ~isempty(xfix)
        N_set = setdiff(N_set, xfix(:, 1), 'stable');
    end
    bin_var = x_lb(N_set);    
    viol = sortrows([N_set, min(bin_var, 1-bin_var)], 2, 'ascend'); % pairs (id_x, viol_bin)
    %disp(viol)
            
    result.time_fixing = 0;
    result.sdp_fixing = 0;
    result.xfix = [];
    
    tol = params.fixing_tol; % initial fixing tolerance
    while true
        k = length(find(viol(:, 2) <= tol));
        if k == 0
            break;
        end
        fprintf('Trying to fix %d variables...\n', k);
        fix0 = [];
        fix1 = [];
        for i=1:k
            % build fix0 e fix1
            id_i = viol(i, 1);
            if x_lb(id_i) <= 0.5
                fix0 = [fix0; id_i];
            else
                fix1 = [fix1; id_i];
            end   
        end
        params.sdp_verbose = 0;
        out = solve_boxqp_multiple_fixing_sdpnal(Q, c, fix0, fix1, xfix, xin, B_TRI, l_TRI, params);
        result.sdp_fixing = result.sdp_fixing + 1;
        result.time_fixing = result.time_fixing + out.timeSDP + out.timeLP;
        v = out.lb + const;
        fprintf('\nBound Fixing: %4.2f\nGub: %4.2f\n', v, gub); 
        if v >= gub
            n0 = size(fix0, 1);
            n1 = size(fix1, 1);
            fprintf('\nFixed %d variables with tol <= %4.2f\n', n0+n1, tol);
            result.xfix = [result.xfix; fix0, zeros(n0, 1); fix1, ones(n1, 1)];
            break;
        else
            if tol == 0.01
                break;
            end
            tol = tol - 0.01;
            fprintf('\nReducing tol to %4.2f\n', tol);
        end
        
    end

end