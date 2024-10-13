function index = find_branching_index_P(P_set, xfix, Q, x, Xshr)

    original_n = size(Q, 1);
    shr_n = size(Xshr, 1);
    
    if (shr_n < original_n) && (~isempty(xfix))
        shr_idx = setdiff(1:original_n, xfix(:, 1)', 'stable');
        X = zeros(original_n, original_n);
        X(shr_idx, shr_idx) = Xshr;
    else
        X = Xshr;
    end
    
    p_vars = length(P_set);
    values = zeros(p_vars, 1);
    for k=1:p_vars
        i = P_set(k);
        sum = 0;
        for j=1:original_n
            %sum = sum + abs(Q(i, j)*(X(i, j)-x(i)*x(j)));
            sum = sum + Q(i, j)*(x(j)*x(i) - X(i, j));
        end
        values(k) = sum;
    end
    
    pairs = sortrows([P_set, values], 2, 'descend'); % pairs (id_x, viol)
    %disp(pairs)
    value = pairs(1, 2);
    if value <= 1e-4
        index = -1;
    else
        index = pairs(1, 1);
    end

end
