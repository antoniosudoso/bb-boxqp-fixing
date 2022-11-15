function [or_B_TRI, or_l_TRI] = recover_triangle(xfix, n, shr_B_TRI, shr_l_TRI)

    n_ineq = size(shr_B_TRI, 2);
    or_B_TRI = cell(1, n_ineq);
    or_l_TRI = zeros(n_ineq, 1);
    
    x_id = setdiff(1:(n+1), (xfix(:, 1)+1)', 'stable');
    
    %disp(x_id)
    
    counter = 1;
    
    for c=1:n_ineq
        
        [id_i, id_j, v] = find(shr_B_TRI{c});        
        or_B_TRI{counter} = sparse(x_id(id_i), x_id(id_j), v, n+1, n+1);
        or_l_TRI(counter) = shr_l_TRI(c);
        counter = counter + 1;
        
    end

end