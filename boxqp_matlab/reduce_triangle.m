function [new_B_TRI, new_l_TRI] = reduce_triangle(xfix, original_n, parent_B_TRI, parent_l_TRI)

    n_ineq = size(parent_B_TRI, 2);
    new_B_TRI = cell(1, n_ineq);
    new_l_TRI = zeros(n_ineq, 1);
    new_n = original_n - size(xfix, 1);
    
    x_id = setdiff(1:(original_n+1), (xfix(:, 1)+1)', 'stable');
    
    %disp(x_id)
    
    counter = 1;
    
    for c=1:n_ineq
        
        [id_i, id_j, v] = find(parent_B_TRI{c});
        [~, new_id_i] = ismember(id_i, x_id);
        [~, new_id_j] = ismember(id_j, x_id);
        
        if any(new_id_i == 0) || any(new_id_j == 0)
            continue;
        end
        
        new_B_TRI{counter} = sparse(new_id_i, new_id_j, v, new_n+1, new_n+1);
        new_l_TRI(counter) = parent_l_TRI(c);
        counter = counter + 1;
        
    end
    
    n_ineq = counter - 1;
    new_B_TRI = new_B_TRI(1:n_ineq);
    new_l_TRI = new_l_TRI(1:n_ineq);

end