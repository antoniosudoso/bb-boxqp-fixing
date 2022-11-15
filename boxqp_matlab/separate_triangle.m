function [Bcell, l, viol_vect] = separate_triangle(Y, eps, max_ineq, perc_ineq)

    n = size(Y, 1)-1;
    
    Bcell = cell(1, max_ineq);
    l = zeros(max_ineq, 1);
    viol_vect = zeros(max_ineq, 1);
    c = 1;
    
    stop = false;
    for i=2:n+1
        for j=i+1:n+1
            for t=j+1:n+1
                 if (i ~= t)
                     % -Y1i-Y1j-Y1t+Yij+Yit+Yjt >= -1
                     viol = -Y(1,i)-Y(1,j)-Y(1,t)+Y(i,j)+Y(i,t)+Y(j,t) + 1;
                     if viol <= -eps
                         Bcell{c} = sparse([1, i, 1, j, 1, t, i, j, i, t, j, t], ...
                             [i, 1, j, 1, t, 1, j, i, t, i, t, j], ...
                             [-0.5, -0.5, -0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5], n+1, n+1);
                         l(c) = -1;
                         viol_vect(c) = viol;
                         c = c + 1;
                     end
                     % -Yij-Yit+Y1i+Yjt >= 0
                     viol = -Y(i,j)-Y(i,t)+Y(1,i)+Y(j,t);
                     if viol <= -eps
                         Bcell{c} = sparse([i, j, i, t, 1, i, j, t], ...
                             [j, i, t, i, i, 1, t, j], ...
                             [-0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5], n+1, n+1);
                         l(c) = 0;
                         viol_vect(c) = viol;
                         c = c + 1;
                     end
                     % -Yij-Yjt+Y1j+Yit >= 0
                     viol = -Y(i,j)-Y(j,t)+Y(1,j)+Y(i,t);
                     if viol <= -eps
                         Bcell{c} = sparse([i, j, j, t, 1, j, i, t], ...
                             [j, i, t, j, j, 1, t, i], ...
                             [-0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5], n+1, n+1);
                         l(c) = 0;
                         viol_vect(c) = viol;
                         c = c + 1;
                     end
                     % -Yit-Yjt+Y1t+Yij >= 0
                     viol = -Y(i,t)-Y(j,t)+Y(1,t)+Y(i,j);
                     if viol <= -eps
                         Bcell{c} = sparse([i, t, j, t, 1, t, i, j], ...
                             [t, i, t, j, t, 1, j, i], ...
                             [-0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5], n+1, n+1);
                         l(c) = 0;
                         viol_vect(c) = viol;
                         c = c + 1;
                     end
                     
                     if c >= max_ineq
                         stop=true;
                         break;
                     end
                     
                 end
            end
            
            if stop
                break;
            end
            
        end
        
        if stop
            break;
        end
        
    end
    
    n_ineq = c-1;
    selected_ineq = max_ineq * perc_ineq;
    if n_ineq <= selected_ineq
        Bcell = Bcell(1:n_ineq);
        l = l(1:n_ineq);
    else    
        viol_vect = abs(viol_vect(1:n_ineq));
        [~, id_sorted] = sort(viol_vect, 'descend');
        Bcell = Bcell(id_sorted);
        l = l(id_sorted);
        n_ineq = floor(selected_ineq);
        Bcell = Bcell(1:n_ineq);
        l = l(1:n_ineq);
    end
    
end