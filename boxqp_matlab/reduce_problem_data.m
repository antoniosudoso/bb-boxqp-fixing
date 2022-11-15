function [newQ, newc, const] = reduce_problem_data(Q, c, xfix)

    if isempty(xfix)
        newQ = Q;
        newc = c;
        const = 0;
        return
    end
    
    n = size(Q, 1);
    fix0 = xfix(xfix(:, 2) == 0);
    fix1 = xfix(xfix(:, 2) == 1);
    nfix1 = size(fix1, 1);
    
    xblock = setdiff(1:n, [fix0; fix1]', 'stable');
    tempQ = Q([xblock, fix1'], [xblock, fix1']);
    tempc = c([xblock, fix1']);
    tempn = size(tempQ, 1);
    k = tempn-nfix1;
    newQ = tempQ(1:k, 1:k);
    Q22 = tempQ(k+1:tempn, k+1:tempn);
    Q12 = tempQ(1:k, k+1:tempn);
    c1 = tempc(1:k);
    c2 = tempc(k+1:tempn);
    
    newc = Q12*ones(nfix1, 1)+c1;
    const = 0.5*sum(Q22(:))+sum(c2);
    
end