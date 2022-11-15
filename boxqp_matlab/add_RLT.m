function [Bcell, l] = add_RLT(n)

    n_ineq = (nchoosek(n, 2)+n)*3;
    Bcell = cell(1, n_ineq);
    l = zeros(n_ineq, 1);
    c = 1;
    for i=2:n+1
        for j=i:n+1
            Bcell{c} = sparse([1, i, i, j], [i, 1, j, i], [0.5, 0.5, -0.5, -0.5], n+1, n+1);
            l(c) = 0;
            c = c + 1;
            Bcell{c} = sparse([i, j, 1, i, 1, j], [j, i, i, 1, j, 1], [0.5, 0.5, -0.5, -0.5, -0.5, -0.5], n+1, n+1);
            l(c) = -1;
            c = c + 1;
            Bcell{c} = sparse([i, j, 1, j], [j, i, j, 1], [-0.5, -0.5, 0.5, 0.5], n+1, n+1);
            l(c) = 0;
            c = c + 1;
        end
    end

end