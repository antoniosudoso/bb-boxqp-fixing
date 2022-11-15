function f = qp_obj(x, H, f)
    f = 0.5*x'*H*x+f'*x;
end