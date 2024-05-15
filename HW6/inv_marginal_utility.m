function c = inv_marginal_utility(v, par)
    if par.gamma == 1
        c = 1 ./ v;
    else
        c = v.^(-1/par.gamma);
    end
end