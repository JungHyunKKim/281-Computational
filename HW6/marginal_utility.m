function du = marginal_utility(c, par)
    if par.gamma == 1
        du = 1 ./ c;
    else
        du = c.^(-par.gamma);
    end
end