function u = utility(c, par)

    if par.gamma == 1
        u = log(c);
    else
        u = c.^(1 - par.gamma) / (1 - par.gamma);
    end

end