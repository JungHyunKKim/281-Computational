function profit = profitfunc_good(k, param)
    
    profit = param.Agood * max(k-param.kappa, 0).^param.alpha - param.r * k;

end