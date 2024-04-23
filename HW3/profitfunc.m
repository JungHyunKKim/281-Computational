function profit = profitfunc(k, param)
    
    profit = max(param.Agood * max(k-param.kappa, 0).^param.alpha, param.Abad * k.^param.alpha) - param.r * k;

end