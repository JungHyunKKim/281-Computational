function profit = profitfunc_bad(k, param)
    
    profit = param.Abad * k.^param.alpha - param.r * k;

end