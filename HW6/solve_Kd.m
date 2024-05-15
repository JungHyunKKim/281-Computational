function [Kd] = solve_Kd(r, par, Ld)

    Kd = ( (par.alpha * exp(par.logTFP)) / (r + par.delta) )^( 1 / (1-par.alpha) ) * Ld ; 

end