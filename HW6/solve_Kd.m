function [Kd] = solve_Kd(r, par, Ld)

    Kd = ( (par.alpha * par.z) / (r + par.delta) )^( 1 / (1-par.alpha) ) * Ld ; 

end