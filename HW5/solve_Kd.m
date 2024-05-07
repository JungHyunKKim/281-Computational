function [Kd] = solve_Kd(r, par, grids, g)

    Kd = ( (par.alpha * par.z) / (r + par.delta) )^( 1 / (1-par.alpha) ) ...
         * sum(sum(par.e.*g.*grids.da)) ; 

end