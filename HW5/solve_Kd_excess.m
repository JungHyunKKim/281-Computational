function [Kd_excess] = solve_Kd_excess(r, par, num, grids, eqtype)
    
    switch eqtype
        case 'GE'
            w  = (1 - par.alpha) * par.z^( 1 / (1-par.alpha) ) ...
                 *( par.alpha / (r + par.delta) )^( par.alpha/(1-par.alpha) ); 
        case 'PE'
            w  = par.w0; 
        otherwise 
            error(['Invalid input specified.' ...
                  'Use "GE" (General Equilibrium) or "PE" (Partial Equilibrium).'])
    end

    [Ks, g] = solve_Ad(r, w, par, num, grids); 
    [Kd]    = solve_Kd(r, par, grids, g); 

    Kd_excess = Kd - Ks; 

end