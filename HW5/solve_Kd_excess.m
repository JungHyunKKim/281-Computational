function [Kd_excess, Kd, Ks, g, w] = solve_Kd_excess(r, par, num, grids, eqtype)
    
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

    [~, ~, Ks, Ls, g] = solve_HH(r, w, par, num, grids); 
    [Kd]              = solve_Kd(r, par, Ls); 

    Kd_excess = Kd - Ks; 

end