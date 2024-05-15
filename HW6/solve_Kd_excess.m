function [Kd_excess,A,u,c,v_new,g,Va_Upwind,Vaf,Vab,If,Ib,I0,aprimedot,Ks,Ls,Kd,w,Aswitch] = ...
         solve_Kd_excess(r, par, num, grids, eqtype)
    
    switch eqtype
        case 'GE'
            w  = (1 - par.alpha) * exp(par.logTFP)^( 1 / (1-par.alpha) ) ...
                 *( par.alpha / (r + par.delta) )^( par.alpha/(1-par.alpha) ); 
        case 'PE'
            w  = par.w0; 
        otherwise 
            error(['Invalid input specified.' ...
                  'Use "GE" (General Equilibrium) or "PE" (Partial Equilibrium).'])
    end
    
    [A,u,c,v_new,g,Va_Upwind,Vaf,Vab,If,Ib,I0,aprimedot,Ks,Ls,Aswitch] ...
              = solve_HH(r, w, par, num, grids); 

    [Kd]      = solve_Kd(r, par, Ls); 
    Kd_excess = Kd - Ks; 

end