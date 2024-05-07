% -------------------------------------------------------------------------
% Solve Aiyagari Model
% Minimize the square of excess capital demand using fminsearch
% Functions: solve_Kd_excess
%            solve_Kd
%            solve_Ad
% Note: if r_vec - rho > 0, accumulate infinite amount of asset
%                           no stationary distribution
%                           ratchet effect (beta*(1+r)>=1 in discrete case)
% -------------------------------------------------------------------------

function [rstar, Kstar, g] = solve_Aiyagari(par, num, grids, eqtype, r0)

tic

% Define square of excess capital demand
Kd_excess_eq = @(r) solve_Kd_excess(r, par, num, grids, eqtype)^2; 

% Solve interest rate that minimizes the squared of excess demand
[rstar, excess_star] = fminsearch(Kd_excess_eq, r0); 

% Get the optimal capital
switch eqtype
    case 'GE'
        w  = (1 - par.alpha) * par.z^( 1 / (1-par.alpha) ) ...
             *( par.alpha / (rstar + par.delta) )^( par.alpha/(1-par.alpha) ); 
    case 'PE'
        w  = par.w0; 
end

[Kstar, g] = solve_Ad(rstar, w, par, num, grids); 

toc
fprintf('Error (excess demand squared) is %.4f \n', excess_star)

end
