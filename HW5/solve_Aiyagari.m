% -------------------------------------------------------------------------
% Solve Aiyagari Model
% Minimize the square of excess capital demand using fminsearch
% Functions: solve_Kd_excess
%            solve_Kd
%            solve_Ad
% Note: 
%  - if r_vec - rho > 0, (1) accumulate infinite amount of asset
%                        (2) no stationary distribution
%                        (3) ratchet effect (beta*(1+r)>=1 in discrete case)
%
%  - Can be solved using root-finding methods (ex. Newton Raphson, bisection)
% -------------------------------------------------------------------------

function [rstar, wstar, Kstar, gstar] = solve_Aiyagari(par, num, grids, eqtype, r0)

% -------------------------------------------------------------------------
% Minimization
% -------------------------------------------------------------------------

tic

% Variables to hold final results
Ks_opt = [];
g_opt  = [];
w_opt  = []; 

% Redefine the objective function to return squared excess demand and capture g and Ks
function [error, Ks, g, w] = objective(r)

    [Kd_excess, ~, Ks, g, w]  = solve_Kd_excess(r, par, num, grids, eqtype);
    error                     = Kd_excess^2;
    
end

% Set the value of the Outputfcn field of the options structure to a function handle to outfun.
opts = optimset('Display','iter','TolFun',1e-08,'OutputFcn', @outfun);

% Solve for the interest rate that minimizes the squared excess demand
[rstar, error, ~, ~] = fminsearch(@objective, r0, opts);
Kstar               = Ks_opt;
gstar               = g_opt;
wstar               = w_opt;

toc

fprintf('Error (excess demand squared) is %.4f \n', error)

% -------------------------------------------------------------------------
% Define output function to save information in last iteration
% https://www.mathworks.com/help/matlab/math/output-functions.html
% -------------------------------------------------------------------------

function stop = outfun(r, ~, state)
        stop = false;
        if isequal(state, 'iter')
            % Update the persistent storage only if necessary
            % This block can be used to log or monitor progress
        elseif isequal(state, 'done')  % Check if the optimization is done
            % Capture the last values when done
            [~, Ks, g, w] = objective(r);
            Ks_opt        = Ks;
            g_opt         = g;
            w_opt         = w; 
            
        end
end

end
