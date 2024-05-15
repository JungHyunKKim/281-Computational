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

function [rSS,wSS,KSS,ASS,uSS,cSS,VSS,gSS,Va_UpwindSS,VafSS,VabSS,IfSS,IbSS,I0SS,AswitchSS] = ...
         solve_Aiyagari(par, num, grids, eqtype, r0)

% -------------------------------------------------------------------------
% Minimization
% -------------------------------------------------------------------------

tic

% Variables to hold final results
w_opt = []; 
Ks_opt = []; 
A_opt  = []; 
u_opt = []; 
c_opt = []; 
V_opt  = []; 
g_opt = []; 
Va_Upwind_opt = []; 
Vaf_opt = []; 
Vab_opt = []; 
If_opt = []; 
Ib_opt = []; 
I0_opt = []; 
Aswitch_opt = []; 

% Redefine the objective function to return squared excess demand and capture g and Ks
function [error,A,u,c,v_new,g,Va_Upwind,Vaf,Vab,If,Ib,I0,aprimedot,Ks,Ls,Kd,w,Aswitch] ...
          = objective(r)
     
    [Kd_excess,A,u,c,v_new,g,Va_Upwind,Vaf,Vab,If,Ib,I0,aprimedot,Ks,Ls,Kd,w,Aswitch] ...
          = solve_Kd_excess(r, par, num, grids, eqtype);

    error = Kd_excess^2;
    
end

% Set the value of the Outputfcn field of the options structure to a function handle to outfun.
opts = optimset('Display','iter','TolFun',1e-08,'OutputFcn', @outfun);

% Solve for the interest rate that minimizes the squared excess demand
[rSS, error, ~, ~] = fminsearch(@objective, r0, opts);

wSS = w_opt; 
KSS = Ks_opt; 
ASS = A_opt; 
uSS = u_opt; 
cSS = c_opt; 
VSS = V_opt; 
gSS = g_opt; 
Va_UpwindSS = Va_Upwind_opt; 
VafSS = Vaf_opt; 
VabSS = Vab_opt; 
IfSS = If_opt; 
IbSS = Ib_opt; 
I0SS = I0_opt; 
AswitchSS = Aswitch_opt; 

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
            [~,A,u,c,v_new,g,Va_Upwind,Vaf,Vab,If,Ib,I0,~,Ks,~,~,w,Aswitch] = objective(r);
            
            w_opt = w; 
            Ks_opt = Ks; 
            A_opt  = A; 
            u_opt = u; 
            c_opt = c; 
            V_opt  = v_new; 
            g_opt = g; 
            Va_Upwind_opt = Va_Upwind; 
            Vaf_opt = Vaf; 
            Vab_opt = Vab; 
            If_opt = If; 
            Ib_opt = Ib; 
            I0_opt = I0; 
            Aswitch_opt = Aswitch; 

        end
end

end
