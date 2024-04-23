% -------------------------------------------------------------------------
% 281-Computational Homework 3
% Jung Hyun Kim
% HJB for Consumption-Saving Problem Without Uncertainty
% ODE Solution Method : Finite-Difference Method
%                       (1) Explicit (vector form) : explicit.m
%                       (2) Explicit (matrix form) : explicitM.m
%                       (3) Implicit (matrix form) : implicit.m
% Important Functions           
%   1. vfi_iteration  : Iteration using explicit.m, explicitM.m, or implicit.m
%   2. vp_upwind      : Derivative of value function using upwind method
%   3. create_A       : Creates the intensity matrix 
%   4. solve_golden   : Maximum solver using golden section method
% -------------------------------------------------------------------------

clear; close all; clc;

% -------------------------------------------------------------------------
% Question 0 : Compare explicit, explicit (matrix form), and implicit at 
%              step-size of Delta = 0.01 (baseline)
% 
% Answers    : Given Delta = 0.01, the number of iterations for all methods
%              are similar. The elapsed time is longest for the implicit 
%              method, followed by explicit method (matrix form) and 
%              explicit method (vector form)
% Code Output     : 
% Explicit method (vector form) converged in 81755 iterations. Elapsed time is 1.173298 seconds.
% Explicit method (matrix form) converged in 81755 iterations. Elapsed time is 3.506651 seconds.
% Implicit method converged in 81764 iterations. Elapsed time is 4.478139 seconds.
% -------------------------------------------------------------------------

% create structure with structural parameters
call_parameters;

% create structure with numerical parameters
numerical_parameters;

% create structure with grids and initial guesses for v
grids = create_grid(param,num);

% Explicit method (vector form)
explicit; 

v_explicit = v_new; 
c_explicit = c; 
adot_explicit = adot;

% Explicit method (Matrix form)
explicitM; 

v_explicitM = v_new; 
c_explicitM = c; 
adot_explicitM = adot;

% Implicit method (Matrix form)
implicit; 

v_implicit = v_new; 
c_implicit = c; 
adot_implicit = adot;

% Figures
figure_Q0; 

% -------------------------------------------------------------------------
% Question 1 : Compare explicit, explicit (matrix form), and implicit at 
%              step-size of Delta = 0.5
% 
% Answers    : Given Delta = 0.5 and max_iter=100000, the explicit methods 
%              do not converge, while the implicit method converges within
%              much shorter iterations and elapsed time (lower than all 
%              methods under Delta=0.01). As shown in the figure, the
%              approximation is also accurate. 
% Code output : 
% Explicit method (vector form) did not converge in 100000 iterations. Elapsed time is 4.810668 seconds.
% Explicit method (matrix form) did not converge in 100000 iterations. Elapsed time is 8.410006 seconds.
% Implicit method converged in 2421 iterations. Elapsed time is 0.129145 seconds.
% -------------------------------------------------------------------------

% Larger step-size 
num.Delta = 0.5; 

% Explicit method (vector form)
explicit; 

v_explicit_smallD = v_new; 
c_explicit_smallD = c; 
adot_explicit_smallD = adot;

% Explicit method (Matrix form)
explicitM; 

v_explicitM_smallD = v_new; 
c_explicitM_smallD = c; 
adot_explicitM_smallD = adot;

% Implicit method (Matrix form)
implicit; 

v_implicit_smallD = v_new; 
c_implicit_smallD = c; 
adot_implicit_smallD = adot;

% Figures
figure_Q1; 

% -------------------------------------------------------------------------
% Question 2 : Compare solutions using discrete and continuous time methods
% 
% Answers    : The solutions between discrete and continuous time
%              frameworks coincide if the discount rate is set in the right
%              way, i.e. e^(-rho*t) = beta^t <-> rho = -log(beta) and if 
%              r = rho. 
% -------------------------------------------------------------------------

%%%%%% Set r = rho %%%%%%

% Implicit method (Matrix form)
implicit; 

v_implicit = v_new; 
c_implicit = c; 

% Discrete-time true solution :
c_discrete = param.y + param.r * grids.a; 
v_discrete = utility(c_discrete) / (1 - param.beta); 

% Figures
figure_Q2;
saveas(gcf,'Q2_baseline.png')


%%%%%% Set r = 0 %%%%%%

param.r = 0; 

% Implicit method (Matrix form)
implicit; 

v_implicit = v_new; 
c_implicit = c; 
adot_implicit = adot;

% Discrete-time true solution :
c_discrete = param.y + param.r * grids.a; 
v_discrete = utility(c_discrete) / (1 - param.beta); 

% Figures
figure_Q2;
saveas(gcf,'Q2_r0.png')
delete('Q2.png')


% -------------------------------------------------------------------------
% Question 3 : Solve two-stage problem
%              Stage 1 : Profit maximization given r
%              Stage 2 : Utility maximization given r, maximized profit
%
% Answers    : See Q3.png (Comparison between endowment & production
%                          economy)
% -------------------------------------------------------------------------

% create structure with structural parameters
call_parameters;

% create structure with numerical parameters
numerical_parameters;
num.Delta = 0.5; 

% create structure with grids and initial guesses for v
grids = create_grid(param,num);

% Stage 1 : profit maximization given r
k_max = 1000; 
[kstar, pistar]  = solve_golden('profitfunc_bad', 0, k_max, param); 

% Stage 2 : HH maximization given r
param.y = pistar; 
num.max_iter = 1000000; 
% Implicit method (Matrix form)
implicit; 
v_firm = v_new; 
c_firm = c; 
adot_firm = adot;

% Figures
figure_Q3;


% -------------------------------------------------------------------------
% Question 4 : Choice of technology + 
%              Non-convex production technologies + financial frictions
%
%              Stage 1 : Profit maximization given r
%              Stage 2 : Utility maximization given r, maximized profit
%
% Note       : kappa almost negligible -> use instead [0 0.5 2]
% Answers    : See Q4.png
% -------------------------------------------------------------------------

% create structure with structural parameters
call_parameters;

% create structure with numerical parameters
numerical_parameters;
num.Delta = 0.5; 
num.max_iter  = 1000000; 
k_max         = grids.a; 

% create structure with grids and initial guesses for v
grids = create_grid(param,num);

kappa_vec     = [0 0.5 2]; 
vec_size      = length(kappa_vec); 
v_mat         = NaN(num.a_n,vec_size); 
c_mat         = NaN(num.a_n,vec_size); 
adot_mat      = NaN(num.a_n,vec_size); 
pistar_mat    = NaN(num.a_n,vec_size); 
kstar_mat     = NaN(num.a_n,vec_size); 
good_tech_mat = NaN(num.a_n,vec_size); 

for kk = 1:vec_size

    param.kappa = kappa_vec(kk);

    % Stage 1 : profit maximization given r, a 
    [kstar, pistar]                = solve_golden('profitfunc', 0, k_max, param); 
    pistar_mat(:,kk)               = pistar; 
    kstar_mat(:,kk)                = kstar; 
    good_tech_mat(:,kk)            = profitfunc_good(grids.a, param) >= profitfunc_bad(grids.a, param); 

    % Stage 2 : HH maximization given r
    param.y        = pistar_mat(:,kk); 
    implicit; 
    v_mat(:,kk)    = v_new; 
    c_mat(:,kk)    = c; 
    adot_mat(:,kk) = adot;

end

% Figures
figure_Q4;
