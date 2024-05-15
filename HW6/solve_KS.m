% -------------------------------------------------------------------------
% Solves the Krusell and Smith (1998)
% -------------------------------------------------------------------------

function [vAggregateTFP, vAggregateOutput, vAggregateConsumption, vAggregateInvestment, vAggregateTFP_reduced, ...
          vAggregateOutput_reduced, vAggregateConsumption_reduced, vAggregateInvestment_reduced, vTime] ...
         = solve_KS(par, num, grids, T, N, vAggregateShock)

tstart = tic;

%% Step 1: Solve for Steady State
tStart = tic;

% Spline parameters
global n_splined
n_knots = 12;
c_power = 7;
x = grids.a(:,1)';
n_post = 2;	% This is from two income states
n_prior = 1;
n_splined = n_prior*n_knots*n_post;

fprintf('Computing steady state...\n')

r0 = par.rho*0.8;

[rSS,wSS,KSS,ASS,uSS,cSS,VSS,gSS,dVUSS,VafSS,VabSS,IfSS,IbSS,I0SS,Aswitch] = ...
    solve_Aiyagari(par, num, grids, 'GE', r0); 

fprintf('Time to compute steady state: %.3g seconds\n\n\n',toc(tStart));

% Store steady state values in column
varsSS = zeros(num.nVars,1);
varsSS(1:2*num.a_n,1) = reshape(VSS,2*num.a_n,1);
ggSS = reshape(gSS,2*num.a_n,1);
varsSS(2*num.a_n+1:4*num.a_n-1,1) = ggSS(1:2*num.a_n-1);
varsSS(4*num.a_n,1) = par.logTFP;
varsSS(4*num.a_n+1,1) = KSS;
varsSS(4*num.a_n+2,1) = rSS;
varsSS(4*num.a_n+3,1) = wSS;
varsSS(4*num.a_n+4,1) = (KSS ^ par.alpha) * (par.zAvg ^ (1 - par.alpha));
CSS = sum(cSS(:) .* gSS(:) * grids.da);
varsSS(4*num.a_n+5,1) = CSS;
varsSS(4*num.a_n+6,1) = par.delta * KSS;

%% Step 2: Linearize Model Equations
% For computing derivatives, the codes written for solving for the
%    steady-state can be used almost verbatim using automatic
%    differentiation toolbox as long as only the functions supported by
%    automatic differentation are used. For list of supported functions and
%    documentation of relevant syntax check <<https://github.com/sehyoun/MATLABAutoDiff>>
fprintf('Taking derivatives of equilibrium conditions...\n')
t0 = tic;

% Prepare automatic differentiation
vars = zeros(2*num.nVars+num.nEErrors+1,1);
vars = myAD(vars); %dependency of different entries w.r.t. other entries
% type

% Evaluate derivatives
derivativesIntermediate = equilibrium_conditions(vars, Aswitch, IfSS, IbSS, I0SS, varsSS, num, par, grids);

% Extract out derivative values
derivs = getderivs(derivativesIntermediate);

% Unpackage derivatives
mVarsDerivs = derivs(:,1:num.nVars);
mVarsDotDerivs = derivs(:,num.nVars+1:2*num.nVars);
mEErrorsDerivs = derivs(:,2*num.nVars+1:2*num.nVars+num.nEErrors);
mShocksDerivs = derivs(:,2*num.nVars+num.nEErrors+1);

% rename derivatives to match notation in paper
g0 = mVarsDotDerivs;
g1 = -mVarsDerivs;
c = sparse(num.nVars,1);
psi = -mShocksDerivs;
pi = -mEErrorsDerivs;


[state_red,inv_state_red,g0,g1,c,pi,psi] = clean_G0_sparse(g0,g1,c,pi,psi);
% clear out variables that are linear transformation of others
n_g_red = num.n_g;


%% Step 4: Solve Linear System
t0 = tic;
fprintf('Solving reduced linear system...\n')

[G1,~,impact,eu,F] = schur_solver(g0,g1,c,psi,pi,1,1,1); %as in gensys / dynare 

fprintf('...Done!\n')
fprintf('Existence and uniqueness? %2.0f and %2.0f\n',eu);
fprintf('Time to solve linear system: %2.4f seconds\n\n\n',toc(t0))

%% Step 5: Simulate Impulse Response Functions
fprintf('Simulating Model...\n')
t0 = tic;

trans_mat = inv_state_red;
[simulated,vTime] = simulate(G1,impact,T,N,vAggregateShock,'implicit',trans_mat);

fprintf('...Done!\n')
fprintf('Time to simulate model: %2.4f seconds\n\n\n',toc(t0))

% Add state-states back in to get values in levels
%varsSS_small = varsSS(4*num.a_n:4*num.a_n+6,1);
vAggregateTFP = simulated(400,:) + varsSS(400);
vAggregateOutput = simulated(404,:) + varsSS(404);
vAggregateConsumption = simulated(405,:) + varsSS(405);
vAggregateInvestment = simulated(406,:) + varsSS(406);

% Compute log differences for plotting
vAggregateTFP_reduced = vAggregateTFP - varsSS(400);
vAggregateOutput_reduced = log(vAggregateOutput) - log(varsSS(404));
vAggregateConsumption_reduced = log(vAggregateConsumption) - log(varsSS(405));
vAggregateInvestment_reduced = log(vAggregateInvestment) - log(varsSS(406));

toc(tstart)


