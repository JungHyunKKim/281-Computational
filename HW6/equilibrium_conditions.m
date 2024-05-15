% Inputs:  (1) vars: vector which contains the variables in the system, the time derivatives of 
%			         those variables, the expectational errors, and the shocks
%
% Outputs: (1) vResduals: residuals of equilibrium conditions, evaluated at vars

function vResidual = equilibrium_conditions(vars, Aswitch, IfSS, IbSS, I0SS, varsSS, num, par, grids)

%----------------------------------------------------------------
% Housekeeping
%----------------------------------------------------------------
	
% Unpack vars that you want to be accurate (need to keep track of the order)
V=vars(1:2*num.a_n) + varsSS(1:2*num.a_n);
g=vars(2*num.a_n+1:4*num.a_n-1) + varsSS(2*num.a_n+1:4*num.a_n-1);	% vector of distribution, removing last point		
g_end=1/grids.da-sum(g);		% ensures that distribution integrates to 1
logAggregateTFP = vars(4*num.a_n);
KHat = vars(4*num.a_n+1) + varsSS(4*num.a_n+1);
rHat = vars(4*num.a_n+2) + varsSS(4*num.a_n+2);
wHat = vars(4*num.a_n+3) + varsSS(4*num.a_n+3);
output = vars(4*num.a_n+4) + varsSS(4*num.a_n+4);
C = vars(4*num.a_n+5) + varsSS(4*num.a_n+5);
investment = vars(4*num.a_n+6) + varsSS(4*num.a_n+6);


V = reshape(V,num.a_n,2);

VDot = vars(num.nVars+1:num.nVars+2*num.a_n);
gDot = vars(num.nVars+2*num.a_n+1:num.nVars+4*num.a_n-1);
logAggregateTFPDot = vars(num.nVars+4*num.a_n);

VEErrors = vars(2*num.nVars+1:2*num.nVars+2*num.a_n); %keep track of expentational errors

aggregateTFPShock = vars(2*num.nVars+num.nEErrors+1);

% Initialize other variables, using vars to ensure everything is a dual number
dVf = V;
dVb = V;

K = sum(grids.aaa .* [g;g_end] * grids.da);
r = exp(logAggregateTFP) * par.alpha * (KHat ^ (par.alpha - 1)) * (par.zAvg ^ (1 - par.alpha)) - par.delta;
w = exp(logAggregateTFP) * (1 - par.alpha) * (KHat ^ par.alpha) * (par.zAvg ^ (-par.alpha)); 

%----------------------------------------------------------------
% Compute one iteration of HJB Equation
%----------------------------------------------------------------

c0 = w * ((1 - par.tau) * par.evec + par.mu * (1 - par.evec)) + r * grids.a;

% Compute forward difference
dVf(1:num.a_n-1,:) = (V(2:num.a_n,:)-V(1:num.a_n-1,:))/grids.da;
dVf(num.a_n,:) = marginal_utility(c0(num.a_n,:), par); %will never be used, but impose state constraint a<=amax just in case

% Compute backward difference
dVb(2:num.a_n,:) = (V(2:num.a_n,:)-V(1:num.a_n-1,:))/grids.da;
dVb(1,:) = marginal_utility(c0(1,:), par); %state constraint boundary condition

% Compute consumption and savings with forward difference
cf = inv_marginal_utility(dVf, par);
ssf = c0 - cf;

% Compute consumption and savings with backward difference
cb = inv_marginal_utility(dVb, par);
ssb = c0 - cb;

% Compute consumption and derivative of value function for no drift
dV0 = marginal_utility(c0, par);

% Compute upwind difference
dV_Upwind = dVf.*IfSS + dVb.*IbSS + dV0.*I0SS;
c = inv_marginal_utility(dV_Upwind, par);
u = utility(c, par);
savings = c0 - c;

% Construct A matrix
X = -ssb.*IbSS/grids.da;
Y = -ssf.*IfSS/grids.da + ssb.*IbSS/grids.da;
Z = ssf.*IfSS/grids.da;

X(1,:)=0;
lowdiag=reshape(X,2*num.a_n,1);
Z(num.a_n,:)=0;

A = spdiags(reshape(Y,2*num.a_n,1),0,2*num.a_n,2*num.a_n)...
    +spdiags(lowdiag(2:2*num.a_n),-1,2*num.a_n,2*num.a_n)...
    +spdiags([0,reshape(Z,1,2*num.a_n)]',1,2*num.a_n,2*num.a_n)...
    +Aswitch;

%----------------------------------------------------------------
% Compute equilibrium conditions
%----------------------------------------------------------------

% HJB Equation
hjbResidual = reshape(u,2*num.a_n,1) + A * reshape(V,2*num.a_n,1) + VDot + VEErrors - par.rho * reshape(V,2*num.a_n,1);

% KFE 
gIntermediate = A' * [g;g_end];
gResidual = gDot - gIntermediate(1:2*num.a_n-1,1);

% Aggregates
kResidual = K - KHat;
rResidual = r - rHat;
wResidual = w - wHat;
yResidual = output - exp(logAggregateTFP) * (sum(grids.aaa .* [g;g_end] * grids.da) ^ par.alpha) * (par.zAvg ^ (1 - par.alpha));
cResidual = C - sum(c(:) .* [g;g_end] * grids.da);
iResidual = investment - sum((savings(:) + par.delta * grids.aaa) .* [g;g_end] * grids.da); 

% Law of motion for aggregate shocks
tfpResidual = logAggregateTFPDot + (1 - par.rhoTFP) * logAggregateTFP - par.sigmaTFP * aggregateTFPShock;

% Order should be same as how you unpacked vars
vResidual = [hjbResidual;gResidual;tfpResidual;kResidual;rResidual;wResidual;yResidual;cResidual;iResidual];
