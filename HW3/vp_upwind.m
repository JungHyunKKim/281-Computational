function [vprime_upwind, sb, sf] = vp_upwind(value,param,num,grid)

% Unpack the initial guess for value and call it V
V = value; 

% Initialize forward and backwards differences with zeros Vaf Vab
Vaf = zeros(num.a_n,1); 
Vab = zeros(num.a_n,1); 

% Use V to compute backward and forward differences
Vaf(1:end-1) = (V(2:end) - V(1:end-1))/grid.da; 
Vaf(end)     = 0; 

Vab(2:end)   = (V(2:end) - V(1:end-1))/grid.da; 
Vab(1)       = (param.r*grid.a(1) + param.y(1))^(-1); 

% Impose the following boundary conditions
Vaf(end)     = 0; 
Vab(1)       = (param.r*grid.a(1) + param.y(1))^(-1); %state constraint boundary condition    

% Consumption and savings with forward difference
cf = max(Vaf,1e-08).^(-1);
sf =  param.r*grid.a + param.y - cf;

% Consumption and savings with backward difference
cb = max(Vab,1e-08).^(-1);
sb = param.r*grid.a + param.y - cb;

% Consumption and derivative of value function at steady state
c0 = param.r*grid.a + param.y;
Va0 = c0.^(-1);

% Compute indicator functions that capture the upwind scheme.
If = sf > 0; 
Ib = sb < 0; 
I0 = (1 - If - Ib); 

% Compute the upwind scheme
vprime_upwind = Vaf.*If + Vab.*Ib + Va0.*I0; %important to include third term

end