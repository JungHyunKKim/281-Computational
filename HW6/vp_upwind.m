function [Va_Upwind, Vaf, Vab, sf, sb, If, Ib, I0] = vp_upwind(v0,r,w,par,num,grids)

V = v0 ;

% Initialize forward and backwards differences with zeros: Vaf Vab
Vaf = zeros(num.a_n,2) ;
Vab = zeros(num.a_n,2) ;

% Use v0 to compute backward and forward differences
Vaf(1:end-1,:) = (V(2:end,:) - V(1:end-1,:))/grids.da ;
Vaf(end,:) = 0; 

Vab(2:end,:) = (V(2:end,:)-V(1:end-1,:))./grids.da;
Vab(1,:) = marginal_utility(par.rshock.*r.*grids.a(1,:) + w*((1 - par.tau) * par.e + par.mu * (1 - par.e)), par); %state constraint boundary condition    

% Consumption and savings with forward difference
cf = inv_marginal_utility(max(Vaf,1e-08), par);
sf = w*((1 - par.tau) * par.e + par.mu * (1 - par.e)) + par.rshock.*r.*grids.a - cf ;

% Consumption and savings with backward difference
cb = inv_marginal_utility(max(Vab,1e-08), par);
sb = w*((1 - par.tau) * par.e + par.mu * (1 - par.e)) + par.rshock.*r.*grids.a - cb ;

% Consumption and derivative of value function at steady state
c0 = w*((1 - par.tau) * par.e + par.mu * (1 - par.e)) + par.rshock.*r.*grids.a ;
Va0 = marginal_utility(max(c0,1e-08), par);

If = sf > 0; % Positive drift --> forward difference
Ib = sb < 0; % Negative drift --> backward difference
I0 = (1-If-Ib); % At steady state
Va_Upwind = Vaf.*If + Vab.*Ib + Va0.*I0; % Important to include third term

end