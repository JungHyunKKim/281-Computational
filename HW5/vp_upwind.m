function [Va_Upwind,sb,sf] = vp_upwind(v0,r,w,par,num,grids)

V = v0 ;

% Initialize forward and backwards differences with zeros: Vaf Vab
Vaf = zeros(num.a_n,2) ;
Vab = zeros(num.a_n,2) ;

% Use v0 to compute backward and forward differences
Vaf(1:end-1,:) = (V(2:end,:) - V(1:end-1,:))/grids.da ;
Vaf(end,:) = 0; 

Vab(2:end,:) = (V(2:end,:)-V(1:end-1,:))./grids.da;
Vab(1,:) = (par.rshock.*r.*grids.a(1,:) + w.*par.e).^(-1); %state constraint boundary condition    

% Consumption and savings with forward difference
cf = max(Vaf,1e-08).^(-1);
sf = w.*par.e + par.rshock.*r.*grids.a - cf ;

% Consumption and savings with backward difference
cb = max(Vab,1e-08).^(-1);
sb = w.*par.e + par.rshock.*r.*grids.a - cb ;

% Consumption and derivative of value function at steady state
c0 = w.*par.e + par.rshock.*r.*grids.a ;
Va0 = c0.^(-1);

If = sf > 0; % Positive drift --> forward difference
Ib = sb < 0; % Negative drift --> backward difference
I0 = (1-If-Ib); % At steady state
Va_Upwind = Vaf.*If + Vab.*Ib + Va0.*I0; % Important to include third term

end