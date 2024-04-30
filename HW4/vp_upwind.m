function [Va_Upwind,sb,sf] = vp_upwind(v0,r,param,num,grids)

V = v0 ;

% Initialize forward and backwards differences with zeros: Vaf Vab
Vaf = zeros(num.a_n,2) ;
Vab = zeros(num.a_n,2) ;

% Use v0 to compute backward and forward differences
Vaf(1:end-1,:) = (V(2:end,:) - V(1:end-1,:))/grids.da ;
Vaf(end,:) = 0; 

Vab(2:end,:) = (V(2:end,:)-V(1:end-1,:))./grids.da;
Vab(1,:) = (r*grids.a(1,:) + param.y).^(-1); %state constraint boundary condition    

% Consumption and savings with forward difference
cf = max(Vaf,1e-08).^(-1);
sf = param.y + r.*grids.a - cf ;

% Consumption and savings with backward difference
cb = max(Vab,1e-08).^(-1);
sb = param.y + r.*grids.a - cb ;

% Consumption and derivative of value function at steady state
c0 = param.y + r.*grids.a ;
Va0 = c0.^(-1);

If = sf > 0; % Positive drift --> forward difference
Ib = sb < 0; % Negative drift --> backward difference
I0 = (1-If-Ib); % At steady state
Va_Upwind = Vaf.*If + Vab.*Ib + Va0.*I0; % Important to include third term

end