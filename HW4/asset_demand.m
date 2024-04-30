function [asset_D] = asset_demand(rho, lambdaL, r, param, num, grids)

% -------------------------------------------------------------------------
% Value Function Iteration
% -------------------------------------------------------------------------
 
v_old = grids.v0 ;

dist = 1 ;
while dist > num.tol 
    [v_new,~,~] = vfi_iteration(v_old,rho,lambdaL,r,param,num,grids) ;
    dist = max(abs((v_new(:) - v_old(:)))) ;
    v_old = v_new ;
end

[v_new,c,A] = vfi_iteration(v_new,rho,lambdaL,r,param,num,grids) ;

% -------------------------------------------------------------------------
% Komogorov Forward Equation: 
% -------------------------------------------------------------------------

[gg] = kf_equation(A,grids,num) ;
g    = [gg(1:num.a_n) , gg(num.a_n+1:2*num.a_n)]; % restack to Ix2 matrix

% -------------------------------------------------------------------------
% Aggregation
% -------------------------------------------------------------------------

asset_D = g(:,1)'*grids.a(:,1)*grids.da + g(:,2)'*grids.a(:,2)*grids.da; 

end