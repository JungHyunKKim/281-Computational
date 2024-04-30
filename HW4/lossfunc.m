function [loss] = lossfunc(target_vec, weight_vec, rho, lambdaL, r, param, num, grids)

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

% Average wealth-to-income ratio
a2y = g(:,1)' * ( ( 1./param.y(:,1) ).*grids.a(:,1) ) * grids.da + ...
      g(:,2)' * ( ( 1./param.y(:,2) ).*grids.a(:,2) ) * grids.da; 

% Ratio btw avg high-income wealth to avg low-income wealth
ahigh2alow = (g(:,2)' * ( ( 1./param.y(:,2) ).*grids.a(:,2)) * grids.da) / ...
             (g(:,1)' * ( ( 1./param.y(:,1) ).*grids.a(:,1)) * grids.da) ; 

% Variance of wealth
exp_asq = g(:,1)' * grids.a(:,1).^2 * grids.da + ...
          g(:,2)' * grids.a(:,2).^2 * grids.da; 

exp_as = g(:,1)' * grids.a(:,1) * grids.da + ...
         g(:,2)' * grids.a(:,2) * grids.da; 

var_a = exp_asq - exp_as^2; 

% Average length of unemployment
exp_uespell = g(:,1)' * ( ( 1/lambdaL )*ones(num.a_n,1) * grids.da) + ...
              g(:,2)' * ( ( 1/lambdaL )*ones(num.a_n,1) * grids.da ); 

% Weighted loss
loss = ([a2y, ahigh2alow, var_a, exp_uespell] - target_vec).^2 * weight_vec'; 

end