function [Ad, g] = solve_Ad(r, w, par, num, grids)

% -------------------------------------------------------------------------
% Value Function Iteration
% -------------------------------------------------------------------------
 
v_old = grids.v0 ;

dist = 1 ;
while dist > num.tol 
    [v_new,~,~] = vfi_iteration(v_old,r,w,par,num,grids) ;
    dist = max(abs((v_new(:) - v_old(:)))) ;
    v_old = v_new ;
end

[~,~,A] = vfi_iteration(v_new,r,w,par,num,grids) ;

% -------------------------------------------------------------------------
% Komogorov Forward Equation: 
% -------------------------------------------------------------------------

[gg] = kf_equation(A,grids,num) ;
g    = [gg(1:num.a_n) , gg(num.a_n+1:2*num.a_n)]; % restack to Ix2 matrix

% -------------------------------------------------------------------------
% Aggregation
% -------------------------------------------------------------------------

Ad   = g(:,1)'*grids.a(:,1)*grids.da + g(:,2)'*grids.a(:,2)*grids.da; 

end