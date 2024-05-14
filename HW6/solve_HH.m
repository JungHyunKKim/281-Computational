function [c, aprimedot, Ad, Ls, g] = solve_HH(r, w, par, num, grids)

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

[~,c,A] = vfi_iteration(v_new,r,w,par,num,grids) ;
aprimedot = r*grids.a + w*repmat(par.e, num.a_n, 1) - c; 

% -------------------------------------------------------------------------
% Komogorov Forward Equation
% -------------------------------------------------------------------------

[gg] = kf_equation(A,grids,num) ;
g    = [gg(1:num.a_n) , gg(num.a_n+1:2*num.a_n)]; % restack to Ix2 matrix

% -------------------------------------------------------------------------
% Aggregation
% -------------------------------------------------------------------------

Ad   = sum(sum(grids.a.*g.*grids.da)); 
Ls   = sum(sum(par.e.*g.*grids.da)); 

end