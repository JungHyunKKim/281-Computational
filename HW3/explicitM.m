% -------------------------------------------------------------------------
% Explicit Method: Matrix Form
% -------------------------------------------------------------------------

v_old         = grids.v0;
dist          = 1; 
upwind_method = 'explicitM'; 
converged     = false; 

tic;
for iter = 1:num.max_iter
    [v_new,~]  = vfi_iteration(v_old, param, num, grids, upwind_method); 
    dist       = max(abs((v_new - v_old))); 

    if dist < num.tol
        converged = true; 
        break
    end

    v_old      = v_new; 
end
    [v_new, c] = vfi_iteration(v_old, param, num, grids, upwind_method); 
    adot = param.r*grids.a + param.y - c;
disp('---------------------------------------------------------')
if converged
    fprintf('Explicit method (matrix form) converged in %d iterations. ', iter)
else
    fprintf('Explicit method (matrix form) did not converge in %d iterations. ', num.max_iter)
end
toc;

v_explicitM = v_new; 
c_explicitM = c; 
adot_explicitM = adot;
