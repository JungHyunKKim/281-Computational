function [v_new,c] = vfi_iteration(v_old, param, num, grid, method)

[vprime_upwind, sb, sf] = vp_upwind(v_old, param, num, grid);

% Infer consumption from Va_Upwind
c = vprime_upwind.^(-1);
u = utility(c);

switch method

    case 'explicit'
        % Update the value function with a step parameter Delta
        vchange_step = u + vprime_upwind.*(param.y + param.r*grid.a - c) - param.rho*v_old; 
        v_new        = v_old + num.Delta*vchange_step;

    case 'explicitM' % HJB in matrix form
        [A]          = create_A(sb, sf, grid, num); 
        % Update the value function with a step parameter Delta
        vchange_step = u + A*v_old - param.rho*v_old; 
        v_new        = v_old + num.Delta*vchange_step; 

    case 'implicit'
        [A]          = create_A(sb, sf, grid, num); 
        B            = (param.rho + 1/num.Delta)*speye(num.a_n) - A; 
        b            = u + 1/num.Delta * v_old; 
        v_new        = B\b; 

    otherwise
        error('Invalid method'); 
end

end