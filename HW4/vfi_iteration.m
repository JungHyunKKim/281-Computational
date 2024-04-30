function [v_new,c,A] = vfi_iteration(v0,rho,lambdaL,r,param,num,grids)

I = num.a_n ;

% -------------------------------------------------------------------------
% Solve HJB in matrix form using implicit method
% See Achdou et al. 2022 Numerical Appendix page 6
% 1. Construct Matrix A of size (2xI, 2xI)
%    1.1. Aswitch  : Construct matrix with jump components (lambdas) only
%                    (stochastic components)
%    1.2. vp_upwind: Find forward and backward difference approximations of the
%                    derivative of the value function and get v' 
%    1.3. A1, A2   : Construct matrix with only drift components for each lambda
%                    (deterministic components)
%    1.4. A        : Create block matrix A using A1, A2 and adding Aswitch
% 2. Solve for value function (Equation 16 of page 6
% -------------------------------------------------------------------------

% 1.1
% Construct matrix with jump components
Aswitch = [-speye(I)*lambdaL, speye(I)*param.lambda(1); ...
            speye(I)*lambdaL, -speye(I)*param.lambda(2)];


% 1.2
% Find forward and backward difference and get approximated v' (Ix2) matrix
[Va_Upwind,ssb,ssf] = vp_upwind(v0,r,param,num,grids) ;

c = Va_Upwind.^(-1);
u = utility(c);

% 1.3
% Construct diagonals of drift components in columns
X = - min(ssb,0)/grids.da;
Y = - max(ssf,0)/grids.da + min(ssb,0)/grids.da;
Z = max(ssf,0)/grids.da;
    
% Construct matrix with drift components
A1 = spdiags(Y(:,1), 0, I, I) + ...
     spdiags(X(2:I,1), -1, I, I) + ... % can directly plug in truncated column for regions below diagonal
     spdiags([0; Z(1:I-1,1)], 1, I, I); % for regions above diagonal, need to plug in column of size as large as the diagonal
     
A2 = spdiags(Y(:,2), 0, I, I) + ...
     spdiags(X(2:I,2), -1, I, I) + ... 
     spdiags([0; Z(1:I-1,2)], 1, I, I);

% 1.4
% Create block matrix A using A1, A2 and adding Aswitch
A = [A1, sparse(I,I); sparse(I,I), A2] + Aswitch;

% 2. Solve for value function after unstacking to (Ix2,1) column
B = (1/num.Delta + rho)*speye(2*I) - A;

u_stacked = [u(:,1); u(:,2)];
V_stacked = [v0(:,1); v0(:,2)];

b = u_stacked + V_stacked/num.Delta;
V_stacked = B\b; 

v_new = [V_stacked(1:I), V_stacked(I+1:2*I)]; % restack

end