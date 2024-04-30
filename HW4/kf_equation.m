function [gg] = kf_equation(A,grids,num)

% Transpose A 
AT = A';

% Stack state space grid into one column (Ix2,1)
b  = zeros(2*num.a_n,1);

% Need to fix one value; otherwise matrix is singular
% -------------------------------------------------------------------------
% (1) 0 = -d[s_{i}(a)g_{i}(a)]/da - lambda_i g_{i}(a) + lambda_j g_{j}(a)
% (2) 1 = sum_i g_{i,1} delta a + sum_i g_{i,1} g_{i,2} delta a
%
% Solve (1) by solving eigenvalue problem: A'g = b = 0 while imposing (2): 
%    1. Fix g_{i,j} = 0.1 (or any number) : 
%         - replace corresponding entry of zero vector in RHS of (1) by 0.1
%         - replace corresponding row of A' by rows of zeros and 1 at diagonal
%    2. Get \tilde_g and renormalize: 
%       g_{i,j} = \tilde_g_{i,j} / 
%                 sum_i \tilde_g_{i,1} delta a + sum_i \tilde_g_{i,1} g_{i,2} delta a
% 
% -------------------------------------------------------------------------

i_fix       = 1; % pick any index
b(i_fix)    = .1; % replace the zero verctor to 0.1 at the index
row         = [zeros(1,i_fix-1), 1, zeros(1,2*num.a_n-i_fix)]; % replacement for A' [empty, 1, zeros]
AT(i_fix,:) = row;

gg    = AT\b; % solve system
g_sum = gg'*ones(2*num.a_n,1)*grids.da; % denominator for renormalization
gg    = gg./g_sum; % renormalize g

end