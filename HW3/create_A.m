function [A] = create_A(sb, sf, grid, num)

% -------------------------------------------------------------------------
% Compute sparse matrix A
% Matlab function: spdiags
% (1) [Bout, id] = spdiags(A)  : extract diagonals from A
% (2) S = spdiags(Bin, d, m,n) : Create m x n sparse matrix, placing
%                                d columns of bins as diagonals
% When m>=n, 
% (1) Diagonals below the main diagonal take elements from the tops of the columns first.
% (2) Diagonals above the main diagonal take elements from the bottoms of columns first.
% When m<n, 
% (1) Diagonals above the main diagonal take elements from the tops of the columns first.
% (2) Diagonals below the main diagonal take elements from the bottoms of columns first.
% 
% ** In our case, the diagonal below main should take elements from the
%    bottoms of columns first, and vice versa --> transpose the matrix
% -------------------------------------------------------------------------

% Get columns of flows used as diagonals for A
inflow_lag  = -min(sb, 0) ./ grid.da; 
inflow_lead = max(sf, 0) ./ grid.da;
outflow     = -inflow_lag - inflow_lead; 

% Create sparse matrix A
Aprime = spdiags([inflow_lead outflow inflow_lag], -1:1, num.a_n, num.a_n); 
A      = Aprime.'; 
%spy(A); 

% % Check if diagols are correctly fed in: 
% Afull = full(A);
% Afull(1:12,1:12)
% inflow_lag(1:14)
% sum(spdiags(A) - [ [inflow_lag(2:end)', 0]', outflow, [0, inflow_lead(1:end-1)']' ])

% % Moll's easier code: 
% A = spdiags(inflow_lag(2:end), -1, num.a_n, num.a_n) + ...
%     spdiags(outflow, 1, num.a_n, num.a_n) + ...
%     spdiags(inflow_lead(1:end-1), 1, num.a_n, num.a_n)

end
