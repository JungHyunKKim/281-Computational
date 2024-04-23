function [grids] = create_grid(param,num)

% Uses structures of structural parameters and numerical parameters to
% create a grid and initial guesses for v.
grids.a = linspace(num.a_min,num.a_max,num.a_n)' ;
grids.v0 = utility(0.1*grids.a) ;
grids.da = grids.a(2) - grids.a(1) ;

end