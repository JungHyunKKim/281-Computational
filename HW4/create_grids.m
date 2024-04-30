function [grids] = create_grids(param,num)
% Uses structures of structural parameters and numerical parameters to
% create a grids and initial guesses for v.
grids.a = repmat(linspace(num.a_min,num.a_max,num.a_n)',[1 2]) ;
grids.v0 = utility(0.1*grids.a + param.y) ;
grids.da = grids.a(2) - grids.a(1) ;

end