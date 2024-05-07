function [grids] = create_grids(par,num)
% Uses structures of structural pareters and numerical pareters to
% create a grids and initial guesses for v.
grids.a = repmat(linspace(num.a_min,num.a_max,num.a_n)',[1 2]) ;
grids.v0 = utility(0.01*grids.a + par.e*1.5) ;
grids.da = grids.a(2) - grids.a(1) ;

end