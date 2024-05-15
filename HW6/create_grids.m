function [grids] = create_grids(par,num)

% Uses structures of structural pareters and numerical parameters to
% create grids and initial guesses for v.
grids.a = repmat(linspace(num.a_min,num.a_max,num.a_n)',[1 2]) ;

w0  = (1 - par.alpha) * exp(par.logTFP)^( 1 / (1-par.alpha) ) ...
     *( par.alpha / (par.rho*0.8 + par.delta) )^( par.alpha/(1-par.alpha) ); 

grids.v0 = utility(par.rho*0.8*grids.a + w0*((1 - par.tau) * par.e + par.mu * (1 - par.e)), par); 
% Note: doesn't converge any initial is -inf: ex) utility(0.01*grids.a + par.e*1.5);

grids.da = grids.a(2) - grids.a(1) ;

grids.aaa = reshape(grids.a, 2*num.a_n, 1); % stack to column

end