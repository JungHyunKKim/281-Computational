par.e      = [0.1 0.2] ;
par.rho    = 0.03 ;
par.lambda = [0.1 0.1] ;
par.delta  = 0.05 ;
par.alpha  = 0.35 ;
par.z      = 1 ; 
par.w0     = (1 - par.alpha) * par.z^( 1 / (1-par.alpha) ) ... % initialize wage
                 *( par.alpha / (par.rho + par.delta) )^( par.alpha/(1-par.alpha) ); 
par.rshock = [1 1]; 