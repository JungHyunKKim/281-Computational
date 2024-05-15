% Preferences
par.gamma  = 2; 
par.rho    = 0.01 ;

% Production function
par.logTFP = log(1); 
par.delta  = 0.025 ;
par.alpha  = 1/3 ;

% Aggregate shock
par.sigmaTFP   = 0.007 ; 
par.rhoTFP     = 0.95; 

% Idiosyncratic shocks
zz1        = 0; 
zz2        = 1; 
par.e      = [zz1 zz2] ;

par.evec   = ones(num.a_n,1) * par.e;
par.estack = reshape(par.evec,2*num.a_n,1);

% Transition probabilities
llambda1 = 1 / 2;  % expected duration of unemployment is 2 quarters
llambda2 = (llambda1 / (zz2 * .93 - zz1))*(zz2 - zz2 * .93); % unemployment rate 7%
par.lambda = [llambda1,llambda2];

% Tax system
par.mu = .15;        % UI replacement rate 15%
par.tau = (par.mu / zz2) * (par.lambda(2) / par.lambda(1));	     % labor income tax

% Labor supply (used as parameter not variable; also solved in solve_Aiyagari.m)
par.zAvg = (par.lambda(1) * par.e(2) + par.lambda(2) * par.e(1)) / (par.lambda(1) + par.lambda(2));

% Initialize wage
par.w0     = (1 - par.alpha) * exp(par.logTFP)^( 1 / (1-par.alpha) ) ... 
              *( par.alpha / (par.rho + par.delta) )^( par.alpha/(1-par.alpha) ); 

% Heterogeneous interest
par.rshock = [1 1]; 