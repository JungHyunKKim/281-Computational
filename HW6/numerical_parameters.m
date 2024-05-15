num.Delta = 1e4; %100; %0.5 ;
num.a_n = 100 ;
num.a_min = 1e-5; % setting to 0 can be unstable (especially for higher gamma & 0 income)
num.a_max = 100; 
num.tol = 1e-7 ;

% Number of variables in the system
num.n_v = 2 * num.a_n; % value function for each asset income level
num.n_g = 2 * num.a_n-1 + 1; % distribution for each asset and income level
num.n_p = 6; % par.logTFP, KSS, rSS, wSS, ySS, cSS, iSS
num.nVars = num.n_v + num.n_g + num.n_p;
num.nEErrors = 2 * num.a_n;

% Note: 
% varsSS = zeros(nVars,1);
% varsSS(1:2*I,1) = reshape(VSS,2*I,1);
% ggSS = reshape(gSS,2*I,1);
% varsSS(2*I+1:4*I-1,1) = ggSS(1:2*I-1);
% varsSS(4*I,1) = par.logTFP;
% varsSS(4*I+1,1) = KSS;
% varsSS(4*I+2,1) = rSS;
% varsSS(4*I+3,1) = wSS;
% varsSS(4*I+4,1) = (KSS ^ par.alpha) * (zAvg ^ (1 - par.alpha));
% CSS = sum(cSS(:) .* gSS(:) * grids.da);
% varsSS(4*I+5,1) = CSS;
% varsSS(4*I+6,1) = par.delta * KSS;
