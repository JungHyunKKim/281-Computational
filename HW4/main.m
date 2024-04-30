% -------------------------------------------------------------------------
% 281-Computational Homework 4
% Jung Hyun Kim
%
% Part 1. Indirect Inference
% Part 2. Homework
%
% Important functions  :
%       lossfunc.m     : function of rho, lambdaL, r ; 
%                        weighted sum of squared difference btw moments of
%                        interest and the empirical target
%       asset_demand.m : compute aggregate asset demand given r
%
% Functions from class : create_grids.m, kf_equation.m, utility.m, 
%                        vfi_iteration.m, vp_upwind.m
% -------------------------------------------------------------------------
clear ; close all ; clc ;

% Create structure with structural parameters
call_parameters;

% Create structure with numerical parameters
numerical_parameters ;

% Create structure with grids and initial guesses for v
grids = create_grids(param,num) ;

% -------------------------------------------------------------------------
% Part 1.1) Calibration of rho to target avg wealth-to-income ratio of 3
% -------------------------------------------------------------------------
rho0       = param.rho; 
lambdaL    = param.lambda(1); 
r          = param.r; 
target_vec = [3, 2, 0.02, 3]; 
weight_vec = [1, 0, 0, 0]; 

% Define loss
loss       = @(rho) lossfunc(target_vec, weight_vec, rho, lambdaL, r, param, num, grids);  

% Minimize loss
tic; 
[rhostar, minloss]    = fminsearch(loss,rho0); 
toc; 

fprintf(['rho is calibrated as %f to target the average wealth-to-income ratio of 3 \n' ...
         'Minimized squared difference between simulated and empirical moments is %f \n'], ...
         [rhostar,minloss]);

% Code output: rho is calibrated as 0.038086 to target the average wealth-to-income ratio of 3 
%              Minimized squared difference between simulated and empirical moments is 0.000029 

% -------------------------------------------------------------------------
% Part 1.2) Calibration of lambdaL to target ratio btw avg high-income
%           wealth to avg low-income wealth of 2
% -------------------------------------------------------------------------
lambdaL0   = param.lambda(1); 
rho        = param.rho; 
r          = param.r; 
target_vec = [3, 2, 0.02, 3]; 
weight_vec = [0, 1, 0, 0]; 

% Define loss
loss       = @(lambdaL) lossfunc(target_vec, weight_vec, rho, lambdaL, r, param, num, grids);  

% Minimize loss
tic; 
[lambda_star, minloss] = fminsearch(loss,lambdaL0); 
toc; 

fprintf(['lambdaL is calibrated as %f to target the ratio between average ' ...
         'high-income wealth to average low-income wealth of 2 \n' ...
         'Minimized squared difference between simulated and empirical moments is %f \n'], ...
         [lambda_star,minloss]);

% Code output: lambdaL is calibrated as 0.351660 to target the ratio between average high-income wealth to average low-income wealth of 2 
%              Minimized squared difference between simulated and empirical moments is 0.000000 

% -------------------------------------------------------------------------
% Part 1.3) Calibration of lambdaL to target variance of wealth of 0.02
% -------------------------------------------------------------------------
lambdaL0   = param.lambda(1); 
rho        = param.rho; 
r          = param.r; 
target_vec = [3, 2, 0.02, 3]; 
weight_vec = [0, 0, 1, 0]; 

% Define loss
loss       = @(lambdaL) lossfunc(target_vec, weight_vec, rho, lambdaL, r, param, num, grids);  

% Minimize loss
tic; 
[lambda_star, minloss] = fminsearch(loss,lambdaL0); 
toc; 

fprintf(['lambdaL is calibrated as %f to target variance of wealth of 0.02 \n' ...
         'Minimized squared difference between simulated and empirical moments is %f \n'], ...
         [lambda_star,minloss]);

% Code output: lambdaL is calibrated as 0.542480 to target variance of wealth of 0.02 
%              Minimized squared difference between simulated and empirical moments is 0.000000 


% -------------------------------------------------------------------------
% Part 2.1) Calibration of lambdaL to target avg length of unemployment of 3 months
% Note: - y_1 is the number of U->E incidences per month 
%       - y_1 follows possion dstrn with parameter lambdaL (job-finding rate)
%       - Let w be the waiting time until U->E. Then, it follows
%         exponential distribution with mean 1/lambdaL
%       - Hence, we should get lambdaL = 1/3 as our answer
% -------------------------------------------------------------------------
lambdaL0   = param.lambda(1); 
rho        = param.rho; 
r          = param.r; 
target_vec = [3, 2, 0.02, 3]; 
weight_vec = [0, 0, 0, 1]; 

% Define loss
loss       = @(lambdaL) lossfunc(target_vec, weight_vec, rho, lambdaL, r, param, num, grids);  

% Minimize loss
tic; 
[lambda_star, minloss] = fminsearch(loss,lambdaL0); 
toc; 

fprintf(['lambdaL is calibrated as %f to target avg length of unemployment of 3 months \n' ...
         'Minimized squared difference between simulated and empirical moments is %f \n'], ...
         [lambda_star,minloss]);

% Code output: lambdaL is calibrated as 0.333301 to target avg length of unemployment of 3 months 
%              Minimized squared difference between simulated and empirical moments is 0.000000 


% -------------------------------------------------------------------------
% Part 2.2) Aggregate savings as function of interest rate
% -------------------------------------------------------------------------
lambdaL   = param.lambda(1); 
rho       = param.rho; 
r_vec     = 0:0.001:0.04; 
asset_D_vec = NaN(length(r_vec),1); 

for ii = 1:length(r_vec)
    r  = r_vec(ii);
    asset_D_vec(ii,1) = asset_demand(rho, lambdaL, r, param, num, grids);
end

% Plot asset demand function
set(groot, 'DefaultAxesFontSize', 14);
set(groot, 'DefaultAxesLabelFontSizeMultiplier', 1.1); 
figure()
hold on
plot(asset_D_vec,r_vec,'LineWidth',3,'Color',[178/255,34/255,34/255])
grid on
ylabel('Interest Rate')
xlabel('Asset Demand')
ylim([r_vec(1), r_vec(end)])
title('Aggregate Asset Demand Function', 'FontSize', 17)
hold off
saveas(gcf,'asset_demand.png')

% see asset_demand.png 
