% -------------------------------------------------------------------------
% 281-Computational Homework 5
% Jung Hyun Kim
%
% Part 0. Solve for stationary equilibrium under baseline parameters
% Part 1. Tracing supply and demand using instruments
% Part 2. Homework: Compare income shocks with homogeneous vs.heterogeneous returns
%
% Important functions  :
%      solve_Aiyagari.m  - Solve Aiyagari model by minimizing squared
%                          excess capital demand
%                        - Need to specify eqtype: 'GE' or 'PE'
%                             
%      solve_Kd_excess.m - Compute excess capital demand given r and eqtype
%
%      solve_HH.m        - Compute asset demand (capital supply) given r and w
%                          Solve HJB using VFI using implicit method
%                          Solve law of motion for distribution using kf_equation.m
%      solve_Kd.m        - Compute capital demand given r and aggregate Ld
%                        - FOC conditions imposed
%
% Functions from class : create_grids.m, kf_equation.m, utility.m, 
%                        vfi_iteration.m, vp_upwind.m
%
% Note: There is no need to solve to clear the labor market due to CRS
%       technology (can analytically derive wage as a function of the
%       rental rate of capital) 
% -------------------------------------------------------------------------

clear; 
close all; 
clc;

% Create structure with structural parameters
call_parameters;

% Create structure with numerical parameters
numerical_parameters;

% Create structure with grids and initial guesses for v
grids = create_grids(par, num);

% Default figure parameters
set(groot, 'DefaultAxesFontSize', 14);
set(groot, 'DefaultAxesLabelFontSizeMultiplier', 1.1); 

% -------------------------------------------------------------------------
% Part 0. Solve for Stationary Equilibrium
%
% (1) Get capital demand and supply functions to plot them
% (2) Solve for stationary equilibrium to check if it is at the
%     intersection of the plotted capital demand and supply curves
% -------------------------------------------------------------------------

% ----------------------------------------
% (1) Get capital demand and supply functions

r_vec  = 0:0.001:par.rho; 
Ks_vec = NaN(length(r_vec),1); 
Kd_vec = NaN(length(r_vec),1); 

for ii = 1:length(r_vec)

    r  = r_vec(ii);
    w  = (1 - par.alpha) * par.z^( 1 / (1-par.alpha) ) ...
          *( par.alpha / (r + par.delta) )^( par.alpha/(1-par.alpha) ); 
    
    [~, ~, Ad, Ls, g] = solve_HH(r, w, par, num, grids); 
    Ks_vec(ii,1) = Ad;
    Kd_vec(ii,1) = solve_Kd(r, par, Ls); 

end

% ----------------------------------------
% (2) Solve for Stationary Equilibrium

r0 = 0.025; %par.rho - 0.005; 

[r_seq, w_seq, K_seq, g_seq] = solve_Aiyagari(par, num, grids, 'GE', r0); 

% ----------------------------------------

% Plot capital demand and supply functions 
figure()
hold on
plot(Ks_vec,r_vec,'LineWidth',3,'Color',[178/255,34/255,34/255])
plot(Kd_vec,r_vec,'--','LineWidth',3,'Color',[0,0.4470,0.7410,1])
plot(K_seq, r_seq, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'LineWidth',2);
legend('Capital Supply', 'Capital Demand', 'Stationary Eq.', 'Location', 'southeast')
grid on
ylabel('Interest Rate')
xlabel('Capital')
ylim([r_vec(1), r_vec(end)])
title('Aggregate Capital', 'FontSize', 17)
set(gcf,'Position',[200 0 600 400]) 
saveas(gcf,'figures/stationary_capital.png')

% Plot stationary distribution
figure()
hold on
plot(grids.a(:,1), g_seq(:,1)/sum(g_seq(:,1)),'LineWidth',3,'Color',[178/255,34/255,34/255])
plot(grids.a(:,1), g_seq(:,2)/sum(g_seq(:,2)),'--','LineWidth',3,'Color',[0,0.4470,0.7410,1])
legend('Low Income', 'High Income')
hold off
grid on
ylabel('Density')
xlabel('Asset')
title('Stationary Distribution Conditional on Income', 'FontSize', 17)
set(gcf,'Position',[200 0 600 400]) 
saveas(gcf,'figures/stationary_distribution.png')

% -------------------------------------------------------------------------
% Part 1.1. Trace Capital Supply Curve Using Different Levels of TFP
%           -> Shifting capital demand curve
% (i)  Let wages adjust (capital supply in general eq)
% (ii) Fix wage to initial stationary eq (capital supply in partial eq)
% -------------------------------------------------------------------------

% Fix w0 to equilibrium wage (used for PE)
par.w0 = w_seq; 

% Save baseline z
z_baseline = par.z; 

% Create vector of TFP's
z_vec = 0.1:0.02:1; 

% Initalize vectors of interest rate and capital (for GE and PE)
zloop.r_vec_ge = NaN(length(z_vec),1); 
zloop.K_vec_ge = NaN(length(z_vec),1); 
zloop.r_vec_pe = NaN(length(z_vec),1); 
zloop.K_vec_pe = NaN(length(z_vec),1); 

% Initialize interest rate
r0_ge = r_seq; 
r0_pe = r_seq; 

% Loop over z_vec to trace capital supply curve
for ii = 1:length(z_vec)

    par.z = z_vec(ii); 
    
    % Solve GE
    [rstar_ge, ~, Kstar_ge, ~] = solve_Aiyagari(par, num, grids, 'GE', r0_ge); 
    zloop.r_vec_ge(ii)      = rstar_ge; 
    zloop.K_vec_ge(ii)      = Kstar_ge; 

    r0_ge = rstar_ge; 
    
    % Solve PE
    [rstar_pe, ~, Kstar_pe, ~] = solve_Aiyagari(par, num, grids, 'PE', r0_pe); 
    zloop.r_vec_pe(ii)      = rstar_pe; 
    zloop.K_vec_pe(ii)      = Kstar_pe; 

    r0_pe = rstar_pe; 

end

% Set z to baseline
par.z   = z_baseline; 

% -------------------------------------------------------------------------
% Part 1.2. Trace Capital Demand Curve Using Different Levels of Discount Rates
%           -> Shifting capital supply curve
% Note: Wage is not a function of the discount rate, 
%       so PE captial demand curve = GE capital demand curve
% -------------------------------------------------------------------------

% Fix w0 to equilibrium wage (used for PE)
par.w0 = w_seq; 

% Save baseline rho
rho_baseline = par.rho; 

% Create vector of TFP's
rho_vec = 0.02:0.02:0.1; 

% Initalize vectors of interest rate and capital
rholoop.r_vec_ge = NaN(length(rho_vec),1); 
rholoop.K_vec_ge = NaN(length(rho_vec),1); 

% Initialize interest rate
r0_ge = r_seq; 

% Loop over rho_vec to trace capital demand curve
for ii = 1:length(rho_vec)

    par.rho = rho_vec(ii); 

    [rstar_ge, ~, Kstar_ge, ~] = solve_Aiyagari(par, num, grids, 'GE', r0_ge); 
    rholoop.r_vec_ge(ii)    = rstar_ge; 
    rholoop.K_vec_ge(ii)    = Kstar_ge; 

    r0_ge = rstar_ge; 

end

% Set rho to baseline
par.rho = rho_baseline; 

% -------------------------------------------------------------------------
% Part 1.3. Plot the traced capital demand and supply curves together
% -------------------------------------------------------------------------

figure()
hold on
plot(rholoop.K_vec_ge, rholoop.r_vec_ge,'LineWidth',3,'Color',[0.9290,0.6940,0.1250,1])
plot(zloop.K_vec_pe, zloop.r_vec_pe,'LineWidth',3,'Color',[178/255,34/255,34/255])
plot(zloop.K_vec_ge, zloop.r_vec_ge,'LineWidth',3,'Color',[0,0.4470,0.7410])
grid on
ylabel('Interest Rate')
xlabel('Capital')
title('Aggregate Capital', 'FontSize', 17)
legend('Capital Demand','Capital Supply (PE)','Capital Supply (GE)','Location','northeast','FontSize',14)
hold off
saveas(gcf,'figures/capital_curves.png')

% -------------------------------------------------------------------------
% Part 2. Idiosyncratic Income Shocks with Heterogeneous Returns
% -------------------------------------------------------------------------

% Specify heterogeneous returns by income shocks
par.rshock = [0.95 1.05]; 

% Get capital demand and supply functions
r_vec          = 0:0.001:par.rho; 
Ks_vec_heteroR = NaN(length(r_vec),1); 
Kd_vec_heteroR = NaN(length(r_vec),1); 

for ii = 1:length(r_vec)

    r  = r_vec(ii);
    w  = (1 - par.alpha) * par.z^( 1 / (1-par.alpha) ) ...
          *( par.alpha / (r + par.delta) )^( par.alpha/(1-par.alpha) ); 
    
    [~, ~, Ad, Ls, g]     = solve_HH(r, w, par, num, grids); 
    Ks_vec_heteroR(ii,1) = Ad;
    Kd_vec_heteroR(ii,1) = solve_Kd(r, par, Ls); 
    
end

% Find equilibrium capital demand and supply
r0 = 0; 
[r_seq_heteroR, w_seq_heteroR, K_seq_heteroR, g_heteroR] = solve_Aiyagari(par, num, grids, 'GE', r0); 

% Plot capital demand and supply functions and the stationary equilibrium

figure()
sgtitle('Income Shock with Homogeneous vs. Heterogeneous Returns', 'FontSize', 17, 'FontWeight', 'bold');
set(gcf,'Position',[200 0 2000 400]) 

subplot(1,3,1)
hold on
plot(Ks_vec,r_vec,'--','LineWidth',3,'Color',[178/255,34/255,34/255])
plot(Ks_vec_heteroR,r_vec,'LineWidth',3,'Color',[178/255,34/255,34/255])
plot(Kd_vec,r_vec,'--','LineWidth',3,'Color',[0,0.4470,0.7410,1])
plot(Kd_vec_heteroR,r_vec,'LineWidth',3,'Color',[0,0.4470,0.7410,1])
plot(K_seq, r_seq, 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'LineWidth',2);
plot(K_seq_heteroR, r_seq_heteroR, 'x', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'LineWidth',2);
lgd = legend('Ks (Baseline)', 'Ks (Hetero returns)', 'Kd (Baseline)' , 'Kd (Hetero returns)', ...
             'SE (Baseline)','SE (Hetero returns)', 'Location', 'northeast'); 
grid on
ylabel('Interest Rate')
xlabel('Capital')
ylim([r_vec(1), r_vec(end)])
title('Capital Curves', 'FontSize', 14,'FontWeight','Normal')
pause(0.01)
legendPos = lgd.Position; % Adjust the legend position slightly downwards
legendPos(2) = legendPos(2) - 0.03;  % Adjust this value to
lgd.Position = legendPos;

subplot(1,3,2)
hold on
plot(grids.a(:,1), g_seq(:,1)/sum(g_seq(:,1)),'--','LineWidth',3,'Color',[178/255,34/255,34/255])
plot(grids.a(:,1), g_heteroR(:,1)/sum(g_seq(:,1)),'LineWidth',3,'Color',[0,0.4470,0.7410,1])
legend('Low Income (Baseline)', 'Low Income (Hetero returns)')
hold off
grid on
ylabel('Density')
xlabel('Asset')
title('Stationary Distribution (Low Income)', 'FontSize', 14 ,'FontWeight','Normal')

subplot(1,3,3)
hold on
plot(grids.a(:,1), g_seq(:,1)/sum(g_seq(:,1)),'--','LineWidth',3,'Color',[178/255,34/255,34/255])
plot(grids.a(:,1), g_heteroR(:,1)/sum(g_seq(:,1)),'LineWidth',3,'Color',[0,0.4470,0.7410,1])
legend('High Income (Baseline)', 'High Income (Hetero returns)')
hold off
grid on
ylabel('Density')
xlabel('Asset')
title('Stationary Distribution (High Income)', 'FontSize', 14 ,'FontWeight','Normal')

saveas(gcf,'figures/compare_shocktype.png')
