% -------------------------------------------------------------------------
% 281-Computational Homework 5
% Jung Hyun Kim
%
% Homework: (1) IRFs when log(TFP) = 0 
%           (2) IRFs when log(TFP) = -0.1
%
% - "solve_KS.m" solves Krusell and Smith (1998) 
% - See Homework 4 for description of the function "solve_Aiyagari.m"
% - Note: model unstable for cases with 0 income, higher gamma, and small a_max(<=7)
% -------------------------------------------------------------------------

clear; clc; close all; 

% -------------------------------------------------------------------------
% Setup the toolbox and set parameters
% -------------------------------------------------------------------------

% Just need to include folders containing the files in the path.
addpath('C:\Users\JungHyun\Documents\GitHub\MATLABAutoDiff');
addpath('C:\Users\JungHyun\Documents\GitHub\phact');

% Initialize shocks for simulation
T = 200;
N = 2000;
vAggregateShock = zeros(1,N);
vAggregateShock(1,1) = 1;

% Set parameters
numerical_parameters; 
call_parameters; 
grids = create_grids(par,num);

% Solve K&S at log(TFP)=0 
[TFP, Y, C, I,TFP_reduced, Y_reduced, C_reduced, I_reduced, vTime] ...
= solve_KS(par, num, grids, T, N, vAggregateShock); 

% Solve K&S at log(TFP)=-0.1
par.logTFP = -0.1; 

[TFP2, Y2, C2, I2, TFP_reduced2, Y_reduced2, C_reduced2, I_reduced2, vTime2] ...
= solve_KS(par, num, grids, T, N, vAggregateShock); 

% -------------------------------------------------------------------------
% Plot impulse response functions
% -------------------------------------------------------------------------

figure
subplot(2,2,1)
hold on
plot(vTime,100 * TFP,'linewidth',1.5)
plot(vTime,100 * TFP2,'linewidth',1.5)
set(gcf,'color','w')
xlim([vTime(1) vTime(end)])
title('TFP Level','interpreter','latex','fontsize',14)
xlabel('Quarters','interpreter','latex')
ylabel('$\%$','interpreter','latex')
legend('logTFP = 0', 'logTFP = -0.1','Location','southeast')
hold off

subplot(2,2,2)
hold on
plot(vTime,100 * Y_reduced,'linewidth',1.5)
plot(vTime,100 * Y_reduced2,'linewidth',1.5)
set(gcf,'color','w')
xlim([vTime(1) vTime(end)])
title('Output','interpreter','latex','fontsize',14)
xlabel('Quarters','interpreter','latex')
ylabel('$\%$ deviation from s.s.','interpreter','latex')
legend('logTFP = 0', 'logTFP = -0.1')
hold off

subplot(2,2,3)
hold on
plot(vTime,100 * C_reduced,'linewidth',1.5)
plot(vTime,100 * C_reduced2,'linewidth',1.5)
set(gcf,'color','w')
xlim([vTime(1) vTime(end)])
title('Consumption','interpreter','latex','fontsize',14)
ylabel('$\%$ deviation from s.s.','interpreter','latex')
xlabel('Quarters','interpreter','latex')
legend('logTFP = 0', 'logTFP = -0.1')
hold off

subplot(2,2,4)
hold on
plot(vTime,100 * I_reduced,'linewidth',1.5)
plot(vTime,100 * I_reduced2,'linewidth',1.5)
set(gcf,'color','w')
xlim([vTime(1) vTime(end)])
title('Investment','interpreter','latex','fontsize',14)
xlabel('Quarters','interpreter','latex')
legend('logTFP = 0', 'logTFP = -0.1')
hold off

set(gcf,'Position',[200 0 1000 800]) 
saveas(gcf,'IRF_TFP.png')