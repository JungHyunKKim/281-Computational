% -------------------------------------------------------------------------
% Plot figures for Question 4
% -------------------------------------------------------------------------

% Set default font size for x-axis and y-axis tick labels
set(groot, 'DefaultAxesFontSize', 14);

% Set default font size for x-label and y-label
% Multiplies the axes font size for the label font size
set(groot, 'DefaultAxesLabelFontSizeMultiplier', 1.1); 

figure()

% Choice of Technology
subplot(2,3,1)
hold on
plot(grids.a,good_tech_mat(:,1),'LineWidth',3,'Color',[178/255,34/255,34/255])
plot(grids.a,good_tech_mat(:,2),'--','LineWidth',3,'Color',[0,0,0,1])
plot(grids.a,good_tech_mat(:,3), ':','LineWidth',3,'Color',[0.9290,0.6940,0.1250,1])
grid on
xlabel('a')
ylabel('1\{Good Technology\}')
ylim([-0.75,1.75])
legend('\kappa = 0','\kappa = 0.5','\kappa = 2','Location','northwest','FontSize',12, 'interpreter', 'tex'); 
xlim([num.a_min num.a_max])
title('Choice of technology', 'FontSize', 17)
hold off

% Capital Demand
subplot(2,3,2)
hold on
plot(grids.a,kstar_mat(:,1),'LineWidth',3,'Color',[178/255,34/255,34/255])
plot(grids.a,kstar_mat(:,2),'--','LineWidth',3,'Color',[0,0,0,1])
plot(grids.a,kstar_mat(:,3), ':','LineWidth',3,'Color',[0.9290,0.6940,0.1250,1])
grid on
xlabel('a')
ylabel('k(a)')
legend('\kappa = 0','\kappa = 0.5','\kappa = 2','Location','northwest','FontSize',12, 'interpreter', 'tex'); 
xlim([num.a_min num.a_max])
title('Capital Demand', 'FontSize', 17)
hold off

% Maximized Profit
subplot(2,3,3)
hold on
plot(grids.a,pistar_mat(:,1),'LineWidth',3,'Color',[178/255,34/255,34/255])
plot(grids.a,pistar_mat(:,2),'--','LineWidth',3,'Color',[0,0,0,1])
plot(grids.a,pistar_mat(:,3), ':','LineWidth',3,'Color',[0.9290,0.6940,0.1250,1])
grid on
xlabel('a')
ylabel('\pi(a)')
legend('\kappa = 0','\kappa = 0.5','\kappa = 2','Location','northwest','FontSize',12, 'interpreter', 'tex'); 
xlim([num.a_min num.a_max])
title('Maximized Profit', 'FontSize', 17)
hold off

% Value function
subplot(2,3,4)
hold on
plot(grids.a,v_mat(:,1),'LineWidth',3,'Color',[178/255,34/255,34/255])
plot(grids.a,v_mat(:,2),'--','LineWidth',3,'Color',[0,0,0,1])
plot(grids.a,v_mat(:,3), ':','LineWidth',3,'Color',[0.9290,0.6940,0.1250,1])
grid on
xlabel('a')
ylabel('V(a)')
legend('\kappa = 0','\kappa = 0.5','\kappa = 2','Location','southeast','FontSize',12, 'interpreter', 'tex'); 
xlim([num.a_min num.a_max])
title('Value Function', 'FontSize', 17)
hold off

% Instantaneous Saving Rate
subplot(2,3,5)
hold on
plot(grids.a,adot_mat(:,1),'LineWidth',3,'Color',[178/255,34/255,34/255])
plot(grids.a,adot_mat(:,2),'--','LineWidth',3,'Color',[0,0,0,1])
plot(grids.a,adot_mat(:,3),':','LineWidth',3,'Color',[0.9290,0.6940,0.1250,1])
grid on
xlabel('a')
ylabel('$\dot{a}$', 'Interpreter', 'latex', 'FontSize',20)
legend('\kappa = 0','\kappa = 0.5','\kappa = 2','Location','southeast','FontSize',12, 'interpreter', 'tex'); 
xlim([num.a_min num.a_max])
title('Instantaneous Saving Rate', 'FontSize', 17)
hold off

% Consumption policy function
subplot(2,3,6)
hold on
plot(grids.a,c_mat(:,1),'LineWidth',3,'Color',[178/255,34/255,34/255])
plot(grids.a,c_mat(:,2),'--','LineWidth',3,'Color',[0,0,0,1])
plot(grids.a,c_mat(:,3),':','LineWidth',3,'Color',[0.9290,0.6940,0.1250,1])
grid on
xlabel('a')
ylabel('c(a)')
legend('\kappa = 0','\kappa = 0.5','\kappa = 2','Location','northwest','FontSize',12, 'interpreter', 'tex'); 
xlim([num.a_min num.a_max])
title('Consumption Policy Function', 'FontSize', 17)
hold off

set(gcf,'Position',[200 0 1800 900]) 
saveas(gcf,'Q4.png')