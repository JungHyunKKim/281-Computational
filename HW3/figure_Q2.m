% -------------------------------------------------------------------------
% Plot figures for Question 2
% -------------------------------------------------------------------------

% Set default font size for x-axis and y-axis tick labels
set(groot, 'DefaultAxesFontSize', 14);

% Set default font size for x-label and y-label
% Multiplies the axes font size for the label font size
set(groot, 'DefaultAxesLabelFontSizeMultiplier', 1.1); 

% Value function
figure()

subplot(1,2,1)
hold on
plot(grids.a,v_implicit,'LineWidth',3,'Color',[178/255,34/255,34/255])
plot(grids.a,v_discrete,'--','LineWidth',3,'Color',[0,0,0,1])
grid on
xlabel('a')
ylabel('V(a)')
legend('Continuous','Discrete','Location','northwest','FontSize',12); 
xlim([num.a_min num.a_max])
title('Value Function', 'FontSize', 17)
hold off

% Consumption policy function
subplot(1,2,2)
hold on
plot(grids.a,c_implicit,'LineWidth',3,'Color',[178/255,34/255,34/255])
plot(grids.a,c_discrete,'--','LineWidth',3,'Color',[0,0,0,1])
grid on
xlabel('a')
ylabel('c(a)')
legend('Continuous','Discrete','Location','northwest','FontSize',12); 
xlim([num.a_min num.a_max])
title('Consumption Policy Function', 'FontSize', 17)
hold off

set(gcf,'Position',[200 0 1000 400]) 
saveas(gcf,'Q2.png')