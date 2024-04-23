% -------------------------------------------------------------------------
% Plot figures for Question 1
% -------------------------------------------------------------------------

% Set default font size for x-axis and y-axis tick labels
set(groot, 'DefaultAxesFontSize', 14);

% Set default font size for x-label and y-label
% Multiplies the axes font size for the label font size
set(groot, 'DefaultAxesLabelFontSizeMultiplier', 1.1); 

% Value function
figure()

subplot(1,3,1)
hold on
plot(grids.a,v_implicit,'LineWidth',3,'Color',[178/255,34/255,34/255])
plot(grids.a,v_implicit_smallD,'--','LineWidth',3,'Color',[0,0,0,1])
grid on
xlabel('a')
ylabel('V(a)')
legend('Implicit (\Delta=0.01)','Implicit (\Delta=0.5)','Location','northwest','FontSize',12, 'Interpreter', 'tex'); 
xlim([num.a_min num.a_max])
title('Value Function', 'FontSize', 17)
hold off

% Instantaneous Saving Rate
subplot(1,3,2)
hold on
plot(grids.a,adot_implicit,'LineWidth',3,'Color',[178/255,34/255,34/255])
plot(grids.a,adot_implicit_smallD,'--','LineWidth',3,'Color',[0,0,0,1])
grid on
xlabel('a')
ylabel('$\dot{a}$', 'Interpreter', 'latex', 'FontSize',20)
legend('Implicit (\Delta=0.01)','Implicit (\Delta=0.5)','Location','northwest','FontSize',12, 'Interpreter', 'tex'); 
xlim([num.a_min num.a_max])
title('Instantaneous Saving Rate', 'FontSize', 17)
hold off

% Consumption policy function
subplot(1,3,3)
hold on
plot(grids.a,c_implicit,'LineWidth',3,'Color',[178/255,34/255,34/255])
plot(grids.a,c_implicit_smallD,'--','LineWidth',3,'Color',[0,0,0,1])
grid on
xlabel('a')
ylabel('c(a)')
legend('Implicit (\Delta=0.01)','Implicit (\Delta=0.5)','Location','northwest','FontSize',12, 'Interpreter', 'tex'); 
xlim([num.a_min num.a_max])
title('Consumption Policy Function', 'FontSize', 17)
hold off

set(gcf,'Position',[200 0 1600 400]) 
saveas(gcf,'Q1.png')