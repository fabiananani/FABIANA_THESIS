clear all
clc

load('plot_nanimodel.mat');

% Data for both models
values_efficiency = [annualEffReduction1, New_Lifetime1, increase_mflow1];
values_Taylor = [annualEffReduction2, New_Lifetime2, increase_mflow2];
% x_positions = [1, 4, 7];
% Define categories
labels = {'Efficiency reduction', 'Predicted Lifetime', 'H_{2} consumption'};

% Define colormap
cmap = crameri('batlow',6); 

% figure
% b = bar( [values_efficiency; values_Taylor]', 'FaceColor', 'flat', 'EdgeColor', 'k','FaceAlpha',0.3, 'LineWidth', 0.7,'BarWidth',1);
% b(1).CData = repmat(cmap(1, :), length(labels), 1); % Color for Original Model
% b(2).CData = repmat(cmap(4, :), length(labels), 1); % Color for Simplified Model
% 
% 
% % set(gca, 'XTick', x_positions); 
% set(gca, 'XTickLabel', labels, 'FontSize', 9, ...
%          'XTickLabelRotation', 0, 'FontName', 'Times New Roman','FontWeight', 'bold');
% 
% ylim([0 max([values_efficiency, values_Taylor]) * 1.2]); 
% box on; 
% set(gca, 'LineWidth',0.7);
% legend('Original Model', 'Simplified Model', 'FontSize', 9, 'Location', 'bestoutside','FontWeight', 'bold');

% rel error [%]
relative_error = abs(eff_deg_FC1 - eff_deg_FC2) ./ abs(eff_deg_FC1) * 100;

% subplots
figure;

% 1 subplot
subplot(2,1,1)
hold on;
plot(Time, eff_deg_FC1, '--','LineWidth', 2, 'Color',cmap(1,:),'DisplayName', 'Original Model');
plot(Time, eff_deg_FC2, 'LineWidth', 2, 'Color', cmap(5,:),'DisplayName', 'Simplified Model');
fill([Time; flipud(Time)], [eff_deg_FC1; flipud(eff_deg_FC2)], [0.8, 0.8, 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Area

xlabel('Time [h]', 'FontSize', 120, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
ylabel('Efficiency [-]', 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
box on;
legend('Case 1', 'Case 2', 'FontSize', 10,'FontWeight', 'bold');
set(gca, 'FontSize', 10, 'LineWidth', 0.7,'FontName', 'Times New Roman');

% 2subplot
subplot(2,1,2)
plot(Time, relative_error, 'k', 'LineWidth', 2.5);
xlabel('Time [h]', 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
ylabel('Relative Error [%]', 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
box on;
set(gca, 'FontSize', 10, 'LineWidth', 0.7,'FontName', 'Times New Roman');

set(gcf, 'Position', [100, 100, 800, 600]); 


