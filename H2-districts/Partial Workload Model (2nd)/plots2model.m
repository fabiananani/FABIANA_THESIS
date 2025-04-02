clear all
clc

load('plot_shehzadmodel.mat');

% Data for both models
values_original = [annualEffReduction, New_Lifetime, increase_mflowH2];
values_simplified = [annualEffReduction_simp, New_Lifetime_simp, increase_mflowH2_simp];
% x_positions = [1, 4, 7];
% Define categories
labels = {'Efficiency reduction', 'Predicted Lifetime', 'H_{2} consumption'};

% Define colormap
cmap = crameri('batlow', 5); 

figure
b = bar( [values_original; values_simplified]', 'FaceColor', 'flat', 'EdgeColor', 'k','FaceAlpha',0.3, 'LineWidth', 0.7,'BarWidth',1);
b(1).CData = repmat(cmap(1, :), length(labels), 1); % Color for Original Model
b(2).CData = repmat(cmap(4, :), length(labels), 1); % Color for Simplified Model


% set(gca, 'XTick', x_positions); 
set(gca, 'XTickLabel', labels, 'FontSize', 9, ...
         'XTickLabelRotation', 0, 'FontName', 'Times New Roman','FontWeight', 'bold');

ylim([0 max([values_original, values_simplified]) * 1.2]); 
box on; 
set(gca, 'LineWidth',0.7);
legend('Original Model', 'Simplified Model', 'FontSize', 9, 'Location', 'bestoutside','FontWeight', 'bold');

% rel error [%]
relative_error = abs(eff_FC_deg1 - eff_FC_deg1_simp) ./ abs(eff_FC_deg1) * 100;

% subplots
figure;

% 1 subplot
subplot(2,1,1)
hold on;
plot(Time, eff_FC_deg1, '-.','LineWidth', 2.5, 'Color',cmap(3,:),'DisplayName', 'Original Model');
plot(Time, eff_FC_deg1_simp, 'LineWidth', 2.5, 'Color', cmap(1,:),'DisplayName', 'Simplified Model');
fill([Time; flipud(Time)], [eff_FC_deg1; flipud(eff_FC_deg1_simp)], [0.8, 0.8, 0.8], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Area

xlabel('Time [h]', 'FontSize', 120, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
ylabel('Efficiency [-]', 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
box on;
legend('Original Model', 'Simplified Model', 'FontSize', 10);
set(gca, 'FontSize', 12, 'LineWidth', 0.7,'FontName', 'Times New Roman');

% 2subplot
subplot(2,1,2)
plot(Time, relative_error, 'k', 'LineWidth', 2.5);
xlabel('Time [h]', 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
ylabel('Relative Error [%]', 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
box on;
set(gca, 'FontSize', 12, 'LineWidth', 0.7,'FontName', 'Times New Roman');

set(gcf, 'Position', [100, 100, 800, 600]); 

%% 

% absolute error
Error_abs = abs(eff_FC_deg1 - eff_FC_deg1_simp);

% Compute cumulative error over time
Cumulative_Error = cumsum(Error_abs); % Sum of absolute errors over time


% Plot Cumulative Error
figure;
plot(Time, Cumulative_Error, 'LineWidth', 2);
xlabel('Time [h]');
ylabel('Cumulative Error');
title('Cumulative Error Over Time');
grid on;
