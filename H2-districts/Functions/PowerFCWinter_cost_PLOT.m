function [] = PowerFCWinter_cost_PLOT(linew,font,numCosts,Time,start,finish,P_FC_opt_matrix,S_FC_matrix,h2_cost)
%   this function enables to plot the trend of the power and the optimal
%   size of the Fuel Cell with respect the different costs of H2 over a
%   selected week 

year=365*2023+126;
start_date=start/24+year;
finish_date=finish/24+year;
tt=datetime(start_date:(1/24):finish_date, 'ConvertFrom', 'datenum');


% Define colors for each hydrogen price plot
cmap = crameri('batlow', numCosts);
colors = cmap;  
% Define a dark red color for size points
sizeColor = [0.5, 0, 0.1]; 

%% Plot with 2 y-axes
figure;
hold on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', font, 'LineWidth', 1.2); 
% Create primary y-axis for FC power
yyaxis left;
legendEntries = [];
for i = 1:numCosts
    areaPatch = fill([tt, fliplr(tt)], [zeros(1, length(tt)), fliplr(P_FC_opt_matrix(start:finish, i)')], ...
        colors(i, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
    plot(tt, P_FC_opt_matrix(start:finish, i), 'Color', colors(i,:), ...
        'LineWidth', 1.2, 'LineStyle', '-');  

    legendEntries = [legendEntries, areaPatch]; 
end
ylabel('Power [kW]', 'FontSize', font, 'FontWeight', 'bold');
set(gca, 'YColor', 'k'); 
xlim([min(tt) max(tt)]);

% Create secondary y-axis for S_FC (Fuel Cell Size)
yyaxis right;

numPoints = numCosts;
x_positions = linspace(tt(1), tt(end), numPoints); 
x_positions = flip(x_positions); 

% Plot points and connecting line
sizeLegendEntry = plot(x_positions, S_FC_matrix(start, :), '-o', 'Color', sizeColor, ...
    'LineWidth', 1, 'MarkerSize', 6, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', sizeColor); % Empty circles with a line

ylabel('Size [kW]', 'FontSize', font, 'FontWeight', 'bold','Color', sizeColor);
set(gca, 'YColor', sizeColor); 

borderWidth = 1;
ax = gca;
set(ax, 'Box', 'on', 'XColor', 'k','LineWidth', 0.7);
% LEGEND
legendLabels = arrayfun(@(x) sprintf('c_{H2} = %.0f CHF/kg', x), h2_cost, 'UniformOutput', false);
legendLabels{end+1} = 'S_{FC}'; 

legend([legendEntries, sizeLegendEntry], legendLabels, 'Location', 'eastoutside');

set(gca, 'FontSize', font - 4); 
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0, 0, 9 , 4.5]);
set(findall(gcf, 'Type', 'line'), 'LineWidth', 1.2); 
hold off;

end