clear all
clc

% Load cost data 
data1 = load('scenariocost1.mat'); % Load Scenario 1 (MILP)
data2 = load('scenariocost2.mat'); % Load Scenario 2 (MILP+DEG)

% Extract cost data from loaded structures
costs_MILP = [data1.cost_inst_1, data1.cost_maint_1, data1.cost_imp_1, -data1.cost_exp_1]; % MILP ScenariO
costs_MILP_DEG = [data2.cost_inst_2, data2.cost_maint_2, data2.cost_imp_2, -data2.cost_exp_2]; % MILP+DEG Scenario

costs = [costs_MILP; costs_MILP_DEG]';


% Define cost category labels (Ensure all five categories are included)
cost_labels = {'C_{inst}', 'C_{maint}', 'C_{imp}', 'C_{exp}'};

% Define x-axis labels for the two scenarios
x_labels = {'MILP', 'MILP+DEG'};

% Create the stacked bar plot
figure;
b = bar(costs', 'stacked','BarWidth', 0.5);
num_categories = size(costs, 1); 
cmap = crameri('batlow', num_categories); 


for i = 1:num_categories
    b(i).FaceColor = cmap(i, :); 
    b(i).FaceAlpha = 0.4;
end

set(gca, 'XTick', 1:2); % Due tick per MILP e MILP+DEG
set(gca, 'XTickLabel', x_labels); 
set(gca, 'FontSize', 12);
set(gca, 'FontName', 'Times New Roman','FontWeight','bold'); 
ylabel('Cost [kCHF/year]', 'FontName', 'Times New Roman');

xtickangle(0); 
xlim([0.5, 2.5]); 


% Calcolare la somma delle barre (altezza totale per il posizionamento dei punti)
y_positions = [data1.cost_tot_1, data2.cost_tot_2];
x_positions = [1, 2]; 

% Aggiungere i punti sopra le colonne
burgundy = [128/255, 0, 32/255]; 
hold on;
plot(x_positions, y_positions, '-o', 'Color', burgundy, 'LineWidth', 1.5); 
% scatter(x_positions, y_positions, 100, burgundy, 'filled', 'MarkerEdgeColor', burgundy); 
hold off;
legendHandle = legend([cost_labels,'C_{tot}'], 'Location', 'northeastoutside');
set(gca, 'LineWidth', 0.7, 'Box', 'on'); 
