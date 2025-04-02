clc
data1 = load('MILP_MILPDEG1.mat'); % Load Scenario 1 (MILP)
data2 = load('MILP_MILPDEG2.mat'); % Load Scenario 2 (MILP+DEG)
% Extract cost data from loaded structures
Pfc_MILP = data1.P_FC_opt_1; % MILP Scenario
Pfc_MILP_DEG = data2.P_FC_opt_2; % MILP+DEG Scenario

% for a specific week between start and finish
startS=9*26*24; % August
finishS=startS+24*7+1;
startW=13*26*24; % December
finishW=startW+24*7+1;

%% SUMMER
year=365*2023+126;
start_dateS=startS/24+year;
finish_dateS=finishS/24+year;
ttS=datetime(start_dateS:(1/24):finish_dateS, 'ConvertFrom', 'datenum');

% Colors for differentiation
colorMILP = [0 0.447 0.741]; % Blue
colorMILP_DEG = [0.85 0.325 0.098]; % Red

figure('Position', [100, 100,600, 400]);
h1=area(ttS,Pfc_MILP(startS:finishS), 'FaceColor', colorMILP, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on
plot(ttS,Pfc_MILP(startS:finishS),'LineWidth',2,'Color',colorMILP);
hold on
h2=area(ttS,Pfc_MILP_DEG(startS:finishS), 'FaceColor', colorMILP_DEG, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on
plot(ttS,Pfc_MILP_DEG(startS:finishS),'LineWidth',2,'Color',colorMILP_DEG);

xlim([min(ttS) max(ttS)])
ylim ([0 20])
ylabel ('Power [kW]','fontweight','bold','FontName','Times New Roman');
legend([h1, h2],'MILP','MILP+DEG','location','eastoutside','FontName','Times New Roman')

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0, 0, 7, 4.5]);
set(gca, 'FontSize', font-8,'FontName','Times New Roman');
set(findall(gcf, 'Type', 'line'), 'LineWidth', 1.2); 
set(gca, 'LineWidth', 0.7, 'Box', 'on', 'XColor', 'k', 'YColor', 'k'); 


%% WINTER
start_dateW=startW/24+year;
finish_dateW=finishW/24+year;
ttW=datetime(start_dateW:(1/24):finish_dateW, 'ConvertFrom', 'datenum');

figure('Position', [100, 100,600, 400]);
h3=area(ttW,Pfc_MILP(startW:finishW), 'FaceColor', colorMILP, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on
plot(ttW,Pfc_MILP(startW:finishW),'LineWidth',2,'Color',colorMILP);
hold on
h4=area(ttW,Pfc_MILP_DEG(startW:finishW), 'FaceColor', colorMILP_DEG, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on
plot(ttW,Pfc_MILP_DEG(startW:finishW),'LineWidth',2,'Color',colorMILP_DEG);

xlim([min(ttW) max(ttW)])
ylim ([0 20])
ylabel ('Power [kW]','fontweight','bold','FontName','Times New Roman');
legend([h3, h4],'MILP','MILP+DEG','location','eastoutside','FontName','Times New Roman')

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0, 0, 7, 4.5]);
set(gca, 'FontSize', font-8,'FontName','Times New Roman');
set(findall(gcf, 'Type', 'line'), 'LineWidth', 1.2); 
set(gca, 'LineWidth', 0.7, 'Box', 'on', 'XColor', 'k', 'YColor', 'k'); 


