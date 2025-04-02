clear all
close all
clc

%% Set up Gurobi: Adding path to optimizer

addpath('C:\gurobi1103\win64\examples\matlab')
gurobi_setup
%% Input data

% Define the main path manually
    path_main = 'C:\Users\nanf\Desktop\H2-districts';
% path for input files
    path_input = fullfile(path_main, 'Input');
% Define the path for functions
    path_functions = fullfile(path_main, 'Functions');

% Reading an Excel file from the Input directory
    excelFilePath1 = fullfile(path_input, 'Inputs_BoundaryLoads_2022.xlsx');
    excelFilePath2 = fullfile(path_input, 'nestel_demand.csv');
    excelFilePath3 = fullfile(path_input,'PriceImpExp_ele.xlsx');
    Input1         = readtable(excelFilePath1);
    Input2         = readtable(excelFilePath2);
    Input3         = readtable(excelFilePath3);
    
%% Pre-processing

    irradiance     = Input1.G;                     % Hourly solar irradiance    [W/m^2]
    % T_amb          = Input1.T_amb;                 % Hourly temperature         [°C]
    nHours         = 8760;                        % Number of hours simulated (8760h)    
    idxHr2ToEnd    = (2:nHours)';                 % Hours until the end
    Time           = (1:nHours)';                 % Time vector
    days           = nHours/24;                   % Number of days simulated
    weeks          = days/7;                      % Number of weeks simulated
    clear t; 
    linew          = 1.5;             
    font           = 18;
    
%% INPUT PARAMETERS

 % general inputs

    bigM            = 600;      % Large number for big-M constraints
    MaxSimTime      = 900;     % Maximum time for MILP solver                 [s]      
    deltat          = 3600;    % time step di 1h                              [s]       

 % Electrical district demand
        
    P_load         = Input2.P_el_nest_h;       % NEST electrical demand      [kW]

 % Fixed technical parameters  

    HHV         = 141.9*10^3;           % (141.9 MJ/kg) Hydrogen higher heating value [kJ/kg]  
    LHV         = 120.1*10^3;           % (120.1 MJ/kg) Hydrogen lower heating value [kJ/kg]


    
 % Flag for electricity cost (true/false)
   
    use_Empa_prices = false; % change with "false" to use DSO prices
    
    if use_Empa_prices

        c_gridimp   = Input3.PriceimpEmpa;                     % [CHF/kWh]
        c_gridexp   = Input3.PriceexpEmpa;                     % [CHF/kWh]
    else 
                c_gridimp   = Input3.Priceimpgreen;                  % [CHF/kWh]
                c_gridexp   = Input3.PriceexpPV;                     % [CHF/kWh]
    end 

    % Efficencies and max lifetime from datasheet

    eff_ch       = 0.95;                                  % charging efficiency of the battery 
    eff_disch    = 0.95;                                  % discharging efficency of the battery
    eff_PV       = 0.17;                                  % constant efficiency of the PV system
    eff_FC       = 0.51;                                  % constant electrical efficency of FC  
    max_LT       = 40000;                                 % [h]

% H2 cost
    c_h2             = 3; 
% % Sensitivity analysis for H2 prices
%     c_h2            = [2 3 6 10];                          
%     numCosts        = numel(c_h2);
%     P_FC_matrix = zeros(nHours, numCosts);                 
%     S_FC_matrix = zeros(nHours, numCosts);
%     A_PV_matrix =  zeros(nHours, numCosts);
%     C_b_matrix  =  zeros(nHours, numCosts);
%     c_tot_matrix  =  zeros(nHours, numCosts);

% Fuel Cell (FC)
    S_FC_min     = 0;                            % Minimum size FC             [kW]
    S_FC_max     = 60;                           % Maximum size FC             [kW]

% PV   
    P_PV_peak_ref  = 100;                 % Reference value peak power PV      [kW]
    Area_PV_ref = P_PV_peak_ref/eff_PV;   % Reference area PV                  [m2]
    Area_PV_min = 0;                      % Minimum area PV                    [m2]
    Area_PV_max = Area_PV_ref*2;          % Maximum area PV                    [m2]

% Battery

    C_b_max     = 96 * 3.6 * 10^3;     % Battery capacity                    [kJ]
    C_b_min     = C_b_max / 3;         % MIN Battery capacity                [kJ]
    P_b_max     = C_b_max / 3600;      % Battery capacity (1C rate)          [kWh]

%% unit prices and lifetime components

    d           = 0.05;                     % Discount rate
    ann         = d / (1 - (1 + d)^(-20));  % annuity factor calculated with plant lifetime        
   
    UP_PV       = 1300;               % Unit price PV                                    [CHF/kW_p]
    life_PV     = 25;                 % Lifetime PV                                      [years]
    maint_PV    = 0.01;               % Annual cost maintenance PV, frac UP_PV
    ann_PV      = d / (1 - (1 + d)^(-life_PV));
   
    UP_b        = 1000;          % Unit price battery [CHF/kWh]
    life_b      = 12;%10;        % Lifetime battery   [years]
    maint_b     = 0.02;          % frac UP_b 
    ann_b       = d / (1 - (1 + d)^(-life_b));
    
    UP_FC        = 950;          % Unit price FC [CHF/kW]
    life_FC      = 10; %7.5         % Lifetime FC   [years]
    maint_FC     = 0.024;        % Annual cost maintenance FC, frac UP_FC
    ann_FC       = d / (1 - (1 + d)^(-life_FC));
    
    % %MILP
    % for i = 1:numCosts
    % current_c_h2 = c_h2(i);
        
%% Define the optimization problem and the optimization variables

    sizingprob = optimproblem;

%% DECISION VARIABLES

% SIZING decision variables

    % Battery capacity in [kJ]
    C_b            = optimvar('C_b','LowerBound',C_b_min,'UpperBound',C_b_max);  % C_b actual installed capacity,C_b_max= max capacity we can have 
    % Fuel Cell size in [kW]
    S_FC           = optimvar('S_FC','LowerBound',S_FC_min,'UpperBound',S_FC_max);
    % PV area in [m2]
    Area_PV        = optimvar('Area_PV','LowerBound',Area_PV_min,'UpperBound',Area_PV_max);

% OPERATIONAL decision variables
    
    % fuel cell generated power in [kW]
    P_FC           = optimvar('P_FC',nHours,'LowerBound',0,'UpperBound',S_FC_max);
    % additional variables for FC operation 
    P_FC_On        = optimvar('P_FC_On',nHours,'LowerBound',0,'UpperBound',S_FC_max);
    delta_On       = optimvar('delta_On',nHours,'Type','integer','LowerBound',0,'UpperBound',1);
    % imported power from the grid in [kW]
    P_imp          = optimvar('P_imp',nHours,'LowerBound',0);
    % exported power to the grid in [kW]
    P_exp          = optimvar('P_exp',nHours,'LowerBound',0);
        
    % Battery energy content in [kWh]
    E_b            = optimvar('E_b',nHours,'LowerBound',0,'UpperBound',C_b_max);   
    % Battery charging power in [kW]
    P_b_ch         = optimvar('P_b_ch',nHours,'LowerBound',0,'UpperBound',P_b_max);
    % Battery discharging power in [kW]
    P_b_disch      = optimvar('P_b_disch',nHours,'LowerBound',0,'UpperBound',P_b_max);
    
    %% Derived variables

    % PV generated power
    P_PV      = irradiance.*eff_PV*Area_PV/1000;                           % [kW]
    P_PV_peak = 1000*eff_PV*Area_PV/1000;                                  % [kW]
    
    % C-rate of battery [kWh]
    P_b_lim     = C_b / 3600;                                              % new max capacity in [kWh], 1-C rate means that the battery charge/dsicharge in 1 hour

    % mass flow of H2
    m_flow_H2   = (P_FC/eff_FC/LHV)*3600;                                  % Consumed mass flow rate of H2 [kg/h]
   
    
    %% CONSTRAINTS

    % energy balance

    sizingprob.Constraints.EnBalance    = (P_FC + P_PV + P_b_disch + P_imp) == (P_b_ch + P_load + P_exp);
   
    % BATTERY

    % sizingprob.Constraints.NoSimultaneousChDisch  = discharging_on + charging_on <= ones(nHours,1);
    % sizingprob.Constraints.NoSimultaneousChDisch  = discharging_on + charging_on <= ones(nHours,1);
    sizingprob.Constraints.PowerBatt_ch_0         = P_b_ch(1) == 0; 
    sizingprob.Constraints.PowerBatt_disch_0      = P_b_disch(1) == 0; 

    %sizingprob.Constraints.Ch_on1                = P_b_ch <= P_b_max * charging_on;
    sizingprob.Constraints.Ch_on2                 = P_b_ch <= P_b_lim;
    %sizingprob.Constraints.Disch_on1             = P_b_disch <= P_b_max * discharging_on;
    sizingprob.Constraints.Disch_on2              = P_b_disch <= P_b_lim;
    sizingprob.Constraints.E_b                    = E_b(idxHr2ToEnd) - E_b(idxHr2ToEnd-1) == eff_ch*P_b_ch(idxHr2ToEnd)*deltat - P_b_disch(idxHr2ToEnd)/eff_disch*deltat;
    sizingprob.Constraints.E_b_cont               = E_b(1) == E_b(end); 
    sizingprob.Constraints.E_b_max                = E_b <= C_b;
    % sizingprob.Constraints.E_b_min              = E_b >= 0.2 * C_b;
    
    % FC    
    sizingprob.Constraints.MaxPowerFC             = P_FC <= P_FC_On;
    % sizingprob.Constraints.MinPowerFC             = P_FC >= 0.2 * P_FC_On;
    sizingprob.Constraints.PowerFC1               = P_FC_On <= S_FC_max * delta_On;
    sizingprob.Constraints.PowerFC2               = P_FC_On <= S_FC;
    sizingprob.Constraints.PowerFC3               = P_FC_On >= S_FC - S_FC_max.*(ones(nHours,1) - delta_On);
    sizingprob.Constraints.PowerFC4               = P_FC <= S_FC_max * delta_On;

   
%% OBJECTIVE FUNCTION

cost_inst     = (P_PV_peak * UP_PV * ann_PV + S_FC * UP_FC * ann_FC + C_b/3600 * UP_b * ann_b)/1000;    % [kCHF/y]
cost_imp     = sum(c_h2*m_flow_H2)/1000 + sum(P_imp.*c_gridimp)/1000;                                   % [kCHF/y]
% cost_imp      = sum(current_c_h2 .* m_flow_H2) / 1000 + sum(P_imp .* c_gridimp) / 1000;               % [kCHF/y]
cost_exp      = sum(P_exp .* c_gridexp) / 1000;                                                         % [kCHF/y]
cost_maint    = (maint_PV * P_PV_peak * UP_PV+ maint_FC * S_FC * UP_FC + maint_b * UP_b * C_b/3600)/1000;  % [kCHF/y]

cost = cost_inst + cost_imp - cost_exp + cost_maint;

% set objective
sizingprob.Objective = cost;

%% Solve optimization problem
% intcon = [];
% options = optimoptions('intlinprog','MaxTime',MaxSimTime);
% [solution,fval,reasonSolverStopped] = solve(sizingprob,'Options',options);
[solution,fval,reasonSolverStopped] = solve(sizingprob);
% show problem
% show(sizingprob);

% % Store the results for H2 costs
% P_FC_matrix(:, i)   = solution.P_FC; 
% S_FC_matrix(:, i)   = solution.S_FC;
% A_PV_matrix(:, i)   = solution.Area_PV;
% C_b_matrix(:, i)    = solution.C_b;  
% c_tot_matrix(:, i)  = evaluate(cost, solution);   
% end
%% Post-processing and results overview

cost_tot_1      = evaluate(cost, solution);                                      % optimization problem result
cost_inst_1       = evaluate(cost_inst, solution);                                 % kCHF/y
cost_imp_1        = evaluate(cost_imp, solution);                                  % kCHF/y
cost_maint_1      = evaluate(cost_maint, solution);                                % kCHF/y
cost_exp_1       = evaluate(cost_exp,solution);

Area_PV_opt     = solution.Area_PV;                                              % [m2]
P_PV_opt        = irradiance.*eff_PV.*solution.Area_PV./1000;                    % [kW]
P_FC_opt_1        = solution.P_FC;                                               % [kW]
P_b_disch_opt   = solution.P_b_disch;                                            % [kW]
P_b_ch_opt      = solution.P_b_ch;                                               % [kW]
P_imp_opt       = solution.P_imp;                                                % [kW]
P_exp_opt       = solution.P_exp;                                                % [kW]
SOC_opt         = solution.E_b/solution.C_b;                                     % [-]
S_FC_opt        = solution.S_FC;                                                 % [kW]
C_b_opt         = solution.C_b;                                                  % [kJ]

E_b_opt         = solution.E_b;
E_ch            = sum(solution.P_b_ch);
E_FC            = sum(solution.P_FC);                                            %[kWh]
E_exp           = sum(solution.P_exp);
E_PV            = sum(P_PV_opt);
E_disch         = sum(solution.P_b_disch);
E_imp           = sum(solution.P_imp);
E_load          = sum(P_load); 
E_consumed      = E_ch + E_exp + E_load;
E_supplied      = E_PV + E_disch + E_imp + E_FC;

delta_On_opt    = evaluate(delta_On,solution); 
tot_hours_opt   = sum(delta_On_opt);           % [h]
m_flow_H2_opt   = P_FC_opt_1./eff_FC/LHV*3600; % [kg/h]

save('scenariocost1.mat', 'cost_inst_1', 'cost_imp_1', 'cost_exp_1', 'cost_maint_1', 'cost_tot_1');
save('MILP_MILPDEG1.mat', 'P_FC_opt_1');
%% INDICATORS 
%for ch2=3
CF= sum(P_FC_opt_1)/(S_FC_opt*nHours)*100;% [%]
% energy indicators
U_FC   = E_FC/(E_load+E_exp+E_ch)*100;  % [%] 
U_PV   = E_PV/(E_load+E_exp+E_ch)*100;  % [%]
U_grid = E_imp/(E_load+E_exp+E_ch)*100; % [%]

% %% Start-Stop cycle per hour calculation
% state_changes = diff(delta_On_opt);
% % Identify start transitions (0 → 1) and stop transitions (1 → 0)
% starts = state_changes == 1; % Logical array for start transitions
% stops = state_changes == -1; % Logical array for stop transitions
% % Total number of transitions (both starts and stops)
% total_transitions = sum(starts) + sum(stops);
% % Each start + stop pair constitutes one cycle
% total_cycles = total_transitions / 2;
% 
% %% Time duration of low-power operation condition In Hours 
% % Low power threshold (5% of maximum power)
% low_power_threshold = 0.05 * S_FC_opt;
% low_power_hours = (0 < P_FC_opt_1)&(P_FC_opt_1 <= low_power_threshold);
% % total duration in hours
% low_power_duration = sum(low_power_hours); %[h]
% %% Time duration of high-power operation condition in Hours
% % High power threshold (90% of maximum power)
% high_power_threshold = 0.9 * S_FC_opt;
% high_power_hours = P_FC_opt_1 >= high_power_threshold;
% % total duration in hours
% high_power_duration = sum(high_power_hours); %[h]
% %% Time duration of large load variations
% % Load variations threshold (10% of maximum power)
% load_var_threshold = 0.10 * S_FC_opt;
% load_var_rate      = diff(P_FC_opt_1)/1; % [kW]
% load_cycles        = abs(load_var_rate) >= load_var_threshold;
% num_load_cycles    = sum(load_cycles); % Count the number of cycles

%%  plots
% % Bar plot for costs
% costs = [cost_inst, cost_imp,-cost_exp,cost_maint,cost_total];
% labels = {'Installation', 'Import', 'Export', 'Maintenance','Total'};
% cmap = crameri('batlow', length(costs));
% figure;
% b = bar(costs);
% b.FaceColor = 'flat';
% b.FaceAlpha = 0.5;
% for i = 1:length(costs)
%     b.CData(i, :) = cmap(i, :);
% end
% grid on
% set(gca, 'xticklabel', labels);
% set(gca, 'FontSize', font-5);
% % ylim([-cost_total cost_total+3])
% ylabel('Costs (kEUR/y)');

% Input data plots
cmap = crameri('bamako', 5);
colors = cmap;  

figure('Position', [100, 100, 900, 1000]); 
subplot(4,1,1)
plot(Time, P_load,'Color',colors(3,:),'LineWidth',0.6)
xlabel('Time [h]','FontWeight', 'bold','FontName','Times New Roman');
xlim([Time(1) Time(end)]);
ylabel('Power [kW]','FontWeight', 'bold','FontName','Times New Roman');
title ('Power demand','FontWeight', 'bold','FontName','Times New Roman');

subplot(4,1,2)
plot(Time, irradiance/1000,'Color',colors(4,:),'LineWidth',0.6)
xlabel('Time [h]','FontWeight', 'bold','FontName','Times New Roman');
xlim([Time(1) Time(end)]);
ylabel('G [kW/m^2]','FontWeight', 'bold','FontName','Times New Roman');
% title ('Irradiance','FontWeight', 'bold','FontName','Times New Roman');

subplot(4,1,3)
plot(Time(1:168),c_gridimp(1:168),'Color',colors(1,:),'LineWidth',linew)
xlabel('Time [h]','FontWeight', 'bold','FontName','Times New Roman');
xlim([Time(1) Time(168)]);

ylabel('Cost [CHF/kWh]','FontWeight', 'bold','FontName','Times New Roman');
title('Cost of imported electricity','FontWeight', 'bold','FontName','Times New Roman')

subplot(4,1,4)
plot(Time,c_gridexp,'Color',colors(2,:),'LineWidth',linew)
xlabel('Time [h]','FontWeight', 'bold','FontName','Times New Roman');
xlim([Time(1) Time(end)]);

ylabel('Cost [CHF/kWh]','FontWeight', 'bold','FontName','Times New Roman');
title('Cost of exported electricity','FontWeight', 'bold','FontName','Times New Roman')

% %% Figures - single week plots
% % functions directory
% addpath(path_functions);
% 
% % for a specific week between start and finish
% startS=9*26*24; % August
% finishS=startS+24*7+1;
% startW=13*26*24; % December
% finishW=startW+24*7+1;
% 
% SelectedWeekSummer_PLOT(linew,font,Time,startS,finishS,P_PV_opt,P_FC_opt_1,P_imp_opt,P_exp_opt,S_FC,SOC_opt,P_load,delta_On_opt)
% SelectedWeekWinter_PLOT(linew,font,Time,startW,finishW,P_PV_opt,P_FC_opt_1,P_imp_opt,P_exp_opt,S_FC,SOC_opt,P_load,delta_On_opt)
% PowerFCSummer_cost_PLOT(linew,font,numCosts,Time,startS,finishS,P_FC_matrix,S_FC_matrix,c_h2)
% movegui('center');
% PowerFCWinter_cost_PLOT(linew,font,numCosts,Time,startW,finishW,P_FC_matrix,S_FC_matrix,c_h2)
% movegui('center');
% Energy_SizeIndicators_PLOT(linew,font,Time,nHours,S_FC_matrix,P_FC_matrix,P_load,eff_FC,E_load,E_ch,E_exp,numCosts,c_h2)