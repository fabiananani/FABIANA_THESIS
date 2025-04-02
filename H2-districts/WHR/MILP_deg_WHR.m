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
    excelFilePath3 = fullfile(path_input, 'PriceImpExp_ele.xlsx');
    excelFilePath4 = fullfile(path_input, 'DegradationModel.xlsx');

    Input1         = readtable(excelFilePath1);
    Input2         = readtable(excelFilePath2);
    Input3         = readtable(excelFilePath3);
    Input4         = readtable(excelFilePath4);
%% Pre-processing

    
    hh = 2; 
    nHours         = 8760*hh;                        % Number of hours simulated   [8760h]
    idxHr2ToEnd    = (2:nHours)';                 % Hours until the end
    Time           = (1:nHours)';                 % Time vector
    days           = nHours/24;                   % Number of days simulated
    weeks          = days/7;                      % Number of weeks simulated
    clear t; 
    linew          = 1.5;             
    font           = 18;
    irradiance     = Input1.G;                    % Hourly solar irradiance    [W/m^2]
%% INPUT PARAMETERS

 % general inputs

    bigM            = 600;     % Large number for big-M constraints
    MaxSimTime      = 1000;    % Maximum time for MILP solver                 [s]      
    deltat          = 3600;    % time step di 1h                              [s]       

 % Electrical district demand
        
    P_load            = Input2.P_el_nest_h;       % NEST electrical demand      [kW]
    P_th_load         = Input1.Heating;           % Heating Demand [kW] -low temperature
 
    % Fixed technical parameters  

    % HHV         = 141.9*10^3;           % (141.9 MJ/kg) Hydrogen higher heating value [kJ/kg]  
     LHV         = 120.1*10^3;           % (120.1 MJ/kg) Hydrogen lower heating value [kJ/kg]

 % Costs and Efficencies 

    c_h2         = 3;                                     % cost of hydrogen [CHF/kg]
    c_gridimp    = Input3.Priceimpgreen;                  % [CHF/kWh]
    c_gridexp    = Input3.PriceexpPV;                     % [CHF/kWh]
 
 % export price of heat (weighting factors from Baldini et al.)

    weight_winter = [0.432; 0.240]; % low temperature is 2nd factor 
    weight_mid    = [0.306; 0.107];
    weight_summer = [0.137; -0.025];
    
    c_exp_heat = c_gridimp;
   
    
    winter_1 = 31*24 + 28*24 + 21*24;
    mid_1    = winter_1 + 10*24 + 30*24 + 31*24 + 21*24;
    summer   = mid_1 + 9*24 + 31*24 + 31*24 + 21*24;
    mid_2    = summer + 9*24 + 31*24 + 30*24 + 21*24;
    winter_2 = mid_2 + 10*24;
    
    c_exp_heat(1:winter_1)        = weight_winter(2).*c_exp_heat(1:winter_1);
    c_exp_heat(winter_1+1:mid_1)  = weight_mid(2).*c_exp_heat(winter_1+1:mid_1);
    c_exp_heat(mid_1+1:summer)    = weight_summer(2).*c_exp_heat(mid_1+1:summer);
    c_exp_heat(summer+1:mid_2)    = weight_mid(2).*c_exp_heat(summer+1:mid_2);
    c_exp_heat(mid_2+1:winter_2)  = weight_winter(2).*c_exp_heat(mid_2+1:winter_2);
    

    % efficiency 
    eff_ch       = 0.95;                                  % charging efficiency of the battery 
    eff_disch    = 0.95;                                  % discharging efficency of the battery
    eff_PV       = 0.17;                                  % constant efficiency of the PV system
    eff_FC       = 0.51;
    eff_HEX      = 0.95;                                  % Heat exchanger efficiency

    % Lifetime in hours for PEMFC from Technology RoadMap-Hyydrogen and Fuel Cells,IEA(2015) 
    nH_lifetime  =  40000;                                   % [h]
    % PEM fuel cell lifetime was defined by the time elapsed until 10% of the initial voltage/performance is lost
    Delta_eff_FC = 0.1*eff_FC;                               % [-]
    % Fuel Cell degradation rate as a constant value in Nani-Model
    Const1 = Delta_eff_FC/nH_lifetime;                       % [1/h]
    % further constant for MILP model from Taylor expansion
    Const2 = Const1/LHV/eff_FC^2; %[kg/h/kJ]
    

    % Fuel Cell (FC)
        S_FC_min     = 0;                   % Minimum size FC                      [kW]
        S_FC_max     = 60;                  % Maximum size FC                      [kW]

    % PV   
        P_PV_peak_ref  = 100;                 % Reference value peak power PV      [kW]
        Area_PV_ref = P_PV_peak_ref/eff_PV;   % Reference area PV                  [m2]
        Area_PV_min = 0;                      % Minimum area PV                    [m2]
        Area_PV_max = Area_PV_ref*2;          % Maximum area PV                    [m2]
    
    % Battery
    
        C_b_max     = 96 * 3.6 * 10^3;     % Battery capacity                    [kJ]
        C_b_min     = C_b_max / 3;         % MIN Battery capacity                [kJ]
        P_b_max     = C_b_max / 3600;      % Battery capacity (1C rate)          [kWh]
    
    % HEX

    % Max thermal power recoverable from the fuel cell

    P_th_FC_max = (1 - eff_FC) * S_FC_max/eff_FC;       % Maximum heat from FC in [kW]
    U_HEX     = 2;                                   % Overall heat transfer coeff HEX           [kW/m2/K]
  
    T_hot_in   = 51; % Hot side inlet [°C]
    T_hot_out  = 41; % Hot side outlet [°C]
    T_cold_in  = 26; % Cold side inlet [°C]
    T_cold_out = 36; % Cold side outlet [°C]

    Delta_T1 = T_hot_out - T_cold_in;
    Delta_T2 = T_hot_in - T_cold_out;
    % LMTD = (Delta_T1 - Delta_T2) / log(Delta_T1 / Delta_T2);   % Logarithmic Mean Temperature Difference [K]
    LMTD = Delta_T1; % DeltaT1=DeltaT2
    A_HEX_max = P_th_FC_max/(U_HEX*LMTD);   % HEX max size [m2]
   

    % unit prices and lifetime components

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
    life_FC      = 10;          % Lifetime FC   [years]
    maint_FC     = 0.024;        % Annual cost maintenance FC, frac UP_FC
    ann_FC       = d / (1 - (1 + d)^(-life_FC));
    
    UP_HEX      = 77.79;    % Unit price heat exchanger                        [CHF/m2]         
    Fixed_HEX   = 0;        % Fixed cost for heat exchanger                    [CHF]    
    life_HEX    = 20;       % Lifetime heat exchanger                          [years]           
    maint_HEX   = 0.01;     % Annual cost maintenance HEX, frac total ann cost                   
    ann_HEX     = d / (1 - (1 + d)^(-life_HEX));
    

%% Assumed value for the power

P_FC_assumed = 17.3; % [kW]
 %% Input variables for more years
    irradiance_years = repmat(irradiance,hh,1);
    P_load_years     = repmat(P_load,hh,1);
    P_th_load_years     = repmat(P_th_load,hh,1);
    c_gridimp_years  = repmat(c_gridimp,hh,1);
    c_gridexp_years  = repmat(c_gridexp,hh,1);
    c_exp_heat_years  = repmat(c_exp_heat,hh,1);

        
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
    % Heat Exchanger size in m2
    A_HEX          = optimvar('A_HEX','LowerBound',0,'UpperBound', A_HEX_max);
% OPERATIONAL decision variables
    
    % fuel cell generated power in [kW]
    P_FC           = optimvar('P_FC',nHours,'LowerBound',0,'UpperBound',S_FC_max);
    % additional variables for FC operation implementation
    P_FC_On        = optimvar('P_FC_On',nHours,'LowerBound',0,'UpperBound',S_FC_max);
    delta_On       = optimvar('delta_On',nHours,'Type','integer','LowerBound',0,'UpperBound',1);
    % cumulative sum of working hours
    cum_hours        = optimvar('cum_hours',nHours,'LowerBound',0,'UpperBound',8760); %cumulative summation of FC_On
    cum_hours_active = optimvar('cum_hours_active', nHours, 'LowerBound', 0, 'UpperBound', 8760); % variable to truck the effective operation of fuel cell for cumulative sum 
    % mass flow rate of hydrogen 
    m_flow_H2        = optimvar('m_flow_H2',nHours,'LowerBound',0,'UpperBound',S_FC_max/eff_FC/LHV.*3600*2); %[kg/h]
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
    % Heat exchanger [kW]
    P_th_HEX_out = optimvar('P_th_HEX_out', nHours, 'LowerBound', 0, 'UpperBound', P_th_FC_max);
    %% Derived variables

    % PV generated power
    P_PV      = irradiance_years.*eff_PV*Area_PV/1000;                     % [kW]
    P_PV_peak = 1000*eff_PV*Area_PV/1000;                            % [kW]
    
    % C-rate of battery [kWh]
    P_b_lim   = C_b / 3600;                                             % new max capacity in [kWh]
    %m_flow_H2 = P_FC./eff_FC./HHV.*3600 + P_FC*Const2.*w_hours*3600;    % [kg/h]
    
  
    %% CONSTRAINTS

    % energy balance

    sizingprob.Constraints.EnBalance              = (P_FC + P_PV + P_b_disch + P_imp) == (P_b_ch + P_load_years + P_exp);
 
   
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

   
    % Constraints for cumulative summation
    sizingprob.Constraints.cum_hours_0                  = cum_hours(1) == delta_On(1);
    sizingprob.Constraints.cum_hours_update             = cum_hours(idxHr2ToEnd) == cum_hours(idxHr2ToEnd-1) + delta_On(idxHr2ToEnd);
    % ensure that cum_hours_active follows the cum_hours
    sizingprob.Constraints.cum_hours_active_ub       = cum_hours_active <= cum_hours; 
    sizingprob.Constraints.cum_hours_active_lb       = cum_hours_active >= cum_hours - (1 - delta_On) * nHours; 
    
    % sizingprob.Constraints.massflow = m_flow_H2 == P_FC./eta_FC./HHV.*deltat + deg_rate.*cum_hours;
     sizingprob.Constraints.massflow                  = m_flow_H2 == P_FC./eff_FC./LHV.*3600 + Const2*P_FC_assumed.* cum_hours_active*3600; 

   % sizingprob.Constraints.HEX_InputLimit = P_th_HEX_out <= eff_HEX* (1 - eff_FC) * P_FC/eff_FC;
 sizingprob.Constraints.HEX_InputLimit  = P_th_HEX_out <= eff_HEX* (m_flow_H2*LHV/3600-P_FC); % [kW]
 sizingprob.Constraints.HeatBalance     = P_th_HEX_out <= P_th_load_years; % [kW]
 sizingprob.Constraints.HEXsize         = P_th_HEX_out / (U_HEX * LMTD) <= A_HEX;
%% OBJECTIVE FUNCTION

cost_inst     = (P_PV_peak * UP_PV * ann_PV + S_FC * UP_FC * ann_FC + C_b/3600 * UP_b * ann_b + A_HEX*UP_HEX * ann_HEX)/1000; % [kCHF/y]
cost_imp      = (sum(c_h2 * m_flow_H2) / 1000  + sum(P_imp .* c_gridimp_years) / 1000)/hh;                                   % [kCHF/y]
cost_exp      = (sum(P_exp .* c_gridexp_years) / 1000 + sum( P_th_HEX_out.*c_exp_heat_years)/1000)/hh;                                              % [kCHF/y]
cost_maint    = (maint_PV * P_PV_peak * UP_PV + maint_FC * S_FC * UP_FC + maint_b * UP_b * C_b/3600 + maint_HEX *(A_HEX * UP_HEX)* ann_HEX)/1000 ;  % [kCHF/y]

cost = cost_inst + cost_imp - cost_exp + cost_maint;

% set objective 
sizingprob.Objective = cost;


%% Solve optimization problem

% intcon = [];
% options = optimoptions('intlinprog','MaxTime',MaxSimTime);
% [solution,fval,reasonSolverStopped] = solve(sizingprob,'Options',options);
[solution,fval,reasonSolverStopped] = solve(sizingprob);

%% show problem
% show(sizingprob);

%% Post-processing and results overview
cost_total      = evaluate(cost, solution);                                      % optimization problem result
cost_inst       = evaluate(cost_inst, solution);                                 % kCHF/y
cost_imp        = evaluate(cost_imp, solution);                                  % kCHF/y
cost_maint      = evaluate(cost_maint, solution);                                % kCHF/y
cost_exp        = evaluate(cost_exp,solution);

Area_PV_opt     = solution.Area_PV;                                              % [m2]
P_PV_opt        = irradiance_years.*eff_PV.*solution.Area_PV./1000;                    % [kW]
P_FC_opt        = solution.P_FC;                                                 % [kW]
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
E_load          = sum(P_load_years); 
E_consumed      = E_ch + E_exp + E_load;
E_supplied      = E_PV + E_disch + E_imp + E_FC;

P_th_HEX_out_opt   = solution.P_th_HEX_out;                                      % [kW]
A_HEX_opt = solution.A_HEX;                                                      % [m2]
Q_HEX = sum(solution.P_th_HEX_out);
Q_load = sum(P_th_load_years);
Q_amb = Q_load-Q_HEX;

m_flow_H2_opt   = solution.m_flow_H2;
m_flow_H2_id    = P_FC_opt./eff_FC./LHV.*3600; 

etadeg_opt      = P_FC_opt./m_flow_H2_opt./LHV*3600; 

P_th_amb = P_th_load_years-P_th_HEX_out_opt;
P_th_FC_opt  = m_flow_H2_opt*LHV/3600-P_FC_opt;

delta_On_opt    = evaluate(delta_On,solution);
hours_year      = sum(delta_On_opt); 
cum_hours_opt   = solution.cum_hours; 


%% further computations 
deltamass_flowh2    = abs(sum(m_flow_H2_opt)-sum(m_flow_H2_id))*100/sum(m_flow_H2_id); %[%]
annualEffReduction  = abs(etadeg_opt(end-3)-etadeg_opt(1))*100/etadeg_opt(1); % efficency reduction [%]
New_Lifetime        = 10*1/annualEffReduction;

CF                  = sum(P_FC_opt)/(S_FC_opt*nHours)*100;    % [%]
U_FC                = E_FC/(E_load)*100; % [%] 
U_PV                = E_PV/(E_load)*100; % [%]
U_grid              = E_imp/(E_load)*100; % [%]

% %% Figures - single week plots
% % functions directory
% addpath(path_functions);
% 
% % for a specific week between start and finish
% startS=9*26*24; % August
% finishS=startS+24*7+1;
% startW=13*26*24; % December
% finishW=startW+24*7+1;
% startA=11*26*24; % October
% finishA=startA+24*7+1;
% startSp=5*26*24; % May
% finishSp=startSp+24*7+1;
% 
% SelectedWeekSummer_QPLOT(linew,font,Time,startS,finishS,P_th_FC_opt,P_th_HEX_out_opt,P_th_load,P_th_amb)
% SelectedWeekWinter_QPLOT(linew,font,Time,startW,finishW,P_th_FC_opt,P_th_HEX_out_opt,P_th_load,P_th_amb)
% SelectedWeekAutumn_QPLOT(linew,font,Time,startA,finishA,P_th_FC_opt,P_th_HEX_out_opt,P_th_load,P_th_amb)
% SelectedWeekSpring_QPLOT(linew,font,Time,startSp,finishSp,P_th_FC_opt,P_th_HEX_out_opt,P_th_load,P_th_amb)