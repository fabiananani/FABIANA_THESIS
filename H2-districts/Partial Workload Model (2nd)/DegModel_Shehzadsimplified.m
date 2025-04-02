clear all
close all
clc
%% Input data

% Define the main path manually
    path_main = 'C:\Users\nanf\Desktop\H2-districts';
% path for input files
    path_input = fullfile(path_main, 'Input');
% Reading an Excel file from the Input directory
    excelFilePath = fullfile(path_input, 'DegradationModel');
    Input         = readtable(excelFilePath);
% Input data from paper
%d_f            = 0.02;               % degradation rate of the fuel cell at maximum output power [-]
delta_eff_max   = 0.1;                % maximum allowable degradation of efficency        
P_FC_max       = 15.9555;                % [kW]
eff_FC         = 0.51;                 % Nominal efficency
% Input vector from optimization problem
P_FC_opt          = Input.P_FC;                   % optimal fuel cell power [kW]
delta_On_opt      = Input.delta_On;               % integer variable for fuel cell operation

% Time vector of simulation of 1 year
nHours         = 8760; 
Time           = (1: nHours)';
nH_lifetime  =  40000; %[h]

LHV         = 33.3 * 3.6 * 10^3;           % Hydrogen higher heating value [kJ/kg]
rho_H2      = 0.0898765;                    % H2 density @ T=0°C,1 bar [kg/m3]
LHV1        = 33.33*rho_H2;                 % [kWh/Nm3]

% Initialization of degradation rate
fuel_cons        = zeros(nHours, 1);
fuel_cons(1)     = eff_FC*LHV1;  % Initial fuel consumption value [kWh/Nm^3] 

% Degradation calculation loop
for i = 2:nHours
    fuel_cons(i) = (1 - ((delta_eff_max /nH_lifetime) * delta_On_opt(i-1))) * fuel_cons(i-1); % [kWh/Nm3]
end
% Mass flow rate considering degradation
flow_H2                   = (P_FC_opt./ fuel_cons)* rho_H2;         % [kg/h]
flow_H2_nom               = (P_FC_opt./eff_FC/LHV1)* rho_H2;        % [kg/h]
tot_m_flow_H2             = sum(flow_H2);     
tot_m_flow_H2_nom         = sum(flow_H2_nom);
increase_mflowH2_simp     = (tot_m_flow_H2-tot_m_flow_H2_nom)*100/tot_m_flow_H2_nom; % [%]
% eff_FC_deg               = P_FC_opt.*FC_On_opt*3600./flow_H2./HHV;
eff_FC_deg1_simp          = fuel_cons./ LHV1;
annualEffReduction_simp   = abs(eff_FC_deg1_simp(end)-eff_FC_deg1_simp(1))*100/eff_FC_deg1_simp(1); % percentage annual efficency reduction [%]
New_Lifetime_simp         = 10/annualEffReduction_simp;

% Plot of degradation rate over time
figure;
plot(Time, eff_FC_deg1_simp, 'LineWidth', 2);
xlabel('Time (hours)');
ylabel('Efficnency Degradation Rate');
title('Degradation efficency rate with respect Operational Hours');


save('plot_shehzadmodel.mat', 'annualEffReduction_simp', 'New_Lifetime_simp','increase_mflowH2_simp' ,'eff_FC_deg1_simp','-append');
