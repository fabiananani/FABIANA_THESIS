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

%% INPUT DATA 
% From paper of M.F. Shehzad et al. (2019)
%d_f            = 0.02;               % degradation rate of the fuel cell at maximum output power [-]
delta_eff_max   = 0.1;                % maximum allowable degradation of efficency 
% Datasheet NEXA 1200 
nH_lifetime    = 40000;   % Number of life hours of the fuel cell [h]

% From optimization problem
P_FC_max       = 15.955;                      % [kW]
eff_FC         = 0.51;                        % Nominal efficency
P_FC_opt       = repmat(Input.P_FC, 8, 1);    % optimal fuel cell power [kW]
delta_On_opt   = repmat(Input.delta_On, 8, 1);  % integer variable for fuel cell operation

LHV         = 33.3 * 3.6 * 10^3;           % Hydrogen higher heating value [kJ/kg]
rho_H2      = 0.0898765;                    % H2 density @ T=0°C,1 bar [kg/m3]
LHV1        = 33.33*rho_H2;                 % [kWh/Nm3]HHV         = 39.39 * 3.6 * 10^3;           % Hydrogen higher heating value [kJ/kg]

% Time vector of simulation of 1 year
nHours         = 8760*8; 
Time           = (1:nHours)';

% Initialization of degradation rate
fuel_cons        = zeros(nHours, 1);
fuel_cons(1)    = eff_FC*LHV1;  % Initial fuel consumption value [kWh/Nm^3] 

% Degradation of fuel consumption
for i = 2:nHours
    fuel_cons(i) = (1 - ((delta_eff_max / (P_FC_max * nH_lifetime)) * P_FC_opt(i-1) * delta_On_opt(i-1))) * fuel_cons(i-1); % [kWh/Nm3]
end
% Mass flow rate considering degradation 
flow_H2     = (P_FC_opt./ fuel_cons)* rho_H2;         % [kg/h]
% massflowrate without degradation
flow_H2_nom = (P_FC_opt./eff_FC/LHV1)* rho_H2;        % [kg/h]
%% No sense 
tot_m_flow_H2            = sum(flow_H2);     
tot_m_flow_H2_nom        = sum(flow_H2_nom);
increase_mflowH2         =(tot_m_flow_H2-tot_m_flow_H2_nom)*100/tot_m_flow_H2_nom; % [%]
% eff_FC_deg               = P_FC_opt.*3600./flow_H2./HHV;
eff_FC_deg1              = fuel_cons./ LHV1;
totEffReduction          = abs(eff_FC_deg1(end)-eff_FC_deg1(1))*100/eff_FC_deg1(1); % [%]

%% Plot of degradation rate over time
cmap = crameri('batlow',2);
colors = cmap;  

figure('Position', [100, 100, 600, 400]);
plot(Time, fuel_cons, 'LineWidth', 5,'Color',colors(1,:));
xlabel('Time [h]','FontWeight','bold');
xlim([Time(1) Time(end)]);
ylabel('Fuel consumption [kWh/Nm3]','FontWeight','bold');
title('Yearly profile of fuel consumption');

%% plot efficency
figure('Position', [100, 100, 600, 400]);
plot(Time,eff_FC_deg1,'LineWidth', 5,'Color',colors(1,:))
xlabel('Time [h]','FontWeight','bold');
xlim([Time(1) Time(end)]);
ylabel('Efficency[-]','FontWeight','bold');


set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0, 0, 7, 4.5]); % Width=3.5in, Height=2.5in
set(gca, 'FontSize', 10); % Set axis font size and font name
set(findall(gcf, 'Type', 'line'), 'LineWidth', 1.2); % Set line width


%% PLOT MASS FLOW RATE 
figure
plot(Time(43600:43800),flow_H2_nom(43600:43800))
hold on
plot(Time(43600:43800),flow_H2(43600:43800))
legend('no deg','deg')