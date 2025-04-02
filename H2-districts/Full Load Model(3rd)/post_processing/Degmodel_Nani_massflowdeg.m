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

% Input vector from optimization problem
P_FC_max       = 15.955;                       % [kW]
P_FC_opt       = Input.P_FC;                   % optimal fuel cell power [kW]
delta_On_opt      = Input.delta_On;            % integer variable for fuel cell operation

% Time vector of simulation of 1 year
font           = 18;
linew          = 1.5; 
nHours         = 8760; 
Time           = (1: nHours)';
eff_FC         = 0.51;                          % Nominal electrical efficency of FC 
LHV            = 120.1 * 10^3;                  % Hydrogen lower heating value [kJ/kg]

% Lifetime in hours for PEMFC from paper of Shehzad et al. 
nH_lifetime  =  40000;                                  % [h]
% PEM fuel cell lifetime was defined by the time elapsed until 10% of the
% initial voltage/performance is lost fro STAYERS project
Delta_eff_FC = 0.1*0.5;                                 % [-]
% Fuel Cell degradation rate as a constant value in Nani-Model
Const = Delta_eff_FC / nH_lifetime;                     % [1/h]

%% Dynamic calculation
% Cumulative summation of the operating hours (FC_On=1 if FC is ON, FC_On=0
% if FC is OFF) 
cumsum = zeros(nHours,1);
cumsum(1) = delta_On_opt(1);

% H2 consumed mass flow rate with respect a degradated efficency 

m_flow_H2 = (P_FC_opt./eff_FC/LHV)*3600;  % "Ideal value": mass flow rate vector in no-degradation model 
m_flow_H2_MILP = zeros(nHours,1); 
m_flow_H2_MILP(1)= m_flow_H2(1)+Const.*cumsum(1);   % Equation used in MILP: mH2(t) = mH2,id(t)+ CONST*cumsum(delta_ON)

for i=2:nHours
    cumsum(i)               = cumsum(i-1)+delta_On_opt(i);
    m_flow_H2_MILP(i)       = m_flow_H2(i)+Const.*cumsum(i); % [kg/h]
    
end 
 eff_deg_FC          = P_FC_opt./m_flow_H2_MILP./LHV*3600; % efficency behavior considering the equation for massflowrate implemented in MILP
 
%% Post-processing
annualEffReduction1          = abs(eff_deg_FC(end)-eff_deg_FC(1))*100/eff_deg_FC(1); % efficency reduction [%]
New_Lifetime1                = 10*1/annualEffReduction1;
tot_m_flow_H2_nom            = sum(m_flow_H2);     
tot_m_flow_H2_new1           = sum(m_flow_H2_MILP);
increase_mflow1              = abs(tot_m_flow_H2_nom-tot_m_flow_H2_new1)*100/tot_m_flow_H2_nom; % [%]

