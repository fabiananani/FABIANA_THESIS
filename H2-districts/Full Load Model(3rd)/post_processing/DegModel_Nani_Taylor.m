clear all
close all
clc

%% Input data

% Define the main path manually
    path_main = 'C:\Users\nanf\Desktop\H2-districts';
% path for input files
    path_input = fullfile(path_main, 'Input');
% Reading an Excel file from the Input directory
    %excelFilePath = fullfile(path_input, 'DegradationModel.xlsx');
    excelFilePath = fullfile(path_input, 'Error.xlsx');
    Input         = readtable(excelFilePath);

% Input vector from optimization problem
% P_FC_max       = 15.9555;                      % [kW]
P_FC_max       = 15.5475;                        % [kW]
P_FC_opt       = Input.P_FC;                   % optimal fuel cell power [kW]
delta_On_opt      = Input.delta_On;            % integer variable for fuel cell operation
% P_FC_years = repmat(P_FC_opt,5,1);
% delta_ON_years = repmat(delta_On_opt,5,1);
 
% Time vector of simulation of 1 year
font           = 18;
linew          = 1.5; 
nHours         = 8760; 
Time           = (1: nHours)';
eff_FC         = 0.51;                          % Nominal electrical efficency of FC 
LHV            = 120.1*10^3;                    % Hydrogen lowe heating value [kJ/kg]
% Lifetime in hours for PEMFC from paper of Shehzad et al. 
nH_lifetime  =  40000;                                  % [h]
% PEM fuel cell lifetime was defined by the time elapsed until 10% of the
% initial voltage/performance is lost fro STAYERS project
Delta_eff_FC = 0.1*eff_FC;                                 % [-]
% Fuel Cell degradation rate as a constant value in Full-Load-Model
Const1 = Delta_eff_FC / nH_lifetime;                     % [1/h]
% further constant for MILP model from Taylor expansion
Const3 = Const1/LHV/eff_FC^2; % [kg/h/kJ]
Const2 = Const3*P_FC_opt; % [kg/h/s]
%% Dynamic calculation
% Cumulative summation of the operating hours (FC_On=1 if FC is ON, FC_On=0
% if FC is OFF) 
cumsum = zeros(nHours,1);
cumsum(1) = delta_On_opt(1);

% H2 consumed mass flow rate with respect a degradated efficency 

m_flow_H2 = (P_FC_opt./eff_FC/LHV)*3600;  % mass flow rate vector in no-degradation model 
m_flow_H2_MILP = zeros(nHours,1); 
m_flow_H2_MILP(1)= m_flow_H2(1)+(Const2(1).*cumsum(1))*3600;   % Equation used in MILP: mH2(t) = mH2,id(t)+ CONST*cumsum(delta_ON)

for i=2:nHours
    cumsum(i)               =cumsum(i-1)+delta_On_opt(i);
    m_flow_H2_MILP(i)       = m_flow_H2(i)+Const2(i).*cumsum(i)*3600; % [kg/h]
end 
 eff_deg_FC2          = P_FC_opt./m_flow_H2_MILP./LHV*3600; % efficency behavior considering the equation for massflowrate implemented in MILP
%% Post-processing
annualEffReduction2          = abs(eff_deg_FC2(end)-eff_deg_FC2(1))*100/eff_deg_FC2(1); % efficency reduction [%]
New_Lifetime2                = 10*1/annualEffReduction2;
tot_m_flow_H2_nom            = sum(m_flow_H2);     
tot_m_flow_H2_new            = sum(m_flow_H2_MILP);
increase_mflow2              = abs(tot_m_flow_H2_nom-tot_m_flow_H2_new)*100/tot_m_flow_H2_nom; % [%]

%% plot
% cmap = crameri('batlow',5);
% colors = cmap;
% 
% %1
% figure('Position', [100, 100, 800, 500]);
% plot(Time,eff_deg_FC,'LineWidth',3,'Color',colors(1,:))
% ylabel('FC efficiency [-]','FontWeight','bold')
% xlabel('Time [h]','FontWeight','bold')
% 
% %2
% %period 
% startS=6*26*24;
% finishS=startS+24*7+1;
% year=365*2023+126;
% start_date=startS/24+year;
% finish_date=finishS/24+year;
% tt=datetime(start_date:(1/24):finish_date, 'ConvertFrom', 'datenum');
% 
% %figure
% figure('Position', [100, 100, 1800, 300]);
% plot(tt,m_flow_H2_MILP(startS:finishS),'LineWidth',linew,'Color',colors(4,:))
% hold on 
% plot(tt,m_flow_H2(startS:finishS),'LineWidth',linew,'Color',colors(2,:))
% hold on
% legend('H2 Mass flow rate with degradation','H2 Mass flow rate without degradation')
% ylabel('massflow rate [kg/h]','FontWeight','bold')
% xlim([tt(1) tt(end)])
% set(gca, 'FontSize', font-5);
% leg = legend('Location', 'eastoutside');
% 
% save('plot_nanimodel.mat','annualEffReduction2','New_Lifetime2','increase_mflow2','eff_deg_FC2','-append');