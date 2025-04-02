close all
clear all
clc
 %% INPUT
font   = 18;
linew  = 1.5; 
nHours = 8760; 
Time =(1:nHours)';
%% PAPER MODEL 
% operating parameters
k_r  = 1; %=1.72                    % Correction coefficient of real-world driving cycle
n11  = 2058/8760;                   % number of load cycles [cycles/h]
n12  = 367/8760;                    % number of start/stop cycles [cycles/h]
t11  = 44/8760;                     % time duration of low-power operation [h/h]
t12  = 2773/8760;                   % time duration of high-power operation [h/h]

% Coefficient in degradation rate from manufacturer's datasheet of a PEMFC ( PEI et al. paper-2008)
k1    = 0.0000593;     % [%/cycle] Degradation coefficient for load changing
k2    = 0.00196;       % [%/cycle] Degradation coefficient for start/stop cycling (PEI 2008)
k3    = 0.00126;       % [%/h]     Degradation coefficient for low-power operation (idle time)
k4PEI = 0.00147;       % [%/h]     Degradation coefficient for high-power operation (PEI 2008)
k4DES = 0.00103;       % [%/h]     Degradation coefficient for high-power operation (DESANTES 2022)

%% Fuel Cell Degradation Rate

D_fc1PEI  = (k1 * n11 + k2 * n12 + k3 * t11 + k4PEI *t12); % Degradation rate [%/h]
D_fc1DES  = (k1 * n11 + k2 * n12 + k3 * t11 + k4DES *t12); % Degradation rate [%/h]

D_lowpow1         = k3*t11;        %[%/h]
D_highpow1PEI     = k4PEI*t12;     %[%/h]
D_startstop1PEI   = k2*n12;        %[%/h]
D_loading1        = k1*n11;        %[%/h]
D_highpow1DES     = k4DES*t12;     %[%/h]

% Lifetime predicted
Newlifetime1PEI = 10/D_fc1PEI; %[h]
Newlifetime1DES = 10/D_fc1DES; %[h]


%% plot
% PEI
degcoeffsPEI = [D_lowpow1, D_highpow1PEI, D_loading1, D_startstop1PEI];
labelsPEI = {'Low Power Operation', 'High Power Operation', 'Loading Operation', 'Start/Stop'};

cmap = crameri('batlow', length(degcoeffsPEI));

figure;
bPEI = bar(degcoeffsPEI, 'FaceColor', 'flat', 'EdgeColor', 'k','FaceAlpha',0.4, 'LineWidth', 0.7); 

for i = 1:length(degcoeffsPEI)
    bPEI.CData(i, :) = cmap(i, :);
end
set(gca, 'XTickLabel', labelsPEI, 'FontSize', font-7, 'FontWeight', 'bold', ...
         'XTickLabelRotation', 25, 'FontName', 'Times New Roman');
ylabel('Degradation rate [%/h]', 'FontSize', font-7, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
ylim([0 5e-4]);
% DES
degcoeffsDES = [D_lowpow1, D_highpow1DES, D_loading1, D_startstop1PEI];
labelsDES = {'Low Power Operation', 'High Power Operation', 'Loading Operation', 'Start/Stop'};

cmap = crameri('batlow', length(degcoeffsDES));

figure;
bDES = bar(degcoeffsDES, 'FaceColor', 'flat', 'EdgeColor', 'k','FaceAlpha',0.4, 'LineWidth', 0.7); 

for i = 1:length(degcoeffsDES)
    bDES.CData(i, :) = cmap(i, :);
end
set(gca, 'XTickLabel', labelsDES, 'FontSize', font-7, 'FontWeight', 'bold', ...
         'XTickLabelRotation', 25, 'FontName', 'Times New Roman');
ylabel('Degradation rate [%/h]', 'FontSize', font-7, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
ylim([0 5e-4]);
box on;
set(gca, 'LineWidth', 0.7);