clear all
close all
clc
% PAPER: "A quick evaluating method for automotive fuel cell lifetime" (2008)
% Coefficient in degradation rate from manufacturer's datasheet of a PEMFC (paper)
k1 = 0.0000593;     % [%/cycle] Degradation coefficient for load changing 
k2 = 0.00196;       % [%/cycle] Degradation coefficient for start/stop operation 
k3 = 0.00126;       % [%/h]     Degradation coefficient for low-power operation (idle time)
k4 = 0.00147;       % [%/h]     Degradation coefficient for high-power operation

% Operating parameters per hour
n21 = 56;    % cycle/h - average load change cycles per hour
n22 = 0.99;  % cycle/h -start/stop cycle
t21 = 13/60; % h/h - average idling time per hour
t22 = 14/60; % h/h - average high power load operation time per hour
% DEGRADATION RATE 
D_fc1 = (n21*k1+n22*k2+t21*k3+t22*k4); % [%/h]
Newlifetime2  = 10/D_fc1;
D_lowpow2    = k3*t21;  %[%/h]
D_highpow2   = k4*t22;  %[%/h]
D_startstop2 = k2*n22;  %[%/h]
D_loading2   = k1*n21;  %[%/h]
%% plot
degcoeffs= [D_lowpow2, D_highpow2,D_loading2,D_startstop2];
labels = {'Low Power Operation', 'High Power Operation', 'Loading Operation', 'Start/Stop cycling'};
cmap = crameri('batlow', length(degcoeffs));
figure;
b = bar(degcoeffs);
b.FaceColor = 'flat';
b.FaceAlpha = 0.8;
for i = 1:length(degcoeffs)
    b.CData(i, :) = cmap(i, :);
end

set(gca, 'xticklabel', labels);
ylabel('Degradation rate (%)','FontWeight','bold');
ylim([0 3.5e-3])
