function [] = Energy_SizeIndicators_PLOT(linew,font,Time,nHours,S_FC_matrix,P_FC_matrix,P_load,eff_FC,E_load,E_ch,E_exp,numCosts,c_h2);
%% INDICATORS 
% Design indicator
% for ch2 vector costs 
%R_FC =zeros(numCosts,1);
CF = zeros(numCosts,1);
U_FC =zeros(numCosts,1);
for j=1:numCosts
    % design indicators 
    % R_FC(j)   = S_FC_matrix(1,j)*eff_FC/mean(P_load); % [-] fuel cell ratio
    % Capacity factor 
    CF(j) = sum(P_FC_matrix(:,j))/(S_FC_matrix(1,j)*nHours)*100; 
    % energy indicators
    U_FC(j)   = (sum(P_FC_matrix(:,j)*1)/(E_load+E_ch+E_exp))*100; % [%] FC utilization
end 

colors = [1 44 86 129 172 214];  
cmap = crameri('batlow', max(colors)); 

figure('Position', [100, 100,700, 300]);
plot(c_h2, U_FC, '-d', 'Color', cmap(colors(1), :), 'LineWidth', linew); 
ylabel ('FC utilization [%]', 'FontWeight', 'bold', 'FontName', 'Times New Roman');
xlabel('Cost of hydrogen [CHF/kg]', 'FontWeight', 'bold', 'FontName', 'Times New Roman');

figure('Position', [100, 100, 700, 300]);
plot(c_h2, CF, '-d', 'Color', cmap(colors(3), :), 'LineWidth', linew);
ylabel ('Capacity Factor [%]', 'FontWeight', 'bold', 'FontName', 'Times New Roman');
xlabel('Cost of hydrogen [CHF/kg]', 'FontWeight', 'bold', 'FontName', 'Times New Roman');


set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0, 0, 7, 4.5]); 
set(gca, 'FontSize', font-8, 'FontName', 'Times New Roman'); 
set(findall(gcf, 'Type', 'line'), 'LineWidth', 1.2); 

% Black box for clear visualization
set(gca, 'LineWidth', 0.7, 'Box', 'on', 'XColor', 'k', 'YColor', 'k'); 

end
