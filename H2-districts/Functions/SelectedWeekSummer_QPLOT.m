function [] =SelectedWeekSummer_QPLOT(linew,font,Time,start,finish,P_th_FC_opt,P_th_HEX_out_opt,P_th_load,P_th_amb)

year=365*2023+126;
start_date=start/24+year;
finish_date=finish/24+year;
tt=datetime(start_date:(1/24):finish_date, 'ConvertFrom', 'datenum');


% colors = [1 44 86 129 172 214];
% cmap = crameri('batlow');
cmap = colormap(hot(256));
colors = round(linspace(50, 200, 8)); 

figure('Position', [100, 100, 1700, 300]);
h1=area(tt,P_th_HEX_out_opt(start:finish), 'FaceColor', cmap(colors(1), :), 'FaceAlpha', 0.4, 'EdgeColor', 'none');
hold on
plot(tt,P_th_HEX_out_opt(start:finish),'LineWidth',linew,'Color',cmap(colors(1), :));
hold on
h2=area(tt, P_th_FC_opt(start:finish), 'FaceColor', cmap(colors(6), :), 'FaceAlpha', 0.4, 'EdgeColor', 'none');
hold on
plot(tt,P_th_FC_opt(start:finish),'LineWidth',linew,'Color',cmap(colors(6), :));
hold on
h3=area(tt,- P_th_load(start:finish), 'FaceColor',cmap(colors(7), :), 'FaceAlpha', 0.4, 'EdgeColor', 'none');
hold on
plot(tt,-P_th_load(start:finish),'LineWidth',linew,'Color',cmap(colors(7), :));
hold on 


xlim([min(tt) max(tt)])
ylim([-10 20])
% ylim([-180 150])
% ylim([-1.2*max([max(P_imp_opt) max(P_e_opt)]) 1.2*max(P_PV_opt) ])
ylabel ('Power [kW]','fontweight','bold','FontName','Times New Roman');
% xlabel('Time','fontweight','bold');
legend([h1, h2, h3],'Q_{HEX}','Q_{FC}','Q_{load}','location','eastoutside','FontName','Times New Roman')

set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0, 0, 7, 4.5]); % Width=3.5in, Height=2.5in
set(gca, 'FontSize', font-8,'FontName','Times New Roman'); 
set(findall(gcf, 'Type', 'line'), 'LineWidth', 1.2); 
set(gca, 'LineWidth', 0.7, 'Box', 'on', 'XColor', 'k', 'YColor', 'k'); 

end