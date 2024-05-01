function showrate_robust(h,Err1,Err2,c1,c2,name1,name2)
% Copyright (C) Andres Rubiano.

tiledlayout(1,2)

Lc1 = length(c1);
Lc2 = length(c2);
cMap = [[0.6350 0.0780 0.1840];
        [0.4660 0.6740 0.1880];
        [0.4940 0.1840 0.5560];
        [0.8500 0.3250 0.0980];
        [0      0.4470 0.7410]];

% Plot the different  curves.

nexttile;

for k = 1 : Lc1
    thisColor = cMap(k,:);
    loglog(1./(h.^2), Err1(:,k), 'Color',thisColor, 'LineWidth', 4,'LineStyle',':','Marker','o','MarkerSize',10,'MarkerFaceColor',thisColor);
    hold on;
    tickLabels1{k} = num2str(c1(k,1), '%.3e');
end

clim([1, Lc1]);
grid off;
colormap(cMap);
h1 = colorbar('Ticks', 1:Lc1, 'TickLabels', tickLabels1, 'Location', 'westoutside', 'TickLabelInterpreter','latex');
h1.Label.String = name1;
h1.Label.VerticalAlignment = 'middle';
h1.Label.Rotation = 0;
h1.Label.Interpreter = 'latex';
h1.Label.Position = [0.5 Lc1+0.2];

xlabel('Number of elements');
ylabel('$||(\mathbf{u}-\mathbf{u}_h,\tilde{p}-\tilde{p}_h)||_{\mathbf{V}_1 \times Q_1}$');
ylim([min([min(Err1),min(Err2)]), max([max(Err1),max(Err2)])]);
xticks(1./(h.^2));
set(gca,'Linewidth',2);
set(gca,'Fontsize',20);

nexttile;

for k = 1 : Lc2
    thisColor = cMap(k,:);
    loglog(1./(h.^2), Err2(:,k), 'Color',thisColor, 'LineWidth', 4,'LineStyle',':','Marker','o','MarkerSize',10,'MarkerFaceColor',thisColor);
    hold on;
    tickLabels2{k} = num2str(c2(k,1), '%.3e');
end

clim([1, Lc2]);
grid off;
colormap(cMap);
h2 = colorbar('Ticks', 1:Lc2, 'TickLabels', tickLabels2, 'Location', 'eastoutside', 'TickLabelInterpreter','latex');
h2.Label.String = name2;
h2.Label.VerticalAlignment = 'middle';
h2.Label.Rotation = 0;
h2.Label.Interpreter = 'latex';
h2.Label.Position = [0.5 Lc2+0.2];

xlabel('Number of elements');
ylabel('$||(\mbox{\boldmath$\zeta$}-\mbox{\boldmath$\zeta$}_h,\varphi-\varphi_h)||_{\mathbf{V}_2 \times Q_2}$');
ylim([min([min(Err1),min(Err2)]), max([max(Err1),max(Err2)])]);
xticks(1./(h.^2));
set(gca,'Linewidth',2);
set(gca,'Fontsize',20);



