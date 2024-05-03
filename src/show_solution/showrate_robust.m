function showrate_robust(h,Err,c1,c2,name1,name2,pk)
% Copyright (C) Andres Rubiano.

tiledlayout(1,1)

Lc = length(c1);


cMap = [[0.6350 0.0780 0.1840];
        [0.4660 0.6740 0.1880];
        [0.4940 0.1840 0.5560];
        [0.8500 0.3250 0.0980];
        [0      0.4470 0.7410]];

if contains(name1, 'mu')
    c1 = c1.*(1/2);
end

% Plot the different  curves.

nexttile;

for k = 1 : Lc
    thisColor = cMap(k,:);
    loglog(1./(h.^2), Err(:,k), 'Color',thisColor, 'LineWidth', 4,'LineStyle',':','Marker','o','MarkerSize',10,'MarkerFaceColor',thisColor);
    hold on;
    tickLabels{k} = strcat(num2str(c1(k,1), '%.e'), ',  ',num2str(c2(k,1), '%.3e'));
end

clim([1, Lc]);
grid off;
colormap(cMap);
h1 = colorbar('Ticks', 1:Lc, 'TickLabels', tickLabels, 'Location', 'westoutside', 'TickLabelInterpreter','latex');
h1.Label.String = strcat(name1, ',  ', name2);
h1.Label.VerticalAlignment = 'middle';
h1.Label.Rotation = 0;
h1.Label.Interpreter = 'latex';
h1.Label.Position = [0.5 Lc+0.05];

xlabel('Number of elements');
ylabel('$\overline{\textnormal{e}}_*$');
xticks(1./(h.^2));
set(gca,'Linewidth',2);
set(gca,'Fontsize',28);
legend('off');