function showsolution_battery(node,elem,u,name)
%Showsolution displays the solution corresponding to a mesh given by [node,elem] in 2-D.
% Example:
%
%  showsolution(node,elem,u, 'facecolor',[0.5 0.9 0.45], 'FaceAlpha',0.6, 'linewidth',1)
%
% Copyright (C) Terence Yu.

% change color here
colormap Turbo

data = [node,u];
if ~iscell(elem)
    h = patch('Faces', elem,...
        'Vertices', data,...
        'FaceColor', 'interp',...
        'CData', u);
else
    max_n_vertices = max(cellfun(@length, elem));
    padding_func = @(vertex_ind) [vertex_ind,...
        NaN(1,max_n_vertices-length(vertex_ind))];  % function to pad the vacancies
    tpad = cellfun(padding_func, elem, 'UniformOutput', false);
    tpad = vertcat(tpad{:});
    h = patch('Faces', tpad,...
        'Vertices', data,...
        'FaceColor', 'interp',...
        'CData', u);
end
%axis square; 
sh = 0.5;
xlim([min(node(:,1)) - sh, max(node(:,1)) + sh])
ylim([min(node(:,2)) - sh, max(node(:,2)) + sh])
%zlim([min(u) - sh, max(u) + sh])
xlabel('$x$'); ylabel('$y$'); %zlabel('u');
xticks([min(node(:,1)) -5 -1 1 5 max(node(:,1))]);
yticks([min(node(:,2)) -5 -1 1 5 max(node(:,2))]);
title({name},'Interpreter','latex');
colorbar;
x0=0.15;
y0=0.15;
width=0.5;
height=0.35;
set(gcf,'position',[x0,y0,width,height])
set(gca,'Linewidth',2);
set(gca,'GridAlpha',0.1);
%set(gca,'Markersize',15);
set(gca,'Fontsize',18);
set(gca,'TickLabelInterpreter','LaTex')
set(gca,'Linewidth',2);

view(2); %3 for view(150,30);
% 
% if nargin>4
%     set(h,varargin{:});
% end

hold on 
t = 0:pi/20:10*pi;
lnh2 = plot3(5*sin(t), 5*cos(t),1.e10+0*sin(t),'-k','LineWidth',3);
lnh1 = plot3(sin(t), cos(t),1.e10+0*sin(t),'-k','LineWidth',3);

uistack(lnh1,'top')
uistack(lnh2,'top')
hold off