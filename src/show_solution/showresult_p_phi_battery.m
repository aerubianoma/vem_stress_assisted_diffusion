function showresult_p_phi_battery(node,elem,ph,phih)

%showresult displays the mesh and the solution
%
% Copyright (C) Terence Yu.

clf;  % clear figures

set(gcf,'Units','normal');
%set(gcf,'Position',[0.25,0.25,0.7,0.25]);

%% Plot mesh

%subplot(1,3,1);
%showmesh(node,elem); 
%hold on;
%plot(node(:,1),node(:,2),'k.', 'MarkerSize', 4);

tiledlayout(1,2)

%% Plot exact solution p

%nexttile;
%p_exact = p(node); p_exact = p_exact(:,1);
%showsolution(node,elem,p_exact(1:size(node,1)),'$\tilde{p}$');

%% Plot numerical solution

nexttile;
showsolution_battery(node,elem,ph(1:size(node,1),1),'$\tilde{p}_h$');
hold on 
fplot(@(t) 5*sin(t), @(t) 5*cos(t),'--k','LineWidth',3);
fplot(@(t) sin(t), @(t) cos(t),'--k','LineWidth',3);

%% Plot exact solution phi

%nexttile;
%phi_exact = phi(node); phi_exact = phi_exact(:,1);
%showsolution(node,elem,phi_exact(1:size(node,1)),'$\varphi$');

%% Plot numerical solution
nexttile;
showsolution_battery(node,elem,phih(1:size(node,1),1),'$\varphi_h$');
hold on 
fplot(@(t) 5*sin(t), @(t) 5*cos(t),'--k','LineWidth',3);
fplot(@(t) sin(t), @(t) cos(t),'--k','LineWidth',3);


