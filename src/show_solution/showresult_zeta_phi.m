function showresult_zeta_phi(nodezetaphi,elemzetaphi,zeta,zetah,phi,phih)

%showresult displays the mesh and the solution
%
% Copyright (C) Terence Yu.

clf;  % clear figures

set(gcf,'Units','normal');
%set(gcf,'Position',[0.25,0.25,0.7,0.25]);

tiledlayout(1,2)

%% Plot exact solution first component

% nexttile;
% zeta_exact = zeta(node); zeta_exact = zeta_exact(:,1);
% showsolution(nodezetaphi,elemzetaphi,vecnorm(zeta_exact,2,2),'$|\mathbf{\zeta}_{h}|$');


%% Plot numerical solution

nexttile;
showsolution(nodezetaphi,elemzetaphi,vecnorm(zetah,2,2),'$|\mbox{\boldmath$\zeta$}_{h}|$');

%% Plot exact solution phi

%nexttile;
%phi_exact = phi(node); phi_exact = phi_exact(:,1);
%showsolution(node,elem,phi_exact(1:size(node,1)),'$\varphi$');

%% Plot numerical solution phi
nexttile;
showsolution(nodezetaphi,elemzetaphi,phih,'$\varphi_h$');