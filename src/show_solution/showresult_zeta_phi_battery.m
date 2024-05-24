function showresult_zeta_phi_battery(nodezetaphi,elemzetaphi,zetah,phih)
%showresult displays the mesh and the solution
%
% Copyright (C) Terence Yu.

clf;  % clear figures

set(gcf,'Units','normal');
%set(gcf,'Position',[0.25,0.25,0.7,0.25]);

%% Plot numerical solution

nexttile;
showsolution(nodezetaphi,elemzetaphi,vecnorm(zetah,2,2),'$|\mbox{\boldmath$\zeta$}_{h}|$');

%% Plot numerical solution phi
nexttile;
showsolution(nodezetaphi,elemzetaphi,phih,'$\varphi_h$');