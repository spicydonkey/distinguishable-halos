function [HFIG]=plotCorrAngle(FIGNUM,CORRDATA,CORRCONFIG,INDEX)
% Plot correlation function in angular configuration

HFIG=figure(FIGNUM);

% create domain grid for surf
[drad,dtheta]=meshgrid(CORRDATA.bCent{INDEX}{1},CORRDATA.bCent{INDEX}{2});

subplot(1,3,1);
surf(drad',dtheta',CORRDATA.G2shot{INDEX},'edgecolor','none');
title_str=['(',num2str(CORRCONFIG.type{INDEX}.comp),') halos,',CORRCONFIG.type{INDEX}.coord,',shots'];
title(title_str);
xlabel('$\delta k$'); ylabel('$\delta\theta$'); zlabel(['$G^{(2)}_{(',num2str(CORRCONFIG.type{INDEX}.comp),')}$']);
axis tight;
shading interp;

subplot(1,3,2);
surf(drad',dtheta',CORRDATA.G2all{INDEX},'edgecolor','none');
title_str=['(',num2str(CORRCONFIG.type{INDEX}.comp),') halos,',CORRCONFIG.type{INDEX}.coord,',all'];
title(title_str);
xlabel('$\delta k$'); ylabel('$\delta\theta$'); zlabel(['$G^{(2)}_{(',num2str(CORRCONFIG.type{INDEX}.comp),')}$']);
axis tight;
shading interp;

subplot(1,3,3);
surf(drad',dtheta',CORRDATA.g2{INDEX},'edgecolor','none');
title_str=['(',num2str(CORRCONFIG.type{INDEX}.comp),') halos,',CORRCONFIG.type{INDEX}.coord,',normalised'];
title(title_str);
xlabel('$\delta k$'); ylabel('$\delta\theta$'); zlabel(['$g^{(2)}_{(',num2str(CORRCONFIG.type{INDEX}.comp),')}$']);
axis tight;
shading interp;