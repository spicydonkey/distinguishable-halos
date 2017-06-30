function plotCorrAngle(CORRDATA,CORRCONFIG)
% Plot correlation function in angular configuration

% create domain grid for surf
[drad,dtheta]=meshgrid(CORRDATA.bCent{1},CORRDATA.bCent{2});

subplot(1,3,1);
surf(drad',dtheta',CORRDATA.G2_corr,'edgecolor','none');
title_str=['(',num2str(CORRCONFIG.type.comp),') halos,',CORRCONFIG.type.coord,',shots'];
title(title_str);
xlabel('$\delta k$'); ylabel('$\delta\theta$'); zlabel(['$G^{(2)}_{(',num2str(CORRCONFIG.type.comp),')}$']);
axis tight;
shading interp;

subplot(1,3,2);
surf(drad',dtheta',CORRDATA.G2_uncorr,'edgecolor','none');
title_str=['(',num2str(CORRCONFIG.type.comp),') halos,',CORRCONFIG.type.coord,',all'];
title(title_str);
xlabel('$\delta k$'); ylabel('$\delta\theta$'); zlabel(['$G^{(2)}_{(',num2str(CORRCONFIG.type.comp),')}$']);
axis tight;
shading interp;

subplot(1,3,3);
surf(drad',dtheta',CORRDATA.g2,'edgecolor','none');
title_str=['(',num2str(CORRCONFIG.type.comp),') halos,',CORRCONFIG.type.coord,',normalised'];
title(title_str);
xlabel('$\delta k$'); ylabel('$\delta\theta$'); zlabel(['$g^{(2)}_{(',num2str(CORRCONFIG.type.comp),')}$']);
axis tight;
shading interp;