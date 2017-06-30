function plotCorrCart(CORRDATA,CORRCONFIG)
% Plot correlation function in Cartesian coord

% create domain grid for surf
% TODO - only plots Z axis atm
[dX,dY]=meshgrid(CORRDATA.bCent{2},CORRDATA.bCent{3});

% TODO for asymetric binning? - is this reliable?
ind_zero=round((CORRCONFIG.nBin+1)/2); % zero-cent'd bin index for sampling 3D-g2

subplot(1,3,1);
surf(dX',dY',squeeze(CORRDATA.G2_corr(ind_zero(1),:,:)),'edgecolor','none');

title_str=['(',num2str(CORRCONFIG.type.comp),') halos,',CORRCONFIG.type.coord,',shots'];
title(title_str);
xlabel('$\delta k_i$'); ylabel('$\delta k_j$'); 
zlabel(['$G^{(2)}_{',CORRCONFIG.type.opt,'(',num2str(CORRCONFIG.type.comp),')}$']);
axis tight;
shading interp;

subplot(1,3,2);
surf(dX',dY',squeeze(CORRDATA.G2_uncorr(ind_zero(1),:,:)),'edgecolor','none');

title_str=['(',num2str(CORRCONFIG.type.comp),') halos,',CORRCONFIG.type.coord,',all'];
title(title_str);
xlabel('$\delta k_i$'); ylabel('$\delta k_j$'); 
zlabel(['$G^{(2)}_{',CORRCONFIG.type.opt,'(',num2str(CORRCONFIG.type.comp),')}$']);
axis tight;
shading interp;

subplot(1,3,3);
surf(dX',dY',squeeze(CORRDATA.g2(ind_zero(1),:,:)),'edgecolor','none');

title_str=['(',num2str(CORRCONFIG.type.comp),') halos,',CORRCONFIG.type.coord,',normalised'];
title(title_str);
xlabel('$\delta k_i$'); ylabel('$\delta k_j$'); 
zlabel(['$g^{(2)}_{',CORRCONFIG.type.opt,'(',num2str(CORRCONFIG.type.comp),')}$']);
axis tight;
shading interp;