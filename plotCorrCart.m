function [HFIG]=plotCorrCart(FIGNUM,CORRDATA,CORRCONFIG,INDEX)
% Plot correlation function in Cartesian coord

HFIG=figure(FIGNUM);

% create domain grid for surf
% TODO - only plots Z axis atm
[dX,dY]=meshgrid(CORRDATA.bCent{INDEX}{2},CORRDATA.bCent{INDEX}{3});

% TODO for asymetric binning? - is this reliable?
ind_zero=round((CORRCONFIG.nBin{INDEX}+1)/2); % zero-cent'd bin index for sampling 3D-g2

subplot(1,3,1);
surf(dX',dY',squeeze(CORRDATA.G2shot{INDEX}(ind_zero(1),:,:)),'edgecolor','none');

title_str=['(',num2str(CORRCONFIG.type{INDEX}.comp),') halos,',CORRCONFIG.type{INDEX}.coord,',shots'];
title(title_str);
xlabel('$\delta k_i$'); ylabel('$\delta k_j$'); 
zlabel(['$G^{(2)}_{BB(',num2str(CORRCONFIG.type{INDEX}.comp),')}$']);
axis tight;
shading interp;

subplot(1,3,2);
surf(dX',dY',squeeze(CORRDATA.G2all{INDEX}(ind_zero(1),:,:)),'edgecolor','none');

title_str=['(',num2str(CORRCONFIG.type{INDEX}.comp),') halos,',CORRCONFIG.type{INDEX}.coord,',all'];
title(title_str);
xlabel('$\delta k_i$'); ylabel('$\delta k_j$'); 
zlabel(['$G^{(2)}_{BB(',num2str(CORRCONFIG.type{INDEX}.comp),')}$']);
axis tight;
shading interp;

subplot(1,3,3);
surf(dX',dY',squeeze(CORRDATA.g2{INDEX}(ind_zero(1),:,:)),'edgecolor','none');

title_str=['(',num2str(CORRCONFIG.type{INDEX}.comp),') halos,',CORRCONFIG.type{INDEX}.coord,',normalised'];
title(title_str);
xlabel('$\delta k_i$'); ylabel('$\delta k_j$'); 
zlabel(['$g^{(2)}_{BB(',num2str(CORRCONFIG.type{INDEX}.comp),')}$']);
axis tight;
shading interp;