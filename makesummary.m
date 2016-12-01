% MAKE summary
close all; 

dotsize=1;
colors='br';

S=load(configs.files.saveddata);
halo=S.halo;
plot_zxy(1,halo.zxy,dotsize,colors);
plot_zxy(2,halo.zxy0,dotsize,colors);