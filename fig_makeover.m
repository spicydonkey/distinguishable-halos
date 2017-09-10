%% generate pretty figures
%% g2 plots
f=gcf;
ax=gca;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% USER configure %%%%%%%%%%%%%%%%%%%%%%%%%
savefig=1;
path_save='C:\Users\HE BEC\Desktop\dist_halo_summary\g2_plots';

datetimestr=datestr(datetime,'yyyymmdd_HHMMSS');    % timestamp when function called

% mrk_size=7;
fontsize=10;
fontsize_small=7;
% linewidth=1.8;

x_lim=[-0.22,0.22];
y_lim=[-5,110];

x_tick=[-0.3:0.1:0.3];
y_tick=[0:50:150];

plotboxaspectratio=[1,1,1];
boxlinewidth=1;
ticklength=[0.02,0.025];

papersize=[8,8];
paperposition=[0,0,papersize];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% configure the figure
ax.Title.String='';     % turn title off

f.Units='centimeters';
f.PaperUnits='centimeters';
f.PaperPositionMode='manual';
f.PaperSize=papersize;
f.PaperPosition=paperposition;

set(ax,'Units','normalized',...
    'XLim',x_lim,...
    'YLim',y_lim,...
    'XTick',x_tick,...
    'YTick',y_tick,...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',fontsize,...
    'PlotBoxAspectRatio',plotboxaspectratio,...
    'LineWidth',boxlinewidth,...
    'TickLength',ticklength);

% save
if savefig
    figname=sprintf('temp_g2_%s',datetimestr);
    print(f,fullfile(path_save,figname),'-dpdf');
end

%% rotating view of S-G separated halos
% savemov=0;
% path_save='C:\Users\HE BEC\Desktop\dist_halo_summary\SG_separated';
% 
% % have figure ready
% % f=gcf;
% ax=gca;
% 
% %%%%%%%%%%%%%%%%%%
% Nframe=20;
% 
% set(ax,'nextplot','replacechildren','visible','off');
% fr=getframe;
% [im,map]=rgb2ind(fr.cdata,256,'nodither');
% im(1,1,1,Nframe)=0;
% 
% az=linspace(0,2*pi,Nframe);
% el=10;
% 
% for k=1:Nframe
%     % do a transform
%     view([az(k),el]);
%     
%     % next frame
%     fr=getframe;
%     im(:,:,1,k)=rgb2ind(fr.cdata,256,'nodither');
% end
% if savemov
%     figname=sprintf('temp_mov_%s',datetimestr);
%     imwrite(im,map,[fullfile(path_save,figname),'.gif'],'DelayTime',0,'LoopCount',inf);
% end