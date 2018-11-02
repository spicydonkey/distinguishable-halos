% SUMMARY OF BB g2 amplitude and scattering halo mode occupancy
%
% DKS
% 2018-06-04

flag_save_plot=0;


%% Figure properties presets
dir_temp=fullfile(userpath,'temp');

path_base=fileparts(mfilename('fullpath'));

figname = 'src_g2_nsc';
path_save=fullfile(path_base,'out');

datetimestr=datestr(datetime,'yyyymmdd_HHMMSS');    % timestamp when function called
figname=sprintf('%s_%s',figname,datetimestr);   % attach timestamp to name

fontsize_normal=12;
fontsize=fontsize_normal;
fontsize_small=10;
fontsize_large=13;
linewidth = 1.2;
markersize = 5;


paperunits = 'centimeters';
papersize=[8 8];
paperposition=[0,0,papersize];

plotboxaspectratio=[1,1,1];
boxlinewidth=1;
ticklength=[0.02,0.025];

%%% MISC
% gray_col=0.5*ones(1,3);         % gray data points
namearray={'LineWidth','MarkerFaceColor','Color'};      % error bar graphics properties
% valarray={linewidth,'w','k'};                 % 90 deg (normal) data
% valarray_30={linewidth,gray_col,gray_col};    % 30 deg (slow) data

% [c1,c2]=palette(2);
% valarray={linewidth,c2(1,:),c1(1,:)};
% valarray_30={linewidth,c2(2,:),c1(2,:)};

[c1,c2]=palette(10);        % 8th color is is purplish
col_idx1=8;     
col_idx2=8;
valarray_90={linewidth,c2(col_idx1,:),c1(col_idx1,:)};
valarray_30={linewidth,c1(col_idx2,:),c1(col_idx2,:)};

patch_col=0.9*ones(1,3);


%% CONFIGS
% plotting limits
n_lim=[2e-3,1e0];       % axes lims for mode occupancy
g2_lim=[8e-1,1e2];      % axes lims for g2-1

% condition for violoation of Bell's test
g2_bell=(2*sqrt(2)+3)-1;            % Wasak, Chwedenczuk - limit g2 - 1

%% DATA
%%% 90 deg data
% n_exp=[0.0085, 0.0126, 0.0575, 0.604, 0.0887, 0.0102, 0.273 0.052];
% g2_exp=[92.5, 23.3, 31.5, 3.6, 20.3, 97.3, 8 26.6];      % mean fit amplitude
% 
% n_unc=[0.68, 1.1, 0.38, 0.24, 0.39, 0.65, 0.32 0.84];       % relative error in ratio to mean
% g2_err=[9.1, 1.2, 2.1, 0.15, 0.5, 8.5, 0.43 1.5];           % error in SE


path_dir='C:\Users\HE BEC\Documents\lab\bell_momentumspin\pairsource\90';
file_name='pairsource_collated_del_run2.mat';

vars2load={'g2_amp','g2_amp_sdev','n_mocc','n_mocc_err'};

% load from revised analysis
% S=load('C:\Users\HE BEC\bell_data\pairsource\90\pairsource_collated_del_run2.mat',vars2load{:});
file_path=fullfile(path_dir,file_name);
S=load(file_path,vars2load{:});


g2_exp=S.g2_amp;
g2_err=S.g2_amp_sdev;
n_exp=S.n_mocc;
n_unc=S.n_mocc_err;


%%% 30 deg data
n_exp_30=[0.196 0.462 0.774 0.299];
g2_exp_30=[7.6 3.6 2.73 5.9];
    
n_unc_30=[0.44 0.38 0.24 0.4];
g2_err_30=[1.0 0.16 0.06 0.41];


% Use SD for g2 uncertainty
% g2_err=g2_err*sqrt(5);
g2_err_30=g2_err_30*sqrt(5);



%% EVALUATE error limits
% relative error to abs
n_err=n_exp.*n_unc;      % calculated
n_err_30=n_exp_30.*n_unc_30;

% error bounds
ndata=numel(n_exp);         % number of data points
n_err_bnd=cell(1,2);
g2_err_bnd=cell(1,2);

for ii=1:2
    % build min/max bound
    n_err_bnd{ii}=n_exp+(-1)^(ii)*n_err;
    g2_err_bnd{ii}=g2_exp+(-1)^(ii)*g2_err;
    g2_err_bnd{ii}=g2_err_bnd{ii}-1;            % shift to g2-1
end

% cut bound at axes lims (MIN only for log plot)
for ii=1:ndata
    if n_err_bnd{1}(ii)<n_lim(1)
        n_err_bnd{1}(ii)=n_lim(1);
    end
    if g2_err_bnd{1}(ii)<g2_lim(1)
        g2_err_bnd{1}(ii)=g2_lim(1);
    end
end


%% Theory curve
%   Model: spontaneous pair source

r=3;
n=linspace(min(n_exp)/r,r*max(n_exp),1000);       % logs
% n=linspace(0,1.1*max(n_exp),1000);       % linear
% g2_theory=1+1./n;
g2_theory=2+1./n;

%% Plot (g2-1)
% plot mode occupancy vs g2(0) - 1
fig=figure('Name',figname,...
    'Units','centimeters',...
    'PaperUnits',paperunits,...
    'PaperPositionMode','manual',...
    'PaperSize',papersize,...
    'PaperPosition',paperposition);
hold on;
fig.Renderer='painters';
htheory=plot(n,g2_theory-1,'k--','LineWidth',1.5,'DisplayName','Theory');       % theory curve

%%% 30 deg data
hdata_30=ploterr(n_exp_30,g2_exp_30-1,n_err_30,g2_err_30,'d','logxy','hhxy',0);
% set(hdata_30(1),namearray,valarray_30,'MarkerSize',markersize,'DisplayName','Slow collision ($30^\circ$)');              % DATAPOINT
set(hdata_30(1),namearray,valarray_30,'MarkerSize',markersize,'DisplayName','$30^\circ$');              % DATAPOINT
set(hdata_30(2),namearray,valarray_30,'DisplayName','');                  % Y-err
set(hdata_30(3),namearray,valarray_30,'DisplayName','');                  % X-err

%%% 90 deg data
hdata=ploterr(n_exp,g2_exp-1,n_err_bnd,g2_err_bnd,'o','logxy','hhxy',0);
% set(hdata(1),namearray,valarray,'MarkerSize',markersize,'DisplayName','Fast collision ($90^\circ$)');              % DATAPOINT
set(hdata(1),namearray,valarray_90,'MarkerSize',markersize,'DisplayName','$90^\circ$');              % DATAPOINT
set(hdata(2),namearray,valarray_90,'DisplayName','');                  % Y-err
set(hdata(3),namearray,valarray_90,'DisplayName','');                  % X-err

%%% Bell's test regime
% build patch corners
x_corn=[n_lim(1),n_lim(2),n_lim(2),n_lim(1)];
% y_corn=[g2_lim(1),g2_lim(1),g2_bell,g2_bell];     % below condition
y_corn=[g2_bell,g2_bell,g2_lim(2),g2_lim(2)];       % above condition

p_bell=patch(x_corn,y_corn,...
    patch_col,'EdgeColor','none');
%     'b','FaceAlpha',0.1,'EdgeColor','none');

% label regions
% text(5e-3,3,'$\mathcal{B} < \frac{1}{\sqrt{2}}$','FontSize',fontsize_normal);
% text(5e-3,7,'$\textrm{max}\,\mathcal{B} > \frac{1}{\sqrt{2}}$','FontSize',fontsize_normal);
text(3e-3,3.5,'Local (QM)','FontSize',fontsize_small);
text(3e-3,6,'Nonlocal (QM)','FontSize',fontsize_small);

uistack(p_bell,'bottom');   % move the patch object to bottom


% annotate plot
box on;
% oleg=legend([hdata(1),hdata_30(1),htheory]);
oleg=legend([hdata(1),hdata_30(1)]);
title(legend,'Beam angle');
set(oleg,'FontSize',fontsize_small);

ax=gca;
ax.XScale='log';
ax.YScale='log';
ax.XTick=10.^[-2,-1,0];
ax.YTick=10.^[0,1,2];

set(gca,...
    'FontSize',fontsize,...
    'LineWidth',boxlinewidth,...
    'TickLength',ticklength,...
    'PlotBoxAspectRatio',plotboxaspectratio);

set(gca,'Layer','Top');     % graphics axes should be always on top

xlim(n_lim);
ylim(g2_lim);

% xlabel('$n_\textrm{sc}$');          % average atom occupation in an arbitrary (check) scattering mode in halo
xlabel('Mode occupancy $n$');
% ylabel('$g^{(2)}_{\mathrm{BB},\uparrow\downarrow}-1$');
ylabel('$g^{(2)}_{\uparrow\downarrow}-1$');


%%% FIG postprocess
set(ax,'FontSize',fontsize_normal);

fig.Units = paperunits;
fig.Position = paperposition;

fig.PaperSize = papersize;
fig.PaperUnits = paperunits;
fig.PaperPosition = paperposition;

% save fig to temp directory
% saveas(fig, [figname, '.eps'], 'psc2');     
% saveas(fig, [figname, '.pdf']);

% print to temp directory
% print(fig,fullfile(dir_temp,strcat(figname,'.pdf')),'-dpdf');


if flag_save_plot
    if ~exist(path_save,'dir')
        warning('Output direcotry %s does not exist. Creating it!',path_save);
        mkdir(path_save);
    end
    print(fig,fullfile(path_save,figname),'-dpdf');
end