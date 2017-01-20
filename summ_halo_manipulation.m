%% Analyse halo manipulation study
% DK Shin
% 17/01/2017

clear all; close all; clc;

% Directory to data - should contain o_XXX.mat files and configuration file
dir_data='C:\Users\HE BEC\Documents\lab\halo_analysis\data\dist_halo\halo_manipulation\aniso_scaling\aniso_scaling_190117_1';

%% Load manipulation data
dir_orig=pwd;   % original directory
cd(dir_data);

% create dir to write output
mkdir('output');

% load configuration settings used
% TODO: set_config file curr req'd for 'g2 bin' properties - which should
% be in a separate config file?
config_S=dir('set_config*');    % get config file from filename token
run(config_S.name);         % load config vars

% get dimension of search parameters
config_search=dir('config_search_*');   % get config file for param search
run(config_search.name);    % load config vars: parlim_zxy_scale, pardiv_zxy_scale

% TODO: this should be set in a config file: currently "mani_search.m"
% script incorporates var setting which is NOT GOOD!
% parlim_zxy_scale=[0.9,1.1];     % param search limits for zxy-scaling
% pardiv_zxy_scale=7;             % linear divisions to search within lims

zxy_scale=linspace(parlim_zxy_scale(1),parlim_zxy_scale(2),pardiv_zxy_scale);   % TODO should be set in config
% param_dim=[7,7,7];   % dims: e_z,e_x,e_y  (TODO: length(zxy_scale*[1,1,1]) )
param_dim=pardiv_zxy_scale*ones(1,3);   % assuming all 3-param search space is equal in size

%% Load common config params
iCorr=1;    % corr analysis number performed

bin_dim=length(analysis.corr.lim{iCorr});  % get binning dims
bin_edge_tmp=cell(1,bin_dim);
bin_cent_tmp=cell(1,bin_dim);
for i=1:bin_dim
    % make bin edge and centre vectors
    bin_edge_tmp{i}=linspace(analysis.corr.lim{iCorr}{i}(1),...
        analysis.corr.lim{iCorr}{i}(2),analysis.corr.nBin{iCorr}(i)+1);
    bin_cent_tmp{i}=0.5*(bin_edge_tmp{i}(1:end-1)+bin_edge_tmp{i}(2:end));
end
drad=bin_cent_tmp{1};
dtheta=bin_cent_tmp{2};     % g2 delta theta vector
[DRAD,DTHETA]=meshgrid(drad,dtheta);    % create mesh grid for 3D-plots

%% Collate runs
% get data files
data_out=dir('o_*.mat');    % not sorted in iter num
iterN=length(data_out);

% Measures
g2_max=zeros(param_dim);

for i=1:iterN
    data_tmp=load(data_out(i).name);    % SLOW!
    
    % get param set for this run
    % Can use "token N" and param_dims OR use the scale_vect variable to place in correct bin
    %%% 1. Filename token N to determine bin location
    thisN=str2double(data_out(i).name(regexp(data_out(i).name,'\d')));
    kk=mod(thisN-1,pardiv_zxy_scale)+1;
    jj=mod(floor((thisN-1)/pardiv_zxy_scale),pardiv_zxy_scale)+1;
    ii=mod(floor((thisN-1)/(pardiv_zxy_scale^2)),pardiv_zxy_scale)+1;
    
    % evaluate normalised g2 in angular-radial
    g2=data_tmp.nshot*data_tmp.G2_single./data_tmp.G2_all;  % this is noisy, don't understand well why
    
    % get radially integrated g2 (angular g2) (avg along dim-1)
    G2_single_ang=(sum(data_tmp.G2_single,1)/size(data_tmp.G2_single,1));
    G2_all_ang=(sum(data_tmp.G2_all,1)/size(data_tmp.G2_all,1));
    g2_ang=data_tmp.nshot*G2_single_ang./G2_all_ang;    % get angular-g2 (1D)
    
    % store ang-g2 amplitude
    g2_max(ii,jj,kk)=max(g2_ang);   % max of integrated

    %% plot angular-g2 over all params searched
    fig_g2_1D=figure(11);
    hold on;
    plot(dtheta,g2_ang);
    if i==1
        title('Integrated $g^2$ from manipulated halo');
        xlabel('$\Delta\theta$');
        ylabel('$g^2$');
        xlim([min(dtheta),max(dtheta)]);
        box on;
        set(gca,'LooseInset',get(gca,'TightInset'));    % set figure tight
    end

    %% plot g2 measures
    fig_g2=figure(21);
    % G2 in each shot
    subplot(1,3,1);
    surf(DRAD',DTHETA',data_tmp.G2_single,'edgecolor','none');
    xlabel('$\delta k$'); ylabel('$\delta\theta$'); zlabel('$G^{(2)}_{corr}$');
    axis tight;
    shading interp;
    
    % G2 across uncorrelated runs (scaled)
    subplot(1,3,2);
    surf(DRAD',DTHETA',data_tmp.G2_all/data_tmp.nshot,'edgecolor','none');
    xlabel('$\delta k$'); ylabel('$\delta\theta$'); zlabel('$G^{(2)}_{uncorr}$');
    axis tight;
    shading interp;
    
    % normalised g2
    subplot(1,3,3);
    surf(DRAD',DTHETA',g2,'edgecolor','none');
%     title(strcat('ZXY-stretch = ',num2str([zxy_scale(ii),zxy_scale(jj),zxy_scale(kk)])));
    xlabel('$\delta k$'); ylabel('$\delta\theta$'); zlabel('$g^{(2)}$');
    axis tight;
    shading interp;
    
    % title
    suplabel(strcat('ZXY-stretch = ',num2str([zxy_scale(ii),zxy_scale(jj),zxy_scale(kk)])),'t');
    
    drawnow;
    
    saveas(fig_g2,['output\g2_',num2str(i),'.png']);
end

% Save figure
saveas(fig_g2_1D,'output\g2_angular.png');

%% Summary from collated
% Max g2 in eps_x and eps_y
eps_x=zxy_scale-1;
eps_y=eps_x;    % symmetric param space
[EPS_X,EPS_Y]=meshgrid(eps_x',eps_y');

for ii_plot=1:size(g2_max,1)
    H=figure(1);
    surf(EPS_X,EPS_Y,squeeze(g2_max(ii_plot,:,:)));
    title(sprintf('Max g2 landscape at Z-stretch = %.2f',zxy_scale(ii_plot)));
    xlabel('$\epsilon_x$');
    ylabel('$\epsilon_y$');
    zlabel('max($\tilde g^2$)');
    xlim([min(eps_x),max(eps_x)]);
    ylim([min(eps_y),max(eps_y)]);
    
    drawnow;
    saveas(H,['output\zz_',num2str(ii_plot),'.png']);
end