%% Iterate manipulation parameter space for g2 CORRELATION analysis
% DK Shin

clear all; close all; clc;

set_config;
configs=usrconfigs;     % copy to protect usrconfigs

% Manipulation search space
zxy_scale=linspace(0.9,1.1,3);    % test manip param space - scaling
% zxy_scale=linspace(1,1,1); % test

%% load halo data
S_temp=load(configs.files.saveddata,'zxy');
halok=S_temp.zxy;   % load the corr-ready data (isoscaled to unit rad)
clear S_temp;

%% CORR ANALYSIS
% Set up bins
iCorr=1;

bin_dim=length(analysis.corr.lim{iCorr});  % get binning dims
bin_edge_tmp=cell(1,bin_dim);
bin_cent_tmp=cell(1,bin_dim);
for i=1:bin_dim
    % make bin edge and centre vectors
    bin_edge_tmp{i}=linspace(analysis.corr.lim{iCorr}{i}(1),...
        analysis.corr.lim{iCorr}{i}(2),analysis.corr.nBin{iCorr}(i)+1);
    bin_cent_tmp{i}=0.5*(bin_edge_tmp{i}(1:end-1)+bin_edge_tmp{i}(2:end));
end

halo_manip=halok;       % manipulated halo cell
nshot=size(halok,1);    % num shots
N=0;
for ii_1=1:length(zxy_scale)
    for ii_2=1:length(zxy_scale)
        for ii_3=1:length(zxy_scale)
            % iter N
            N=N+1
            
            % get manip params for this iteration
            scale_vect=[zxy_scale(ii_1),zxy_scale(ii_2),zxy_scale(ii_3)];   % zxy scaling
            
            % manipulate halo 2 (mj=1)
            for jj=1:nshot
                halo_manip{jj,2}=bsxfun(@times,halok{jj,2},scale_vect);  % anisotropic scaling    
            end
            
            [G2_single,G2_all]=G2_angular(halo_manip,bin_edge_tmp,verbose);
            g2=nshot*G2_single./G2_all;
            
            %% g2 plot
%             [drad,dtheta]=meshgrid(bin_cent_tmp{1},bin_cent_tmp{2});
%             figure(1);
%             surf(drad',dtheta',g2,'edgecolor','none');
%             title_str=['Iter num: ',num2str(N),', scale=',num2str(scale_vect)];
%             title(title_str);
%             xlabel('$\delta k$'); ylabel('$\delta\theta$'); zlabel('$g^{(2)}_{(1,2)}$');
%             axis tight;
%             shading interp;
%             
%             saveas(gcf,[configs.files.archive,'\o_',num2str(N),'.png']);
% %             saveas(gcf,[configs.files.archive,'\o_test',num2str(N),'.png']);    % test
            
            %% save data
            save([configs.files.archive,'\o_',num2str(N),'.mat'],'G2_single','G2_all','g2','scale_vect');
%             save([configs.files.archive,'\o_test',num2str(N),'.mat'],'G2_single','G2_all','g2','scale_vect');    % test
        end
    end
end