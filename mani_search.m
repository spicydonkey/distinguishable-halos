%% Iterate manipulation parameter space for g2 CORRELATION analysis
% DK Shin

clear all; close all; clc;

set_config;
configs=usrconfigs;     % copy to protect usrconfigs

% Manipulation search space
zxy_scale=linspace(0.9,1.1,7);    % test manip param space - scaling
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
            N=N+1;
            fprintf('N = %6d\t(%3d,%3d,%3d)\n',N,ii_1,ii_2,ii_3);
            
            % get manip params for this iteration
            scale_vect=[zxy_scale(ii_1),zxy_scale(ii_2),zxy_scale(ii_3)];   % zxy scaling
            
            % manipulate halo 2 (mj=1)
            for jj=1:nshot
                halo_manip{jj,2}=bsxfun(@times,halok{jj,2},scale_vect);  % anisotropic scaling    
            end
            
            [G2_single,G2_all]=G2_angular(halo_manip,bin_edge_tmp,verbose);
%             g2=nshot*G2_single./G2_all;    % don't really need g2 (since
%             it can be calculated from unnormalised G2's?
            
            %% save data
            % TODO: shouldn't need nshot - should be collected from some
            % kind of summary of original halo data (since manipulation
            % keeps all shots)
            % TODO: for same reason, wouldn't need scale_vect as well since
            % it will be determined from file number token "N"
            save([configs.files.archive,'\o_',num2str(N),'.mat'],'G2_single','G2_all','scale_vect','nshot');
        end
    end
end