function corr_out=halo_g2_manager(zxy0,configs,verbose)
% Correlation task manager

% input check
if ~exist('verbose','var')
    warning('verbose is not provided - setting to quiet (0)');
    verbose=0;
end

if verbose>0, fprintf('Starting correlation analysis...\n'), end;
%% MAIN
t_fun_start=tic;

% get paths
dir_output=configs.files.dirout;

%% Do tasked correlation analysis
for iCorr=1:length(configs.corr.type)
    % Set up bins
    bin_dim=length(configs.corr.lim{iCorr});  % get binning dims
    bin_edge_tmp=cell(1,bin_dim);
    bin_cent_tmp=cell(1,bin_dim);
    for i=1:bin_dim
        % make bin edge and centre vectors
        bin_edge_tmp{i}=linspace(configs.corr.lim{iCorr}{i}(1),...
            configs.corr.lim{iCorr}{i}(2),configs.corr.nBin{iCorr}(i)+1);
        bin_cent_tmp{i}=0.5*(bin_edge_tmp{i}(1:end-1)+bin_edge_tmp{i}(2:end));
    end
    
    % Evaluate G2 correlation
    % TODO - insert debug code for G2 and zxy manipulator
    [G2_shot_tmp,G2_all_tmp]=G2_caller(zxy0(:,configs.corr.type{iCorr}.comp),...
        bin_edge_tmp,configs.corr.type{iCorr}.coord,configs.corr.type{iCorr}.opt,verbose);
    g2_tmp=size(zxy0,1)*G2_shot_tmp./G2_all_tmp;      % normalised g2
    
    % Get results
    corr_out.bEdge{iCorr}=bin_edge_tmp;
    corr_out.bCent{iCorr}=bin_cent_tmp;
    corr_out.G2shot{iCorr}=G2_shot_tmp;
    corr_out.G2all{iCorr}=G2_all_tmp;
    corr_out.g2{iCorr}=g2_tmp;
end
% clear workspace
clear bin_dim bin_edge_tmp bin_cent_tmp G2_shot_tmp G2_all_tmp g2_tmp;

%% Plot the original g2 correlation function
if configs.flags.graphics
    for iCorr=1:length(configs.corr.type)
        nfig_tmp=10+iCorr;  % g2 figures start from figure 11
        hfig=plotCorr(nfig_tmp,corr_out,configs.corr,iCorr);
        
        % save figs
        fname_str=['corr_',num2str(iCorr)];
        saveas(hfig,fullfile(dir_output,[fname_str,'.fig']));
        saveas(hfig,fullfile(dir_output,[fname_str,'.png']));
    end
end

%% 1D correlation profile and Gaussian fit
for iCorr=1:length(configs.corr.type)
    if configs.flags.graphics   % plotting
        nfig_tmp=20+iCorr;  % 1D starts from fig 21
        hfig=figure(nfig_tmp);
        ax=gca;
    end
    
    this_corr_type=configs.corr.type{iCorr}.coord;
    % Get integrated or sliced 1D correlation profile
    if isequal(this_corr_type,'angular')
        % TODO: correlations - dk integrated (may reduce peak height)
        g2_1d_tmp=size(zxy0,1)*sum(corr_out.G2shot{iCorr},1)./sum(corr_out.G2all{iCorr},1);
        
        % Gaussian fit (angular G2 allows fitting both BB,CL)
        % BB fit
        param0=[4,pi,0.1,1];     % fit estimate [amp,mu,sigma,offset]
        [fitparam_tmp,fit_g2_tmp]=gaussfit(corr_out.bCent{iCorr}{2},g2_1d_tmp,param0,verbose);
        
        % CL fit
        param0=[2,0,0.1,1];
        [fitparam_tmpCL,fit_g2_tmpCL]=gaussfit(corr_out.bCent{iCorr}{2},g2_1d_tmp,param0,verbose);
        
        % Plot
        if configs.flags.graphics   % plotting
            % calculated 1D correlation profile
            plot(corr_out.bCent{iCorr}{2},g2_1d_tmp,'*');
            
            hold(ax,'on');
            title_str=['$\Delta k$-integrated, ',...
                '(',num2str(configs.corr.type{iCorr}.comp),') halos'];
            title(title_str);
            xlabel('$\Delta\theta$'); ylabel('$\bar{g}^{(2)}$');
            xlim([0,pi]); ylim auto;
            
            % plot the fits
            plot(ax,fit_g2_tmp.x,fit_g2_tmp.y,'r');     
            plot(ax,fit_g2_tmpCL.x,fit_g2_tmpCL.y,'b--');
            hold(ax,'off');
        end
        
    elseif isequal(this_corr_type,'cart')
        % TODO - do in X/Y?
        % Get line through Z-axis
        ind_zero_tmp=round((configs.corr.nBin{iCorr}+1)/2);    % zero-cent'd bin index for sampling 3D-g2
        
        % Gaussian fit
        % Set initial params
        if strcmp(configs.corr.type{iCorr}.opt,'BB')
            param0=[4,0,0.1,1];     % fit estimate [amp,mu,sigma,offset]
        else
            param0=[2,0,0.1,1];     % CL peaks at 2
        end
        [fitparam_tmp,fit_g2_tmp]=gaussfit(corr_out.bCent{iCorr}{1},...
            corr_out.g2{iCorr}(:,ind_zero_tmp(2),ind_zero_tmp(3)),param0,verbose);
        
        % Plot
        if configs.flags.graphics   % plotting
            plot(corr_out.bCent{iCorr}{1},...
                corr_out.g2{iCorr}(:,ind_zero_tmp(2),ind_zero_tmp(3)),'*');
            hold(ax,'on');
            
            title_str=['(',num2str(configs.corr.type{iCorr}.comp),') halos, ',...
                configs.corr.type{iCorr}.opt,', ','$Z$-axis'];
            title(title_str);
            xlabel('$\Delta K_z$'); ylabel(['$g^{(2)}_{',configs.corr.type{iCorr}.opt,'}$']);
            
            plot(ax,fit_g2_tmp.x,fit_g2_tmp.y,'r');     % plot the fit
            hold(ax,'off');
        end
    else
        warning('SOMETHING IS WRONG!');
    end
    
    % Get fit params
    corr_out.fit{iCorr}=fitparam_tmp;
    
    % Save figs
    if configs.flags.graphics   % plotting
        fname_str=['corr1d_',num2str(iCorr)];
        saveas(hfig,fullfile(dir_output,[fname_str,'.fig']));
        saveas(hfig,fullfile(dir_output,[fname_str,'.png']));
    end
end
clear g2_1d_tmp param0 fitparam_tmp fit_g2_tmp ax this_corr_type;

%% END
t_fun_end=toc(t_fun_start);   % end of code
if verbose>0
    disp('-----------------------------------------------');
    fprintf('Total elapsed time for %s (s): %7.1f\n','corrTaskManager',t_fun_end);
    disp('-----------------------------------------------');
end