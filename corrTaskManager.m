function [result]=corrTaskManager(analysis,configs,VERBOSE)
% Correlation task manager

% input check
if ~exist('VERBOSE','var')
    warning('VERBOSE is not provided - setting to quiet (0)');
    VERBOSE=0;
end

vars_save={'analysis','result'};

if VERBOSE>0, fprintf('Starting correlation analysis...\n'), end;
%% MAIN
t_fun_start=tic;

% get paths
dir_output=configs.files.dirout;

% load counts for correlation analysis
S_temp=load(configs.files.saveddata,'zxy');
zxy=cell(size(S_temp.zxy));    % preallocate size for MATLAB (not sure if this helps)
zxy=S_temp.zxy;   % load the corr-ready data (transformed)
clear S_temp;

% "result" should not be defined previously in data file (will be overwritten)

%% Do tasked correlation analysis
for iCorr=1:length(analysis.corr.type)
    % Set up bins
    bin_dim=length(analysis.corr.lim{iCorr});  % get binning dims
    bin_edge_tmp=cell(1,bin_dim);
    bin_cent_tmp=cell(1,bin_dim);
    for i=1:bin_dim
        % make bin edge and centre vectors
        bin_edge_tmp{i}=linspace(analysis.corr.lim{iCorr}{i}(1),...
            analysis.corr.lim{iCorr}{i}(2),analysis.corr.nBin{iCorr}(i)+1);
        bin_cent_tmp{i}=0.5*(bin_edge_tmp{i}(1:end-1)+bin_edge_tmp{i}(2:end));
    end
    
    % Evaluate G2 correlation
    % TODO - insert debug code for G2 and zxy manipulator
    [G2_shot_tmp,G2_all_tmp]=G2_caller(zxy(:,analysis.corr.type{iCorr}.comp),...
        bin_edge_tmp,analysis.corr.type{iCorr}.coord,analysis.corr.type{iCorr}.opt,VERBOSE);
    g2_tmp=size(zxy,1)*G2_shot_tmp./G2_all_tmp;      % normalised g2
    
    % Get results
    result.corr.bEdge{iCorr}=bin_edge_tmp;
    result.corr.bCent{iCorr}=bin_cent_tmp;
    result.corr.G2shot{iCorr}=G2_shot_tmp;
    result.corr.G2all{iCorr}=G2_all_tmp;
    result.corr.g2{iCorr}=g2_tmp;
end
% clear workspace
clear bin_dim bin_edge_tmp bin_cent_tmp G2_shot_tmp G2_all_tmp g2_tmp;

%% Plot the original g2 correlation function
for iCorr=1:length(analysis.corr.type)
    nfig_tmp=10+iCorr;  % g2 figures start from figure 11
    hfig=plotCorr(nfig_tmp,result.corr,analysis.corr,iCorr);
    
    % save figs
    fname_str=['corr_',num2str(iCorr)];
    saveas(hfig,[dir_output,fname_str,'.fig']);
    saveas(hfig,[dir_output,fname_str,'.png']);
end

%% 1D correlation profile and Gaussian fit
for iCorr=1:length(analysis.corr.type)
    nfig_tmp=20+iCorr;  % 1D starts from fig 21
    hfig=figure(nfig_tmp);
    ax=gca;
    
    this_corr_type=analysis.corr.type{iCorr}.coord;
    % Get integrated or sliced 1D correlation profile
    if isequal(this_corr_type,'angular')
        % integrate dk
        g2_1d_tmp=size(zxy,1)*sum(result.corr.G2shot{iCorr},1)./sum(result.corr.G2all{iCorr},1);
        
        plot(result.corr.bCent{iCorr}{2},g2_1d_tmp,'*');
        
        hold(ax,'on');
        title_str=['$\Delta k$-integrated, ',...
            '(',num2str(analysis.corr.type{iCorr}.comp),') halos'];
        title(title_str);
        xlabel('$\Delta\theta$'); ylabel('$\bar{g}^{(2)}$');
        xlim([0,pi]); ylim auto;
        
        % Gaussian fit
        param0=[4,pi,0.1,1];     % fit estimate [amp,mu,sigma,offset]
        [fitparam_tmp,fit_g2_tmp]=gaussfit(result.corr.bCent{iCorr}{2},g2_1d_tmp,param0,VERBOSE);
        plot(ax,fit_g2_tmp.x,fit_g2_tmp.y,'r');     % plot the fit
        hold(ax,'off');
        
    elseif isequal(this_corr_type,'cart')
        % TODO - do in X/Y?
        % Get line through Z-axis
        ind_zero_tmp=round((analysis.corr.nBin{iCorr}+1)/2);    % zero-cent'd bin index for sampling 3D-g2
        plot(result.corr.bCent{iCorr}{1},result.corr.g2{iCorr}(:,ind_zero_tmp(2),ind_zero_tmp(3)),'*');
        
        hold(ax,'on');
        title_str=['(',num2str(analysis.corr.type{iCorr}.comp),') halos, ',...
            analysis.corr.type{iCorr}.opt,', ','$Z$-axis'];
        title(title_str);
        xlabel('$\Delta K_z$'); ylabel(['$g^{(2)}_{',analysis.corr.type{iCorr}.opt,'}$']);
        
        % Gaussian fit
        param0=[4,0,0.1,1];     % fit estimate [amp,mu,sigma,offset]
        [fitparam_tmp,fit_g2_tmp]=gaussfit(result.corr.bCent{iCorr}{1},result.corr.g2{iCorr}(:,ind_zero_tmp(2),ind_zero_tmp(3)),param0,VERBOSE);
        plot(ax,fit_g2_tmp.x,fit_g2_tmp.y,'r');     % plot the fit
        hold(ax,'off');
        
    else
        warning('SOMETHING IS WRONG!');
    end
    
    % Get fit params
    result.corr.fit{iCorr}=fitparam_tmp;
    
    % Save figs
    fname_str=['corr1d_',num2str(iCorr)];
    saveas(hfig,[dir_output,fname_str,'.fig']);
    saveas(hfig,[dir_output,fname_str,'.png']);
end
clear g2_1d_tmp param0 fitparam_tmp fit_g2_tmp ax this_corr_type;

%% Save processed data
% Append to existing data file
if VERBOSE>0,fprintf('Saving data...\n'); end;
for i = 1:length(vars_save)
    if ~exist(vars_save{i},'var')
        warning(['Variable "',vars_save{i},'" does not exist.']);
        continue;
    end
    save(configs.files.saveddata,vars_save{i},'-v6','-append');
end

%% END
t_fun_end=toc(t_fun_start);   % end of code
if VERBOSE>0
    disp('-----------------------------------------------');
    fprintf('Total elapsed time for %s (s): %7.1f\n','corrTaskManager',t_fun_end);
    disp('-----------------------------------------------');
end