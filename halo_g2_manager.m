function corr_out=halo_g2_manager(zxy,configs,verbose)
% Correlation task manager

% input check
if ~exist('verbose','var')
    warning('verbose is not provided - setting to quiet (0)');
    verbose=0;
end

if verbose>0, fprintf('Starting correlation analysis...\n'), end;
%% MAIN
t_fun_start=tic;

% misc graphics
markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
colors = distinguishable_colors(3);
coords = {'Z','X','Y'};

% get paths
dir_output=configs.files.dirout;

n_corr_analysis=length(configs.corr);
corr_out=cell(n_corr_analysis,1);
%% Do tasked correlation analysis
for iCorr=1:n_corr_analysis
    corr_this=configs.corr{iCorr};
    
    % Set up bins
    bin_dim=length(corr_this.lim);  % get binning dims
    bin_edge_tmp=cell(1,bin_dim);
    bin_cent_tmp=cell(1,bin_dim);
    for i=1:bin_dim
        % make bin edge and centre vectors
        bin_edge_tmp{i}=linspace(corr_this.lim(i,1),...
            corr_this.lim(i,2),corr_this.nBin(i)+1);
        bin_cent_tmp{i}=0.5*(bin_edge_tmp{i}(1:end-1)+bin_edge_tmp{i}(2:end));
    end
    
    % Evaluate G2 correlation
    % TODO - insert debug code for G2 and zxy manipulator
    [G2_shot_tmp,G2_all_tmp]=G2_caller(zxy(:,corr_this.type.comp),...
        bin_edge_tmp,corr_this.type.coord,corr_this.type.opt,verbose);
    g2_tmp=size(zxy,1)*G2_shot_tmp./G2_all_tmp;      % normalised g2
    
    % Get results
    corr_out_this.bEdge=bin_edge_tmp;
    corr_out_this.bCent=bin_cent_tmp;
    corr_out_this.G2shot=G2_shot_tmp;
    corr_out_this.G2all=G2_all_tmp;
    corr_out_this.g2=g2_tmp;
    
    corr_out{iCorr}=corr_out_this;
    
    %% Plot g2 result
    if configs.flags.graphics
        hfig_g2_this=figure();
        plotCorr(corr_out_this,configs.corr{iCorr});
        
        drawnow;
        
        % save figs
        fname_str=['corr_',num2str(iCorr)];
        saveas(hfig_g2_this,fullfile(dir_output,[fname_str,'.fig']));
        saveas(hfig_g2_this,fullfile(dir_output,[fname_str,'.png']));
    end
end
% clear workspace
clear bin_dim bin_edge_tmp bin_cent_tmp G2_shot_tmp G2_all_tmp g2_tmp;

%% 1D correlation profile and Gaussian fit
for iCorr=1:n_corr_analysis
    corr_this=configs.corr{iCorr};
    corr_out_this=corr_out{iCorr};
    
    % prepare graphics
    if configs.flags.graphics   
        hfig_g2_1d_this=figure();
        ax=gca;
    end
    
    this_corr_type=corr_this.type.coord;
    % Get integrated or sliced 1D correlation profile
    if isequal(this_corr_type,'angular')
        % TODO: correlations - dk integrated (may reduce peak height)
        g2_1d_tmp=size(zxy,1)*sum(corr_out_this.G2shot,1)./sum(corr_out_this.G2all,1);
        
        % Gaussian fit (angular G2 allows fitting both BB,CL)
        % BB fit
        param0=[4,pi,0.1,1];     % fit estimate [amp,mu,sigma,offset]
        [fitparam_tmp,fit_g2_tmp]=gaussfit(corr_out_this.bCent{2},g2_1d_tmp,param0,0);
        
        % CL fit
        param0=[2,0,0.1,1];
        [fitparam_tmpCL,fit_g2_tmpCL]=gaussfit(corr_out_this.bCent{2},g2_1d_tmp,param0,0);
        
        % Plot
        if configs.flags.graphics   % plotting
            % calculated 1D correlation profile
            plot(corr_out_this.bCent{2},g2_1d_tmp,'*');
            
            hold(ax,'on');
            title_str=['$\Delta k$-integrated, ',...
                '(',num2str(corr_this.type.comp),') halos'];
            title(title_str);
            xlabel('$\Delta\theta$'); ylabel('$\bar{g}^{(2)}$');
            xlim([0,pi]); ylim auto;
            
            % plot the fits
            plot(ax,fit_g2_tmp.x,fit_g2_tmp.y,'r');     
            plot(ax,fit_g2_tmpCL.x,fit_g2_tmpCL.y,'b--');
            hold(ax,'off');
            
            legend({'Data','Gaussian fit'});            
        end
        
    elseif isequal(this_corr_type,'cart')
        % TODO - do in X/Y?
        %%% Get line through Z-axis @ dkx=dky=0
        ind_zero_tmp=round((corr_this.nBin+1)/2);    % zero-cent'd bin index for sampling 3D-g2
        
        % Gaussian fit
        % Set initial params
        if strcmp(corr_this.type.opt,'BB')
            param0=[4,0,0.1,1];     % fit estimate [amp,mu,sigma,offset]
        else
            param0=[2,0,0.1,1];     % CL peaks at 2
        end
            
        perm_vect=[2,3,1];      % 1,2,3-->2,3,1 cyclic map
        g2=corr_out_this.g2;    % 3D g2 in ZXY coord
        for jj=1:3
            % take line profile thru the bin center (d~=0)
            [fitparam_tmp{jj},fit_g2_tmp{jj}]=gaussfit(corr_out_this.bCent{jj},...
                g2(:,ind_zero_tmp(2),ind_zero_tmp(3)),param0,0);
            
            % Plot
            if configs.flags.graphics
                hold on;
                
                % plot data
                plot(corr_out_this.bCent{jj},g2(:,ind_zero_tmp(2),ind_zero_tmp(3)),...
                    markers{jj},'MarkerEdgeColor',colors(jj,:),...
                    'DisplayName',coords{jj});
                
                % plot gaussian fit
                plot(ax,fit_g2_tmp{jj}.x,fit_g2_tmp{jj}.y,...
                    'Color',colors(jj,:),'LineWidth',1.5,...
                    'HandleVisibility','off');
            end
            
            % cyclically permute the 3D vectors: CAREFUL! don't use after
            % this loop
            g2=permute(g2,perm_vect);
            ind_zero_tmp=ind_zero_tmp(perm_vect);
        end
        % figure annotation
        if configs.flags.graphics
            title_str=['(',num2str(corr_this.type.comp),') halos, ',...
                corr_this.type.opt,', ','single-axis'];
            title(title_str);
            xlabel('$\Delta K_i$'); ylabel(['$g^{(2)}_{',corr_this.type.opt,'}$']);
            legend('show');
            box on;
        end
        
    else
        warning('SOMETHING IS WRONG!');
    end
    drawnow;
    
    % Get fit params
    corr_out{iCorr}.fit=fitparam_tmp;
    
    % Save figs
    if configs.flags.graphics   % plotting
        fname_str=['corr1d_',num2str(iCorr)];
        saveas(hfig_g2_1d_this,fullfile(dir_output,[fname_str,'.fig']));
        saveas(hfig_g2_1d_this,fullfile(dir_output,[fname_str,'.png']));
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