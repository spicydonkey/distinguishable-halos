function [corr_out,hfig]=halo_g2_manager(zxy,configs,verbose)
% Correlation task manager

% input check
if ~exist('verbose','var')
    warning('verbose is not provided - setting to quiet (0)');
    verbose=0;
end

if verbose>0, fprintf('Starting correlation analysis...\n'), end;
%% MAIN
t_fun_start=tic;
hfig={};

% misc graphics
markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
colors = distinguishable_colors(3);
coords = {'Z','X','Y'};

n_shots=size(zxy,1);

n_corr_analysis=length(configs.corr);
corr_out=cell(n_corr_analysis,1);

%% Do tasked correlation analysis
for iCorr=1:n_corr_analysis
    corr_this=configs.corr{iCorr};
    
    % Set up bins
    bin_dim=size(corr_this.lim,1);      % get binning dims
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
    [G2_corr,G2_uncorr]=G2_caller(zxy(:,corr_this.type.comp),...
        bin_edge_tmp,corr_this.type.coord,corr_this.type.opt,verbose);
    
%     % normalisation factor
%     n_counts=cell2mat(cellfun(@(C) size(C,1),zxy,'UniformOutput',false));   % number of counts per shot
%     nn1=n_counts(:,1);      % counts in species 1
%     nn2=n_counts(:,end);    % counts in species 2
%     % normalisation factor for g2 correlation depends only on number of
%     % components
%     if length(corr_this.type.comp)==2
%         % X-state
%         N_pair_uncorr=sum(nn1)*sum(nn2)-sum(nn1.*nn2);
%         N_pair_corr=sum(nn1.*nn2);
%     else
%         % Single-halo
%         N_pair_corr=sum(nn1.*(nn1-1)/2);
%         N_pair_uncorr=(sum(nn1)^2-sum(nn1.^2))/2;
%     end
%     norm_factor=N_pair_uncorr/N_pair_corr;
    
    % TODO - try smoothing G2 before calculating g2
%     g2=norm_factor*G2_corr./G2_uncorr;      % normalised g2 RAW
    g2=G2_corr./G2_uncorr;
    
    % Get results
    corr_out_this.bEdge=bin_edge_tmp;
    corr_out_this.bCent=bin_cent_tmp;
    corr_out_this.G2_corr=G2_corr;
    corr_out_this.G2_uncorr=G2_uncorr;
    corr_out_this.g2=g2;
    
    corr_out{iCorr}=corr_out_this;
    
    %% Plot g2 result
    if configs.flags.graphics
        hfig_g2_this=figure();
        hfig{length(hfig)+1}=gcf;
        
        plotCorr(corr_out_this,configs.corr{iCorr});
        
        drawnow;
    end
end
% clear workspace
clear bin_dim bin_edge_tmp bin_cent_tmp G2_corr G2_uncorr g2;

%% 1D correlation profile and Gaussian fit
for iCorr=1:n_corr_analysis
    corr_this=configs.corr{iCorr};
    corr_out_this=corr_out{iCorr};
    
    % prepare graphics
    if configs.flags.graphics   
        hfig_g2_1d=figure();
        hfig{length(hfig)+1}=gcf;
        ax=gca;
    end
    
    this_corr_type=corr_this.type.coord;
    % Get integrated or sliced 1D correlation profile
    if isequal(this_corr_type,'angular')
        % TODO - correlations - dk integrated (may reduce peak height)
        g2_1d=sum(corr_out_this.G2_corr,1)./sum(corr_out_this.G2_uncorr,1);     % average dk

        % Gaussian fit (angular G2 allows fitting both BB,CL)
        % BB fit (mu=pi; c=1)
        param0=[4,0.03];            % fit estimate [amp,sigma]
        parameq={[],pi,[],1};       % fix params
%         [fitparam{1},fit_g2{1}]=gaussfit2(corr_out_this.bCent{2},g2_1d,param0,0);
        [fitparam{1},fit_g2{1}]=fit_gauss_1d(corr_out_this.bCent{2},g2_1d,param0,parameq);

        % CL fit (mu=0; c=1)
        param0=[2,0.03];            % fit estimate [amp,sigma]
        parameq={[],0,[],1};        % fix params
%         [fitparam{2},fit_g2{2}]=gaussfit2(corr_out_this.bCent{2},g2_1d,param0,0);
        [fitparam{2},fit_g2{2}]=fit_gauss_1d(corr_out_this.bCent{2},g2_1d,param0,parameq);
        
        % Plot
        if configs.flags.graphics   % plotting
            % calculated 1D correlation profile
            plot(corr_out_this.bCent{2},g2_1d,'*');
            
            hold(ax,'on');
            title_str=['$\Delta k$-integrated, ',...
                '(',num2str(corr_this.type.comp),') halos'];
            title(title_str);
            xlabel('$\Delta\theta$'); ylabel('$\bar{g}^{(2)}$');
            xlim([0,pi]); ylim auto;
            
            % plot the fits
            plot(ax,fit_g2{1}.x,fit_g2{1}.y,'r','LineWidth',2);     
            plot(ax,fit_g2{2}.x,fit_g2{2}.y,'b--','LineWidth',2);
            hold(ax,'off');
            
            legend({'Data','Gaussian fit (BB)','Gaussian fit (CL)'});            
        end
        
    elseif isequal(this_corr_type,'cart')
        %%% Get line through Z-axis @ dkx=dky=0
        [~,I_zero]=cellfun(@(x) min(abs(x)),corr_out_this.bCent);   % index to bin centre nearest 0
        
        % Gaussian fit
        % Set initial params; fix background/constant offset to 1
        if strcmp(corr_this.type.opt,'BB')
            param0=[4,0,0.03];      % BB [amp,mu,sigma]
            parameq={[],[],[],1};   % fix params
        else
            % CL --> symmetric around 0
            param0=[2,0.03];        % CL [amp,sigma]
            parameq={[],0,[],1};
        end
            
        perm_vect=[2,3,1];      % cyclic map (1,2,3-->2,3,1) to take 1D profile in Z,X,Y
        g2=corr_out_this.g2;    % 3D g2 in ZXY coord
        for jj=1:3
            % Gaussian fit to line profile thru zero!
            [fitparam{jj},fit_g2{jj}]=fit_gauss_1d(corr_out_this.bCent{jj},...
                g2(:,I_zero(2),I_zero(3)),param0,parameq);
            
            % Plot
            if configs.flags.graphics
                hold on;
                
                % plot data
                plot(corr_out_this.bCent{jj},g2(:,I_zero(2),I_zero(3)),...
                    markers{jj},'MarkerEdgeColor',colors(jj,:),...
                    'DisplayName',coords{jj});
                
                % plot gaussian fit
                plot(ax,fit_g2{jj}.x,fit_g2{jj}.y,...
                    'Color',colors(jj,:),'LineWidth',1.5,...
                    'HandleVisibility','off');
            end
            
            % cyclically permute the 3D vectors: CAREFUL! don't use after
            % this loop
            g2=permute(g2,perm_vect);
            I_zero=I_zero(perm_vect);
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
    corr_out{iCorr}.fit=fitparam;
    
    clear fitparam fit_g2;    % clean temp looping data
end

%% END
t_fun_end=toc(t_fun_start);   % end of code
if verbose>0
    disp('-----------------------------------------------');
    fprintf('Total elapsed time for %s (s): %7.1f\n',mfilename,t_fun_end);
    disp('-----------------------------------------------');
end