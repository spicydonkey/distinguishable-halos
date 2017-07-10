% plot g2 with error bars
%   requires: corr, g2_cluster_se

% plot configs
coords = {'Z','X','Y'};
colors={'b','r','g'};
markers = {'s','o','^','.','x','+','d','*','v','>','<','p','h'};
lstyle={'-','--',':','-.'};

namearray={'LineWidth','MarkerFaceColor'};
valarray={1.5,'w'};

hfig_g2=cell(ncorr,1);
for ii=1:ncorr
    corr_this=corr{ii};
    corrcfg=configs.corr{ii};
    
    g2=corr_this.g2;
    g2_se=g2_cluster_se{ii};
    
    if isequal(corrcfg.type.coord,'cart')
        % plot g2 line profiles thru Z,X,Y axis
        [~,I_zero]=cellfun(@(x) min(abs(x)),corr_this.bCent);   % index to bin centre nearest 0
        
        perm_vect=[2,3,1];      % cyclic map (1,2,3-->2,3,1) to take 1D profile in Z,X,Y
        
        hfig_g2{ii}=figure();
        hold on;
        for jj=1:3
            % get data to plot
            X=corr_this.bCent{jj};               % bin centers along each axis (transverse = 0)
            Y=g2(:,I_zero(2),I_zero(3));        % g2 along this axis
            Yerr=g2_se(:,I_zero(2),I_zero(3));  % SE of g2
            
            % plot
            hdata=errorbar(X,Y,Yerr,[colors{jj},markers{jj}],...
                'DisplayName',coords{jj},namearray,valarray);  
%             alpha(0.5);
            
            % cyclically permute the 3D vectors: CAREFUL! don't use after
            % this loop
            g2=permute(g2,perm_vect);
            g2_se=permute(g2_se,perm_vect);
            I_zero=I_zero(perm_vect);
        end
        %% annotate
        legend('show');
        box on;
        
        % set ylim from 0 -- auto max
        axis auto;
        ylim_auto=get(gca,'YLim');
        set(gca,'YLim',[0,ylim_auto(2)]);   % set min to 0
        
        % clip data to axis lims
        set(gca,'Clipping','on');
        
        % labelling
        title_str=['(',num2str(corrcfg.type.comp),') halos, ',...
            corrcfg.type.opt,', ','single-axis'];
        title(title_str);
        xlabel('$\Delta K_i$'); ylabel(['$g^{(2)}_{',corrcfg.type.opt,'}$']);
        
    else
        % TODO for angular g2
        X=corr_this.bCent{2};        % theta bins
        
        % radially integrated g2 (may have reduced peak height)
        Y=sum(corr_this.G2_corr,1)./sum(corr_this.G2_uncorr,1);
        
        %% evaluate error
        % METHOD 1 - radial SD
        Yerr=std(g2,[],1);      % take SD at each theta along dk
        
        %% plot 
        hfig_g2{ii}=figure();
        hdata=errorbar(X,Y,Yerr,'o',...
            'DisplayName','Data',namearray,valarray);
        
        %% annotate
        legend('show');
        box on;
        
        % set ylim from 0 -- auto max / xlim 0 -- pi
        axis auto;
        ylim_auto=get(gca,'YLim');
        set(gca,'XLim',[0,pi]);
        set(gca,'YLim',[0,ylim_auto(2)]);   % set min to 0
        
        % clip data to axis lims
        set(gca,'Clipping','on');
        
        % labelling
        title_str=['$\Delta k$-integrated, ',...
            '(',num2str(corrcfg.type.comp),') halos'];
        title(title_str);
        xlabel('$\Delta\theta$'); ylabel('$\bar{g}^{(2)}$');
    end
end