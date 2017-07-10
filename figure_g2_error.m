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
    corrthis=corr{ii};
    corrcfg=configs.corr{ii};
    
    g2=corrthis.g2;
    g2_se=g2_cluster_se{ii};
    
    if isequal(corrcfg.type.coord,'cart')
        % plot g2 line profiles thru Z,X,Y axis
        [~,I_zero]=cellfun(@(x) min(abs(x)),corrthis.bCent);   % index to bin centre nearest 0
        
        perm_vect=[2,3,1];      % cyclic map (1,2,3-->2,3,1) to take 1D profile in Z,X,Y
        
        hfig_g2{ii}=figure();
        hold on;
        for jj=1:3
            % get data to plot
            X=corrthis.bCent{jj};               % bin centers along each axis (transverse = 0)
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
        % annotate
        legend('show');
        box on;
        axis auto;
        ylim_auto=get(gca,'YLim');
        axis auto;
        set(gca,'YLim',[0,ylim_auto(2)]);
        set(gca,'Clipping','on');
        title_str=['(',num2str(corrcfg.type.comp),') halos, ',...
            corrcfg.type.opt,', ','single-axis'];
        title(title_str);
        xlabel('$\Delta K_i$'); ylabel(['$g^{(2)}_{',corrcfg.type.opt,'}$']);
        
    else
        % TODO for angular g2
        continue;
    end
end