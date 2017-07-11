% SUMMARY OF BB g2 amplitude and scattering halo mode occupancy

%% CONFIGS
% plotting limits
n_lim=[2e-3,1e0];       % axes lims for mode occupancy
g2_lim=[1e0,2e2];      % axes lims for g2-1

% condition for violoation of Bell's test
g2_bell=(2*sqrt(2)+3)-1;            % Wasak, Chwedenczuk - limit g2 - 1


%% DATA
%%% 90 deg data
n_exp=[0.0085, 0.0126, 0.0575, 0.604, 0.0887, 0.0102, 0.273];
g2_exp=[92.5, 23.3, 31.5, 3.6, 20.3, 97.3, 8];      % mean fit amplitude

% uncertainties
n_unc=[0.68, 1.1, 0.38, 0.24, 0.39, 0.65, 0.32];      % relative error in ratio to mean

g2_err=[9.1, 1.2, 2.1, 0.15, 0.5, 8.5, 0.43];   % error in SE


%%% 30 deg data
n_exp_30=[0.224 0.577 0.82];
g2_exp_30=[7.5 3.7 2.6];
    
n_unc_30=[0.46 0.31 0.28];
g2_err_30=[1.0 0.27 0.1];


% Use SD for g2 uncertainty
g2_err=g2_err*sqrt(5);
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
r=3;
n=linspace(min(n_exp)/r,r*max(n_exp),1000);       % logs
% n=linspace(0,1.1*max(n_exp),1000);       % linear
g2_theory=1+1./n;


%% Plot
% plot configs
gray_col=0.5*ones(1,3);
namearray={'LineWidth','MarkerFaceColor','Color'};
valarray={1.5,'w','k'};
valarray_30={1.5,gray_col,gray_col};

% plot mode occupancy vs g2(0) - 1
figure();
hold on;
htheory=plot(n,g2_theory-1,'k--','LineWidth',1.5,'DisplayName','Theory');       % theory curve

% 30 deg data
hdata_30=ploterr(n_exp_30,g2_exp_30-1,n_err_30,g2_err_30,'^','logxy','hhxy',0);
set(hdata_30(1),namearray,valarray_30,'DisplayName','Slow collision');              % DATAPOINT
set(hdata_30(2),namearray,valarray_30,'DisplayName','');                  % Y-err
set(hdata_30(3),namearray,valarray_30,'DisplayName','');                  % X-err

% 90 deg data
hdata=ploterr(n_exp,g2_exp-1,n_err_bnd,g2_err_bnd,'o','logxy','hhxy',0);
set(hdata(1),namearray,valarray,'DisplayName','Fast collision');              % DATAPOINT
set(hdata(2),namearray,valarray,'DisplayName','');                  % Y-err
set(hdata(3),namearray,valarray,'DisplayName','');                  % X-err

% annotate plot
box on;
grid on;
legend([hdata(1),hdata_30(1),htheory]);

ax=gca;
ax.XScale='log';
ax.YScale='log';

xlim(n_lim);
ylim(g2_lim);

xlabel('Mode occupancy');
ylabel('$g^{(2)}_{BB}(0)-1$');

axis square;        % square plot

%% Bell's test regime
% build patch corners
x_corn=[n_lim(1),n_lim(2),n_lim(2),n_lim(1)];
% y_corn=[g2_lim(1),g2_lim(1),g2_bell,g2_bell];     % below condition
y_corn=[g2_bell,g2_bell,g2_lim(2),g2_lim(2)];       % above condition

p_bell=patch(x_corn,y_corn,...
    'b','FaceAlpha',0.1,'EdgeColor','none');

uistack(p_bell,'bottom');   % move the patch object to bottom