% SUMMARY OF BB g2 amplitude and scattering halo mode occupancy

%% CONFIGS
% plotting
n_lim=[2e-3,2e0];       % axes lims for mode occupancy
g2_lim=[3e-1,3e2];      % axes lims for g2-1

%%% DATA
%% 90 deg data
n_exp=[0.0080, 0.0156, 0.0389, 0.604, 0.086, 0.011, 0.273];
g2_exp=[87, 22, 32, 3.6, 20, 95, 8];

% UNCERTAINTY
n_unc=[0.68, 1.1, 0.38, 0.24, 0.39, 0.65, 0.32];      % relative error in ratio to mean

% g2_err
g2_unc=0.75;     % simple

% %% 30 deg data
% n_exp=[n_exp 0.224];
% n_unc=[n_unc 0.46];
% 
% g2_exp=[g2_exp 7.5];
%     

%% EVALUATE error limits
% relative error
n_err=n_exp.*n_unc;      % calculated
g2_err=g2_exp*g2_unc;

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

% cut bound at axes lims (MIN only)
for ii=1:ndata
    if n_err_bnd{1}(ii)<n_lim(1)
        n_err_bnd{1}(ii)=n_lim(1);
    end
    if g2_err_bnd{1}(ii)<g2_lim(1)
        g2_err_bnd{1}(ii)=g2_lim(1);
    end
end

%% Theory
r=3;
n=linspace(min(n_exp)/r,r*max(n_exp),1000);       % logs
% n=linspace(0,1.1*max(n_exp),1000);       % linear
g2_theory=1+1./n;


%% Plot
% plot configs
namearray={'LineWidth','MarkerFaceColor'};
valarray={1.5,'w'};

% plot mode occupancy vs g2(0) - 1
figure();
hold on;
htheory=plot(n,g2_theory-1,'r--','LineWidth',1.5,'DisplayName','Theory');       % theory curve

% hdata=ploterr(n_exp,g2_exp-1,n_err,g2_err,'ko','logxy','hhxy',0);
hdata=ploterr(n_exp,g2_exp-1,n_err_bnd,g2_err_bnd,'ko','logxy','hhxy',0);
set(hdata(1),namearray,valarray,'DisplayName','Data');              % DATAPOINT
set(hdata(2),namearray,valarray,'DisplayName','');                  % Y-err
set(hdata(3),namearray,valarray,'DisplayName','');                  % X-err

% annotate plot
legend([htheory,hdata(1)]);

ax=gca;
ax.XScale='log';
ax.YScale='log';

xlim(n_lim);
ylim(g2_lim);

axis tight;
box on;
grid on;
xlabel('$n$');
ylabel('$g^{(2)}(0)-1$');
axis tight;