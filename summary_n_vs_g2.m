% SUMMARY OF BB g2 amplitude and scattering halo mode occupancy

%% DATA
n_exp=[0.0080, 0.0156, 0.0389, 0.86, 0.086, 0.011, 0.3];
g2_exp=[87, 22, 32, 3.1, 20, 95, 7];

% UNCERTAINTY
n_unc=[0.68, 0.8, 0.38, 0.8, 0.39, 0.65, 0.8];      % simple ±80% uncertainty
n_err=n_exp.*n_unc;      % calculated

% g2_err
g2_unc=0.5;     % simple
g2_err=g2_exp*g2_unc;

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

hdata=ploterr(n_exp,g2_exp-1,n_err,g2_err,'ko','logxy','hhxy',0);
set(hdata(1),namearray,valarray,'DisplayName','Data');              % DATAPOINT
set(hdata(2),namearray,valarray,'DisplayName','');                  % Y-err
set(hdata(3),namearray,valarray,'DisplayName','');                  % X-err

% annotate
legend([htheory,hdata(1)]);

ax=gca;
ax.XScale='log';
ax.YScale='log';

axis tight;
box on;
grid on;
xlabel('$n$');
ylabel('$g^{(2)}(0)-1$');
axis tight;