% SUMMARY OF BB g2 amplitude and scattering halo mode occupancy

n_exp=[0.017,0.015,0.04,0.86,0.09,2.3e-3,0.3];
g2_exp=[85,22,32,3.1,16,140,7];

r=3;
n=linspace(min(n_exp)/r,r*max(n_exp),1000);       % logs
% n=linspace(0,1.1*max(n_exp),1000);       % linear
g2_theory=1+1./n;

% Plot
% plot mode occupancy vs g2(0) - 1
figure();
plot(n_exp,g2_exp-1,'*');
hold on;
plot(n,g2_theory-1);

ax=gca;

% ax.XScale='lin';
% ax.YScale='lin';
ax.XScale='log';
ax.YScale='log';

axis tight;
grid on;
xlabel('$n$');
ylabel('$g^{(2)}(0)-1$');
axis tight;