function [fit_param,fit_xy]=gaussfit2(x,y,param0,VERBOSE)
% Gaussian fit to the correlation plots and return fit params as estimate(amp,mu,c),SE
% gaussfit2 fixes offset to 1

% INPUT
% x: x-data
% y: y-data
% params0: 1x3 array of initial conditions to fitting algorithm - [amp_0, mu_0, sigma_0]
%
% OUTPUT
% fit_params: summary of fit parameters - [avg_params, sd_params]
% fit_xy: struct with fields x, y - fitted curve

% Input parsing
if ~exist('VERBOSE','var')
    VERBOSE=0;
end

fo = statset('TolFun',10^-10,...
    'TolX',10^-10,...
    'MaxIter',10^6,...
    'UseParallel',0);

fitobject=fitnlm(x,y,...
    'y~amp*exp(-1*(x1-mu)^2/(2*sigma^2))+1',...
    param0,...
    'CoefficientNames',{'amp','mu','sigma'},'Options',fo);

fit_param=[fitobject.Coefficients.Estimate,fitobject.Coefficients.SE];
fit_xy.x=linspace(min(x),max(x),300);
fit_xy.y=feval(fitobject,fit_xy.x);

if VERBOSE>0
    figure();
    plot(x,y,'*');
    hold on;
    plot(fit_xy.x,fit_xy.y);
    legend('Data','Gaussian fit');
    xlabel('x'); ylabel('y');
end
end