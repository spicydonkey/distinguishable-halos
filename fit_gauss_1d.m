function [paramfit,fitval,fitobject]=fit_gauss_1d(x,y,param0,parameq,VERBOSE)
% Gaussian fit to the correlation plots and return fit params as estimate(amp,mu,,sigma,c) and fit param SE

% INPUT
% x: x-data
% y: y-data
% params0: 1x4(max) array of initial conditions to fitting algorithm - [amp_0, mu_0, sigma_0, c_0]
% parameq: 1x4 cell-array for Gaussian parameter constraint (empty to
% unconstrain)

% OUTPUT
% paramfit: summary of fit parameters - [avg_params, sd_params]
% evalfit: struct with fields x, y - fitted curve

% Input parsing
if ~exist('VERBOSE','var')
    VERBOSE=0;
end
if ~exist('parameq','var')
    parameq=cell(1,4);      % all params unconstrained
end

% parse problem
gausscoeffs={'amp','mu','sigma','c'};       % the full param list
funcoeffs=gausscoeffs;      % function coefficients for this particular problem
coeffnames={};              % coefficient names to optimise in this problem
for ii=1:length(parameq)
    if ~isempty(parameq{ii})
        % fix this param
        funcoeffs{ii}=num2str(parameq{ii});     % convert number with default precision
    else
        % this param is free
        coeffnames=[coeffnames,gausscoeffs{ii}];    % concat the cell array
    end
end

% build function to fit
funtofit=sprintf('y~(%s)*exp(-1*(x1-(%s))^2/(2*(%s)^2))+(%s)',funcoeffs{:});

% define optimiser
fo = statset('TolFun',10^-10,...
    'TolX',10^-10,...
    'MaxIter',10^6,...
    'UseParallel',0);

% do the fit!
fitobject=fitnlm(x,y,...
    funtofit,param0,...
    'CoefficientNames',coeffnames,'Options',fo);

% get fit results
paramfit=[fitobject.Coefficients.Estimate,fitobject.Coefficients.SE];
fitval.x=linspace(min(x),max(x),1000);
fitval.y=feval(fitobject,fitval.x);

% plot fit against data
if VERBOSE>0
    figure();
    plot(x,y,'*');
    hold on;
    plot(fitval.x,fitval.y);
    legend({'Data','Gaussian fit'});
    xlabel('x'); ylabel('y');
end
end