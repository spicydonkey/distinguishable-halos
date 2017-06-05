function [Nsc,dk,gfit] = halo_characterise(halos,zcap,verbose)
% Characterises halos:
%   Number of scattered atoms, Nsc
%   Halo thickness, dk
% INPUT
%   halos is a "centered" Nx1 cell of Mx3 ZXY array

% define params
% zcap=config_halo.zcap;  % zcap spherical cap culling fraction to elim bec
det_qe=0.1;     % detector QE


%% Halo thickness
% get COM radial distance to each atom
halo_R=cellfun(@(C) sqrt(sum(C.^2,2)),halos,'UniformOutput',false);
halo_R_all=vertcat(halo_R{:});  % collate all shots

nbins=ceil(length(halo_R_all)/100);      % number of bins to histogram radial distribution (avg 100 counts per bin)
[N_r,r_edge]=histcounts(halo_R_all,nbins);  % do the histogram
N_pdf=N_r./diff(r_edge);                % in probability density function
r_cent=r_edge(1:end-1)+0.5*diff(r_edge);    % get bin centres

% fit Gaussian
gaussfun='y~amp*exp(-1*(x1-mu)^2/(2*sigma^2))';
param0=[100,1,0.1];
fo = statset('TolFun',10^-10,...
    'TolX',10^-10,...
    'MaxIter',10^6,...
    'UseParallel',0);
gfit=fitnlm(r_cent,N_r,...
     gaussfun,param0,...
     'CoefficientNames',{'amp','mu','sigma'},'Options',fo);

dk=gfit.Coefficients.Estimate(3);   % halo thickness in standard deviation
 
fitx=linspace(min(r_cent),max(r_cent),300);
fity=feval(gfit,fitx);

if verbose>0
    figure();
    hold on;
    plot(r_cent,N_r,'o');
    plot(fitx,fity);
    box on;
    
    title('Halo radial distribution');
    xlabel('Radius');
    ylabel('Number in radial bin');
    
    drawnow;
end

%% number of scattered atoms
Nsc=cellfun(@(C) size(C,1),halos);      % number detected
Nsc=Nsc/(zcap*det_qe);      % compensate for culled caps and detector QE

% summary statistics
Nsc_avg=mean(Nsc);
Nsc_std=std(Nsc);
Nsc_se=Nsc_std/sqrt(size(halos,1));


%% Summary
if verbose>0
    fprintf('===================HALO SUMMARY===================\n');
    fprintf('Halo width (rms): %0.2g\n',dk);
    fprintf('Number of scattered atoms, Nsc: %0.3g ± %0.2g\n',Nsc_avg,Nsc_se);
    fprintf('std Nsc: %0.2g\n',Nsc_std);
    fprintf('==================================================\n');
end