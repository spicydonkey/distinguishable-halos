%% g2 normalisation factor for noisy datasets
% this illustrates exact method of calculating the G2 normalisation factor
%   usually large per-shot number and low number fluctuation are assumed to
%   give the approximation: norm_factor=number_of_shots

Nshots=1000;    % number of repeated experimental runs (shots)

% randomly generate atoms counts detected in each experimental run
n_mean=10;
n_sigma=0;
n_counts=normrnd(n_mean,n_sigma,[Nshots,1]);    % generate from normal dist
n_counts=round(n_counts);   % round to nearest integer

%% Combinatorics: unique pair counting
% correlated pairs
Z_corr=sum(0.5*(n_counts.*(n_counts-1)));   % total number of correlated pairs

% uncorrelated pairs
Z_uncr=0.5*((sum(n_counts).^2)-sum(n_counts.^2));

norm_factor=Z_uncr/Z_corr;