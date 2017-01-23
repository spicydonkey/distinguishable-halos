% PLOT ZXY DATA
% Plots Nshot-by-Mspecies cell of ZXY counts
% DKS 

function plot_zxy(ZXY,SIZE,COLORS)
if ~iscell(ZXY)
    warning('ZXY must be a cell array.');
    return;
end
if ~exist('COLORS','var')
    COLORS='brgmcyk';  % default COLORS
end
if ~exist('SIZE','var')
    SIZE=1;     % default scatter dot size
end

hold on;    % hold current figure - to plot on
n_species=size(ZXY,2);
for i=1:n_species
    temp_zxy=vertcat(ZXY{:,i});
    scatter_zxy(temp_zxy,SIZE,COLORS(i));
end

end