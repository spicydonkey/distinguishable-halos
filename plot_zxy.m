% PLOT ZXY DATA
% DKS 

function plot_zxy(FIG,ZXY,SIZE,COLORS)
if ~iscell(ZXY)
    warning('ZXY must be a cell array.');
    return;
end

if ~exist('FIG','var')
    FIG=1;      % default fig number = 1
end
if ~exist('COLORS','var')
    COLORS='brgmcyk';  % default COLORS
end
if ~exist('SIZE','var')
    SIZE=1;     % default scatter dot size
end

figure(FIG); hold on;
n_species=size(ZXY,2);
for i=1:n_species
    temp_zxy=vertcat(ZXY{:,i});
    scatter_zxy(FIG,temp_zxy,SIZE,COLORS(i));
end

end