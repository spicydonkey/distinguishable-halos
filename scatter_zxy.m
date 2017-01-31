% SCATTER PLOT ZXY ARRAY
% Scatter plot for a N-by-3 "ZXY"-array
% DKS 31/10/2016

function scatter_zxy(ZXY_ARRAY,SIZE,COLOR)

if isempty(ZXY_ARRAY)
    warning('Empty array passed as data.');
    return;
end
if ~exist('COLOR','var')
    COLOR='k';  % default COLOR is black
end
if ~exist('SIZE','var')
    SIZE=1;     % default scatter dot size
end

scatter3(ZXY_ARRAY(:,2),ZXY_ARRAY(:,3),ZXY_ARRAY(:,1),SIZE,[COLOR,'.']);
axis vis3d; 
%axis equal;

view(3);    %default view
end