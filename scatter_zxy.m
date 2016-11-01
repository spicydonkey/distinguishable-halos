% SCATTER PLOT ZXY ARRAY
% DKS 31/10/2016

function scatter_zxy(FIG,ZXY_ARRAY,SIZE,COLOR)

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

figure(FIG); hold on;
scatter3(ZXY_ARRAY(:,2),ZXY_ARRAY(:,3),ZXY_ARRAY(:,1),SIZE,[COLOR,'.']);
axis vis3d; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
end