% SCATTER PLOT ZXY ARRAY
% DKS 31/10/2016

function scatter_zxy(FIG,ZXY_ARRAY,COLOR)
if ~exist('COLOR','var')
    COLOR='k';  % default COLOR is black
end

figure(FIG); hold on;
scatter3(ZXY_ARRAY(:,2),ZXY_ARRAY(:,3),ZXY_ARRAY(:,1),1,[COLOR,'.']);
axis vis3d; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
end