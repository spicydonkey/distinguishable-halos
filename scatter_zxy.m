% SCATTER PLOT ZXY ARRAY
% DKS 31/10/2016

function scatter_zxy(FIG,ZXY_ARRAY)
figure(FIG); hold on;
scatter3(ZXY_ARRAY(:,2),ZXY_ARRAY(:,3),ZXY_ARRAY(:,1),1,'k.');
axis vis3d; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
end