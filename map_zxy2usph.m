% ZXY cell-array mapping: ellipsoid --> usphere

function zxy_usph = map_zxy2usph(zxy,param_ellipsoid)
% circ transform ZXY-->XYZ coord sys; tranform to usph by given param; reverse XYZ-->ZXY
zxy_usph=cellfun(@(ZXY) circshift(ell2usph(circshift(ZXY,-1,2),param_ellipsoid),1,2),...
    zxy,'UniformOutput',false);
end