% objective function for mapping halo to a unit sphere from an ellipsoid
% fit: see ell2usph
function width_rms=halo_sph_rms(halo_zxy,ellip_param)
xyz=circshift(halo_zxy,[0,-1]);          % to XYZ coord

xyz_mapped=ell2usph(xyz,ellip_param);

r=sqrt(sum(xyz_mapped.^2,2));           % distances from centre

width_rms=std(r);               % rms width - unit sphere mapped

end