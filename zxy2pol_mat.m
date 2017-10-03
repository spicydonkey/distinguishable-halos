function sph_out = zxy2pol_mat(zxy_in)
% Converts counts in Cartesian coordinate to Spherical coord
%
% I/O
% ZXY_IN: should be a Nx3 array
%
% SPH_OUT: Nx3 array of sph polar coord [NORM,AZIM,ELEV]

rad_tmp=sqrt(sum(zxy_in.^2,2));    % radius
azim_tmp=atan2(zxy_in(:,3),zxy_in(:,2));  % azimuthal ang=atan2(Y/X) - 4-quadrant inverse tangent
elev_tmp=asin(zxy_in(:,1)./rad_tmp);       % elevation ang=asin(Z/NORM)

sph_out=[rad_tmp, azim_tmp, elev_tmp];   % SPHPOL=(R,azim,elev)