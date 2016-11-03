function SPH_OUT = zxy2pol(ZXY_IN)
% Converts counts in Cartesian coordinate to Spherical coord
%
% I/O
% ZXY_IN: should be a Nx1 cell of ZXY-array
%
% SPH_OUT: Nx1 cell of arrays of sph polar coord [NORM,AZIM,ELEV]

nShot=size(ZXY_IN,1);

SPH_OUT=cell(nShot,1);
for i=1:nShot
    zxy_tmp=ZXY_IN{i};
    rad_tmp=sqrt(sum(zxy_tmp.^2,2));    % radius
    azim_tmp=atan(zxy_tmp(:,3)./zxy_tmp(:,2));  % azimuthal ang=atan(Y/X)
    elev_tmp=asin(zxy_tmp(:,1)./rad_tmp);       % elevation ang=asin(Z/NORM)
    
    SPH_OUT{i}=[rad_tmp, azim_tmp, elev_tmp];   % SPHPOL=(R,azim,elev)
end