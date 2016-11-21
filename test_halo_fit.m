% Test script for ellipsoid fit to halo
% DKS 21/11/2016

set_config;
configs=usrconfigs;

DATA=load(configs.files.saveddata);     % load saved data

halo=DATA.halo.zxy0;        % get oscillation compensated halo in cell form

halo=collate_shots(halo);   % collate the halo counts (ZXY)

% Ellipsoid fit to data

for ii=1:2
    % get data points
    x=halo{ii}(:,2);
    y=halo{ii}(:,3);
    z=halo{ii}(:,1);
    
    % do the fit (with principal axis in xyz-directions)
    [center{ii}, radii{ii}, evecs{ii}, v{ii}, chi2{ii} ] = ellipsoid_fit( [ x y z ],'0');
    
    % draw fit
    figure();
    plot3(x,y,z,'.r');
    hold on;
    
    mind = min( [ x y z ] );
    maxd = max( [ x y z ] );
    nsteps = 50;
    
    step = ( maxd - mind ) / nsteps;
    [ x, y, z ] = meshgrid( linspace( mind(1) - step(1), maxd(1) + step(1), nsteps ), linspace( mind(2) - step(2), maxd(2) + step(2), nsteps ), linspace( mind(3) - step(3), maxd(3) + step(3), nsteps ) );
    
    Ellipsoid = v{ii}(1) *x.*x +   v{ii}(2) * y.*y + v{ii}(3) * z.*z + ...
        2*v{ii}(4) *x.*y + 2*v{ii}(5)*x.*z + 2*v{ii}(6) * y.*z + ...
        2*v{ii}(7) *x    + 2*v{ii}(8)*y    + 2*v{ii}(9) * z;
    
    p{ii} = patch( isosurface( x, y, z, Ellipsoid, -v{ii}(10) ));
    hold off;
    set( p{ii}, 'FaceColor', 'g', 'EdgeColor', 'none','FaceAlpha',0.7);
    view( -70, 40 );
    axis vis3d equal;
    camlight;
    lighting phong;
    
    title(sprintf('Ellipsoid fit to halo #%d',ii));
    xlabel('X'); ylabel('Y'); zlabel('Z');
end

halo_orig=DATA.halo.zxy0;
nShot=size(halo_orig,1);
halo_fitted=cell(size(halo_orig));
% Filter well-fit counts
for ii=1:2
    for jj=1:nShot
        % some kind of radial mask around the ellipsoid fit
    end
end