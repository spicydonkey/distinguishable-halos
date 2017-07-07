% mf=(1,1) collision analysis
% Exmaple of single-halo capture

%% set config
config_mf1;

verbose=configs.flags.verbose;


%% load TXY
txy=load_txy(configs.files.path,configs.load.id,...
            configs.load.window,configs.load.mincount,configs.load.maxcount,...
            configs.load.rot_angle,configs.flags.build_txy,verbose,configs.flags.graphics);
        
zxy=txy2zxy(txy);       % zxy


%% Halo capture
% capture BECs at poles
[bec_cent,bool_bec]=capture_bec(zxy,configs.bec.pos,configs.bec.Rmax,verbose);

% pop bec + thermal counts
% bool_bec_combined=cellfun(@(x,y)x|y,bool_bec(:,1),bool_bec(:,2),'UniformOutput',false);
bool_bec_combined=cell_horzcat(bool_bec);       % horizontal concat the boolean list
bool_bec_combined=cellfun(@(x)sum(x,2)>0,bool_bec_combined,'UniformOutput',false);    % row-wise OR

halo_zxy=cellfun(@(ZXY,BOOL)ZXY(~BOOL,:),zxy,bool_bec_combined,'UniformOutput',false);

% halo centre: midpoint of 'BEC poles'
halo_cent=cellfun(@(c1,c2)(c1+c2)/2,bec_cent(:,1),bec_cent(:,2),'UniformOutput',false);

% centre to halo
halo_zxy0=cellfun(@(zxy,c0)zxy-repmat(c0,[size(zxy,1),1]),halo_zxy,halo_cent,'UniformOutput',false);


%% radial capture (spherical shell)
% distance from halo centre
r0=cellfun(@(x)sqrt(sum(x.^2,2)),halo_zxy0,'UniformOutput',false);
rlim=configs.halo.R{1}*(1+configs.halo.dR{1}*[-1,1]);       % radial limits for sph-shell capture

bool_halo=cellfun(@(r)(r<rlim(2))&(r>rlim(1)),r0,'UniformOutput',false);    % determine halo points

% capture halo!
halo_zxy0=cellfun(@(ZXY,BOOL)ZXY(BOOL,:),halo_zxy0,bool_halo,'UniformOutput',false);

% remove caps
dz_poles=cellfun(@(c1,c2)abs(c1(1)-c2(1)),bec_cent(:,1),bec_cent(:,2),'UniformOutput',false);
dz_poles=mean(vertcat(dz_poles{:}));    % pole-pole distance (BEC) in Z

bool_zcap=cellfun(@(zxy)abs(zxy(:,1))>(0.5*dz_poles*configs.halo.zcap),halo_zxy0,'UniformOutput',false);

halo_zxy0=cellfun(@(ZXY,BOOL)ZXY(~BOOL,:),halo_zxy0,bool_zcap,'UniformOutput',false);


%% Cull aliased hits
% cull from centred zxy0 halo data
alias_deadz=100e-9*configs.misc.vel_z;      % dZ in 100 ns
bool_alias=cellfun(@(x) findalias(x,alias_deadz),halo_zxy0,'UniformOutput',false);
halo_zxy0_filt=cellfun(@(x,y) x(~y,:),halo_zxy0,bool_alias,'UniformOutput',false);

% TODO - better var management
halo_zxy0_orig=halo_zxy0;       % store original
halo_zxy0=halo_zxy0_filt;       % alias filtered


%% Elipsoid fit
efit_flag='';
[ecent,erad,evecs,v,echi2]=ellipsoid_fit(circshift(vertcat(halo_zxy0{:}),-1,2),efit_flag);


%% Tranform to unit sphere (k-space)
% Initialise k-vector in XYZ coord
halo_k=cellfun(@(x) circshift(x,-1,2),halo_zxy0,'UniformOutput',false);

% Centre to fitted ellipsoid
halo_k(:)=boost_zxy(halo_k,-ecent');

% Transform to ellipsoid principal axis (Rotation)
M_rot=evecs;   % rotation matrix: X'Y'Z'(principal ellipsoid) --> XYZ
halo_k=cellfun(@(x) (M_rot\x')',halo_k,'UniformOutput',false);     % inverse transform

% Stretch to UNIT sphere: unit in collision wave-vector/momenta
halo_k=cellfun(@(x) x./repmat(erad',size(x,1),1),halo_k,'UniformOutput',false);

% Reverse transform to original/detector axis
halo_k=cellfun(@(x) (M_rot*x')',halo_k,'UniformOutput',false);

% transform to ZXY system
halo_k=cellfun(@(x) circshift(x,1,2),halo_k,'UniformOutput',false);


%% Characterise halo
[Nsc,dk]=halo_characterise(halo_k,configs.halo.zcap,verbose);


%% Halo centering
% boost for best g2 BB centering
halo_k=boost_zxy(halo_k,configs.halo.boost{1});


%% Correlation analysis
corr=halo_g2_manager(halo_k,configs,verbose);