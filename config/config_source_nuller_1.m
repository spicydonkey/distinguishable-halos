% Configuration file for distinguishable_halos.m

%%% FLAGS
configs.flags.do_corr_analysis=0;
    configs.flags.do_corr_err=0;
configs.flags.force_all_stages=0;    % force all the stages to run (useful for debug)
configs.flags.verbose=2;
configs.flags.savedata=0;       % TODO - req'd currently since each stage passes data by save/load to disk
configs.flags.archive_txy=0;        % archives loaded TXY as .mat file for future reuse
configs.flags.graphics=1;       % toggle to control graphics/plotting options
configs.flags.build_txy=1;

%%% MISCELLANEOUS
configs.misc.vel_z=9.8*0.416;    % atom free-fall vert v at detector hit for T-to-Z conversion;
vz=configs.misc.vel_z;
configs.misc.deadtime=300e-9;    %100e-9;    % data acquired after disc adjustment

%% FILES
configs.files.path='C:\Users\HE BEC\bell\source_nuller\1\d';

% WARNING: MODIFYING BELOW DIR SETTINGS ARE NOT RECOMMENDED
configs.files.dir_data=fileparts(configs.files.path);    % fullpath to data directory
configs.files.archive=fullfile(configs.files.dir_data,'archive');   % dir to archive folder
configs.files.dirout=fullfile(configs.files.dir_data,'output');      % output directory (will be time-stamped)


%% LOAD
configs.load.version=1;         % TXY load stage version number

% file ID and simple pass/fail
configs.load.id=[];
configs.load.mincount=0;
configs.load.maxcount=Inf;

% Detector/trap alignment
configs.load.rot_angle=0.61;

% TXY window - region of interest ( [] --> no crop )
configs.load.window{1}=[0.38,0.44];      % T [s]
configs.load.window{2}=[-35e-3,35e-3];    % X [m]
configs.load.window{3}=[-35e-3,35e-3];    % Y [m]

%% HALO
%%% HALO PARAMS: BEC counts + oscillation removal for broad capture of halos
configs.bec.pos{1}=[1.6525,-3.3e-3,4.5e-3];   % approx condensate locations (z,x,y)
configs.bec.Rmax{1}=8e-3;      % max condensate sph radius
configs.bec.dR_tail{1}=1;     % BEC tail radial frac diff

configs.bec.pos{2}=[1.6289,-2.9e-3,3.9e-3];
configs.bec.Rmax{2}=8e-3;
configs.bec.dR_tail{2}=1;

configs.halo.R{1}=26e-3;     % estimated radius of halo
configs.halo.dR{1}=0.4;      % broad radial mask fractional width (in/out)

configs.halo.R{2}=24e-3;
configs.halo.dR{2}=0.4;

configs.halo.zcap=0.8;

configs.halo.boost{1}=zeros(1,3);
configs.halo.boost{2}=zeros(1,3);

