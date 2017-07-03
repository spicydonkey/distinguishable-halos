% Configuration file for distinguishable_halos.m

%%% FLAGS
configs.flags.do_corr_analysis=0;
configs.flags.force_all_stages=0;    % force all the stages to run (useful for debug)
configs.flags.verbose=2;
configs.flags.savedata=0;       % TODO - req'd currently since each stage passes data by save/load to disk
configs.flags.archive_txy=1;        % archives loaded TXY as .mat file for future reuse
configs.flags.graphics=1;       % toggle to control graphics/plotting options
configs.flags.build_txy=0;

%%% MISCELLANEOUS
configs.misc.vel_z=9.8*0.416;    % atom free-fall vert v at detector hit for T-to-Z conversion;
vz=configs.misc.vel_z;


%% FILES
configs.files.path='\\AMPLPC29\He BEC Archive\EXPERIMENT-DATA\xstate_mom_corr\90deg_raman_beams\9_vvlownum\d';

% WARNING: MODIFYING BELOW DIR SETTINGS ARE NOT RECOMMENDED
configs.files.dir_data=fileparts(configs.files.path);    % fullpath to data directory
configs.files.archive=fullfile(configs.files.dir_data,'archive');   % dir to archive folder
configs.files.dirout=fullfile(configs.files.dir_data,'output');      % output directory (will be time-stamped)


%% LOAD
configs.load.version=1;         % TXY load stage version number

% file ID and simple pass/fail
configs.load.id=1:5074;         % file id numbers to use for analysis
configs.load.mincount=700;         % min counts in window - 0 for no min
configs.load.maxcount=1000;          % max counts in window - Inf for no max

% Detector/trap alignment
configs.load.rot_angle=0.61;

% TXY window - region of interest ( [] --> no crop )
configs.load.window{1}=[2.44,2.48];      % T [s]
configs.load.window{2}=[-35e-3,35e-3];    % X [m]
configs.load.window{3}=[-35e-3,35e-3];    % Y [m]

%% HALO
%%% HALO PARAMS: BEC counts + oscillation removal for broad capture of halos
configs.bec.pos{1}=[2.464*vz,-3e-3,4.2e-3];   % approx condensate locations (z,x,y)
configs.bec.Rmax{1}=5e-3;      % max condensate sph radius
configs.bec.dR_tail{1}=1;     % BEC tail radial frac diff
configs.bec.pos{2}=[2.456*vz,-2.1e-3,0.6e-3];
configs.bec.Rmax{2}=5e-3;
configs.bec.dR_tail{2}=1;

configs.halo.R{1}=26e-3;     % estimated radius of halo
configs.halo.dR{1}=0.2;      % broad radial mask fractional width (in/out)
configs.halo.R{2}=24e-3;
configs.halo.dR{2}=0.2;

configs.halo.zcap=0.8;   % z-cutoff (fractional wrt radius)

% TODO - does boost need to be optimised for different g2 analysis?
%   currently SINGLE boost applied to halo2 to obtain best g2_01_BB
configs.halo.boost{1}=zeros(1,3);
% configs.halo.boost{2}=[0.01,-0.025,0.005];
configs.halo.boost{2}=[0.0,0.0,0.0];

%% CORRELATION ANALYSIS
% 1) X-halo Cart BB
configs.corr{1}.type.comp=[1,2];           % components to analysis: cross halo 1,2
configs.corr{1}.type.coord='cart';         % Cartesian (ZXY)
configs.corr{1}.type.opt='BB';             % BB / CL
configs.corr{1}.lim(1,:)=0.2*[-1,1]; % bin limits - Z
configs.corr{1}.lim(2,:)=0.2*[-1,1]; % bin limits - X
configs.corr{1}.lim(3,:)=0.2*[-1,1]; % bin limits - Y
configs.corr{1}.nBin=15*[1,1,1];   % number of bins
