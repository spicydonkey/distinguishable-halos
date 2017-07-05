% Configuration file for distinguishable_halos.m

%%% FLAGS
configs.flags.do_corr_analysis=0;
configs.flags.force_all_stages=0;    % force all the stages to run (useful for debug)
configs.flags.verbose=2;
configs.flags.savedata=0;       % TODO - req'd currently since each stage passes data by save/load to disk
configs.flags.archive_txy=1;        % archives loaded TXY as .mat file for future reuse
configs.flags.graphics=1;       % toggle to control graphics/plotting options
configs.flags.build_txy=1;

%%% MISCELLANEOUS
configs.misc.vel_z=9.8*0.416;    % atom free-fall vert v at detector hit for T-to-Z conversion;
vz=configs.misc.vel_z;


%% FILES
configs.files.path='C:\Users\David\Documents\hebec\halo_analysis\data\dist_halo\mf1_collision\d';

% WARNING: MODIFYING BELOW DIR SETTINGS ARE NOT RECOMMENDED
configs.files.dir_data=fileparts(configs.files.path);    % fullpath to data directory
configs.files.archive=fullfile(configs.files.dir_data,'archive');   % dir to archive folder
configs.files.dirout=fullfile(configs.files.dir_data,'output');      % output directory (will be time-stamped)


%% LOAD
% version 1 - 10ns txy dead time
% version 2 - 100ns txy dead time
configs.load.version=2;         % TXY load stage version number

% file ID and simple pass/fail
% configs.load.id=1:1500;         % file id numbers to use for analysis
% configs.load.mincount=1500;         % min counts in window - 0 for no min
% configs.load.maxcount=1700;          % max counts in window - Inf for no max
configs.load.id=1:350;         % file id numbers to use for analysis
configs.load.mincount=0;         % min counts in window - 0 for no min
configs.load.maxcount=Inf;          % max counts in window - Inf for no max

% Detector/trap alignment
configs.load.rot_angle=0.61;

% TXY window - region of interest ( [] --> no crop )
configs.load.window{1}=[2.44,2.46];      % T [s]
configs.load.window{2}=[-35e-3,35e-3];    % X [m]
configs.load.window{3}=[-35e-3,35e-3];    % Y [m]

%% HALO
%%% HALO PARAMS: BEC counts + oscillation removal for broad capture of halos
configs.bec.pos{1}=[9.9715,-1.9e-3,-0.5e-3];   % approx condensate locations (z,x,y)
configs.bec.Rmax{1}=8e-3;      % max condensate sph radius
configs.bec.dR_tail{1}=1;     % BEC tail radial frac diff
configs.bec.pos{2}=[10.021,-2.2e-3,0.6e-3];
configs.bec.Rmax{2}=8e-3;
configs.bec.dR_tail{2}=1;

configs.halo.R{1}=24e-3;     % estimated radius of halo
configs.halo.dR{1}=0.3;      % broad radial mask fractional width (in/out)
% configs.halo.R{2}=24e-3;
% configs.halo.dR{2}=0.25;

configs.halo.zcap=0.75;   % z-cutoff (fractional wrt radius)

% TODO - does boost need to be optimised for different g2 analysis?
%   currently SINGLE boost applied to halo2 to obtain best g2_01_BB
configs.halo.boost{1}=[-0.005,0.0,-0.033];
% configs.halo.boost{2}=[0.035,-0.002,-0.025];

%% CORRELATION ANALYSIS
% 1) Single-halo Angular - m_J=0
configs.corr{1}.type.comp=1;           
configs.corr{1}.type.coord='angular';
configs.corr{1}.type.opt=[];
configs.corr{1}.lim=[0,0.2;0,pi];
configs.corr{1}.nBin=[21,501];

% 2) Single-halo cart CL - m_J=0
configs.corr{2}.type.comp=1;
configs.corr{2}.type.coord='cart';
configs.corr{2}.type.opt='CL';
configs.corr{2}.lim=0.2*repmat([0,1],[3,1]);
configs.corr{2}.nBin=11*[1,1,1];   % number of bins

% 3) Single-halo cart BB - m_J=0
configs.corr{3}.type.comp=1;
configs.corr{3}.type.coord='cart';
configs.corr{3}.type.opt='BB';
configs.corr{3}.lim=0.2*repmat([-1,1],[3,1]);
configs.corr{3}.nBin=11*[1,1,1];   % number of bins
