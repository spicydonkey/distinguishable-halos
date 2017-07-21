% Configuration file for distinguishable_halos.m

%%% FLAGS
configs.flags.do_corr_analysis=1;
    configs.flags.do_corr_err=1;
configs.flags.force_all_stages=0;    % force all the stages to run (useful for debug)
configs.flags.verbose=2;
configs.flags.savedata=1;       % TODO - req'd currently since each stage passes data by save/load to disk
configs.flags.archive_txy=1;        % archives loaded TXY as .mat file for future reuse
configs.flags.graphics=1;       % toggle to control graphics/plotting options
configs.flags.build_txy=1;

%%% MISCELLANEOUS
configs.misc.vel_z=9.8*0.416;    % atom free-fall vert v at detector hit for T-to-Z conversion;
vz=configs.misc.vel_z;
configs.misc.deadtime=300e-9;    %100e-9;    % data acquired after disc adjustment

%% FILES
configs.files.path='\\AMPLPC29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\xstate_collision_highv\11_atr_set\d';

% WARNING: MODIFYING BELOW DIR SETTINGS ARE NOT RECOMMENDED
configs.files.dir_data=fileparts(configs.files.path);    % fullpath to data directory
configs.files.archive=fullfile(configs.files.dir_data,'archive');   % dir to archive folder
configs.files.dirout=fullfile(configs.files.dir_data,'output');      % output directory (will be time-stamped)


%% LOAD
configs.load.version=1;         % TXY load stage version number

% file ID and simple pass/fail
configs.load.id=1:3045;         % file id numbers to use for analysis
configs.load.mincount=600;         % min counts in window - 0 for no min
configs.load.maxcount=1000;          % max counts in window - Inf for no max

% Detector/trap alignment
configs.load.rot_angle=0.61;

% TXY window - region of interest ( [] --> no crop )
configs.load.window{1}=[0.44,0.48];      % T [s]
configs.load.window{2}=[-35e-3,35e-3];    % X [m]
configs.load.window{3}=[-35e-3,35e-3];    % Y [m]

%% HALO
%%% HALO PARAMS: BEC counts + oscillation removal for broad capture of halos
configs.bec.pos{1}=[1.8926e+00,-3.0e-3,4.1e-3];   % approx condensate locations (z,x,y)
configs.bec.Rmax{1}=6e-3;      % max condensate sph radius
configs.bec.dR_tail{1}=0.8;     % BEC tail radial frac diff
configs.bec.pos{2}=[1.8612e+00,-1.9e-3,1.1e-3];
configs.bec.Rmax{2}=6e-3;
configs.bec.dR_tail{2}=0.8;

configs.halo.R{1}=26e-3;     % estimated radius of halo
configs.halo.dR{1}=0.25;      % broad radial mask fractional width (in/out)
configs.halo.R{2}=24e-3;
configs.halo.dR{2}=0.3;

configs.halo.zcap=0.8;

configs.halo.boost{1}=zeros(1,3);
configs.halo.boost{2}=[0.035,-0.015,-0.02];

%% CORRELATION ANALYSIS
% ERROR ANALYSIS
configs.error.ncluster=5;

% 1) X-halo Cart BB
configs.corr{1}.type.comp=[1,2];           % components to analysis: cross halo 1,2
configs.corr{1}.type.coord='cart';         % Cartesian (ZXY)
configs.corr{1}.type.opt='BB';             % BB / CL
configs.corr{1}.lim=0.2*repmat([-1,1],[3,1]);
configs.corr{1}.nBin=15*[1,1,1];   % number of bins
% 
% % 2) X-halo Cart CL
% configs.corr{2}.type.comp=[1,2];
% configs.corr{2}.type.coord='cart';
% configs.corr{2}.type.opt='CL';
% configs.corr{2}.lim=0.2*repmat([0,1],[3,1]);
% configs.corr{2}.nBin=11*[1,1,1];   % number of bins
% 
% % 3) Single-halo cart CL - m_J=0
% configs.corr{3}.type.comp=1;
% configs.corr{3}.type.coord='cart';
% configs.corr{3}.type.opt='CL';
% configs.corr{3}.lim=0.2*repmat([0,1],[3,1]);
% configs.corr{3}.nBin=11*[1,1,1];   % number of bins
% 
% % 4) Single-halo cart CL - m_J=1
% configs.corr{4}.type.comp=2;
% configs.corr{4}.type.coord='cart';
% configs.corr{4}.type.opt='CL';
% configs.corr{4}.lim=0.2*repmat([0,1],[3,1]);
% configs.corr{4}.nBin=11*[1,1,1];
% 
% % 5) Single-halo cart BB - m_J=0
% configs.corr{5}.type.comp=1;
% configs.corr{5}.type.coord='cart';
% configs.corr{5}.type.opt='BB';
% configs.corr{5}.lim=0.2*repmat([-1,1],[3,1]);
% configs.corr{5}.nBin=11*[1,1,1];   % number of bins
% 
% % 6) Single-halo cart BB - m_J=1
% configs.corr{6}.type.comp=2;
% configs.corr{6}.type.coord='cart';
% configs.corr{6}.type.opt='BB';
% configs.corr{6}.lim=0.2*repmat([-1,1],[3,1]);
% configs.corr{6}.nBin=11*[1,1,1];   % number of bins
% 
% % 7) X-halo Angular
% configs.corr{7}.type.comp=[1,2];
% configs.corr{7}.type.coord='angular';
% configs.corr{7}.type.opt=[];
% configs.corr{7}.lim=[0,0.1;0,pi];
% configs.corr{7}.nBin=[21,501];
% 
% % 8) Single-halo Angular - m_J=0
% configs.corr{8}.type.comp=1;           
% configs.corr{8}.type.coord='angular';
% configs.corr{8}.type.opt=[];
% configs.corr{8}.lim=[0,0.1;0,pi];
% configs.corr{8}.nBin=[21,501];
% 
% % 9) Single-halo Angular - m_J=1
% configs.corr{9}.type.comp=2;           
% configs.corr{9}.type.coord='angular';
% configs.corr{9}.type.opt=[];
% configs.corr{9}.lim=[0,0.1;0,pi];
% configs.corr{9}.nBin=[21,501];