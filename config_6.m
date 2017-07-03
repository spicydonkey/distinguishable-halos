% Configuration file for distinguishable_halos.m

%%% FLAGS
configs.flags.do_corr_analysis=1;
configs.flags.force_all_stages=0;    % force all the stages to run (useful for debug)
configs.flags.verbose=2;
configs.flags.savedata=1;       % TODO - req'd currently since each stage passes data by save/load to disk
configs.flags.archive_txy=1;        % archives loaded TXY as .mat file for future reuse
configs.flags.graphics=1;       % toggle to control graphics/plotting options
configs.flags.build_txy=1;

%%% MISCELLANEOUS
configs.misc.vel_z=9.8*0.416;    % atom free-fall vert v at detector hit for T-to-Z conversion;
vz=configs.misc.vel_z;


%% FILES
configs.files.path='\\AMPLPC29\He BEC Archive\EXPERIMENT-DATA\xstate_mom_corr\90deg_raman_beams\6_mednum\d';

% WARNING: MODIFYING BELOW DIR SETTINGS ARE NOT RECOMMENDED
configs.files.dir_data=fileparts(configs.files.path);    % fullpath to data directory
configs.files.archive=fullfile(configs.files.dir_data,'archive');   % dir to archive folder
configs.files.dirout=fullfile(configs.files.dir_data,'output');      % output directory (will be time-stamped)


%% LOAD
configs.load.version=1;         % TXY load stage version number

% file ID and simple pass/fail
% configs.load.id=1:1500;         % file id numbers to use for analysis
% configs.load.mincount=1500;         % min counts in window - 0 for no min
% configs.load.maxcount=1700;          % max counts in window - Inf for no max
configs.load.id=1:2848;         % file id numbers to use for analysis
configs.load.mincount=1400;         % min counts in window - 0 for no min
configs.load.maxcount=1800;          % max counts in window - Inf for no max

% Detector/trap alignment
configs.load.rot_angle=0.61;

% TXY window - region of interest ( [] --> no crop )
configs.load.window{1}=[2.44,2.48];      % T [s]
configs.load.window{2}=[-35e-3,35e-3];    % X [m]
configs.load.window{3}=[-35e-3,35e-3];    % Y [m]

%% HALO
%%% HALO PARAMS: BEC counts + oscillation removal for broad capture of halos
configs.bec.pos{1}=[2.464*vz,-3.6e-3,4.3e-3];   % approx condensate locations (z,x,y)
configs.bec.Rmax{1}=6e-3;      % max condensate sph radius
configs.bec.dR_tail{1}=1;     % BEC tail radial frac diff
configs.bec.pos{2}=[2.456*vz,-2.2e-3,0.8e-3];
configs.bec.Rmax{2}=6e-3;
configs.bec.dR_tail{2}=1;

configs.halo.R{1}=26e-3;     % estimated radius of halo
configs.halo.dR{1}=0.25;      % broad radial mask fractional width (in/out)
configs.halo.R{2}=24e-3;
configs.halo.dR{2}=0.25;

configs.halo.zcap=0.8;   % z-cutoff (fractional wrt radius)

% TODO - does boost need to be optimised for different g2 analysis?
%   currently SINGLE boost applied to halo2 to obtain best g2_01_BB
configs.halo.boost{1}=zeros(1,3);
configs.halo.boost{2}=[0.035,-0.002,-0.025];

%% CORRELATION ANALYSIS
% 1) X-halo Cart BB
configs.corr{1}.type.comp=[1,2];           % components to analysis: cross halo 1,2
configs.corr{1}.type.coord='cart';         % Cartesian (ZXY)
configs.corr{1}.type.opt='BB';             % BB / CL
configs.corr{1}.lim(1,:)=0.2*[-1,1]; % bin limits - Z
% configs.corr{1}.lim(1,:)=[-0.12,0.17]; % bin limits - Z
configs.corr{1}.lim(2,:)=0.2*[-1,1]; % bin limits - X
% configs.corr{1}.lim(2,:)=[-0.1,0.3]; % bin limits - X
configs.corr{1}.lim(3,:)=0.2*[-1,1]; % bin limits - Y
% configs.corr{1}.lim(3,:)=[-0.3,0.1]; % bin limits - Y
configs.corr{1}.nBin=15*[1,1,1];   % number of bins

% 2) X-halo Cart CL
configs.corr{2}.type.comp=[1,2];
configs.corr{2}.type.coord='cart';
configs.corr{2}.type.opt='CL';
configs.corr{2}.lim=0.3*repmat([0,1],[3,1]);
configs.corr{2}.nBin=7*[1,1,1];   % number of bins

% 3) Single-halo cart CL - m_J=0
configs.corr{3}.type.comp=1;
configs.corr{3}.type.coord='cart';
configs.corr{3}.type.opt='CL';
configs.corr{3}.lim=0.2*repmat([0,1],[3,1]);
configs.corr{3}.nBin=11*[1,1,1];   % number of bins

% 3) Single-halo cart CL - m_J=1
configs.corr{4}.type.comp=2;
configs.corr{4}.type.coord='cart';
configs.corr{4}.type.opt='CL';
configs.corr{4}.lim=0.2*repmat([0,1],[3,1]);
configs.corr{4}.nBin=11*[1,1,1];

% 4) Single-halo cart BB - m_J=0
configs.corr{5}.type.comp=1;
configs.corr{5}.type.coord='cart';
configs.corr{5}.type.opt='BB';
configs.corr{5}.lim=0.2*repmat([-1,1],[3,1]);
configs.corr{5}.nBin=11*[1,1,1];   % number of bins

% 5) Single-halo cart BB - m_J=1
configs.corr{6}.type.comp=2;
configs.corr{6}.type.coord='cart';
configs.corr{6}.type.opt='BB';
configs.corr{6}.lim=0.2*repmat([-1,1],[3,1]);
configs.corr{6}.nBin=11*[1,1,1];   % number of bins


%     % #. <TEMPLATE ANGULAR>
%     analysis.corr.type{#}.comp=[1,2];           % components to analysis: cross halo 1,2
%     analysis.corr.type{#}.coord='angular';      % angular coordinate
%     analysis.corr.type{#}.opt=[];               % 'angular' has no optional feature atm
%     analysis.corr.lim{#}{1}=0.3*[-1,1];  % bin limits - radial separation
%     analysis.corr.lim{#}{2}=[0,pi];      % bin limits - angular separation
%     analysis.corr.nBin{#}=[11,51];          % number of bins
%
%     % #. <TEMPLATE CARTESIAN>
%     analysis.corr.type{#}.comp=[1,2];           % components to analysis: cross halo 1,2
%     analysis.corr.type{#}.coord='cart';         % Cartesian (ZXY)
%     analysis.corr.type{#}.opt='BB';             % BB / CL
%     analysis.corr.lim{#}{1}=0.8*[-1,1]; % bin limits - Z
%     analysis.corr.lim{#}{2}=0.8*[-1,1]; % bin limits - X
%     analysis.corr.lim{#}{3}=0.8*[-1,1]; % bin limits - Y
%     analysis.corr.nBin{#}=[51,13,13];   % number of bins