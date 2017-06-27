% Configuration file for distinguishable_halos.m

%%% FLAGS
configs.flags.do_corr_analysis=0;
configs.flags.force_all_stages=1;    % force all the stages to run (useful for debug)
configs.flags.verbose=2;
configs.flags.savedata=1;       % TODO - req'd currently since each stage passes data by save/load to disk
configs.flags.graphics=1;       % toggle to control graphics/plotting options

%%% MISCELLANEOUS
configs.misc.vel_z=9.8*0.416;    % atom free-fall vert v at detector hit for T-to-Z conversion;
vz=configs.misc.vel_z;

%% LOAD
% files -  data file
% path to unindexed data file (e.g. 'a\b\datadir\$DATA_FNAME_TOKEN$')
configs.load.path='\\AMPLPC29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\xstate_collision_highv\2\d';
configs.load.id=1:5817;         % file id numbers to use for analysis
configs.load.mincount=3000;         % min counts in window - 0 for no min
configs.load.maxcount=5000;          % max counts in window - Inf for no max

%%% DIRECTORIES
configs.files.dir_data=fileparts(configs.load.path);    % fullpath to data directory
configs.files.archive=fullfile(configs.files.dir_data,'archive');   % dir to archive folder
configs.files.dirout=fullfile(configs.files.dir_data,'output');      % output directory (will be time-stamped)

% Detector/trap alignment
configs.load.rot_angle=0.61;

% TXY window - region of interest ( [] --> no crop )
configs.load.window{1}=[0.388,0.42];      % T [s]
configs.load.window{2}=[-35e-3,30e-3];    % X [m]
configs.load.window{3}=[-30e-3,35e-3];    % Y [m]

%% HALO
%%% HALO PARAMS: BEC counts + oscillation removal for broad capture of halos
configs.bec.pos{1}=[0.405*vz,-3.3e-3,4.1e-3];   % approx condensate locations (z,x,y)
configs.bec.Rmax{1}=5e-3;      % max condensate sph radius
configs.bec.dR_tail{1}=2.8;     % BEC tail radial frac diff
configs.bec.pos{2}=[0.403*vz,-3.3e-3,0.0e-3];
configs.bec.Rmax{2}=5e-3;
configs.bec.dR_tail{2}=2.8;

configs.halo.R{1}=25.5e-3;     % estimated radius of halo
configs.halo.dR{1}=0.18;      % broad radial mask fractional width (in/out)
configs.halo.R{2}=21.5e-3;
configs.halo.dR{2}=0.18;

configs.halo.zcap=0.75;   % z-cutoff (fractional wrt radius)

%% CORRELATION ANALYSIS
% 1) X-halo Cart BB
configs.corr{1}.type.comp=[1,2];           % components to analysis: cross halo 1,2
configs.corr{1}.type.coord='cart';         % Cartesian (ZXY)
configs.corr{1}.type.opt='BB';             % BB / CL
configs.corr{1}.lim(1,:)=0.12*[-1,1]; % bin limits - Z
configs.corr{1}.lim(2,:)=0.12*[-1,1]; % bin limits - X
configs.corr{1}.lim(3,:)=0.12*[-1,1]; % bin limits - Y
configs.corr{1}.nBin=[11,11,11];   % number of bins

% % 1. Cross-halo rad/angular correlations
% configs.corr{2}.type.comp=[1,2];           % components to analysis: cross halo 1,2
% configs.corr{2}.type.coord='angular';      % angular coordinate
% configs.corr{2}.type.opt=[];               % 'angular' has no optional feature atm
% configs.corr{2}.lim(1,:)=0.2*[-1,1];  % bin limits - radial separation
% configs.corr{2}.lim(2,:)=pi*[0.7,1];      % bin limits - angular separation
% configs.corr{2}.nBin=[21,41];          % number of bins

% % 2) Single-halo cart CL - m_J=0
% configs.corr{2}.type.comp=1;
% configs.corr{2}.type.coord='cart';
% configs.corr{2}.type.opt='CL';
% configs.corr{2}.lim(1,:)=0.5*[-1,1]; % bin limits - Z
% configs.corr{2}.lim(2,:)=0.5*[-1,1]; % bin limits - X
% configs.corr{2}.lim(3,:)=0.5*[-1,1]; % bin limits - Y
% configs.corr{2}.nBin=[31,31,31];   % number of bins

% % 3) Single-halo cart CL - m_J=1
% configs.corr{3}.type.comp=2;
% configs.corr{3}.type.coord='cart';
% configs.corr{3}.type.opt='CL';
% configs.corr{3}.lim(1,:)=0.5*[-1,1]; % bin limits - Z
% configs.corr{3}.lim(2,:)=0.5*[-1,1]; % bin limits - X
% configs.corr{3}.lim(3,:)=0.5*[-1,1]; % bin limits - Y
% configs.corr{3}.nBin=[31,31,31];   % number of bins

%% correlation analysis
% % 1. Cross-halo rad/angular correlations
% configs.corr.type{1}.comp=[1,2];           % components to analysis: cross halo 1,2
% configs.corr.type{1}.coord='angular';      % angular coordinate
% configs.corr.type{1}.opt=[];               % 'angular' has no optional feature atm
% configs.corr.lim{1}{1}=0.2*[-1,1];  % bin limits - radial separation
% configs.corr.lim{1}{2}=pi*[0.7,1];      % bin limits - angular separation
% configs.corr.nBin{1}=[11,31];          % number of bins

% % 2. Cross-halo cartesian BB-correlations
% configs.corr.type{2}.comp=[1,2];           % components to analysis: cross halo 1,2
% configs.corr.type{2}.coord='cart';         % Cartesian (ZXY)
% configs.corr.type{2}.opt='BB';             % BB / CL
% configs.corr.lim{2}{1}=0.8*[-1,1]; % bin limits - Z
% configs.corr.lim{2}{2}=0.8*[-1,1]; % bin limits - X
% configs.corr.lim{2}{3}=0.8*[-1,1]; % bin limits - Y
% configs.corr.nBin{2}=[51,13,13];   % number of bins
% 
% % 3. Single-halo (1) rad/angular correlations
% configs.corr.type{3}.comp=1;           % single component
% configs.corr.type{3}.coord='angular';
% configs.corr.type{3}.opt=[];
% configs.corr.lim{3}{1}=0.3*[-1,1];  % bin limits - radial separation
% configs.corr.lim{3}{2}=[0,pi];      % bin limits - angular separation
% configs.corr.nBin{3}=[11,101];          % number of bins
% 
% % 4. Single-halo (1) cartesian BB-correlations
% configs.corr.type{4}.comp=1;
% configs.corr.type{4}.coord='cart';
% configs.corr.type{4}.opt='BB';
% configs.corr.lim{4}{1}=0.8*[-1,1]; % bin limits - Z
% configs.corr.lim{4}{2}=0.8*[-1,1]; % bin limits - X
% configs.corr.lim{4}{3}=0.8*[-1,1]; % bin limits - Y
% configs.corr.nBin{4}=[51,13,13];   % number of bins
% 
% % 5. Cross-halo cartesian CL-correlations
% configs.corr.type{5}.comp=[1,2];
% configs.corr.type{5}.coord='cart';
% configs.corr.type{5}.opt='CL';
% configs.corr.lim{5}{1}=0.8*[-1,1]; % bin limits - Z
% configs.corr.lim{5}{2}=0.8*[-1,1]; % bin limits - X
% configs.corr.lim{5}{3}=0.8*[-1,1]; % bin limits - Y
% configs.corr.nBin{5}=[51,13,13];   % number of bins
% 
% % 6. Single-halo (1) cartesian CL-correlations
% configs.corr.type{6}.comp=1;
% configs.corr.type{6}.coord='cart';
% configs.corr.type{6}.opt='CL';
% configs.corr.lim{6}{1}=0.8*[-1,1]; % bin limits - Z
% configs.corr.lim{6}{2}=0.8*[-1,1]; % bin limits - X
% configs.corr.lim{6}{3}=0.8*[-1,1]; % bin limits - Y
% configs.corr.nBin{6}=[51,13,13];   % number of bins

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