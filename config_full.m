% Configuration file for distinguishable_halos.m

%%% GENERAL
do_corr_analysis=1;
force_all_stages=1;    % force all the stages to run (useful for debug)
verbose=3;

%%% FLAGS
configs.flags.savedata=1;       % TODO - req'd currently since each stage passes data by save/load to disk
configs.flags.graphics=0;       % toggle to control graphics/plotting options

%%% Raw data handling
% files -  data file
configs.files.path='C:\Users\HE BEC\Documents\lab\halo_analysis\data\dist_halo\4_separated_lownum\d';    % path to unindexed data file (e.g. 'a\b\datadir\$DATA_FNAME_TOKEN$')
configs.files.id=1:3000;         % file id numbers to use for analysis
configs.files.minCount=100;     % min counts to use for analysis

% Detector/trap alignment
configs.rot_angle=0.61;

% TXY window - region of interest ( [] --> no crop )
configs.window{1}=[4.907,4.918];      % T [s]
configs.window{2}=[-30e-3,20e-3];    % X [m]
configs.window{3}=[-15e-3,20e-3];    % Y [m]

%%% HALO PARAMS: BEC counts + oscillation removal for broad capture of halos
configs.bec.pos{1}=[20.7024,4.74e-3,2.72e-3];   % approx condensate locations (z,x,y)
configs.bec.Rmax{1}=7e-3;    % max condensate sph radius
configs.bec.dR_tail{1}=0.3;	% BEC tail radial frac diff
configs.bec.pos{2}=[20.7005,-7.38e-3,6.55e-3];
configs.bec.Rmax{2}=7e-3;
configs.bec.dR_tail{2}=0.3;

configs.halo.R{1}=11e-3;     % estimated radius of halo
configs.halo.dR{1}=0.15;      % broad radial mask fractional width (in/out)
configs.halo.R{2}=10e-3;
configs.halo.dR{2}=0.15;

configs.halo.zcap=0.7;   % z-cutoff (fractional wrt radius)

%%% MISCELLANEOUS
configs.misc.vel_z=9.8*0.416;    % atom free-fall vert v at detector hit for T-to-Z conversion;

%%% correlation analysis
% 1. Cross-halo rad/angular correlations
analysis.corr.type{1}.comp=[1,2];           % components to analysis: cross halo 1,2
analysis.corr.type{1}.coord='angular';      % angular coordinate
analysis.corr.type{1}.opt=[];               % 'angular' has no optional feature atm
analysis.corr.lim{1}{1}=0.2*[-1,1];  % bin limits - radial separation
analysis.corr.lim{1}{2}=pi*[0.7,1];      % bin limits - angular separation
analysis.corr.nBin{1}=[11,31];          % number of bins

% 2. Cross-halo cartesian BB-correlations
analysis.corr.type{2}.comp=[1,2];           % components to analysis: cross halo 1,2
analysis.corr.type{2}.coord='cart';         % Cartesian (ZXY)
analysis.corr.type{2}.opt='BB';             % BB / CL
analysis.corr.lim{2}{1}=0.8*[-1,1]; % bin limits - Z
analysis.corr.lim{2}{2}=0.8*[-1,1]; % bin limits - X
analysis.corr.lim{2}{3}=0.8*[-1,1]; % bin limits - Y
analysis.corr.nBin{2}=[51,13,13];   % number of bins

% 3. Single-halo (1) rad/angular correlations
analysis.corr.type{3}.comp=1;           % single component
analysis.corr.type{3}.coord='angular';
analysis.corr.type{3}.opt=[];
analysis.corr.lim{3}{1}=0.3*[-1,1];  % bin limits - radial separation
analysis.corr.lim{3}{2}=[0,pi];      % bin limits - angular separation
analysis.corr.nBin{3}=[11,101];          % number of bins

% 4. Single-halo (1) cartesian BB-correlations
analysis.corr.type{4}.comp=1;
analysis.corr.type{4}.coord='cart';
analysis.corr.type{4}.opt='BB';
analysis.corr.lim{4}{1}=0.8*[-1,1]; % bin limits - Z
analysis.corr.lim{4}{2}=0.8*[-1,1]; % bin limits - X
analysis.corr.lim{4}{3}=0.8*[-1,1]; % bin limits - Y
analysis.corr.nBin{4}=[51,13,13];   % number of bins

% 5. Cross-halo cartesian CL-correlations
analysis.corr.type{5}.comp=[1,2];
analysis.corr.type{5}.coord='cart';
analysis.corr.type{5}.opt='CL';
analysis.corr.lim{5}{1}=0.8*[-1,1]; % bin limits - Z
analysis.corr.lim{5}{2}=0.8*[-1,1]; % bin limits - X
analysis.corr.lim{5}{3}=0.8*[-1,1]; % bin limits - Y
analysis.corr.nBin{5}=[51,13,13];   % number of bins

% 6. Single-halo (1) cartesian CL-correlations
analysis.corr.type{6}.comp=1;
analysis.corr.type{6}.coord='cart';
analysis.corr.type{6}.opt='CL';
analysis.corr.lim{6}{1}=0.8*[-1,1]; % bin limits - Z
analysis.corr.lim{6}{2}=0.8*[-1,1]; % bin limits - X
analysis.corr.lim{6}{3}=0.8*[-1,1]; % bin limits - Y
analysis.corr.nBin{6}=[51,13,13];   % number of bins

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

%%% ALGORITHM CONFIGS
% DO NOT ADJUST
configs.files.saveddata=[configs.files.path,'_data.mat'];     % path to store saved data
configs.files.archive=[configs.files.path,'_archive'];   % dir to archive folder
configs.files.dirout=[configs.files.path,'_output'];      % output directory

do_next=force_all_stages;  % flag for executing stages for efficiency (do not change)

%% END