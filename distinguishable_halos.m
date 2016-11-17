% Distinguishable halo
% DKS 16/11/2016

%% CONFIGS
%   <<Group into each stages>>
% [x] raw - raw-data handling
% [-] preproc - broad capture
% [] postproc - manipulations
% [] analysis

%%% GENERAL
force_all_stages=0;    % force all the stages to run (useful for debug)

verbose=2;

%%% Raw data handling
% files -  data file
usrconfigs.files.path='C:\Users\HE BEC\Documents\lab\halo_analysis\data\dist_halo\4_separated_lownum\d';    % path to unindexed data file (e.g. 'a\b\datadir\$DATA_FNAME_TOKEN$')
usrconfigs.files.id=1:100;         % file id numbers to use for analysis
usrconfigs.files.minCount=100;     % min counts to use for analysis

% TXY window - region of interest ( [] --> no crop )
usrconfigs.window{1}=[4.905,4.92];      % T [s]
usrconfigs.window{2}=[-20e-3,18e-3];    % X [m]
usrconfigs.window{3}=[-10e-3,17e-3];    % Y [m]

%%% HALO PARAMS: BEC counts + oscillation removal for broad capture of halos
usrconfigs.bec.pos{1}=[20.7024,4.74e-3,2.72e-3];   % approx condensate locations (z,x,y)
usrconfigs.bec.Rmax{1}=7e-3;    % max condensate sph radius
usrconfigs.bec.dR_tail{1}=0.3;	% BEC tail radial frac diff
usrconfigs.bec.pos{2}=[20.7005,-7.38e-3,6.55e-3];
usrconfigs.bec.Rmax{2}=7e-3;
usrconfigs.bec.dR_tail{2}=0.3;

usrconfigs.halo.R{1}=11e-3;     % estimated radius of halo
usrconfigs.halo.dR{1}=0.3;      % broad radial mask fractional width (in/out)
usrconfigs.halo.R{2}=10e-3;
usrconfigs.halo.dR{2}=0.3;

% MISCELLANEOUS
usrconfigs.misc.vel_z=9.8*0.430;    % atom free-fall vert v at detector hit for T-to-Z conversion;

%%% ALGORITHM CONFIGS
% DO NOT ADJUST
usrconfigs.files.saveddata=[usrconfigs.files.path,'_data.mat'];     % path to store saved data
usrconfigs.files.archive=[usrconfigs.files.path,'_archive'];   % dir to archive folder

do_next=force_all_stages;  % flag for executing stages for efficiency (do not change)

configs=usrconfigs;     % copy to protect usrconfigs

%% Load raw-data DLD/TXY
%   <<Skip if already done and configs are same>>
%   Process rawDLD/TXY archive to:
%     - get counts separated by shot# (ditch low-counts)
%     - save (ZXY)-counts to .mat file: $txy_all
%     - $files contains summary of each TXY file

% Check if can be skipped
if ~do_next
    % check for existing saved data file
    if exist(configs.files.saveddata,'file')
        % same configs?: files, window
        % TODO - will throw error if $configs is not a variable in data
        S_temp=load(configs.files.saveddata,'configs');  % load configs from prev data
        % TODO - will throw error if $files and $window are not fields of $configs
        if ~(isequal(S_temp.configs.files,configs.files)&&isequal(S_temp.configs.window,configs.window))
            warning('Existing data has different configs. Setting do_next=1.');
            do_next=1;
        end
    else
        warning('No existing data. Setting do_next=1.');
        do_next=1;
    end
    clear S_temp;
end

if do_next
    loadExpData(configs,verbose);
end

%% Preprocess: Get halo counts
%   <<Skip if already done AND configs are same>>
%   Process loaded counts to:
%     - locate condensate accurately in each shot
%     - capture all counts in the halo with approx centre location and
%       a broad radial mask
%     - collate all shots with BEC oscillations cancelled (esp. halo)
%     - save to file (think about not saving BEC's - too many counts)


captureHalo(configs,verbose);



%% Fit halos
%   Process the broadly captured halos (ZXY):
%     - fit the point cloud to some shape (sphere/oblate spheroid/tri-axial
%       ellipsoid, etc.)
%     - get "centre", "radii", etc. for each halo
%     - reference to centre of mass (halo centres)
%     - extract well-fitted counts (tight capture)
%     - save to file (fitted halo with shape and fit params)

%% Postprocess: Manipulate halos
%   <<Skip if already done and configs are same>>
%   Manipulate the halos to account for distortions:
%     - e.g. spin in Z-axis, isotropic/anisotropic scaling, etc.
%     - save to file

%% Analysis
%   <<Skip if already done and configs are same>>
%   Correlation analysis (in angular, cartesian):
%     - cross-halo
%     - single-halo