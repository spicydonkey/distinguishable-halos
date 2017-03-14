% Distinguishable halo
% DKS 16/11/2016

%%% CONFIGS (in set_config.m)
% [x] raw - raw-data handling
% [x] preproc - broad capture
% [] fit - fit halo
% [-] postproc - manipulations
% [x] analysis

% TODO - replace saving and loading data files by using functions to return
% nice results

clear all; clc; close all;

%% Load configurations
config_full;     % set configurations


%% MAIN
%%%% 
t_main_start=tic;

% CHECKS
% Output management
if ~isdir(configs.files.dirout)
    warning(['output directory "',configs.files.dirout,'" does not exist. Creating directory...']);
    mkdir(configs.files.dirout);
end


%% Load raw-data DLD/TXY
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
            warning('Raw data: Existing data has different configs. Setting do_next=1.');
            do_next=1;
        end
    else
        warning('No existing data. Setting do_next=1.');
        do_next=1;
    end
    clear S_temp;
end

if do_next
    [txy,files_out]=loadExpData(configs,verbose);
end


%% Preprocess: Get halo counts
%   Process loaded counts to:
%     - [x] locate condensate accurately in each shot
%     - [x] capture all counts in the halo with approx centre location and
%       a broad radial mask
%     - [x] collate all shots with BEC oscillations cancelled (esp. halo)
%     - [x] save to file (think about not saving BEC's - too many counts)

% Check if broad halo capture can be skipped
if ~do_next
    % same configs?: files, window
    % TODO - will throw error if $configs is not a variable in data
    S_temp=load(configs.files.saveddata,'configs');  % load configs from prev data
    % TODO - will throw error if bec,halo,or misc are not fields of $configs
    if ~(isequal(S_temp.configs.bec,configs.bec)&&...
            isequal(S_temp.configs.halo,configs.halo)&&...
            isequal(S_temp.configs.misc,configs.misc))
        warning('Halo capture: Existing data has different configs. Setting do_next=1.');
        do_next=1;
    end
    clear S_temp;
end

if do_next
    [halo,bec,culled,errflag]=captureHalo(configs,verbose);
end


%% Fit halos
%   <<Currently just returns the same halo>>
%   Process the broadly captured halos (ZXY):
%     - [] fit the point cloud to some shape (sphere/oblate spheroid/tri-axial
%       ellipsoid, etc.)
%     - [] get "centre", "radii", etc. for each halo
%     - [] reference to centre of mass (halo centres)
%     - [] extract well-fitted counts (tight capture)
%     - [] save to file (fitted halo with shape and fit params)

% No user settable params yet
if do_next
    [halo_fit]=fitHalo(configs,verbose);
end


%% Postprocess: Transform halos
%   Transform halos to account for distortions:
%     - e.g. spin in Z-axis, isotropic/anisotropic scaling, etc.
%     - save to file
%   * creates the 'zxy' cell array of finalised counts

% no user settable params yet
if do_next
    [zxy_transform]=transformHalo(configs,verbose);
end


%% Analysis
%   Correlation analysis (in angular, cartesian):
%     - cross-halo
%     - single-halo

if ~do_next
    % check if analysis configs has changed
    S_temp=load(configs.files.saveddata,'analysis');    % load prev analysis configs
    if ~isequal(S_temp.analysis.corr,analysis.corr)
        warning('Correlation analysis: prev saved data has different configs. Setting do_next=1.');
        do_next=1;
    end
    clear S_temp;
end

% Correlation analysis
if do_next&&do_corr_analysis
    [result]=corrTaskManager(analysis,configs,verbose);
end


%% end of code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_main_end=toc(t_main_start);   % end of code
disp('-----------------------------------------------');
fprintf('Total elapsed time (s): %7.1f\n',t_main_end);
disp('===================ALL TASKS COMPLETED===================');