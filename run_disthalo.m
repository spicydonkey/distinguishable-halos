% Distinguishable halo

clear all; clc; close all;

%% DEBUG switches
force_all_stages=0;
verbose=2;

%% Load configurations
% loads configuration params in structure "configs"
config_loadtxy_demo;    % load and collate exp data and crop TXY
config_capturehalo_demo;    % capture halos from collated TXY data

%% Set output directories
configs.files.saveddata=[configs.files.path,'_data.mat']; % path to store saved data
configs.files.archive=[configs.files.path,'_archive'];   % dir to archive folder
configs.files.dirout=[configs.files.path,'_output\'];      % output directory

%% MAIN
%%%% Initialise program
t_main_start=tic;
do_next=force_all_stages;

%% Load raw-data DLD/TXY
% Check if appropriate TXY-data is already available
if ~do_next
    % check for existing saved data file
    path_txy_data=[configs.files.path,'_data.mat']; % path to store saved data
    if exist(path_txy_data,'file')
        % check for differences in configs: files, window
        S_txy_saved=load(path_txy_data,'configs');  % load configs from prev data
        if isfield(S_txy_saved,'configs')
            if ~(isfield(S_txy_saved.configs,'files')&&isfield(S_txy_saved.configs,'window'))
                warning('File %s exists, but "configs" is ill-defined. Setting do_next=1.',path_txy_data);
                do_next=1;
            else
                % Check if params are the same
                if ~(isequal(S_txy_saved.configs.files,configs.files)&&...
                        isequal(S_txy_saved.configs.window,configs.window))
                    warning('Existing TXY data has different configs. Setting do_next=1.');
                    do_next=1;
                else
                    fprintf('TXY data exists: %s\n',path_txy_data);
                end
            end
        else
            warning('File %s exists, but "configs" is not defined. Setting do_next=1.',path_txy_data);
            do_next=1;
        end
    % No existing TXY-data file
    else
        warning('No existing TXY data. Setting do_next=1.');
        do_next=1;
    end
    clear S_txy_saved;
end

if do_next
    txy_all=loadExpData(configs,verbose);
else
    % load existing txy data
    fprintf('Loading existing TXY data...\n');
    load(path_txy_data,'txy_all');
end

%% Capture halos
% Check if broad halo capture can be skipped
if ~do_next
    % check for existing saved data file
    path_halo_data=[configs.files.path,'_halo.mat']; % path to store saved data
    
    S_halo=load(path_halo_data,'configs');  % load configs from prev data
    if isfield(S_halo,'configs')
        S_halo=load(path_halo_data,'configs');      % load configs from prev data
        
        % check for differences in configs
        if ~(isfield(S_halo.configs,'bec')&&isfield(S_halo.configs,'halo'))
            warning('File %s exists, but "configs" is ill-defined. Setting do_next=1.',path_txy_data);
            do_next=1;
        else
            if ~(isequal(S_halo.configs.bec,configs.bec)&&...
                    isequal(S_halo.configs.halo,configs.halo))
                warning('Existing data has different configs. Setting do_next=1.');
                do_next=1;
            else
                fprintf('Halo data exists: %s\n',path_halo_data);
            end
        end
    else
        warning('File %s exists, but "configs" is not defined. Setting do_next=1.',path_halo_data);
        do_next=1;
    end
    
    clear S_halo_data;
end

if do_next
    halo=captureHalo(configs,verbose);
else
    % load existing halo data
    fprintf('Loading existing halo data...\n');
    load(path_halo_data,'halo');
end

% 
% %% Fit halos
% %   <<Currently just returns the same halo>>
% %   Process the broadly captured halos (ZXY):
% %     - [] fit the point cloud to some shape (sphere/oblate spheroid/tri-axial
% %       ellipsoid, etc.)
% %     - [] get "centre", "radii", etc. for each halo
% %     - [] reference to centre of mass (halo centres)
% %     - [] extract well-fitted counts (tight capture)
% %     - [] save to file (fitted halo with shape and fit params)
% 
% % No user settable params yet
% if do_next
%     fitHalo(configs,verbose);
% end
% 
% 
% %% Postprocess: Transform halos
% %   Transform halos to account for distortions:
% %     - e.g. spin in Z-axis, isotropic/anisotropic scaling, etc.
% %     - save to file
% %   * creates the 'zxy' cell array of finalised counts
% 
% % no user settable params yet
% if do_next
%     transformHalo(configs,verbose);
% end
% 
% 
% %% Analysis
% %   Correlation analysis (in angular, cartesian):
% %     - cross-halo
% %     - single-halo
% 
% if ~do_next
%     % check if analysis configs has changed
%     S_temp=load(configs.files.saveddata,'analysis');    % load prev analysis configs
%     if ~isequal(S_temp.analysis.corr,analysis.corr)
%         warning('Correlation analysis: prev saved data has different configs. Setting do_next=1.');
%         do_next=1;
%     end
%     clear S_temp;
% end
% 
% % Correlation analysis
% if do_next&&do_corr_analysis
%     corrTaskManager(analysis,configs,verbose);
% end


%% end of code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_main_end=toc(t_main_start);   % end of code
disp('-----------------------------------------------');
fprintf('Total elapsed time (s): %7.1f\n',t_main_end);
disp('===================ALL TASKS COMPLETED===================');