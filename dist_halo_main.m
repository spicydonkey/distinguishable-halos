% Main script for the analysis of distinguishable s-wave scattered halos
% DKS

clear all; close all; clc;

%% USER CONFIG
% GENERAL
use_saved_data=1;   %if false will remake the fully processed data files used in analysis
use_txy=1;          %if false will remake the txy_forc files

verbose=2;

% PLOTS
% TODO

% IN/OUTPUTS
% files -  data file
usrconfigs.files.path='..\data\test\d';    % path to unindexed data file (e.g. 'a\b\datadir\datafile')
usrconfigs.files.id=1:100;         % file id numbers to use for analysis
usrconfigs.files.minCount=100;     % min counts to use for analysis

% windows - captures all the interesting counts in range: [] will not crop
% TODO XY axis is unclear
usrconfigs.window.all_T=[4.905,4.92];       % T [s]
usrconfigs.window.all{1}=[20.678,20.726];   % Z [m] T-window will overide Z
usrconfigs.window.all{2}=[-20e-3,18e-3];    % X [m]
usrconfigs.window.all{3}=[-10e-3,17e-3];    % Y [m]

% DIST HALO PARAMS: params specific to data processing and analysis
usrconfigs.bec.pos{1}=[20.701,5e-3,3e-3];   % approx condensate locations (z,x,y)
usrconfigs.bec.Rmax{1}=7e-3;  % max condensate sph radius
usrconfigs.bec.pos{2}=[20.701,-7.3e-3,6.8e-3];
usrconfigs.bec.Rmax{2}=6e-3;


%% %%%%%%%%%%%%%%%%%%%%%%%%%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

%initalize variables
configs=usrconfigs;    % create an alias to avoid overwriting user's config
missingfiles=zeros(length(configs.files.id),1);     % missing dld/txy file
filestotxy=zeros(length(configs.files.id),1);       % dld files processed to txy
lowcountfiles=zeros(length(configs.files.id),1);    % files with too few counts (skipped in analysis)
importokfiles=zeros(length(configs.files.id),1);    % successfully imported files


%% Load processed data
vars_saved = {'usrconfigs','configs'...
    'halo_data','bec_data',...
    'lowcountfiles','missingfiles','filestotxy','importokfiles',...
    };  % variables important in halo analysis - saved to a .mat file

if use_saved_data     
    if ~fileExists([configs.files.path,'data.mat'])
        warning(['Saved data "', [configs.files.path,'data.mat'], '" does not exist - setting use_saved_data=0 and continue...']);
        use_saved_data=0;
    else
        % Error checks on prev saved file
        S_data = load([configs.files.path,'data.mat']);
        
        is_complete=1;
        for i=1:length(vars_saved)
            is_complete=is_complete&isfield(S_data,vars_saved{i});
        end
        if ~is_complete
            % Completeness of information
            warning('Previously saved file is incomplete - setting use_saved_data=0 and continue...');
            use_saved_data=0;
        elseif ~isequal(usrconfigs,S_data.usrconfigs)
            % Check consistency of requested analysis configs and saved data
            warning('Previously saved data contains settings (windows, files) different to requested - setting use_saved_data=0 and continue...');
            use_saved_data=0;
        else
            % OK load all vars in file
            if verbose>0, disp('Loading saved data...');, end;
            
            % Check saved file is not empty
            if isempty(vertcat(S_data.halo_data{:}))
                error('Saved file is empty - terminating program. Try recreating data file with different configs');
            end
            
            % okay to load everything
            load([configs.files.path,'data.mat']);  
        end
        
        clear S_data;
    end
end


%% Pre-processing
if ~use_saved_data
    %% Prepare TXY files
    % TXY files (usually with suffix _txy_forc and #ID attached to the core file name)
    %   are built by processing raw DLD output files from the TDC
    % TXY files should contain complete result of experiment re. atom positions
    % A number of checks are performed in this stage to filter out missing
    % files and shots of low count, etc.
    if ~use_txy
        % Create TXY from raw DLD files
        if verbose>0, disp('Creating TXY files...');, end;
        for i=1:length(configs.files.id)
            %check for source - raw DLD file
            if ~fileExists([configs.files.path,num2str(configs.files.id(i)),'.txt'])
                missingfiles(i)=1;
                if verbose>0
                    warning(['use_txy=0 but raw DLD data does not exist - file #', num2str(configs.files.id(i))]);
                end
                continue;
            end
            
            % create TXY file from the raw DLD file
            dld_raw_to_txy(configs.files.path,configs.files.id(i),configs.files.id(i));
            filestotxy(i)=1;
            
            if verbose>1
                disp(['TXY file for #',num2str(configs.files.id(i)),' is created.']);
            end
        end
    else
        % Using pre-built TXY files
        if verbose>0, disp('Using pre-built TXY files...');, end;
        % Check for missing files
        for i=1:length(configs.files.id)
            if ~fileExists([configs.files.path,'_txy_forc',num2str(configs.files.id(i)),'.txt'])
                % missing TXY file
                missingfiles(i)=1;
                if verbose>0
                    warning(['TXY file for #',num2str(configs.files.id(i)),' is missing.']);
                end
            end
        end
    end
    
    % Check for low counts
    if verbose>0,disp('Checking for low counts in TXY files...');,end;
    for i=1:length(configs.files.id)
        if fileExists([configs.files.path,'_txy_forc',num2str(configs.files.id(i)),'.txt'])
            txy_temp=txy_importer(configs.files.path,configs.files.id(i));    % import data in this txy to memory
            % Check count in this file
            if size(txy_temp,1)<configs.files.minCount
                lowcountfiles(i)=1;
                if verbose>0
                    warning(['Low count detected in file #',num2str(configs.files.id(i))]);
                end
            end
        end
    end
    
    % Summary of TXY-preprocessing
    importokfiles=~(missingfiles|lowcountfiles);
    if verbose>1
        disp('==============SUMMARY==============');
        disp([int2str(sum(importokfiles)),' imported ok']);
        disp([int2str(sum(lowcountfiles)),' files with low counts']);
        disp([int2str(sum(missingfiles)),' files missing']);
        disp([int2str(sum(filestotxy)),' files converted to txy']);
        disp('===================================');
    end
    
    %% Process raw TXY data
    % This is where the raw TXY data (from DLD) must be processed to do further
    % analysis
    % capture the distinguishable halos
    % find properties of halo, etc
    % SAVE the refined data and results
    % TODO: name vars appropriately and update the vars_saved, etc above
    
    configs.files.idok=configs.files.id(importokfiles);     % ID's for OK files
    
    if verbose>0,disp('Processing TXY files to generate data for analysis...');,end;
    [halo_data,bec_data]=distinguish_halo(configs,verbose);  % all the data processing on raw-TXY to generate halo data
    
    %% Save processed data
    if verbose>0,disp('Saving data..');,end;
    % delete existing data
    if fileExists([configs.files.path,'data.mat'])
        delete([configs.files.path,'data.mat']);    
    end
    % save data
    save([configs.files.path,'data.mat'],vars_saved{1});    % must create *.mat without append option
    for i = 1:length(vars_saved)
        if ~exist(vars_saved{i},'var')
            warning(['Variable "',vars_saved{i},'" does not exist.']);
            continue
        end
        save([configs.files.path,'data.mat'],vars_saved{i},'-append');
    end
end


%% DEBUG
% % halo_collate = vertcat(halo_data{:});       % all halos collated
% figure(99);
% % scatter3(halo_collate(:,2),halo_collate(:,3),halo_collate(:,1),1,'k.');
% rem_collate = vertcat(
% set(gcf,'Color',[1 1 1]);
% axis vis3d;
% axis equal;
% xlabel('X'); ylabel('Y'); zlabel('Z');

%% Correlation analysis
% find correlations