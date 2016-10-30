% Main script for the analysis of distinguishable s-wave scattered halos
% DKS

%% USER CONFIG
% GENERAL
use_saved_data=1;   %if false will remake the fully processed data files used in analysis
use_txy=1;          %if false will remake the txy_forc files

verbose=2;

%IN/OUTPUTS
% files
configs.files.path='..\data\multihalos\multihalos_';    % path to unindexed data file (e.g. 'a\b\datadir\datafile')
configs.files.id=5:10;         % file id numbers to use for analysis
configs.files.minCount=100;     % min counts to use for analysis

%% %%%%%%%%%%%%%%%%%%%%%%%%%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
close all; clc;

%initalize variables
missingfiles=zeros(length(configs.files.id),1);     % missing dld/txy file
filestotxy=zeros(length(configs.files.id),1);       % dld files processed to txy
lowcountfiles=zeros(length(configs.files.id),1);    % files with too few counts (skipped in analysis)
importokfiles=zeros(length(configs.files.id),1);    % successfully imported files


%% Load processed data
vars_saved = {'halo_centered_cells','windows','files','bec_bragg',...
    'halo_radius','counts_plus','counts_minus','counts_zero',...
    'lowcountfiles','missingfiles','filestotxy','bec_or_bragg_zero','importokfiles',...
    'all_points'};

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
        elseif ~isequal(configs,S_data.configs)
            % Check consistency of requested analysis configs and saved data
            warning('Previously saved data contains settings (windows, files) different to requested - setting use_saved_data=0 and continue...');
            use_saved_data=0;
        else
            % OK load all vars in file
            if verbose>0, disp('Loading saved data...');, end;
            
            % Check saved file is not empty
            if isempty(vertcat(S_data.halo_centered_cells{:}))
                error('Saved file is empty - terminating program. Try recreating data file with different configs');
            end
            
            % okay to load everything
            load([configs.files.path,'data.mat']);  
        end
        
        clear S_data;
    end
end

%% Create TXY- from raw DLD files
if ~use_saved_data && ~use_txy
    if verbose>0, disp('Creating TXY files...');, end;
    for i=configs.files.id
        i_str=num2str(i);
        
        %check for source - raw DLD file
        if ~fileExists([configs.files.path,i_str,'.txt'])
            missingfiles(i)=1;
            if verbose>0
                warning(['use_txy=0 but raw DLD data does not exist - file #', i_str]);
            end
            continue;
        end
        
        % create TXY file from the raw DLD file
        dld_raw_to_txy(configs.files.path,i,i);
        filestotxy(n)=1;
        
        if verbose>1
            disp(['TXY file for #',i_str,' is created.']);
        end
    end
end

%% Process raw TXY data
% This is where the raw TXY data (from DLD) must be processed to do further
% analysis
% capture the distinguishable halos
% find properties of halo, etc
% SAVE the refined data and results
% TODO: name vars appropriately and update the vars_saved, etc above

%% Correlation analysis
% find correlations