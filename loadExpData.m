function [txy_all,files]=loadExpData(configs,verbose)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loads raw-data DLD/TXY and get counts in region of interest (after XY rotation) and 
% save to file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: some inputs are not required + DOCUMENTATION
% INPUT
% configs
%   requires fields:    files.id,files.path, files.minCount (not req)
%                       files.dirout (not req)
%                       window - 3x1 cell of Z,X,Y window limits (crop)
%   optional fields:
%
% OUTPUT
% files: flags for processed files
%   has fields:     build_txy,missing,lowcount,id_ok
%
% MISC
% * configs.savedata: saves the data (txy_all, files) to file + configs used
%
%%%%%%%%%%%%%%%%% LOG %%%%%%%%%%%%%%%%%
% 31/01/17 | DKS | XY rotation implemented (before crop)
% TODO 06/03/17 | DKS | some inputs made optional
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('verbose','var')
    verbose=0;  % default verbose is quiet
end

vars_save={'configs','txy_all','files'};  % a list of variables to save to file

%% MAIN
t_fun_start=tic;

f_id=configs.files.id;      % get files id's to process
f_path=configs.files.path;  % get path to file (no id token)

% check for user-specified flag
if ~isfield(configs,'savedata')
    warning('configs.savedata is unspecified. By default data will not be saved.');
    configs.savedata=false;     % default
end
if configs.savedata
    if ~isfield(configs.files,'dirout')
        % check for MATLAB default folder, where startup.m is located
        dir_home=which('startup.m');
        if isequal(dir_home,'')
            warning('configs.files.dirout is unspecified and home directory cannot be found. Setting configs.savedata=false');
            configs.savedata=false;
            configs.files.dirout='';
        else
            dir_home=fileparts(dir_home);   % extract HOME directory
            warning('configs.files.dirout is unspecified. Setting default to %s.',dir_home);
            configs.files.dirout=strcat(dir_home,'\');    % default output path
        end
    end
else
    warning('Resetting configs.files.dirout to NULL.');
    configs.files.dirout='';    % default
end
dir_output=configs.files.dirout;    % output directory

%%% optional fields
if ~isfield(configs.files,'saveddata')
    if configs.savedata
        dir_saveddata=[configs.files.path,'_data.mat'];
        warning('configs.files.saveddata is unspecified. Setting to default %s.',dir_saveddata);
        configs.files.saveddata=dir_saveddata;      % default
    end
end

if ~isfield(configs.files,'archive')
    if configs.savedata
        dir_archive=[configs.files.path,'_archive'];
        warning('configs.files.archive is unspecified. Setting to default %s.',dir_archive);
        configs.files.archive=dir_archive;      % default
    end
end

% minCounts
if ~isfield(configs.files,'minCount')
    warning('configs.files.minCount is unspecified. Setting to default 1000.');
    configs.files.minCount=1000;    % default
end
f_minCount=configs.files.minCount;  % min count in TXY to flag as errorneous


% Initialiase variables
% file output flags
files.build_txy=false(length(f_id),1);      % dld files processed to txy
files.missing=false(length(f_id),1);        % missing dld/txy file
files.lowcount=false(length(f_id),1);       % files with too few counts (skipped in analysis)
files.id_ok=false(length(f_id),1);          % files id's that were successfully processed

% XY rotation
if ~isfield(configs,'rot_angle')    % XY rotation needs to be back-compatible
    warning('rot_angle was not defined. Setting to default (0.61).');
    configs.rot_angle=0.61;     % set param to default value
end
rot_angle=configs.rot_angle;

%% Prepare TXY-files and check for low-counts (possibly errorneous shots)
if verbose>0, fprintf('Preparing TXY-files...\n'), end;
for i=1:length(f_id)
    % TXY-file does not exist
    if ~fileExists([f_path,'_txy_forc',num2str(f_id(i)),'.txt'])
        if verbose>1, warning('Could not find TXY-file #%d.',f_id(i)); end;
        
        % rawDLD source file exists
        if fileExists([f_path,num2str(f_id(i)),'.txt'])
            % Create TXY from DLD
            if verbose>0, warning('Creating TXY-file from raw source #%d.',f_id(i)); end;
            dld_raw_to_txy(f_path,f_id(i),f_id(i));
            files.build_txy(i)=1;
            
        % no source exists
        else
            % error - file # is missing
            warning('Could not load data. Source file #%d is missing.',f_id(i));
            files.missing(i)=1;
            continue;
        end
    end  % TXY-file exists
end

%% Extract a region of interest from TXY-files
% f_id=f_id(~files.missing);      % get ids for existing data files
txy_all=cell(length(f_id),1);   % TXY data cell in window

if verbose>0, fprintf('Getting counts in window from TXY-files...\n'); end;
counter=1;
for i=1:length(f_id)
    % pass for missing files
    if files.missing(i)
        continue
    end
    
    % load the TXY-file to memory
    txy_temp=txy_importer(f_path,f_id(i));
    
%     % Crop to pre-window
%     for dim=1:3
%         if isempty(CONFIGS.window{dim})
%             continue;    % empty window will pass cropping
%         end
%         in_window=((txy_temp(:,dim)>CONFIGS.window{dim}(1))&(txy_temp(:,dim)<CONFIGS.window{dim}(2)));
%         txy_temp=txy_temp(in_window,:);     % filter temporary variable
%     end
    
    % Crop in T to massively cull data
    if isempty(configs.window{1})   % pass crop if empty
        continue;
    end
    in_window=((txy_temp(:,1)>configs.window{1}(1))&...
        (txy_temp(:,1)<configs.window{1}(2)));
    txy_temp=txy_temp(in_window,:);     % pre-filter txy counts in T axis
    
    % Apply rotation in XY plane
    x_temp=txy_temp(:,2);   % temp store x,y vects
    y_temp=txy_temp(:,3);
    txy_temp(:,2)=x_temp*cos(rot_angle)-y_temp*sin(rot_angle);
    txy_temp(:,3)=x_temp*sin(rot_angle)+y_temp*cos(rot_angle);
    
    % Crop in XY plane
    for i_dim=2:3
        if isempty(configs.window{i_dim})   % pass crop if empty
            continue;
        end
        in_window=((txy_temp(:,i_dim)>configs.window{i_dim}(1))&(txy_temp(:,i_dim)<configs.window{i_dim}(2)));
        txy_temp=txy_temp(in_window,:);     % crop counts in XY axis
    end
    
    % Check for errorneous files by low-count in cropped region
    if size(txy_temp,1)<f_minCount
        files.lowcount(i)=1;
        if verbose>0
            warning('ROI Low-count in file #%d. Discarding from further processing.',f_id(i));
        end
        continue
    end
    
    txy_all{counter}=txy_temp;    % save captured counts
    counter=counter+1;
end
txy_all(counter:end)=[];    % delete all empty cells

% evaluate files successfully processed: id_ok
files.id_ok=f_id(~(files.missing|files.lowcount));

%% Plot captured counts (TXY)
if verbose>2
    h_zxy_all=figure();     % create figure
    plot_zxy(txy_all,1,'k');
    title('All counts');
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('T [s]');
    
    % save plot
    if configs.savedata
        fname_str='all_counts';
        saveas(h_zxy_all,[dir_output,fname_str,'.png']);
    end
end

%% Save processed data
if configs.savedata
    % if a file already exists it needs to be replaced
    if verbose>0,fprintf('Saving data...\n'); end;
    if exist(configs.files.saveddata,'file')
        warning('Data file already exists. Moving existing file to archive...');
        % create archive directory
        if ~exist(configs.files.archive,'dir')
            mkdir(configs.files.archive);
        end
        % move existing file to archive
        movefile(configs.files.saveddata,configs.files.archive);  % will overwrite an existing file of same name in the archive dir
    end
    
    % save data
    save(configs.files.saveddata,vars_save{1});    % must create *.mat without append option
    for i = 1:length(vars_save)
        if ~exist(vars_save{i},'var')
            warning(['Variable "',vars_save{i},'" does not exist.']);
            continue;
        end
        save(configs.files.saveddata,vars_save{i},'-v6','-append');     % -v6 version much faster (but no compression)?
    end
end

%% Summary
if verbose>0
    fprintf('===================IMPORT SUMMARY===================\n');
    fprintf('Number of shots successfully loaded: %d\n',length(files.id_ok));
    fprintf('Number of shots with counts below %d: %d\n',f_minCount,sum(files.lowcount));
    fprintf('Number of missing files: %d\n',sum(files.missing));
    fprintf('Number of files converted to TXY: %d\n',sum(files.build_txy));
    fprintf('====================================================\n');
end

%% END
t_fun_end=toc(t_fun_start);   % end of code
if verbose>0
    disp('-----------------------------------------------');
    fprintf('Total elapsed time for %s (s): %7.1f\n','loadExpData',t_fun_end);
    disp('-----------------------------------------------');
end