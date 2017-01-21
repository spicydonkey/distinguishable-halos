function [zxy,files,fout]=load_zxy(config_path,verbose)
% 

if ~exist('verbose','var')
    verbose=0;  % default verbose is quiet
end

if ~exist(config_path,'file')
    error('Given configuration file does not exist.');
    return
end

%% MAIN
t_fun_start=tic;

configs=load(config_path);

f_id=configs.data.id;       % get files id's to process
f_path=configs.data.path;   % get path to file (no id token)
f_minCount=configs.data.minCount;   % min count in TXY to flag as bad file

% Initialise variables
% file output flags
files.build_txy=false(length(f_id),1);      % dld files processed to txy
files.missing=false(length(f_id),1);        % missing dld/txy file
files.lowcount=false(length(f_id),1);       % files with too few counts (skipped in analysis)
files.id_ok=false(length(f_id),1);          % files id's that were successfully processed

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
    
    % Crop to pre-window
    for dim=1:3
        if isempty(configs.data.window{dim})
            continue;    % empty window will pass cropping
        end
        in_window=((txy_temp(:,dim)>configs.data.window{dim}(1))&(txy_temp(:,dim)<configs.data.window{dim}(2)));
        txy_temp=txy_temp(in_window,:);     % filter temporary variable
    end
    
    % Check for errorneous files by low-count
    if size(txy_temp,1)<f_minCount
        files.lowcount(i)=1;
        if verbose>0
            warning('Low-count in file #%d. Discarding from further processing.',f_id(i));
        end
        continue
    end
    
    txy_all{counter}=txy_temp;    % save captured counts
    counter=counter+1;
end
txy_all(counter:end)=[];    % delete all empty cells

% evaluate files successfully processed: id_ok
files.id_ok=f_id(~(files.missing|files.lowcount));

%% T-to-Z conversion
zxy=txy_all;    % create a copy
vel_z=9.8*0.430;    % atom free-fall vert v at detector hit for T-to-Z conversion;
for i=1:length(txy_all)
    zxy{i}(:,1)=zxy{i}(:,1)*vel_z;
end

%% Save processed data
if verbose>0,fprintf('Saving data...\n'); end;
% create processed data directory
dir_output=strcat(configs.proj_path,'\proc_data');
if ~exist(dir_output,'dir')
    mkdir(dir_output);
end
fout=strcat(dir_output,'\zxy_',configs.config.id);

% save data
vars_save={'zxy','files'};  % a list of variables to save to file
save(fout,vars_save{1},'-v6');    % must create *.mat without append option
for i = 1:length(vars_save)
    if ~exist(vars_save{i},'var')
        warning(['Variable "',vars_save{i},'" does not exist.']);
        continue;
    end
    save(fout,vars_save{i},'-v6','-append');     % -v6 version much faster (but no compression)?
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