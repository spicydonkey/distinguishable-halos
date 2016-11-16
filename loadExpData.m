function [txy_all,files]=loadExpData(CONFIGS,VERBOSE)
% Loads raw-data DLD/TXY and get counts in region of interest and save to
% file
%
% INPUT
% CONFIGS
%   has fields:     files.id,.path,.minCount 
%                   window - 3x1 cell of Z,X,Y window limits (crop)
%
% OUTPUT
% files: flags for processed files
%   has fields:     build_txy,missing,lowcount,id_ok
%
% Saves the data to a file
%   txy_all, files
%

if ~exist('VERBOSE','var')
    VERBOSE=0;  % default verbose is quiet
end

vars_save={'txy_all','files'};  % a list of variables to save to file

%% MAIN
t_fun_start=tic;

f_id=CONFIGS.files.id;      % get files id's to process
f_path=CONFIGS.files.path;  % get path to file (no id token)
f_minCount=CONFIGS.files.minCount;  % min count in TXY to flag as errorneous

% Initialiase variables
% file output flags
files.build_txy=false(length(f_id),1);      % dld files processed to txy
files.missing=false(length(f_id),1);        % missing dld/txy file
files.lowcount=false(length(f_id),1);       % files with too few counts (skipped in analysis)
files.id_ok=false(length(f_id),1);          % files id's that were successfully processed

%% Prepare TXY-files and check for low-counts (possibly errorneous shots)
if VERBOSE>0, fprintf('Preparing TXY-files...\n'), end;
for i=1:length(f_id)
    % TXY-file does not exist
    if ~fileExists([f_path,'_txy_forc',num2str(f_id(i)),'.txt'])
        if VERBOSE>1, warning('Could not find TXY-file #%d.',f_id(i));, end;
        
        % rawDLD source file exists
        if fileExists([f_path,num2str(f_id(i)),'.txt'])
            % Create TXY from DLD
            if VERBOSE>0, warning('Creating TXY-file from raw source #%d.',f_id(i));, end;
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

if VERBOSE>0, fprintf('Getting counts in window from TXY-files...\n');, end;
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
        if isempty(CONFIGS.window{dim})
            continue;    % empty window will pass cropping
        end
        in_window=((txy_temp(:,dim)>CONFIGS.window{dim}(1))&(txy_temp(:,dim)<CONFIGS.window{dim}(2)));
        txy_temp=txy_temp(in_window,:);     % filter temporary variable
    end
    
    % Check for errorneous files by low-count
    if size(txy_temp,1)<f_minCount
        files.lowcount(i)=1;
        if VERBOSE>0
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

%% Save processed data
% if a file already exists it needs to be replaced
if VERBOSE>0,fprintf('Saving data...\n');,end;
if exist([f_path,'data.mat'],'file')
    warning('Data file already exists. Moving existing file to archive...');
    % create archive directory
    if ~exist([f_path,'_archive'],'dir')
       mkdir([f_path,'_archive']);
    end
    % move existing file to archive
    movefile([f_path,'data.mat'],[f_path,'_archive']);  % will overwrite an existing file of same name in the archive dir
end

% save data
save([f_path,'data.mat'],vars_save{1});    % must create *.mat without append option
for i = 1:length(vars_save)
    if ~exist(vars_save{i},'var')
        warning(['Variable "',vars_save{i},'" does not exist.']);
        continue;
    end
    save([f_path,'data.mat'],vars_save{i},'-append');
end

%% END
t_fun_end=toc(t_fun_start);   % end of code
if VERBOSE>0
    disp('-----------------------------------------------');
    fprintf('Total elapsed time for %s (s): %7.1f\n','loadExpData',t_fun_end);
    disp('-----------------------------------------------');
end