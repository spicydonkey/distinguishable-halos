%% distinguishable halo analysis
% DKS 15/05/2017

function err=run_dist_halo(configm)
    %% Initialise
    t_main_start=tic;   % for reporting process duration
    datetimestr=datestr(datetime,'yyyymmdd_HHMMSS');    % timestamp when function called
    
    run(configm);   % run the configuration script
    
    % output directory
    if ~isdir(configs.files.dirout)
        warning(['output directory "',configs.files.dirout,'" does not exist. Creating directory...']);
        mkdir(configs.files.dirout);
    end
    
    %% Load TXY data
    if ~do_next
        % check for existing saved file for preloaded data
        mat_list=dir([configs.files.dirout,'/*.mat']);      % list of saved data
        if size(mat_list,1)==0
            do_next=1;  % no previously saved data files
        else
            % look for file with same load configs
            for ii=1:size(mat_list,1)
                this_file=mat_list(ii).name;
                % TODO - will throw error if $configs is not a variable in data
                S_temp=load(this_file,'configs');  % load configs from prev data
                
                % check if relevant configs equal
                if ~isequal(S_temp.configs.load,configs.load)
                    warning('Raw data: Existing data has different configs. Setting do_next=1.');
                    do_next=1;
                end
            end
        end
    end
    
    if do_next
        [txy,files_out]=load_txy(configs.load.path,configs.load.id,...
            configs.load.window,configs.load.minCount,...
            configs.load.rot_angle,verbose,configs.flags.graphics);
    end
end