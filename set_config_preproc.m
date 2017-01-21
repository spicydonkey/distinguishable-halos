%% Configuration setter for TDC preprocessing
% Generates config file required to:
% * Generate TXY from raw DLD files
% * TXY data processing to extract window of interest
%
% DK Shin
% 21/01/2017

proj_path='C:\Users\David\Documents\hebec\halo_analysis\data\dist_halo\4';

%% PROJECT DIRECTORIES
data.path=strcat(proj_path,'\raw_data\d');  % raw data dir
% config directory
config.path=strcat(proj_path,'\config');
if ~isdir(config.path)
    warning(strcat('config directory does not exist. Creating directory: ',config.path));
    mkdir(config.path);
end

%% Raw data handling
data.id=1:100;     % file id numbers to process
data.minCount=100;  % min counts in file to pass

%% TXY window (set to [] for no cropping)
data.window{1}=[4.907,4.918];    % T [s]
data.window{2}=[-20e-3,18e-3];   % X [m]
data.window{3}=[-10e-3,17e-3];   % Y [m]

%% save data
vars_save={'proj_path','data','config'};

% generate random config ID
ran_num=ceil(rand()*1000);
config.id=strcat(date,'_',num2str(ran_num));

config_fname=strcat(config.path,'\config_txy_',config.id,'.mat');
disp(config_fname);

save(config_fname,vars_save{1},'-v6');
for i=1:length(vars_save)
    save(config_fname,vars_save{i},'-v6','-append');
end