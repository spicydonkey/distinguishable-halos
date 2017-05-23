% load and validate user config file and return well structured config
% variable for caller

function [configs,err]=load_config(config_file_path)
% initialise output
configs=[];
err=0;

% open config file
fi_config=fopen(config_file_path,'r');
if fi_config==-1
    warning('Config file %s could not be opened.',config_file_path);
    err=-1;
    return;
end

% run all lines in config file on MATLAB
while true
    this_line=fgetl(fi_config);
    if ~ischar(this_line)
        break;
    end     % enf of file
    eval(this_line);
end     % user defined "config" is built
fclose(fi_config);

% summary
if configs.flags.verbose>0
    fprintf('===================CONFIG SUMMARY===================\n');
    fprintf('config %s loaded successfully.\n',config_file_path);
    fprintf('====================================================\n');
end

end
