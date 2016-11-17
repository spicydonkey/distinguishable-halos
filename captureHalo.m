function [halo,bec,culled,errflag]=captureHalo(CONFIGS,VERBOSE)
% Process TXY-counts to isolate all posssible halo counts and remove
% shot-to-shot oscillations
% DKS 17/11/2016
%

if ~exist('VERBOSE','var')
    warning('VERBOSE is not provided - setting to quiet (0)');
    VERBOSE=0;
end

vars_save={'halo','bec','culled','errflag'};  % a list of variables to save to file

%% MAIN
t_fun_start=tic;
configs=CONFIGS;

% TEST
halo=[];
bec=[];
culled=[];
errflag=[];


%% Save processed data
% Append to existing data file with all counts
if VERBOSE>0,fprintf('Saving data...\n');,end;
for i = 1:length(vars_save)
    if ~exist(vars_save{i},'var')
        warning(['Variable "',vars_save{i},'" does not exist.']);
        continue;
    end
    save(configs.files.saveddata,vars_save{i},'-append');
end

%% END
t_fun_end=toc(t_fun_start);   % end of code
if VERBOSE>0
    disp('-----------------------------------------------');
    fprintf('Total elapsed time for %s (s): %7.1f\n','captureHalo',t_fun_end);
    disp('-----------------------------------------------');
end