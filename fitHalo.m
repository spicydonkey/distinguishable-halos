function [halo_fit]=fitHalo(CONFIGS,VERBOSE)
% FIT the halos
%

% Input check
if ~exist('VERBOSE','var')
    VERBOSE=0;  % default verbose is quiet
end

vars_save={'halo_fit'};  % a list of variables to save to file

if VERBOSE>0, fprintf('Fitting halos...\n'), end;
%% MAIN
t_fun_start=tic;

configs=CONFIGS;

% load captured halo counts
S_temp=load(configs.files.saveddata,'halo');
halo_zxy0=S_temp.halo.zxy0;     % get the oscillation compensated halo counts
clear S_temp;


% Do fitting here
%   get halo_fit: .zxy0, .param, etc.

halo_fit.zxy0=halo_zxy0;     % DUMMY - no fitting


%% Save processed data
% Append to existing data file
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
    fprintf('Total elapsed time for %s (s): %7.1f\n','fitHalo',t_fun_end);
    disp('-----------------------------------------------');
end