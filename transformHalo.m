function [zxy]=transformHalo(CONFIGS,VERBOSE)
% Transform halos
%
% [] iso/aniso scaling
%   [x] normalise
% [] Z-axis spin : but shouldn't have spun much at all
%
% FINAL counts are called: 'zxy'

% input check
if ~exist('VERBOSE','var')
    warning('VERBOSE is not provided - setting to quiet (0)');
    VERBOSE=0;
end

vars_save={'zxy'};

if VERBOSE>0, fprintf('Transforming halos...\n'), end;
%% MAIN
t_fun_start=tic;
configs=CONFIGS;

% load fitted halo counts
S_temp=load(configs.files.saveddata,'halo_fit');
halo_zxy0=S_temp.halo_fit.zxy0;  % get fitted halo counts (real space)
clear S_temp;

% Initialise vars
zxy=cell(size(halo_zxy0));

% Do transforms here
%% k-space transform
% indepdent isotropic scaling to unit radius (normalisation)

% Get mean halo radius for all shots combined
R_halo_temp=zeros(2,1);
for i=1:2
    halo_temp=vertcat(halo_zxy0{:,i});
    R_halo_temp(i)=mean(sqrt(sum(halo_temp.^2,2)));
end
clear halo_temp;

% Isometric scaling (normalisation)
for i=1:size(halo_zxy0,1)
    for j=1:2
        zxy{i,j}=halo_zxy0{i,j}/R_halo_temp(j);
    end
end
clear R_halo_temp;


%% Save processed data
% Append to existing data file
if VERBOSE>0,fprintf('Saving data...\n'); end;
for i = 1:length(vars_save)
    if ~exist(vars_save{i},'var')
        warning(['Variable "',vars_save{i},'" does not exist.']);
        continue;
    end
    save(configs.files.saveddata,vars_save{i},'-v6','-append');
end

%% END
t_fun_end=toc(t_fun_start);   % end of code
if VERBOSE>0
    disp('-----------------------------------------------');
    fprintf('Total elapsed time for %s (s): %7.1f\n','transformHalo',t_fun_end);
    disp('-----------------------------------------------');
end