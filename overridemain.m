% OVERRIDES dist_halo_main.m
clear all; clc; close all;

% PARAMS TO OVERRIDE
OR_configs.halo.R{1}=linspace(9e-3,13e-3,10);
OR_configs.halo.R{2}=linspace(8e-3,11e-3,10);

% create archive directory
OR_dir_archive='C:\Users\HE BEC\Documents\lab\halo_analysis\data\dist_halo\4_separated_lownum\d_archive';
if ~isdir(OR_dir_archive)
    mkdir(OR_dir_archive);
end

% MAIN
for OR_i=1:100
    disp(OR_i); %DEBUG
    
    % Set this override command
    override='testing...';
    
    % call main script with override
    dist_halo_main;
    
    % move all saved files to a new directory
    OR_dir_new=[OR_dir_archive,'\',num2str(OR_i)];
    mkdir(OR_dir_new);
    movefile([dir_output,'*'],OR_dir_new);
    save([OR_dir_new,'\data.mat'],'halo','configs','analysis');
end