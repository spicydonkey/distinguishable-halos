%% <template> CONFIG for loading experimental data

%-----------------------------------------------------------------
% PARAMETERS
%-----------------------------------------------------------------
% configs.files.path        % path to exp data up to file ID token
% configs.files.id          % file ID numbers to load
% configs.files.minCount    % min counts in file to pass as good
% 
% configs.window        % TXY-cropping window
% configs.window{1}         % T [s]
% configs.window{2}         % X [m]
% configs.window{3}         % Y [m]

configs.files.path='C:\Users\HE BEC\Documents\lab\halo_analysis\data\dist_halo\4_separated_lownum\d';
configs.files.id=1:3000;
configs.files.minCount=100;

configs.window{1}=[4.907,4.918];
configs.window{2}=[-20e-3,18e-3];
configs.window{3}=[-10e-3,17e-3];