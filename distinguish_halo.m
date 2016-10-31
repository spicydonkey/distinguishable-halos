function [HALO_ZXY,BEC_ZXY] = distinguish_halo(CONFIGS,VERBOSE)
% Processes the raw TXY data to identify distinguishable halos and their
% properties
%
% INPUT:
% DATA_CONFIGS - struct; with fields 'files' and ...<TODO> that
% configures the analysis
% VERBOSE - default is 0
%
% OUTPUT:
% HALO_DATA - 2D cell of TXY-arrays from each halo
% OUTPUT - optional output report (TODO)
%
% DKS 31/10/16

% T-Z scaling
v_z = 9.8*0.430;    % TODO not sure of exact value to use

% Input check
if ~exist('VERBOSE','var')
    warning('VERBOSE is not provided - setting to quiet (0)');
    VERBOSE=0;
end

% Parse input
f_path = CONFIGS.files.path;
f_id = CONFIGS.files.idok;     % only analyse ok files

% Initialise variables
zxy_proc=cell(length(f_id),1);
HALO_ZXY=cell(length(f_id),2);
BEC_ZXY=cell(length(f_id),2);

bec_cent = cell(length(f_id),2);

%% halo processing
for i=1:length(f_id)
    zxy_shot=txy_importer(f_path,f_id(i));  % import full txy-data to memory
    zxy_shot(:,1)=zxy_shot(:,1)*v_z;    % ToF to Z - scale with z-velocity at detector
    
    %% Crop to pre-filter window
    % calculate z-window to do crop
    if isempty(CONFIGS.window.all{1})
        CONFIGS.window.all{1}=CONFIGS.window.all_T*v_z;
    end

    % apply windows
    for crop_dim=1:3
        if isempty(CONFIGS.window.all{crop_dim})
            continue;    % empty window will pass cropping
        end
        in_window=(zxy_shot(:,crop_dim)>CONFIGS.window.all{crop_dim}(1) & zxy_shot(:,crop_dim)<CONFIGS.window.all{crop_dim}(2));
        zxy_shot=zxy_shot(in_window,:);
    end
    
    %% Locate condensates
    zxy_bec=cell(2,1);  % BEC counts
    for i_cond=1:2
        zxy_temp=zxy_shot-repmat(CONFIGS.bec.pos{i_cond},[length(zxy_shot),1]); % relocate centre to approx BEC position
        zxy_temp=sum(zxy_temp.^2,2);    % radial distance vector
        ind_bec=zxy_temp<(CONFIGS.bec.Rmax{i_cond}^2);  % logical index vector for BEC
        zxy_bec{i_cond}=zxy_shot(ind_bec,:);      % BEC atoms
        zxy_shot=zxy_shot(~ind_bec,:);    % pop BEC out
    end
    
    %% save single-shot to a cell
    BEC_ZXY(i,:)=zxy_bec;  % BEC
    zxy_proc{i}=zxy_shot;   % remaining counts
end


if verbose>1
    
end


HALO_ZXY=zxy_proc;     % TODO temporary soln

% %% DEBUG
% figure(199);
% scatter3(halo_collate(:,2),halo_collate(:,3),halo_collate(:,1),'k.');
% xlabel('X'); ylabel('Y'); zlabel('Z');