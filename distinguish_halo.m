function [HALO,BEC,CULLED] = distinguish_halo(CONFIGS,VERBOSE)
% Processes the raw TXY data to identify distinguishable halos and their
% properties
% Data processing includes:
%   - Crop to pre-processing window (box)
%   - Locate condensates per shot
%   - Clean shot around condensate (thermal fraction)
%   - Capture halos (TODO - improve with better radius estimation)
%   - Shot-to-shot oscillation cancellation (post single-count classification)
%
% INPUT:
% CONFIGS - struct configures the analysis: refer to dist_halo_main
% VERBOSE - default is 0 (TODO - need some print outs?)
%
% OUTPUT:
% HALO,BEC,CULLED are structs with fields zxy, cent, etc. where appropriate 
%   as processed from CONFIGS
%
% DKS 31/10/16

% Input check
if ~exist('VERBOSE','var')
    warning('VERBOSE is not provided - setting to quiet (0)');
    VERBOSE=0;
end

%% Parse input
f_path = CONFIGS.files.path;    % full path + filename token
f_id = CONFIGS.files.idok;      % only analyse ok files

v_z = CONFIGS.misc.vel_z;       % z-velocity for T-Z conversion

R_tail=CONFIGS.bec.dR_tail;     % estimated tail radius from BEC
R_halo=CONFIGS.halo.dR;         % fractional error around estimtaed halo rad for cropping (TODO - sensitivity?)

%% Initialise variables
HALO.zxy=cell(length(f_id),2);  % halos
HALO.cent=cell(length(f_id),2); % centre of halos
BEC.zxy=cell(length(f_id),2);   % BECs
BEC.cent=cell(length(f_id),2);  % centre of BECs
CULLED.tail.zxy=cell(length(f_id),2);   % BEC tails
CULLED.fuzz.zxy=cell(length(f_id),1);   % halo fuzz + background counts - 1D collated since no source to distinguish

%% Halo processing
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
    % Locate condensate by cropping a ball around estimated centres and
    %   radius and AVERAGE count positions - ITERATED until convergence
    zxy_bec=cell(1,2);  % BEC counts
    
    % ITERATE to improve BEC estimation
    % Initialisise
    BEC.cent(i,:)=CONFIGS.bec.pos;
    ball_cent=CONFIGS.bec.pos;      % ball (centre) to capture BEC
    ball_rad=CONFIGS.bec.Rmax;      % radial window to count as BEC
    for i_cond=1:2        
        % iterate until mean position from counts ~= ball centre
        err_cent=Inf;
        n_bec_max=0;
        n_iter=0;
        while err_cent>0.02e-3     
            % get new ball centre for BEC capture
            ball_cent{i_cond}=BEC.cent{i,i_cond};
            
            zxy_temp=zxy_shot-repmat(ball_cent{i_cond},[length(zxy_shot),1]); % centre to ball
            rsq_temp=sum(zxy_temp.^2,2);            % evaluate radial distances
            %TODO - this could be much tighter since single shot shows BEC rad < 4 mm
            ind_bec=rsq_temp<((0.7*ball_rad{i_cond})^2);  % logical index vector for BEC - atoms within Rmax
            zxy_bec{i_cond}=zxy_shot(ind_bec,:);    % collate captured BEC counts
            
            % Evaluate BEC centre
            BEC.cent{i,i_cond}=mean(zxy_bec{i_cond},1);     % approx of BEC centre by mean position
            err_cent=norm(ball_cent{i_cond}-BEC.cent{i,i_cond});    % error this iteration
            
            % total count in this ball - used to tell if BEC well-captured
            n_bec_this=sum(ind_bec);
            if n_bec_this>n_bec_max
                n_bec_max=n_bec_this;
            end

            n_iter=n_iter+1;
        end
        % centre has converged
        
        zxy_shot=zxy_shot(~ind_bec,:);  % pop BEC out
        
        % Summary
        if VERBOSE>2
            disp('-------------------------------------------------');
            disp(['Shot: ',num2str(i),', BEC#: ',num2str(i_cond)]);
            disp(['Iterations: ',num2str(n_iter)]);
            disp(['Total counts in BEC (/max): ',num2str(n_bec_this),' / ',num2str(n_bec_max)]);
            % total deviation from initial guess
            dev_tot=norm(BEC.cent{i,i_cond}-CONFIGS.bec.pos{i_cond});
            disp(['Deviation from initial guess: ',num2str(1e3*dev_tot),' mm']);
        end
    end
    
    
    %% Clean around condensates
    % There is still spherical tail around the condensate capturing sphere which is
    %   much stronger than scattered counts
    % Must be done after locating condensates since the clean-up windows
    %   for two condensates will significantly overlap
    zxy_tail=cell(1,2);
    
    for i_cond=1:2
        zxy_temp=zxy_shot-repmat(CONFIGS.bec.pos{i_cond},[length(zxy_shot),1]); % relocate centre to approx BEC position
        zxy_temp=sum(zxy_temp.^2,2);            % evaluate radial distances
        ind_tail=zxy_temp<(((1+R_tail{i_cond})*CONFIGS.bec.Rmax{i_cond})^2);  % tail counts index
        zxy_tail{i_cond}=zxy_shot(ind_tail,:);  % counts in BEC tail
        zxy_shot=zxy_shot(~ind_tail,:);         % pop tail out
    end
    
    %% Capture halos
    % assuming BEC centre lies on halo's Z-extremity and estimating halo
    % radii, apply radial crop
    zxy_halo=cell(1,2); % temp for halo counts
    
    for i_halo=1:2
        % add/subtract estimated halo radius in Z: HALO 1 should be the
        %   non-magnetic Raman outcoupled atoms - Halo must sit "ABOVE" the
        %   BEC
        HALO.cent{i,i_halo}=BEC.cent{i,i_halo}+((-1)^(i_halo+1))*[1,0,0]*CONFIGS.halo.R{i_halo};
        zxy_temp=zxy_shot-repmat(HALO.cent{i,i_halo},[length(zxy_shot),1]);     % ref about halo G
        zxy_temp=sum(zxy_temp.^2,2);            % evaluate radial distances
        ind_halo=((zxy_temp<(((1+R_halo{i_halo})*CONFIGS.halo.R{i_halo})^2)) & (zxy_temp>(((1-R_halo{i_halo})*CONFIGS.halo.R{i_halo})^2)));  % halo counts index
        zxy_halo{i_halo}=zxy_shot(ind_halo,:);  % counts in halo
        zxy_shot=zxy_shot(~ind_halo,:);         % pop halo out
    end
    
    %% save single-shot to a cell
    BEC.zxy(i,:)=zxy_bec;   % BEC
    HALO.zxy(i,:)=zxy_halo;   % TODO: this is a temporary solution - dump all remaining counts to halo
    CULLED.tail.zxy(i,:)=zxy_tail;  % tail
    CULLED.fuzz.zxy{i}=zxy_shot;    % all remaining counts
end
clear zxy_shot in_window zxy_bec zxy_temp ind_bec zxy_tail ind_tail zxy_halo ind_halo;
clear i crop_dim i_cond i_halo;


%% Shot-to-shot oscillation cancellation
% post-processing: centre BEC-halo pairs to the centre of halo as
%   determined above
if CONFIGS.proc.cancel_oscillation
    for i=1:2
        BEC.zxy(:,i)=zxy_translate(BEC.zxy(:,i),HALO.cent(:,i));
        HALO.zxy(:,i)=zxy_translate(HALO.zxy(:,i),HALO.cent(:,i));
    end
end