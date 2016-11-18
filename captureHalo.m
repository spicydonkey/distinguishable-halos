function [halo,bec,culled,errflag]=captureHalo(CONFIGS,VERBOSE)
% Process TXY-counts to isolate all posssible halo counts and remove
% shot-to-shot oscillations
% DKS 17/11/2016
%

if ~exist('VERBOSE','var')
    warning('VERBOSE is not provided - setting to quiet (0)');
    VERBOSE=0;
end

vars_save={'configs','halo','bec','culled','errflag'};  % a list of variables to save to file
% configs needs to be overwritten to update new bec/halo capture settings

if VERBOSE>0, fprintf('Beginning halo capture...\n'), end;
%% MAIN
t_fun_start=tic;
configs=CONFIGS;

% Parse input
v_z=configs.misc.vel_z;	% z-velocity for T-Z conversion

R_tail=configs.bec.dR_tail;     % estimated tail radius from BEC
R_halo=configs.halo.dR;         % fractional error around estimtaed halo rad for cropping (TODO - sensitivity?)

% load saved file with windowed counts
files=struct([]);   % initialise 'files'
load(configs.files.saveddata,'files');  % load files summary
f_idok=files.id_ok;     % id's of files in txy_all (loaded ok from source)

txy_all=cell(length(f_idok),1); % preallocate txy_all before loading
load(configs.files.saveddata,'txy_all');    % load all counts

% Initialise vars
halo.zxy=cell(length(f_idok),2);    % halo counts
halo.zxy0=cell(length(f_idok),2);   % halo counts (oscillations compensated)
halo.cent=cell(length(f_idok),2);   % halo centres
halo.R=cell(length(f_idok),2);      % halo radii 
bec.zxy=cell(length(f_idok),2);     % BEC counts
bec.cent=cell(length(f_idok),2);    % BEC centres
culled.tail.zxy=cell(length(f_idok),2);   % BEC tails
culled.fuzz.zxy=cell(length(f_idok),1);   % halo fuzz + background counts - 1D collated since no source to distinguish

% Convert T-->Z
for i=1:size(txy_all,1)
    txy_all{i}(:,1)=txy_all{i}(:,1)*v_z;    % ToF to Z
end
ZXY_all=txy_all;    % much easier to handle counts in ZXY
clear txy_all;

%% Process TXY data
% BEC
for i=1:length(f_idok)
    %fprintf('%d: %d\n',i,size(ZXY_all{i},1));
    
    %% Locate condensates
    % Locate condensate by cropping a ball around estimated centres and
    %   radius and AVERAGE count positions - ITERATED until convergence
    zxy_bec=cell(1,2);  % BEC counts
    
    % ITERATE to improve BEC estimation
    % Initialisise
    bec.cent(i,:)=configs.bec.pos;
    ball_cent=configs.bec.pos;      % ball (centre) to capture BEC
    ball_rad=configs.bec.Rmax;      % radial window to count as BEC
    for i_cond=1:2
        % iterate until mean position from counts ~= ball centre
        err_cent=Inf;
        n_bec_max=0;
        n_iter=0;
        while err_cent>0.02e-3
            % get new ball centre for BEC capture
            ball_cent{i_cond}=bec.cent{i,i_cond};
            
            zxy_temp=ZXY_all{i}-repmat(ball_cent{i_cond},[size(ZXY_all{i},1),1]); % centre to ball
            rsq_temp=sum(zxy_temp.^2,2);            % evaluate radial distances
            %TODO - this could be much tighter since single shot shows BEC rad < 4 mm
            ind_bec=rsq_temp<((0.7*ball_rad{i_cond})^2);  % logical index vector for BEC - atoms within Rmax
            zxy_bec{i_cond}=ZXY_all{i}(ind_bec,:);    % collate captured BEC counts
            
            % Evaluate BEC centre
            bec.cent{i,i_cond}=mean(zxy_bec{i_cond},1);     % approx of BEC centre by mean position
            err_cent=norm(ball_cent{i_cond}-bec.cent{i,i_cond});    % error this iteration
            
            % total count in this ball - used to tell if BEC well-captured
            n_bec_this=sum(ind_bec);
            if n_bec_this>n_bec_max
                n_bec_max=n_bec_this;
            end
            
            n_iter=n_iter+1;
        end
        % centre has converged
        
        ZXY_all{i}=ZXY_all{i}(~ind_bec,:);  % pop BEC out of all counts
        
        % Summary
        if n_bec_this/n_bec_max<0.8
            warning('While locating BEC, sudden drop in count observed in shot #%d.',f_idok(i));
        end
        if VERBOSE>2
            disp('-------------------------------------------------');
            fprintf('Shot: %d, BEC#: %d\n',f_idok(i),i_cond);
            disp(['Iterations: ',num2str(n_iter)]);
            disp(['Total counts in BEC (/max): ',num2str(n_bec_this),' / ',num2str(n_bec_max)]);
            % total deviation from initial guess
            dev_tot=norm(bec.cent{i,i_cond}-configs.bec.pos{i_cond});
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
        zxy_temp=ZXY_all{i}-repmat(bec.cent{i,i_cond},[size(ZXY_all{i},1),1]); % relocate centre to approx BEC position
        rsq_temp=sum(zxy_temp.^2,2);            % evaluate radial distances
        ind_tail=rsq_temp<(((1+R_tail{i_cond})*configs.bec.Rmax{i_cond})^2);  % tail counts index
        zxy_tail{i_cond}=ZXY_all{i}(ind_tail,:);  % counts in BEC tail
        ZXY_all{i}=ZXY_all{i}(~ind_tail,:);         % pop tail out
    end
    
    % save single-shot to a cell
    bec.zxy(i,:)=zxy_bec;	% BEC
	culled.tail.zxy(i,:)=zxy_tail;  % tail
end

%% HALO (broad) capture
% assuming BEC centre lies on halo's Z-extremity and estimating halo
% radii, apply radial crop
errflag=false(size(f_idok));  % flag for bad shots
zxy_halo=cell(1,2);     % temp for halo counts
zxy0_halo=cell(1,2);    % temp for centred halo counts
for i=1:length(f_idok)        % TODO: MAJOR BUG: treat f_id properly
    for i_halo=1:2
        % add/subtract estimated halo radius in Z: HALO 1 should be the
        %   non-magnetic Raman outcoupled atoms - Halo must sit "ABOVE" the
        %   BEC
        halo.cent{i,i_halo}=bec.cent{i,i_halo}+((-1)^(i_halo+1))*[1,0,0]*configs.halo.R{i_halo};
        zxy_temp=ZXY_all{i}-repmat(halo.cent{i,i_halo},[size(ZXY_all{i},1),1]);     % ref about halo G
        rsq_temp=sum(zxy_temp.^2,2);            % evaluate radial distances
        ind_halo=((rsq_temp<(((1+R_halo{i_halo})*configs.halo.R{i_halo})^2)) & (rsq_temp>(((1-R_halo{i_halo})*configs.halo.R{i_halo})^2)));  % halo counts index
        
        % Check for no counts captured for halo
        if sum(ind_halo)==0
            warning('Halo capture: shot #%d, halo #%d has no counts.',f_idok(i),i_halo);
            errflag(i)=1;
        end
        
        zxy_halo{i_halo}=ZXY_all{i}(ind_halo,:);       % counts in halo (abs ref)
        zxy0_halo{i_halo}=zxy_temp(ind_halo,:);          % approx halo-centered ref frame
        ZXY_all{i}=ZXY_all{i}(~ind_halo,:);             % pop halo out
    
        % Estimate halo radius: [AVG_RADIUS STD_RADIUS] for counts in halo
        % around the guessed centre
        halo.R{i,i_halo}=[mean(sqrt(rsq_temp(ind_halo))),std(sqrt(rsq_temp(ind_halo)))];
        
        % NaN error occurs for halo radius - Occurs or zero halo count
        if VERBOSE>0
            if isnan(halo.R{i,i_halo}(1))
                warning('Halo capture: shot #%d, halo #%d has returned NaN for halo radius.',f_idok(i),i_halo);
            end
        end
    end
    
    % save single-shot to a cell
    halo.zxy(i,:)=zxy_halo;             % HALO (abs ref frame)
    halo.zxy0(i,:)=zxy0_halo;           % HALO - oscillations compensated
    culled.fuzz.zxy{i}=ZXY_all{i};      % all remaining counts called fuzz
end


%% Remove caps
haloR=configs.halo.R;
zcap=configs.halo.zcap;
for i=1:size(halo.zxy0,1)
    for j=1:2
        % remove caps
        ind_cap_temp=abs(halo.zxy0{i,j}(:,1))>(zcap*haloR{j});
        halo.zxy0{i,j}=halo.zxy0{i,j}(~ind_cap_temp,:);

        % check for zero counts
        if sum(~ind_cap_temp)==0
            warning('Halo capture: shot #%d, halo #%d has no counts after removing cap.',f_idok(i),i_halo);
            errflag(i)=1;
        end
    end
end


%% Handle bad shots
% bad shots will be discarded
nBadshots=sum(errflag);

if nBadshots>0
    warning('%d bad shot(s) in data will be discarded.',nBadshots);
    
    % Cull bad shots for all returned processed data
    halo.zxy=halo.zxy(~errflag,:);
    halo.zxy0=halo.zxy0(~errflag,:);
    halo.R=halo.R(~errflag,:);
    halo.cent=halo.cent(~errflag,:);
    
    bec.zxy=bec.zxy(~errflag,:);
    bec.cent=bec.cent(~errflag,:);
    
    culled.fuzz.zxy=culled.fuzz.zxy(~errflag,:);
    culled.tail.zxy=culled.tail.zxy(~errflag,:);
end


%% Save processed data
% Append to existing data file with all counts
if VERBOSE>0,fprintf('Saving data...\n'); end;
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