function [halo,bec,culled,hfig,errflag]=halo_capture(txy,files,configs,verbose)
% Process TXY-counts to isolate all posssible halo counts and remove
% shot-to-shot oscillations
% DKS 15/05/2017
%

if ~exist('verbose','var')
    warning('verbose is not provided - setting to quiet (0)');
    verbose=0;
end

if verbose>0, fprintf('Beginning halo capture...\n'), end;
%% MAIN
t_fun_start=tic;
hfig={};

% Parse input
R_tail=configs.bec.dR_tail;     % estimated tail radius from BEC
R_halo=configs.halo.dR;         % fractional error around estimtaed halo rad for cropping (TODO - sensitivity?)

shot_id=files.id_ok;        % id's of files in txy (loaded ok from source)
nshot=size(txy,1);          % equivalent to length(shot_id)

% Initialise vars
halo.zxy=cell(nshot,2);    % halo counts
halo.zxy0=cell(nshot,2);   % halo counts (oscillations compensated + cap removed)
halo.cent=cell(nshot,2);   % halo centres
halo.R=cell(nshot,2);      % halo radii 
bec.zxy=cell(nshot,2);     % BEC counts
bec.cent=cell(nshot,2);    % BEC centres
culled.tail.zxy=cell(nshot,2);   % BEC tails
culled.fuzz.zxy=cell(nshot,1);   % halo fuzz + background counts - 1D collated since no source to distinguish

% Convert T-->Z
ZXY_all=txy2zxy(txy);
clear txy;

%% Locate condensates
% get BEC counts
bool_bec=cell(nshot,2);
for ii=1:nshot
    for jj=1:2
        [bec.cent{ii,jj},bool_bec{ii,jj}]=capture_bec(ZXY_all{ii},configs.bec.pos{jj},configs.bec.Rmax{jj});
    end
end

% store BEC counts
for ii=1:numel(configs.bec.pos)
    bec.zxy(:,ii)=cellfun(@(ZXY,BOOL)ZXY(BOOL,:),ZXY_all,bool_bec(:,ii),'UniformOutput',false);
end

% pop BEC
bool_bec_combined=cell_horzcat(bool_bec);       % horizontal concat the boolean list
bool_bec_combined=cellfun(@(x)sum(x,2)>0,bool_bec_combined,'UniformOutput',false);    % row-wise OR

ZXY_all=cellfun(@(ZXY,BOOL)ZXY(~BOOL,:),ZXY_all,bool_bec_combined,'UniformOutput',false);

% pop thermal tail - around condensate
for ii=1:numel(configs.bec.pos)
    % get radial dist from BEC centre
    r=cellfun(@(x,y)zxy2rdist(x,y),ZXY_all,bec.cent(:,ii),'UniformOutput',false);
    
    % get points close to BEC
    bool_tail=cellfun(@(x)x<configs.bec.Rmax{ii}*(1+R_tail{ii}),r,'UniformOutput',false);
    
    % store tail
    culled.tail.zxy(:,ii)=cellfun(@(ZXY,BOOL)ZXY(BOOL,:),ZXY_all,bool_tail,'UniformOutput',false);
    
    % cull tail
    ZXY_all=cellfun(@(ZXY,BOOL)ZXY(~BOOL,:),ZXY_all,bool_tail,'UniformOutput',false);
end

%% HALO (broad) capture
% assuming BEC centre lies on halo's Z-extremity and estimating halo
% radii, apply radial crop
errflag=false(nshot,1);  % flag for bad shots
zxy_halo=cell(1,2);     % temp for halo counts
zxy0_halo=cell(1,2);    % temp for centred halo counts
for i=1:nshot        % TODO: MAJOR BUG: treat f_id properly
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
            warning('Halo capture: shot #%d, halo #%d has no counts.',shot_id(i),i_halo);
            errflag(i)=1;
        end
        
        zxy_halo{i_halo}=ZXY_all{i}(ind_halo,:);       % counts in halo (abs ref)
        zxy0_halo{i_halo}=zxy_temp(ind_halo,:);          % approx halo-centered ref frame
        ZXY_all{i}=ZXY_all{i}(~ind_halo,:);             % pop halo out
    
        % Estimate halo radius: [AVG_RADIUS STD_RADIUS] for counts in halo
        % around the guessed centre
        halo.R{i,i_halo}=[mean(sqrt(rsq_temp(ind_halo))),std(sqrt(rsq_temp(ind_halo)))];
        
        % NaN error occurs for halo radius - Occurs or zero halo count
        if verbose>0
            if isnan(halo.R{i,i_halo}(1))
                warning('Halo capture: shot #%d, halo #%d has returned NaN for halo radius.',shot_id(i),i_halo);
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
            warning('Halo capture: shot #%d, halo #%d has no counts after removing cap.',shot_id(i),i_halo);
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

%% Plot captured halo
if configs.flags.graphics
    % All captured counts (halo only)
    h_haloraw=figure();     % create figure instance for each plot
    hfig{length(hfig)+1}=gcf;
    plot_zxy(halo.zxy,[],1);
    view(3);
    title('captured halos');
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    drawnow;
    
    % Halo (oscillation compensated, etc)
    h_halooc=figure();
    hfig{length(hfig)+1}=gcf;
    plot_zxy(halo.zxy0,[],1);
    view(3);
    axis equal;
    title('halos: oscillation compensated');
    xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
    drawnow;
end

%% Summary
if verbose>0
    bec_num=zeros(size(halo.zxy0));     % BEC nums in shot (same size as halos)
    halo_num=zeros(size(halo.zxy0));    % halo nums in shot
    
    bec_cent=cell(1,2);     % centres of BECs
    bec_cent_mean=cell(1,2);
    bec_cent_std=cell(1,2);
    
    halo_collated=cell(1,2);	% halo counts collated
    
    halo_radius=cell(1,2);      % approx radius of halos
    halo_radius_mean=zeros(1,2);
    halo_radius_std=zeros(1,2);
    for i=1:2
        % get number of counts captured for BEC/halo
        for j=1:size(halo.zxy0,1)
            bec_num(j,i)=size(bec.zxy{j,i},1);
            halo_num(j,i)=size(halo.zxy0{j,i},1);
        end
        
        % collate BEC centres
        bec_cent{i}=vertcat(bec.cent{:,i});
        bec_cent_mean{i}=mean(bec_cent{i},1);
        bec_cent_std{i}=std(bec_cent{i},1);
        
        halo_collated{i}=vertcat(halo.zxy0{:,i});
        % approx halo radius from all counts
        % TODO - this should have already been calculated
        halo_radius{i}=sqrt(sum(halo_collated{i}.^2,2));
        halo_radius_mean(i)=mean(halo_radius{i});
        halo_radius_std(i)=std(halo_radius{i});
    end
    bec_num_mean=mean(bec_num,1);
    halo_num_mean=mean(halo_num,1);
    bec_num_std=std(bec_num,1);
    halo_num_std=std(halo_num,1);
    
    fprintf('====================HALO CAPTURE====================\n');
    fprintf('Number of shots successful: %d\n',sum(~errflag));
    fprintf('Number of shots with zero counts in captured halo: %d\n',sum(errflag));
    fprintf('----------------------------------------------------\n');
    fprintf('BEC counts: [%5.1f,%5.1f]�[%5.1f,%5.1f]\n',bec_num_mean,bec_num_std);
    fprintf('Halo counts: [%3.1f,%3.1f]�[%3.1f,%3.1f]\n',halo_num_mean,halo_num_std);
    fprintf('----------------------------------------------------\n');
    for i=1:2
        fprintf('BEC(%d) centre: [%.5e,%.3e,%.3e]�[%.3e,%.3e,%.3e]\n',i,bec_cent_mean{i},bec_cent_std{i});
    end
    fprintf('----------------------------------------------------\n');
    fprintf('Halo radius: [%.3e,%.3e]�[%.3e,%.3e]\n',halo_radius_mean,halo_radius_std);
    fprintf('====================================================\n');
end


%% END
t_fun_end=toc(t_fun_start);   % end of code
if verbose>0
    disp('-----------------------------------------------');
    fprintf('Total elapsed time for %s (s): %7.1f\n','captureHalo',t_fun_end);
    disp('-----------------------------------------------');
end