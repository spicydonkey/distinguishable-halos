%% distinguishable halo analysis
% DKS 15/05/2017

%% TODO
% Characterise halos
    %   mode occupancy - correlation volume (from CL) from BB (Sean's
    %   paper)
% halo transform scan for g2
    %   center --> g2 peaks at 0 for each x,y,z axis; could do 3D Gaussian
    %   fit
    %   rotation
    %   scaling
    
function [halo_k0,corr,efit,halo,txy,fout,err]=run_dist_halo(config_file)
    %% Initialise
    % configure
    [configs,err_cf]=load_config(config_file);
    if err_cf~=0
        error('config file could not be loaded.');
    end
    
    do_corr_analysis=configs.flags.do_corr_analysis;
    do_next=configs.flags.force_all_stages;
    verbose=configs.flags.verbose;
    
    t_main_start=tic;   % for reporting process duration
    datetimestr=datestr(datetime,'yyyymmdd_HHMMSS');    % timestamp when function called
    configs.files.dirout=[configs.files.dirout,'_',datetimestr];
    HFIG={};
    
    if configs.flags.savedata
        % output directory
        if ~isdir(configs.files.dirout)
            warning(['output directory "',configs.files.dirout,'" does not exist. Creating directory...']);
            mkdir(configs.files.dirout);
        end
    end
    
    % check for archive directory
    if ~isdir(configs.files.archive)
        warning(['archive directory "',configs.files.archive,'"does not exist. Creating directory...']);
        mkdir(configs.files.archive);
    end
    
    
    %% Load TXY data
    % TODO - what if there's multiple TXY files with same load config? e.g.
    % from updating load TXY --> versioning!
    loadconfigs=configs.load;   % abbreviated configs storing only txy load part
    
    if ~do_next
        % check for existing saved file for preloaded data
        arch_mat_list=dir([configs.files.archive,'/txy_*.mat']);      % list of saved data
        
        do_next=1;      % to do all task until saved data found
        if size(arch_mat_list,1)~=0
            % look for file with same load configs
            for ii=1:size(arch_mat_list,1)
                this_file=arch_mat_list(ii).name;
                S_temp=load([configs.files.archive,'/',this_file],'loadconfigs');   % load configs from prev data
                
                if isequal(S_temp.loadconfigs,loadconfigs)
                    % success! this is the data to be loaded
                    warning('Loading TXY: Match found in archive. Loading %s.',[configs.files.archive,'/',this_file]);
                    load([configs.files.archive,'/',this_file],'txy','fout');
                    do_next=0;      % skip full stage
                    break
                end
            end
            
            % no saved txy found
            if do_next
                warning('Loading TXY: No relevant TXY data found in archive. Setting do_next=1.');
            end
            
        end
    end
    clear S_temp this_file arch_mat_list;       % clean workspace
    
    if do_next
        [txy,fout,HFIG{length(HFIG)+1}]=load_txy(configs.files.path,configs.load.id,...
            configs.load.window,configs.load.mincount,configs.load.maxcount,...
            configs.load.rot_angle,configs.flags.build_txy,verbose,configs.flags.graphics);
        
        % save loaded TXY data to archive for fast loading
        if configs.flags.archive_txy&&(~configs.flags.force_all_stages)     % don't save when debugging - creates multiplicity
            % save the loaded txy, output log, and configs used
            fpath_txy_archive=[configs.files.archive,'/txy_',datetimestr,'.mat'];
            warning('Archiving TXY data as .mat file: %s',fpath_txy_archive);
            save(fpath_txy_archive,'txy','fout','loadconfigs');
        end
    end
    
    clear loadconfigs fpath_txy_archive;  % clean workspace
    
    
    %% Density profiles
    % TODO
    %   [] run only when debug flag ON: this is an optional analysis for
    %       checking unprocessed TXY
    %   [*] inspect halo + background density profile: background count density
    %   [] inspect BEC/thermal density profile
    
    % 1. Collision COM frame
    % from TXY crop out a disk perpendicular to the collision axis
    ddiskheight=0.1;   % thickness of disk in halo radius unit
    
    % calculate halo centre in TXY
    cent_halo=configs.bec.pos;      % get rough estimate of BEC position (in ZXY)
    for ii=1:2
        cent_halo{ii}=cent_halo{ii}+(-1)^(ii-1)*configs.halo.R{ii}*[1,0,0];     % estimate of halo centres (TXY)
        cent_halo{ii}(1)=cent_halo{ii}(1)/configs.misc.vel_z;       % Z-->T
    end
    
    % get counts in thin disk
    txy_disk=cell(size(txy,1),2);
    for ii=1:2
        R=configs.halo.R{ii};   % halo capture radius [m]
        txy_disk(:,ii)=cellfun(@(x) cylindercull(x,cent_halo{ii},[2*R,ddiskheight*R/configs.misc.vel_z],1),txy,'UniformOutput',false);  % capture counts in thin disk
    end
    
    % merge shots
    txy_disk_collated=cell(1,2);
    for ii=1:2
        txy_disk_collated{ii}=vertcat(txy_disk{:,ii});
    end
    
    % centre halo
    for ii=1:2
        txy_disk_collated{ii}=txy_disk_collated{ii}-repmat(cent_halo{ii},[size(txy_disk_collated{ii},1),1]);
    end
    
    % plot captured counts in disk
    hfig_halo_disk=figure();
    plot_zxy(txy_disk_collated);
    box on;
    view(2);    % XY plane
    xlabel('X [m]');
    ylabel('Y [m]');
    titlestr='Point cloud through centre of halo';
    title(titlestr);
    legend({'$m_F=0$','$m_F=1$'});
    
    % get radii
    r_disk_com=cell(1,2);
    for ii=1:2
        r_disk_com{ii}=sqrt(sum(txy_disk_collated{ii}(:,2:3).^2,2))/configs.halo.R{ii};     % radii in halo radius unit
    end
    
    % histogram radially
    [n_r_halo,r_edges]=cellfun(@(x) histcounts(x,100),r_disk_com,'UniformOutput',false);
    r_cents=cellfun(@(x) x(1:end-1)+0.5*diff(x),r_edges,'UniformOutput',false);
    
    for ii=1:2
        den_r_halo{ii}=n_r_halo{ii}./(2*pi*r_cents{ii}.*diff(r_edges{ii})*ddiskheight);     % radial density
    end
    
    % plot density profile
    hfig_com_den_profile=figure();
    hold on;
    for ii=1:2
        plot(r_cents{ii},den_r_halo{ii},'*-');
    end
    box on;
    xlabel('$r$ [m]');
    ylabel('$\rho(r)$ [m$^{-3}$]');
    titlestr='Radial density profile about the halo';
    title(titlestr);
    legend({'$m_F=0$','$m_F=1$'});
    
    % 2. BEC-thermal
    
    
    %% Capture halos
    % TODO
    %   check for preexisting saved files:
    %   should have equal LOAD + HALO
    [halo,bec,culled,HFIG{length(HFIG)+1},errflag]=halo_capture(txy,fout,configs,verbose);
    
    
    %% Cull aliased hits
    % TODO - check
    
    % cull from centred zxy0 halo data
    alias_deadz=100e-9*configs.misc.vel_z;      % dZ in 100 ns
    bool_alias=cellfun(@(x) findalias(x,alias_deadz),halo.zxy0,'UniformOutput',false);
    halo_zxy0_filt=cellfun(@(x,y) x(~y,:),halo.zxy0,bool_alias,'UniformOutput',false);
    
    % TODO - better var management
    halo.zxy0_orig=halo.zxy0;       % store original
    halo.zxy0=halo_zxy0_filt;       % alias filtered
    
    %% Characterise halos
    % TODO
    %   num in halo (+ mode occupancy - needs corr volume)
    %   halo thickness + distribution
    %   number distribution
    
    %% Quick ellipsoid fit
    efit_flag='';
    efit=cell(2,1);
    for ii=1:2
        [ecent,eradii,evecs,v,echi2]=ellipsoid_fit(circshift(vertcat(halo.zxy0{:,ii}),-1,2),efit_flag);
        % NOTE: efit params are in XYZ coordinate frame!
        efit{ii}.cent=ecent; 
        efit{ii}.rad=eradii;
        efit{ii}.evecs=evecs;
        efit{ii}.v=v;
        efit{ii}.chi2=echi2;
    end
    
    % plot ellipsoid fit
    if configs.flags.graphics
        COLOR={'b','r'};
        hfig_ellipsoid_fit=figure();
        for ii=1:2
            subplot(1,2,ii);
            hold on;
            
            % draw data
            xyz_temp=circshift(vertcat(halo.zxy0{:,ii}),-1,2);
            plot_zxy(halo.zxy0(:,ii),1e5,5,COLOR{ii});
            
            % draw fit
            mind=min(xyz_temp);
            maxd=max(xyz_temp);
            nsteps=50;
            step=(maxd-mind)/nsteps;
            [x,y,z]=meshgrid( linspace( mind(1) - step(1), maxd(1) + step(1), nsteps ), linspace( mind(2) - step(2), maxd(2) + step(2), nsteps ), linspace( mind(3) - step(3), maxd(3) + step(3), nsteps ) );
            v=efit{ii}.v;
            Ellipsoid = v(1) *x.*x +   v(2) * y.*y + v(3) * z.*z + ...
                2*v(4) *x.*y + 2*v(5)*x.*z + 2*v(6) * y.*z + ...
                2*v(7) *x    + 2*v(8)*y    + 2*v(9) * z;
            p = patch( isosurface( x, y, z, Ellipsoid, -v(10) ) );
            hold off;
            set( p, 'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha',0.5);
            
            view( -70, 40 );
            axis vis3d equal;
            camlight;
            lighting phong;
            xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
        end
        drawnow;
        
        % Update HFIG
        HFIG{length(HFIG)+1}={hfig_ellipsoid_fit};
    end
    
    %% Transform halos to COM frame
    % Initialise k-vector in XYZ coord
    halo_k=cellfun(@(x) circshift(x,-1,2),halo.zxy0,'UniformOutput',false);
    for ii=1:2
        % Centre to fitted ellipsoid
        halo_k(:,ii)=boost_zxy(halo_k(:,ii),-efit{ii}.cent');
        
        % Transform to ellipsoid principal axis (Rotation)
        M_rot=efit{ii}.evecs;   % rotation matrix: X'Y'Z'(principal ellipsoid) --> XYZ
        halo_k(:,ii)=cellfun(@(x) (M_rot\x')',halo_k(:,ii),'UniformOutput',false);     % inverse transform
        
        % Stretch to UNIT sphere: unit in collision wave-vector/momenta
        halo_k(:,ii)=cellfun(@(x) x./repmat(efit{ii}.rad',size(x,1),1),halo_k(:,ii),'UniformOutput',false);
        
        % Reverse transform to original/detector axis
        halo_k(:,ii)=cellfun(@(x) (M_rot*x')',halo_k(:,ii),'UniformOutput',false);
    end
    % transform to ZXY system
    halo_k=cellfun(@(x) circshift(x,1,2),halo_k,'UniformOutput',false);
    
    % Plot k-space mapped halos
    if configs.flags.graphics
        hfig_halo_k=figure();
        plot_zxy(halo_k,1e5);
        
        axis equal;
        axis vis3d;
        view(3);
        xlabel('$K_{x}$ [m]'); ylabel('$K_{y}$ [m]'); zlabel('$K_{z}$ [m]');
        titlestr='k-space mapped halos';
        title(titlestr);
        drawnow;
        
        % Update HFIG
        HFIG{length(HFIG)+1}={hfig_halo_k};
    end
    
    %% Characterise halos - spherically mapped
    Nsc=cell(1,size(halo_k,2));     % total number in halo
    dk=zeros(1,size(halo_k,2));     % halo rms width
%     gfit=cell(1,size(halo_k,2));    % gaussian fit to radial dist
    
    for ii=1:size(halo_k,2)
        [Nsc{ii},dk(ii),~]=halo_characterise(halo_k(:,ii),configs.halo.zcap,verbose);
    end
    
    %% Correct for halo centre
    halo_k0=cell(size(halo_k));     % boosted halo in k-space
    for ii=1:2
        halo_k0(:,ii)=boost_zxy(halo_k(:,ii),configs.halo.boost{ii});   % boost each halo as defined
    end
    
    %% Correlation analysis
    % TODO
    %   check for preexisting saved files passed for above
    if do_corr_analysis
        [corr,HFIG{length(HFIG)+1}]=halo_g2_manager(halo_k0,configs,verbose);
    else
        corr=NaN;       % need to return corr
    end
    
    %% Mode occupancy
    % TODO - uncertainty in mode occupancy
    
    % initialise vars
    wbb=zeros(3,2);     % [g2_BB_rms_width SE]
    I_BB_corr=0;
    sigk=NaN;
    n_mocc=NaN;
    
    % get cart BB correlation task
    for ii=1:numel(configs.corr)
        if isequal(configs.corr{ii}.type.coord,'cart')&&isequal(configs.corr{ii}.type.opt,'BB')
            I_BB_corr=ii;
            break;
        end
    end
    
    % mode occupancy
    if I_BB_corr~=0
        % get cart BB correlation widths
        for ii=1:3      % Z,X,Y 1d g2 fit params
            wbb(ii,:)=abs(corr{I_BB_corr}.fit{ii}(3,:));     % [rms_width,SE]
        end
        % source condensate k width
        sigk=geomean(wbb(:,1))/1.1;
        
        % evaluate mode occupancy
        n_mocc=halo_mocc(1,mean(abs(dk)),mean(vertcat(Nsc{:})),sigk);
    end
    
    %% Save results
    if configs.flags.savedata
        %%% fig
        % TODO - need to be able to save all graphics from each subroutine
        for ii=1:length(HFIG)
            for jj=1:length(HFIG{ii})
                saveas(HFIG{ii}{jj},[configs.files.dirout,'/',sprintf('fig_%d_%d',ii,jj),'.png']);
                saveas(HFIG{ii}{jj},[configs.files.dirout,'/',sprintf('fig_%d_%d',ii,jj),'.fig']);
            end
        end
        clear HFIG;     % clear HFIG graphics handle cell array from workspace
        
        %%% data
        % get all vars in workspace except graphics handles
        allvars=whos;
        tosave=cellfun(@isempty,regexp({allvars.class},'^matlab\.(ui|graphics)\.'));
        
        % doesn't save the whole workspace this way
        save([configs.files.dirout,'/',mfilename,'_data','.mat'],allvars(tosave).name);
    end
    
    %% end of code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t_main_end=toc(t_main_start);   % end of code
    disp('-----------------------------------------------');
    fprintf('Total elapsed time (s): %7.1f\n',t_main_end);
    disp('===================ALL TASKS COMPLETED===================');
    
    err=0;  % everything ok
end
