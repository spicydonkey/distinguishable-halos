%% distinguishable halo analysis
% DKS 15/05/2017
%% TODO
% map to sphere
    
function [halo_k,corr,efit,halo,txy,fout,err]=run_dist_halo(config_file)
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
    HFIG={};
    
    % output directory
    if ~isdir(configs.files.dirout)
        warning(['output directory "',configs.files.dirout,'" does not exist. Creating directory...']);
        mkdir(configs.files.dirout);
    end
    
    %% Load TXY data
    if ~do_next
        % check for existing saved file for preloaded data
        mat_list=dir([configs.files.dirout,'/*.mat']);      % list of saved data
        if size(mat_list,1)==0
            do_next=1;  % no previously saved data files
        else
            % look for file with same load configs
            for ii=1:size(mat_list,1)
                this_file=mat_list(ii).name;
                % TODO - will throw error if $configs is not a variable in data
                S_temp=load(this_file,'configs');  % load configs from prev data
                
                % check if relevant configs equal
                if ~isequal(S_temp.configs.load,configs.load)
                    warning('Raw data: Existing data has different configs. Setting do_next=1.');
                    do_next=1;
                end
            end
            clear S_temp;
        end
    end
    
    if do_next
        [txy,fout,HFIG{length(HFIG)+1}]=load_txy(configs.load.path,configs.load.id,...
            configs.load.window,configs.load.mincount,configs.load.maxcount,...
            configs.load.rot_angle,verbose,configs.flags.graphics);
    end

    %% Capture halos
    % TODO
    %   check for preexisting saved files:
    %   should have equal LOAD + HALO
    [halo,bec,culled,HFIG{length(HFIG)+1},errflag]=halo_capture(txy,fout,configs,verbose);
    

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
    
    %% Correlation analysis
    % TODO
    %   check for preexisting saved files passed for above
    if do_corr_analysis
        corr=halo_g2_manager(halo_k,configs,verbose);
    end
    
    %% Save results
    if configs.flags.savedata
        %%% fig
        % TODO - need to be able to save all graphics from each subroutine
        for ii=1:length(HFIG)
            for jj=1:length(HFIG{ii})
                saveas(HFIG{ii}{jj},[configs.files.dirout,'/',sprintf('fig_%d_%d_',ii,jj),datetimestr,'.fig']);
                saveas(HFIG{ii}{jj},[configs.files.dirout,'/',sprintf('fig_%d_%d_',ii,jj),datetimestr,'.png']);
            end
        end
        clear HFIG;     % clear HFIG graphics handle cell array from workspace
        
        %%% data
        % get all vars in workspace except graphics handles
        allvars=whos;
        tosave=cellfun(@isempty,regexp({allvars.class},'^matlab\.(ui|graphics)\.'));
        
        % doesn't save the whole workspace this way
        save([configs.files.dirout,'/',mfilename,'_output_',datetimestr,'.mat'],allvars(tosave).name);
    end
    
    %% end of code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t_main_end=toc(t_main_start);   % end of code
    disp('-----------------------------------------------');
    fprintf('Total elapsed time (s): %7.1f\n',t_main_end);
    disp('===================ALL TASKS COMPLETED===================');
    
    err=0;  % everything ok
end
