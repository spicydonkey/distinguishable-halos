% Main script for the analysis of distinguishable s-wave scattered halos
% DKS

clear all; close all; clc;

%% USER CONFIG
% GENERAL
use_saved_data=1;   %if false will remake the fully processed data files used in analysis (takes a while)
use_txy=1;          %if false will remake the txy_forc files

verbose=2;

% DEBUG MODE
use_inverted_pairs=0;

% IN/OUTPUTS
% files -  data file
usrconfigs.files.path='C:\Users\HE BEC\Documents\lab\halo_analysis\data\dist_halo\4_separated_lownum\d';    % path to unindexed data file (e.g. 'a\b\datadir\$DATA_FNAME_TOKEN$')
usrconfigs.files.id=1:3000;         % file id numbers to use for analysis
usrconfigs.files.minCount=100;     % min counts to use for analysis

% MISCELLANEOUS
usrconfigs.misc.vel_z=9.8*0.430;    % atom free-fall vert v at detector hit for T-to-Z conversion

% WINDOW - captures only box-cropped counts in files ([]-->no crop)
% TODO - XY axis is unclear
usrconfigs.window.all_T=[4.905,4.92];       % T [s]
usrconfigs.window.all{1}=[20.678,20.726];   % Z [m] T-window will overide Z
usrconfigs.window.all{2}=[-20e-3,18e-3];    % X [m]
usrconfigs.window.all{3}=[-10e-3,17e-3];    % Y [m]

% DISTINGUISHABLE HALO PARAMS: params specific to data processing and analysis
usrconfigs.bec.pos{1}=[20.7024,4.74e-3,2.72e-3];   % approx condensate locations (z,x,y)
usrconfigs.bec.Rmax{1}=7e-3;    % max condensate sph radius
usrconfigs.bec.dR_tail{1}=0.3;	% BEC tail radial frac diff
usrconfigs.bec.pos{2}=[20.7005,-7.38e-3,6.55e-3];
usrconfigs.bec.Rmax{2}=7e-3;
usrconfigs.bec.dR_tail{2}=0.3;

usrconfigs.halo.R{1}=11e-3;     % estimated radius of halo
usrconfigs.halo.dR{1}=0.1;      % halo fractional thickness each dir (in/out)
usrconfigs.halo.R{2}=10e-3;
usrconfigs.halo.dR{2}=0.1;

% POST
usrconfigs.post.removecap=1;    % remove caps on halo (in Z)
    usrconfigs.post.zcap=0.5;   % z-cutoff (kspace;abs)
    
% ANALYSIS
% g2 correlations
analysis.corr.run_g2=1;
    analysis.corr.polar.nBin=[11,51];    % num bins (r,theta) (USE ODD)
        analysis.corr.polar.lim{1}=0.3*[-1,1];  % radial lim
        analysis.corr.polar.lim{2}=[0,pi];      % angular lim
    analysis.corr.cart.nBin=[51,13,13];   % num bins (Z,X,Y) (USE ODD)
        analysis.corr.cart.lim=0.8*[-1,1];    % lims (z,x,y symmetric)
        
%% PLOTS
% 3D real space
doplot.real.all=0;      % real space
doplot.real.ind=[];    % plots the selection of shots

% 3D k-space (normed)
doplot.kspace.all=0;    % k-space
doplot.kspace.ind=[];  % plots the selection of shots

%% CONSTANTS
W_BB_GI=4e-4;   % g2 bb correlation length as determined from RIK's GI

%% Output management
dir_output=[usrconfigs.files.path,'_output\'];
if ~isdir(dir_output)
    warning(['output directory "',dir_output,'" does not exist. Creating directory...']);
    mkdir(dir_output);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%
t_main_start=tic;

%initalize variables
configs=usrconfigs;    % create an alias to avoid overwriting user's config
missingfiles=zeros(length(configs.files.id),1);     % missing dld/txy file
filestotxy=zeros(length(configs.files.id),1);       % dld files processed to txy
lowcountfiles=zeros(length(configs.files.id),1);    % files with too few counts (skipped in analysis)

%% Load processed data
vars_saved = {'usrconfigs','configs'...
    'halo','bec','culled','errflag',...
    'lowcountfiles','missingfiles','filestotxy',...
    };  % variables important in halo analysis - saved to a .mat file

if use_saved_data     
    if ~fileExists([configs.files.path,'data.mat'])
        warning(['Saved data "', [configs.files.path,'data.mat'], '" does not exist - setting use_saved_data=0 and continue...']);
        use_saved_data=0;
    else
        % Error checks on prev saved file
        S_data = load([configs.files.path,'data.mat']);
        
        is_complete=1;
        for i=1:length(vars_saved)
            is_complete=is_complete&isfield(S_data,vars_saved{i});
        end
        if ~is_complete
            % Completeness of information
            warning('Previously saved file is incomplete - setting use_saved_data=0 and continue...');
            use_saved_data=0;
        elseif ~isequal(usrconfigs,S_data.usrconfigs)
            % Check consistency of requested analysis configs and saved data
            warning('Previously saved data contains settings (windows, files) different to requested - setting use_saved_data=0 and continue...');
            use_saved_data=0;
        else
            % OK load all vars in file
            if verbose>0, disp('Loading saved data...');, end;
            
            % Check saved file is not empty
            if isempty(vertcat(S_data.halo.zxy{:}))
                error('Saved file is empty - terminating program. Try recreating data file with different configs');
            end
            
            % okay to load everything
            load([configs.files.path,'data.mat']);  
        end
        
        clear S_data;
    end
end

%% Pre-processing
if ~use_saved_data
    %% Prepare TXY files
    % TXY files (usually with suffix _txy_forc and #ID attached to the core file name)
    %   are built by processing raw DLD output files from the TDC
    % TXY files should contain complete result of experiment re. atom positions
    % A number of checks are performed in this stage to filter out missing
    % files and shots of low count, etc.
    if ~use_txy
        % Create TXY from raw DLD files
        if verbose>0, disp('Creating TXY files...');, end;
        for i=1:length(configs.files.id)
            %check for source - raw DLD file
            if ~fileExists([configs.files.path,num2str(configs.files.id(i)),'.txt'])
                missingfiles(i)=1;
                if verbose>0
                    warning(['use_txy=0 but raw DLD data does not exist - file #', num2str(configs.files.id(i))]);
                end
                continue;
            end
            
            % create TXY file from the raw DLD file
            dld_raw_to_txy(configs.files.path,configs.files.id(i),configs.files.id(i));
            filestotxy(i)=1;
            
            if verbose>1
                disp(['TXY file for #',num2str(configs.files.id(i)),' is created.']);
            end
        end
    else
        % Using pre-built TXY files
        if verbose>0, disp('Using pre-built TXY files...');, end;
        % Check for missing files
        for i=1:length(configs.files.id)
            if ~fileExists([configs.files.path,'_txy_forc',num2str(configs.files.id(i)),'.txt'])
                % missing TXY file
                missingfiles(i)=1;
                if verbose>0
                    warning(['TXY file for #',num2str(configs.files.id(i)),' is missing.']);
                end
            end
        end
    end
    
    % Check for low counts
    if verbose>0,disp('Checking for low counts in TXY files...');,end;
    for i=1:length(configs.files.id)
        if fileExists([configs.files.path,'_txy_forc',num2str(configs.files.id(i)),'.txt'])
            txy_temp=txy_importer(configs.files.path,configs.files.id(i));    % import data in this txy to memory
            % Check count in this file
            if size(txy_temp,1)<configs.files.minCount
                lowcountfiles(i)=1;
                if verbose>0
                    warning(['Low count detected in file #',num2str(configs.files.id(i))]);
                end
            end
        end
    end
    
    %% Process raw TXY data
    % This is where the raw TXY data (from DLD) must be processed to do further
    % analysis
    % capture the distinguishable halos
    % find properties of halo, etc
    % SAVE the refined data and results
    % TODO: name vars appropriately and update the vars_saved, etc above
    
    importokfiles=~(missingfiles|lowcountfiles);
    configs.files.idok=configs.files.id(importokfiles);     % ID's for OK files
    
    if verbose>0,disp('Processing TXY files to generate data for analysis...');,end;
    [halo,bec,culled,errflag]=distinguish_halo(configs,verbose);  % all the data processing on raw-TXY to generate halo data
    
    
    %% Save processed data
    if verbose>0,disp('Saving data..');,end;
    % delete existing data
    if fileExists([configs.files.path,'data.mat'])
        delete([configs.files.path,'data.mat']);    
    end
    % save data
    save([configs.files.path,'data.mat'],vars_saved{1});    % must create *.mat without append option
    for i = 1:length(vars_saved)
        if ~exist(vars_saved{i},'var')
            warning(['Variable "',vars_saved{i},'" does not exist.']);
            continue
        end
        save([configs.files.path,'data.mat'],vars_saved{i},'-append');
    end
end

% SUMMARY
% TXY-preprocessing
importokfiles=~(missingfiles|lowcountfiles);
if verbose>1
    disp('===================IMPORT SUMMARY===================');
    disp([int2str(sum(importokfiles)),' imported ok']);
    disp([int2str(sum(lowcountfiles)),' files with low counts']);
    disp([int2str(sum(missingfiles)),' files missing']);
    disp([int2str(sum(filestotxy)),' files converted to txy']);
    disp('====================================================');
end

% Distinguishable halo processing
if verbose>0
    % Calculated BEC centre
    disp('=================PROCESSING SUMMARY=================');
    disp([num2str(sum(errflag)),' bad shot(s) were discarded during data processing.']);
    for i_mj=1:2
        disp(['BEC ',num2str(i_mj),' centres: [',num2str(mean(vertcat(bec.cent{:,i_mj}))),']±[',num2str(std(vertcat(bec.cent{:,i_mj}))),']']);
        %disp(['HALO ',num2str(i_mj),' centres: [',num2str(mean(vertcat(halo.cent{:,i_mj}))),']±[',num2str(std(vertcat(halo.cent{:,i_mj}))),']']);   % not really useful with this algorithm
        % TODO - with optimisation, report halo radii
        disp(['Halo ',num2str(i_mj),' (R,dR_rms): [',num2str(mean(vertcat(halo.R{:,i_mj}))),']±[',num2str(std(vertcat(halo.R{:,i_mj}))),']']);
    end
    disp('====================================================');
end

% Plot processed counts
if doplot.real.all
    figN=101; dotSize=1;
    hfig=figure(figN);
    scatter_zxy(figN,vertcat(halo.zxy{:,1}),dotSize,'r');
    scatter_zxy(figN,vertcat(halo.zxy{:,2}),dotSize,'b');
    scatter_zxy(figN,vertcat(bec.zxy{:,1}),dotSize,'m');
    scatter_zxy(figN,vertcat(bec.zxy{:,2}),dotSize,'c');
    scatter_zxy(figN,vertcat(culled.tail.zxy{:,1}),dotSize,'k');
    scatter_zxy(figN,vertcat(culled.tail.zxy{:,2}),dotSize,'k');
    scatter_zxy(figN,vertcat(culled.fuzz.zxy{:}),dotSize,'k');
    
    title('All shots (real space)');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    
    saveas(hfig,[dir_output,'1','.fig']);
    saveas(hfig,[dir_output,'1','.png']);
end
if ~isempty(doplot.real.ind)
    figN=102; dotSize=100;
    hfig=figure(figN);
    scatter_zxy(figN,vertcat(halo.zxy{doplot.real.ind,1}),dotSize,'r');
    scatter_zxy(figN,vertcat(halo.zxy{doplot.real.ind,2}),dotSize,'b');
    scatter_zxy(figN,vertcat(bec.zxy{doplot.real.ind,1}),dotSize,'m');
    scatter_zxy(figN,vertcat(bec.zxy{doplot.real.ind,2}),dotSize,'c');
    scatter_zxy(figN,vertcat(culled.tail.zxy{doplot.real.ind,1}),dotSize,'k');
    scatter_zxy(figN,vertcat(culled.tail.zxy{doplot.real.ind,2}),dotSize,'k');
    scatter_zxy(figN,vertcat(culled.fuzz.zxy{doplot.real.ind}),dotSize,'k');
    
    title(['Selected ',num2str(length(doplot.real.ind)),' shots (real space)']);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    
    saveas(hfig,[dir_output,'2','.fig']);
    saveas(hfig,[dir_output,'2','.png']);
end

%% k-space conversion
% transforming to k-space will simplify many analysis for s-wave scattering
% TODO
% Implementation:
%   - translate to set halo centre as the origin
%   - normalisation - isometric scale to let halo radius (i.e. recoil k) to 1

% initialise vars
halo.k=cell(size(halo.zxy));
bec.k=cell(size(bec.zxy));

% Translation
for i=1:2
    halo.k(:,i)=zxy_translate(halo.zxy(:,i),halo.cent(:,i),-1);
    bec.k(:,i)=zxy_translate(bec.zxy(:,i),halo.cent(:,i),-1);
end

% Isometric scaling (normalisation)
for i=1:size(halo.k,1)    % iterate shots
    for j=1:2   % iterate internal states
        halo.k{i,j}=halo.k{i,j}/halo.R{i,j}(1);     % normalise in k-space
        bec.k{i,j}=bec.k{i,j}/halo.R{i,j}(1);
    end
end

% Plot counts in k-space
if doplot.kspace.all
    figN=111; dotSize=1;
    hfig=figure(figN);
    scatter_zxy(figN,vertcat(halo.k{:,1}),dotSize,'r');
    scatter_zxy(figN,vertcat(halo.k{:,2}),dotSize,'b');
    scatter_zxy(figN,vertcat(bec.k{:,1}),dotSize,'m');
    scatter_zxy(figN,vertcat(bec.k{:,2}),dotSize,'c');
    
    title('All halos and BEC (k-space)');
    xlabel('$K_{X}$'); ylabel('$K_{Y}$'); zlabel('$K_{Z}$');
    
    saveas(hfig,[dir_output,'3','.fig']);
    saveas(hfig,[dir_output,'3','.png']);
end
if ~isempty(doplot.kspace.ind)
    figN=112; dotSize=100;
    hfig=figure(figN);
    scatter_zxy(figN,vertcat(halo.k{doplot.kspace.ind,1}),dotSize,'r');
    scatter_zxy(figN,vertcat(halo.k{doplot.kspace.ind,2}),dotSize,'b');
    scatter_zxy(figN,vertcat(bec.k{doplot.kspace.ind,1}),dotSize,'m');
    scatter_zxy(figN,vertcat(bec.k{doplot.kspace.ind,2}),dotSize,'c');
    
    title(['Selected ',num2str(length(doplot.kspace.ind)),' shots (k-space)']);
    xlabel('$K_{X}$'); ylabel('$K_{Y}$'); zlabel('$K_{Z}$');
    
    saveas(hfig,[dir_output,'4','.fig']);
    saveas(hfig,[dir_output,'4','.png']);
end

%% Cull halo caps
if configs.post.removecap
    for i=1:2
        for j=1:size(halo.k,1)
            ind_cap=abs(halo.k{j,i}(:,1))>configs.post.zcap;
            halo.k{j,i}=halo.k{j,i}(~ind_cap,:);
        end
    end
    
    % PLOTS
    if doplot.kspace.all
        figN=121; dotSize=1;
        hfig=figure(figN);
        scatter_zxy(figN,vertcat(halo.k{:,1}),dotSize,'r');
        scatter_zxy(figN,vertcat(halo.k{:,2}),dotSize,'b');
        title('Z-cap removed (k-space)');
        xlabel('$K_{X}$'); ylabel('$K_{Y}$'); zlabel('$K_{Z}$');
        
        saveas(hfig,[dir_output,'5','.png']);
        saveas(hfig,[dir_output,'5','.fig']);
    end
    
    if ~isempty(doplot.kspace.ind)
        figN=122; dotSize=100;
        hfig=figure(figN);
        scatter_zxy(figN,vertcat(halo.k{doplot.kspace.ind,1}),dotSize,'r');
        scatter_zxy(figN,vertcat(halo.k{doplot.kspace.ind,2}),dotSize,'b');
        title(['Z-cap removed Selected ',num2str(length(doplot.kspace.ind)),' shots (k-space)']);
        xlabel('$K_{X}$'); ylabel('$K_{Y}$'); zlabel('$K_{Z}$');
        
        saveas(hfig,[dir_output,'6','.png']);
        saveas(hfig,[dir_output,'6','.fig']);
    end
end

%% Remove ZERO-halo count shots for analysis
zeroShot=false(1,size(halo.k,1));
for i=1:size(halo.k,1)
    if size(halo.k{i,1},1)*size(halo.k{i,2},1)==0
        zeroShot(i)=1;
    end
end
n_zero_shot=sum(zeroShot);
if n_zero_shot>0
    warning([num2str(n_zero_shot),' shots had zero-counts after processing. Removing from analysis...']);
end
halo.k=halo.k(~zeroShot,:);

%% Summarise count numbers in halos for analysis
halo_count=zeros(size(halo.k));
for i=1:size(halo.k,1)
    for j=1:size(halo.k,2)
        halo_count(i,j)=size(halo.k{i,j},1);    % get halo counts
    end
end

% Summary for post processed halos
if verbose>0
    disp('=================POST PROCESSING SUMMARY=================');
    disp(['Total shots suitable for analysis: ',num2str(size(halo.k,1))]);
    disp(['Counts per halo: [',num2str(mean(halo_count,1)),']±[',num2str(std(halo_count,1)),']']);
    disp(['Total counts: ',num2str(sum(sum((halo_count)))),': [',num2str(sum(halo_count,1)),']']);
    disp('=========================================================');
end

%% Create test data to force BB particles
% Mirror the non-magnetic halo (#1) around origin which is more spherical and
% captured better - creates completely "entangled" data, even for
% noises and background
if use_inverted_pairs
    halo_invert=cell(size(halo.k,1),1);
    for i=1:size(halo.k,1)
        halo_invert{i}=-halo.k{i,1};    % mirror (inverted) image of halo#1
    end
    halo_inv_pair=cell(size(halo.k,1),2);
    halo_inv_pair(:,1)=halo.k(:,1);
    halo_inv_pair(:,2)=halo_invert;
    clear halo_invert;
end

%% Cartesian to Spherical polar conversion
% Build k-space counts in the conventional spherical polar system
% Should be simpler to do correlation analysis in sph pol coord system

% Initialise variable
halo.k_pol=cell(size(halo.k));

% Cart-SphPol Conversion
for i=1:2
    halo.k_pol(:,i)=zxy2pol(halo.k(:,i));   % Polar coord of halo in k-space 
end

% DEBUG - manipulating halos
% for i=1:size(halo.k_pol,1)
% %     % scale halo 1
% %     halo.k_pol{i,1}(:,1)=1.1*halo.k_pol{i,1}(:,1);
%     
% %     % rotate halo 1
% %     halo.k_pol{i,1}(:,2)=halo.k_pol{i,1}(:,2)+pi/10;
% end

if use_inverted_pairs
    halo_inv_pair_pol=cell(size(halo.k));
    for i=1:2
        halo_inv_pair_pol(:,i)=zxy2pol(halo_inv_pair(:,i));
    end
end

%% Correlation analysis
if analysis.corr.run_g2
    %% Cross-halo back-to-back: in (dk,dtheta)
    % Set up bins
    for i=1:2
        bin_edge_pol{i}=linspace(analysis.corr.polar.lim{i}(1),analysis.corr.polar.lim{i}(2),analysis.corr.polar.nBin(i)+1);
        bin_cent_pol{i}=0.5*(bin_edge_pol{i}(1:end-1)+bin_edge_pol{i}(2:end));
    end    
    
    % TODO for asymetric binning?
    ind_zero_pol=round((analysis.corr.polar.nBin+1)/2); % zero-cent'd bin index for sampling 2D-g2 (dr-dtheta)
    
    % Evaluate G2 correlation
    if use_inverted_pairs
        % this is to check g2 for ideal, completely B-B paired halos (but
        % possibly with many pair occupations)
        [G2_bb_pol_shot,G2_bb_pol_all]=G2_polar(halo_inv_pair_pol,bin_edge_pol,'BB',2);
    else
        [G2_bb_pol_shot,G2_bb_pol_all]=G2_polar(halo.k_pol,bin_edge_pol,'BB',2);
    end
    g2_bb_pol=size(halo.k,1)*G2_bb_pol_shot./G2_bb_pol_all;  %normalise

    % Plot
    [dR_bin,dtheta_bin]=meshgrid(bin_cent_pol{1},bin_cent_pol{2});  % create xy-grid for surf
    
    hfig=figure(11);
    
    subplot(1,3,1);
    surf(dR_bin',dtheta_bin',G2_bb_pol_shot,'edgecolor','none');
    title('X-halo,BB,$\delta \vec{k}$ (pol),shots');
    xlabel('$\delta k$'); ylabel('$\delta\theta$'); zlabel('$G^{(2)}_{BB(0,1)}$');
    axis tight;
    shading interp;
    
    subplot(1,3,2);
    surf(dR_bin',dtheta_bin',G2_bb_pol_all,'edgecolor','none');
    title('X-halo,BB,$\delta \vec{k}$ (pol),collated');
    xlabel('$\delta k$'); ylabel('$\delta\theta$'); zlabel('$G^{(2)}_{BB(0,1)}$');
    axis tight;
    shading interp;
    
    subplot(1,3,3);
    surf(dR_bin',dtheta_bin',g2_bb_pol,'edgecolor','none');
    title('X-halo,BB,$\delta \vec{k}$ (pol),normalised');
    xlabel('$\delta k$'); ylabel('$\delta\theta$'); zlabel('$g^{(2)}_{BB(0,1)}$');
    axis tight;
    shading interp;
    
    saveas(hfig,[dir_output,'7','.fig']);
    saveas(hfig,[dir_output,'7','.png']);
    
    % dk-integrated g2(dtheta)
    g2_dtheta=size(halo.k_pol,1)*sum(G2_bb_pol_shot,1)./sum(G2_bb_pol_all,1);
    hfig=figure(12);
    plot(bin_cent_pol{2},g2_dtheta,'*');
    
    title('$\Delta k$-integrated X-halo BB correlation');
    xlabel('$\Delta\theta$'); ylabel('$\bar{g}^{(2)}_{BB,(0,1)}$');
    xlim([0,pi]); ylim auto;
    
    saveas(hfig,[dir_output,'7_2','.fig']);
    saveas(hfig,[dir_output,'7_2','.png']);
    
    %% TEST
    [G2test1,G2test2]=G2_angular(halo.k,bin_edge_pol,verbose);
    g2test=size(halo.k,1)*G2test1./G2test2;
    
    hfig=figure(13);
    
    subplot(1,3,1);
    surf(dR_bin',dtheta_bin',G2test1,'edgecolor','none');
    title('X-halo,BB,$\delta \vec{k}$ (pol),shots');
    xlabel('$\delta k$'); ylabel('$\delta\theta$'); zlabel('$G^{(2)}_{BB(0,1)}$');
    axis tight;
    shading interp;
    
    subplot(1,3,2);
    surf(dR_bin',dtheta_bin',G2test2,'edgecolor','none');
    title('X-halo,BB,$\delta \vec{k}$ (pol),collated');
    xlabel('$\delta k$'); ylabel('$\delta\theta$'); zlabel('$G^{(2)}_{BB(0,1)}$');
    axis tight;
    shading interp;
    
    subplot(1,3,3);
    surf(dR_bin',dtheta_bin',g2test,'edgecolor','none');
    title('X-halo,BB,$\delta \vec{k}$ (pol),normalised');
    xlabel('$\delta k$'); ylabel('$\delta\theta$'); zlabel('$g^{(2)}_{BB(0,1)}$');
    axis tight;
    shading interp;
    
    %% Cross-halo back-to-back: in Cartesian delta_k
    % Set up bins
    for i=1:3
        bin_edge_cart{i}=linspace(analysis.corr.cart.lim(1),analysis.corr.cart.lim(2),analysis.corr.cart.nBin(i)+1);
        bin_cent_cart{i}=0.5*(bin_edge_cart{i}(1:end-1)+bin_edge_cart{i}(2:end));
    end
    
    if use_inverted_pairs
        % this is to check g2 for ideal, completely B-B paired halos (but
        % possibly with many pair occupations)
        [G2_bb_cart_shot,G2_bb_cart_all]=G2_cart(halo_inv_pair,bin_edge_cart,'BB',2);
    else
        [G2_bb_cart_shot,G2_bb_cart_all]=G2_cart(halo.k,bin_edge_cart,'BB',2);
    end
    g2_bb_cart=size(halo.k,1)*G2_bb_cart_shot./G2_bb_cart_all;   % normalise

    % Plot
    [dX_bin,dY_bin]=meshgrid(bin_cent_cart{2},bin_cent_cart{3});        % create xy-grid for surf
    
    % TODO for asymetric binning?
    ind_zero_cart=round((analysis.corr.cart.nBin+1)/2); % zero-cent'd bin index for sampling 3D-g2 
    
    hfig=figure(21);
    
    subplot(1,3,1);
    surf(dX_bin',dY_bin',squeeze(G2_bb_cart_shot(ind_zero_cart(1),:,:)),'edgecolor','none');
    title('X-halo,BB,$\delta \vec{k}$ (cart),shots');
    xlabel('$\delta k_i$'); ylabel('$\delta k_j$'); zlabel('$G^{(2)}_{BB(0,1)}$');
    axis tight;
    shading interp;
    
    subplot(1,3,2);
    surf(dX_bin',dY_bin',squeeze(G2_bb_cart_all(ind_zero_cart(1),:,:)),'edgecolor','none');
    title('X-halo,BB,$\delta \vec{k}$ (cart),collated');
    xlabel('$\delta k_i$'); ylabel('$\delta k_j$'); zlabel('$G^{(2)}_{ALL,BB(0,1)}$');
    axis tight;
    shading interp;
    
    subplot(1,3,3);
    surf(dX_bin',dY_bin',squeeze(g2_bb_cart(ind_zero_cart(1),:,:)),'edgecolor','none');
    title('X-halo,BB,$\delta \vec{k}$ (cart),normalised');
    xlabel('$\delta k_i$'); ylabel('$\delta k_j$'); zlabel('$g^{(2)}_{BB(0,1)}$');
    axis tight;
    shading interp;
    
    saveas(hfig,[dir_output,'8','.fig']);
    saveas(hfig,[dir_output,'8','.png']);
    
    % 1-D g2 in Z
    hfig=figure(22);
    plot(bin_cent_cart{1},g2_bb_cart(:,ind_zero_cart(2),ind_zero_cart(3)),'*');
    title('X-halo BB correlations in $Z$-axis');
    xlabel('$\Delta K_z$'); ylabel('$g^{(2)}_{BB,(0,1)}$');
    
    saveas(hfig,[dir_output,'8_1','.fig']);
    saveas(hfig,[dir_output,'8_1','.png']);
    
    %% Solo BB (polar)
    % Back-to-back g2 correlations in the s-wave scattered particles of
    % single collision source
    
    % Halo 1
	[G2_bb_solo_pol_shot,G2_bb_solo_pol_all]=G2_polar(halo.k_pol(:,1),bin_edge_pol,'BB',2);
    g2_bb_solo_pol=size(halo.k_pol,1)*G2_bb_solo_pol_shot./G2_bb_solo_pol_all;  %normalise

    % plot
    hfig=figure(31);
    
    subplot(1,3,1);
    surf(dR_bin',dtheta_bin',G2_bb_solo_pol_shot,'edgecolor','none');
    title('Single-halo,BB,$\delta \vec{k}$ (pol),shots');
    xlabel('$\delta k$'); ylabel('$\delta\theta$'); zlabel('$G^{(2)}_{BB(0,0)}$');
    axis tight;
    shading interp;
    
    subplot(1,3,2);
    surf(dR_bin',dtheta_bin',G2_bb_solo_pol_all,'edgecolor','none');
    title('Single-halo,BB,$\delta \vec{k}$ (pol),collated');
    xlabel('$\delta k$'); ylabel('$\delta\theta$'); zlabel('$G^{(2)}_{BB(0,0)}$');
    axis tight;
    shading interp;
    
    subplot(1,3,3);
    surf(dR_bin',dtheta_bin',g2_bb_solo_pol,'edgecolor','none');
    title('Single-halo,BB,$\delta \vec{k}$ (pol),normalised');
    xlabel('$\delta k$'); ylabel('$\delta\theta$'); zlabel('$g^{(2)}_{BB(0,0)}$');
    axis tight;
    shading interp;
    
    saveas(hfig,[dir_output,'9','.fig']);
    saveas(hfig,[dir_output,'9','.png']);
    
    % dk-integrated g2(dtheta)
    g2_dtheta_solo=size(halo.k_pol,1)*sum(G2_bb_solo_pol_shot,1)./sum(G2_bb_solo_pol_all,1);
    hfig=figure(32);
    plot(bin_cent_pol{2},g2_dtheta_solo,'*');
    
    title('$\Delta k$-integrated Single-halo BB correlation');
    xlabel('$\Delta\theta$'); ylabel('$\bar{g}^{(2)}_{BB,(0,0)}$');
    xlim([0,pi]); ylim auto;
    
    saveas(hfig,[dir_output,'9_2','.fig']);
    saveas(hfig,[dir_output,'9_2','.png']);
    
    %% TEST - single component angular G2
    [G2testsolo1,G2testsolo2]=G2_angular(halo.k(:,1),bin_edge_pol,verbose);
    g2testsolo=size(halo.k,1)*G2testsolo1./G2testsolo2;
    
    hfig=figure(33);
    
    subplot(1,3,1);
    surf(dR_bin',dtheta_bin',G2testsolo1,'edgecolor','none');
    title('Single-halo,BB,$\delta \vec{k}$ (pol),shots');
    xlabel('$\delta k$'); ylabel('$\delta\theta$'); zlabel('$G^{(2)}_{BB(0,0)}$');
    axis tight;
    shading interp;
    
    subplot(1,3,2);
    surf(dR_bin',dtheta_bin',G2testsolo2,'edgecolor','none');
    title('Single-halo,BB,$\delta \vec{k}$ (pol),collated');
    xlabel('$\delta k$'); ylabel('$\delta\theta$'); zlabel('$G^{(2)}_{BB(0,0)}$');
    axis tight;
    shading interp;
    
    subplot(1,3,3);
    surf(dR_bin',dtheta_bin',g2testsolo,'edgecolor','none');
    title('Single-halo,BB,$\delta \vec{k}$ (pol),normalised');
    xlabel('$\delta k$'); ylabel('$\delta\theta$'); zlabel('$g^{(2)}_{BB(0,0)}$');
    axis tight;
    shading interp;
    
    %% Solo BB (cartesian)
    % Back-to-back g2 correlations in the single s-wave scatterer source
    
    % Halo 1
    [G2_bb_solo_cart_shot,G2_bb_solo_cart_all]=G2_cart(halo.k(:,1),bin_edge_cart,'BB',2);
    g2_bb_solo_cart=size(halo.k,1)*G2_bb_solo_cart_shot./G2_bb_solo_cart_all;   % normalise
    
    % plot
    hfig=figure(41);
    
    subplot(1,3,1);
    surf(dX_bin',dY_bin',squeeze(G2_bb_solo_cart_shot(ind_zero_cart(1),:,:)),'edgecolor','none');
    title('Single-halo,BB,$\delta \vec{k}$ (cart),shots');
    xlabel('$\delta k_i$'); ylabel('$\delta k_j$'); zlabel('$G^{(2)}_{BB(0,0)}$');
    axis tight;
    shading interp;
    
    subplot(1,3,2);
    surf(dX_bin',dY_bin',squeeze(G2_bb_solo_cart_all(ind_zero_cart(1),:,:)),'edgecolor','none');
    title('Single-halo,BB,$\delta \vec{k}$ (cart),collated');
    xlabel('$\delta k_i$'); ylabel('$\delta k_j$'); zlabel('$G^{(2)}_{ALL,BB(0,0)}$');
    axis tight;
    shading interp;
    
    subplot(1,3,3);
    surf(dX_bin',dY_bin',squeeze(g2_bb_solo_cart(ind_zero_cart(1),:,:)),'edgecolor','none');
    title('Single-halo,BB,$\delta \vec{k}$ (cart),normalised');
    xlabel('$\delta k_i$'); ylabel('$\delta k_j$'); zlabel('$g^{(2)}_{BB(0,0)}$');
    axis tight;
    shading interp;
    
    saveas(hfig,[dir_output,'10','.fig']);
    saveas(hfig,[dir_output,'10','.png']);
    
   	% 1-D g2 in Z
    hfig=figure(42);
    plot(bin_cent_cart{1},g2_bb_solo_cart(:,ind_zero_cart(2),ind_zero_cart(3)),'*');
    title('Single-halo BB correlations in $Z$-axis');
    xlabel('$\Delta K_z$'); ylabel('$g^{(2)}_{BB,(0,0)}$');
    
    saveas(hfig,[dir_output,'10_1','.fig']);
    saveas(hfig,[dir_output,'10_1','.png']);
    
    %% Compare g2 between scattering partners to non
    hfig=figure(51);
    plot(bin_cent_pol{2},g2_bb_pol(ind_zero_pol(1),:),'-*');
    hold on;
    plot(bin_cent_pol{2},g2_bb_solo_pol(ind_zero_pol(1),:),'-*');
    xlim([0,pi]);   ylim auto;
    title('Correlations in distinguishable $s$-wave scattering');
    xlabel('$\Delta\theta$'); ylabel('$\bar{g}^{(2)}$');
    legend({'(0,1)','(0,0)'},'Location','northwest');
    
    saveas(hfig,[dir_output,'11','.fig']);
    saveas(hfig,[dir_output,'11','.png']);
    
end

%% end of code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_main_end=toc(t_main_start);
disp('-----------------------------------------------');
fprintf('Total elapsed time (s): %7.1f\n',t_main_end);
disp('-------------------COMPLETED-------------------');