% Main script for the analysis of distinguishable s-wave scattered halos
% DKS

clear all; close all; clc;

%% USER CONFIG
% GENERAL
use_saved_data=1;   %if false will remake the fully processed data files used in analysis
use_txy=1;          %if false will remake the txy_forc files

verbose=2;

% IN/OUTPUTS
% files -  data file
usrconfigs.files.path='C:\Users\HE BEC\Documents\lab\halo_analysis\data\dist_halo\5_separated_higherN_short\d';    % path to unindexed data file (e.g. 'a\b\datadir\datafile')
usrconfigs.files.id=1:100;         % file id numbers to use for analysis
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
    analysis.corr.polar.nBin=25*[1,1];    % num bins (r,theta)
        analysis.corr.polar.lim{1}=[-0.25,0.25];  % radial lim
        analysis.corr.polar.lim{2}=[0,pi/4];      % angular lim
    analysis.corr.cart.nBin=10*[1,1,1];   % num bins (Z,X,Y)
        analysis.corr.cart.lim=[-0.25,0.25];    % lims (x,y symmetric)
        
%% PLOTS
% 3D real space
doplot.real.all=1;      % real space
doplot.real.ind=1:3;    % plots the selection of shots

% 3D k-space (normed)   TODO
doplot.kspace.all=1;    % k-space
doplot.kspace.ind=1:3;  % plots the selection of shots

%% Output management
dir_output=[usrconfigs.files.path,'_output\'];
if ~isdir(dir_output)
    warning(['output directory "',dir_output,'" does not exist. Creating directory...']);
    mkdir(dir_output);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

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
        disp(['Halo ',num2str(i_mj),' radius: [',num2str(mean(vertcat(halo.R{:,i_mj}))),']±[',num2str(std(vertcat(halo.R{:,i_mj}))),']']);
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
    
    saveas(hfig,[dir_output,'4','.png']);
end

%% Cull halo caps
if configs.post.removecap
    nShot=size(halo.k,1);
    for i=1:2
        for j=1:nShot
            ind_cap=abs(halo.k{j,i}(:,1))>configs.post.zcap;
            halo.k{j,i}=halo.k{j,i}(~ind_cap,:);
        end
    end
    
    % PLOTS
    figN=121; dotSize=1;
    hfig=figure(figN);
    scatter_zxy(figN,vertcat(halo.k{:,1}),dotSize,'r');
    scatter_zxy(figN,vertcat(halo.k{:,2}),dotSize,'b');
    title('Z-cap removed (k-space)');
    xlabel('$K_{X}$'); ylabel('$K_{Y}$'); zlabel('$K_{Z}$');
    
    figN=122; dotSize=100;
    figure(figN);
    scatter_zxy(figN,vertcat(halo.k{doplot.kspace.ind,1}),dotSize,'r');
    scatter_zxy(figN,vertcat(halo.k{doplot.kspace.ind,2}),dotSize,'b');
    title(['Z-cap removed Selected ',num2str(length(doplot.kspace.ind)),' shots (k-space)']);
    xlabel('$K_{X}$'); ylabel('$K_{Y}$'); zlabel('$K_{Z}$');
    
    saveas(hfig,[dir_output,'5','.png']);
end

%% Cartesian to Spherical polar conversion
% Build k-space counts in the conventional spherical polar system
% Should be simpler to do correlation analysis in sph pol coord system

% Initialise variable
halo.k_pol=cell(size(halo.k));

% Cart-SphPol Conversion
for i=1:2
    halo.k_pol(:,i)=zxy2pol(halo.k(:,1));   % Polar coord of halo in k-space 
end

%% Correlation analysis
if analysis.corr.run_g2
    %% Cross-halo back-to-back: in (dk,dtheta)
    % Set up bins
    for i=1:2
        bin_edge_pol{i}=linspace(analysis.corr.polar.lim{i}(1),analysis.corr.polar.lim{i}(2),analysis.corr.polar.nBin(i)+1);
        bin_cent_pol{i}=0.5*(bin_edge_pol{i}(1:end-1)+bin_edge_pol{i}(2:end));
    end    
    
    % Evaluate G2 correlation
    [G2_bb_pol_shot,G2_bb_pol_all]=G2_polar(halo.k_pol,bin_edge_pol,'BB',2);
    g2_bb_pol=nShot*G2_bb_pol_shot./G2_bb_pol_all;  %normalise

    % Plot
    [dR_bin,dtheta_bin]=meshgrid(bin_cent_pol{1},bin_cent_pol{2});  % create xy-grid for surf
    
    hfig=figure(11);
    
    subplot(1,3,1);
    surf(dR_bin',dtheta_bin',G2_bb_pol_shot,'edgecolor','none');
    title('X-halo,BB,$\delta \vec{k}$ (pol),shots');
    xlabel('$\delta k$'); ylabel('$\delta\theta$'); zlabel('$G^{(2)}_{BB(0,1)}$');
    axis tight;
    
    subplot(1,3,2);
    surf(dR_bin',dtheta_bin',G2_bb_pol_all,'edgecolor','none');
    title('X-halo,BB,$\delta \vec{k}$ (pol),collated');
    xlabel('$\delta k$'); ylabel('$\delta\theta$'); zlabel('$G^{(2)}_{BB(0,1)}$');
    axis tight;
    
    subplot(1,3,3);
    surf(dR_bin',dtheta_bin',g2_bb_pol,'edgecolor','none');
    title('X-halo,BB,$\delta \vec{k}$ (pol),normalised');
    xlabel('$\delta k$'); ylabel('$\delta\theta$'); zlabel('$g^{(2)}_{BB(0,1)}$');
    axis tight;
    
    saveas(hfig,[dir_output,'6','.png']);
    
    
    %% Cross-halo back-to-back: in Cartesian delta_k
    % Set up bins
    for i=1:3
        bin_edge_cart{i}=linspace(analysis.corr.cart.lim(1),analysis.corr.cart.lim(2),analysis.corr.cart.nBin(i)+1);
        bin_cent_cart{i}=0.5*(bin_edge_cart{i}(1:end-1)+bin_edge_cart{i}(2:end));
    end
    
    [G2_bb_cart_shot,G2_bb_cart_all]=G2_cart(halo.k,bin_edge_cart,'BB',2);
    g2_bb_cart=nShot*G2_bb_cart_shot./G2_bb_cart_all;   % normalise
    
    % Plot
    [dX_bin,dY_bin]=meshgrid(bin_cent_cart{2},bin_cent_cart{3});        % create xy-grid for surf
    
    mid_slice=round(analysis.corr.cart.nBin(1)/2);     % TODO: mid-slice is where diff is zero (for symmetric binning)
    
    hfig=figure(31);
    
    subplot(1,3,1);
    surf(dX_bin',dY_bin',squeeze(G2_bb_cart_shot(mid_slice,:,:)),'edgecolor','none');
    title('X-halo,BB,$\delta \vec{k}$ (cart),shots');
    xlabel('$\delta k_i$'); ylabel('$\delta k_j$'); zlabel('$G^{(2)}_{BB(0,1)}$');
    axis tight;
    
    subplot(1,3,2);
    surf(dX_bin',dY_bin',squeeze(G2_bb_cart_all(mid_slice,:,:)),'edgecolor','none');
    title('X-halo,BB,$\delta \vec{k}$ (cart),collated');
    xlabel('$\delta k_i$'); ylabel('$\delta k_j$'); zlabel('$G^{(2)}_{ALL,BB(0,1)}$');
    axis tight;
    
    subplot(1,3,3);
    surf(dX_bin',dY_bin',squeeze(g2_bb_cart(mid_slice,:,:)),'edgecolor','none');
    title('X-halo,BB,$\delta \vec{k}$ (cart),normalised');
    xlabel('$\delta k_i$'); ylabel('$\delta k_j$'); zlabel('$g^{(2)}_{BB(0,1)}$');
    axis tight;
    
    saveas(hfig,[dir_output,'7','.png']);
    
    %% Cross-halo colinear
end

%% end of code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc;