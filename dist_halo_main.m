% Main script for the analysis of distinguishable s-wave scattered halos
% DKS

% TODO - implement clearvars to override usrconfigs from a caller script
% clear all;
clearvars -except OR* override
close all; clc;

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
usrconfigs.halo.dR{1}=0.15;      % halo fractional thickness each dir (in/out)
usrconfigs.halo.R{2}=10e-3;
usrconfigs.halo.dR{2}=0.15;

if exist('override','var')
    % make this more generic
    usrconfigs.halo.R{1}=OR_configs.halo.R{1}(mod(OR_i-1,10)+1);
    usrconfigs.halo.R{2}=OR_configs.halo.R{2}(ceil(OR_i/10));
end

% POST
usrconfigs.post.removecap=1;    % remove caps on halo (in Z)
    usrconfigs.post.zcap=0.5;   % z-cutoff (kspace;abs)
    
%% PLOTS
% 3D real space
do_plot.real.all=0;      % real space
do_plot.real.ind=[1:100];    % plots the selection of shots

% 3D k-space (normed)
do_plot.kspace.all=0;    % k-space
do_plot.kspace.ind=[1:100];  % plots the selection of shots
    
%% ANALYSIS
% correlation analysis
do_corr_analysis=1;
    % 1. Cross-halo rad/angular correlations
    analysis.corr.type{1}.comp=[1,2];           % components to analysis: cross halo 1,2
    analysis.corr.type{1}.coord='angular';      % angular coordinate
    analysis.corr.type{1}.opt=[];               % 'angular' has no optional feature atm
    analysis.corr.lim{1}{1}=0.3*[-1,1];  % bin limits - radial separation
    analysis.corr.lim{1}{2}=[0,pi];      % bin limits - angular separation
    analysis.corr.nBin{1}=[11,51];          % number of bins
    
    % 2. Cross-halo cartesian BB-correlations
    analysis.corr.type{2}.comp=[1,2];           % components to analysis: cross halo 1,2
    analysis.corr.type{2}.coord='cart';         % Cartesian (ZXY)
    analysis.corr.type{2}.opt='BB';             % BB / CL
    analysis.corr.lim{2}{1}=0.8*[-1,1]; % bin limits - Z
    analysis.corr.lim{2}{2}=0.8*[-1,1]; % bin limits - X
    analysis.corr.lim{2}{3}=0.8*[-1,1]; % bin limits - Y
    analysis.corr.nBin{2}=[51,13,13];   % number of bins

    % 3. Single-halo (1) rad/angular correlations
    analysis.corr.type{3}.comp=1;           % single component
    analysis.corr.type{3}.coord='angular';
    analysis.corr.type{3}.opt=[];         
    analysis.corr.lim{3}{1}=0.3*[-1,1];  % bin limits - radial separation
    analysis.corr.lim{3}{2}=[0,pi];      % bin limits - angular separation
    analysis.corr.nBin{3}=[11,51];          % number of bins
    
    % 4. Single-halo (1) cartesian BB-correlations
    analysis.corr.type{4}.comp=1;
    analysis.corr.type{4}.coord='cart';        
    analysis.corr.type{4}.opt='BB';
    analysis.corr.lim{4}{1}=0.8*[-1,1]; % bin limits - Z
    analysis.corr.lim{4}{2}=0.8*[-1,1]; % bin limits - X
    analysis.corr.lim{4}{3}=0.8*[-1,1]; % bin limits - Y
    analysis.corr.nBin{4}=[51,13,13];   % number of bins
    
%     % #. <TEMPLATE ANGULAR>
%     analysis.corr.type{#}.comp=[1,2];           % components to analysis: cross halo 1,2
%     analysis.corr.type{#}.coord='angular';      % angular coordinate
%     analysis.corr.type{#}.opt=[];               % 'angular' has no optional feature atm
%     analysis.corr.lim{#}{1}=0.3*[-1,1];  % bin limits - radial separation
%     analysis.corr.lim{#}{2}=[0,pi];      % bin limits - angular separation
%     analysis.corr.nBin{#}=[11,51];          % number of bins
%    
%     % #. <TEMPLATE CARTESIAN>
%     analysis.corr.type{#}.comp=[1,2];           % components to analysis: cross halo 1,2
%     analysis.corr.type{#}.coord='cart';         % Cartesian (ZXY)
%     analysis.corr.type{#}.opt='BB';             % BB / CL
%     analysis.corr.lim{#}{1}=0.8*[-1,1]; % bin limits - Z
%     analysis.corr.lim{#}{2}=0.8*[-1,1]; % bin limits - X
%     analysis.corr.lim{#}{3}=0.8*[-1,1]; % bin limits - Y
%     analysis.corr.nBin{#}=[51,13,13];   % number of bins
    


%% CONSTANTS
W_BB_GI=4e-4;   % g2 bb correlation length as determined from RIK's GI

%% Output management
dir_output=[usrconfigs.files.path,'_output\'];
if ~isdir(dir_output)
    warning(['output directory "',dir_output,'" does not exist. Creating directory...']);
    mkdir(dir_output);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%
t_main_start=tic;   % start timer for main script

%initalize variables
configs=usrconfigs;    % create an alias to avoid overwriting user's config
files.missing=false(length(configs.files.id),1);        % missing dld/txy file
files.build_txy=false(length(configs.files.id),1);      % dld files processed to txy
files.lowcount=false(length(configs.files.id),1);       % files with too few counts (skipped in analysis)

%% Load processed data
vars_saved = {'usrconfigs','configs'...
    'halo','bec','culled','errflag',...
    'files',...
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
            warning('Previously saved data contains settings (e.g. windows, files) different to requested - setting use_saved_data=0 and continue...');
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
    clear is_complete;
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
                files.missing(i)=1;
                if verbose>0
                    warning(['use_txy=0 but raw DLD data does not exist - file #', num2str(configs.files.id(i))]);
                end
                continue;
            end
            
            % create TXY file from the raw DLD file
            dld_raw_to_txy(configs.files.path,configs.files.id(i),configs.files.id(i));
            files.build_txy(i)=1;
            
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
                files.missing(i)=1;
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
                files.lowcount(i)=1;
                if verbose>0
                    warning(['Low count detected in file #',num2str(configs.files.id(i))]);
                end
            end
        end
    end
    clear txy_temp;
    
    %% Process raw TXY data
    % This is where the raw TXY data (from DLD) must be processed to do further
    % analysis
    % capture the distinguishable halos
    % find properties of halo, etc
    % SAVE the refined data and results
    % TODO: name vars appropriately and update the vars_saved, etc above
    
    importokfiles=~(files.missing|files.lowcount);
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
importokfiles=~(files.missing|files.lowcount);
if verbose>1
    disp('===================IMPORT SUMMARY===================');
    disp([int2str(sum(importokfiles)),' imported ok']);
    disp([int2str(sum(files.lowcount)),' files with low counts']);
    disp([int2str(sum(files.missing)),' files missing']);
    disp([int2str(sum(files.build_txy)),' files converted to txy']);
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
if do_plot.real.all
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
if ~isempty(do_plot.real.ind)
    figN=102; dotSize=100;
    hfig=figure(figN);
    scatter_zxy(figN,vertcat(halo.zxy{do_plot.real.ind,1}),dotSize,'r');
    scatter_zxy(figN,vertcat(halo.zxy{do_plot.real.ind,2}),dotSize,'b');
    scatter_zxy(figN,vertcat(bec.zxy{do_plot.real.ind,1}),dotSize,'m');
    scatter_zxy(figN,vertcat(bec.zxy{do_plot.real.ind,2}),dotSize,'c');
    scatter_zxy(figN,vertcat(culled.tail.zxy{do_plot.real.ind,1}),dotSize,'k');
    scatter_zxy(figN,vertcat(culled.tail.zxy{do_plot.real.ind,2}),dotSize,'k');
    scatter_zxy(figN,vertcat(culled.fuzz.zxy{do_plot.real.ind}),dotSize,'k');
    
    title(['Selected ',num2str(length(do_plot.real.ind)),' shots (real space)']);
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

% Get mean halo radius
for i=1:2
    halo_temp=vertcat(halo.k{:,i});
    R_halo_temp(i)=mean(sqrt(sum(halo_temp.^2,2)));
end
clear halo_temp;

% Isometric scaling (normalisation)
for i=1:size(halo.k,1)    % iterate shots
    for j=1:2   % iterate internal states
        %halo.k{i,j}=halo.k{i,j}/halo.R{i,j}(1);     % TODO: SHOULDN'T SCALE LIKE THIS! normalise in k-space
        %bec.k{i,j}=bec.k{i,j}/halo.R{i,j}(1);
        halo.k{i,j}=halo.k{i,j}/R_halo_temp(j);     % normalise to population mean
        bec.k{i,j}=bec.k{i,j}/R_halo_temp(j);
    end
end

% Plot counts in k-space
if do_plot.kspace.all
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
if ~isempty(do_plot.kspace.ind)
    figN=112; dotSize=100;
    hfig=figure(figN);
    scatter_zxy(figN,vertcat(halo.k{do_plot.kspace.ind,1}),dotSize,'r');
    scatter_zxy(figN,vertcat(halo.k{do_plot.kspace.ind,2}),dotSize,'b');
    scatter_zxy(figN,vertcat(bec.k{do_plot.kspace.ind,1}),dotSize,'m');
    scatter_zxy(figN,vertcat(bec.k{do_plot.kspace.ind,2}),dotSize,'c');
    
    title(['Selected ',num2str(length(do_plot.kspace.ind)),' shots (k-space)']);
    xlabel('$K_{X}$'); ylabel('$K_{Y}$'); zlabel('$K_{Z}$');
    
    saveas(hfig,[dir_output,'4','.fig']);
    saveas(hfig,[dir_output,'4','.png']);
end

%% Cull halo caps
if configs.post.removecap
    for i=1:2
        for j=1:size(halo.k,1)
            ind_cap_tmp=abs(halo.k{j,i}(:,1))>configs.post.zcap;
            halo.k{j,i}=halo.k{j,i}(~ind_cap_tmp,:);
        end
    end
    
    % PLOTS
    if do_plot.kspace.all
        figN=121; dotSize=1;
        hfig=figure(figN);
        scatter_zxy(figN,vertcat(halo.k{:,1}),dotSize,'r');
        scatter_zxy(figN,vertcat(halo.k{:,2}),dotSize,'b');
        title('Z-cap removed (k-space)');
        xlabel('$K_{X}$'); ylabel('$K_{Y}$'); zlabel('$K_{Z}$');
        
        saveas(hfig,[dir_output,'5','.png']);
        saveas(hfig,[dir_output,'5','.fig']);
    end
    
    if ~isempty(do_plot.kspace.ind)
        figN=122; dotSize=100;
        hfig=figure(figN);
        scatter_zxy(figN,vertcat(halo.k{do_plot.kspace.ind,1}),dotSize,'r');
        scatter_zxy(figN,vertcat(halo.k{do_plot.kspace.ind,2}),dotSize,'b');
        title(['Z-cap removed Selected ',num2str(length(do_plot.kspace.ind)),' shots (k-space)']);
        xlabel('$K_{X}$'); ylabel('$K_{Y}$'); zlabel('$K_{Z}$');
        
        saveas(hfig,[dir_output,'6','.png']);
        saveas(hfig,[dir_output,'6','.fig']);
    end
end
clear ind_cap_tmp;

%% Remove ZERO-halo count shots for analysis
zero_final_halo=false(1,size(halo.k,1));
for i=1:size(halo.k,1)
    if size(halo.k{i,1},1)*size(halo.k{i,2},1)==0
        zero_final_halo(i)=1;   % flag fully processed halos with 0 counts
    end
end
nzero_final_tmp=sum(zero_final_halo);
if nzero_final_tmp>0
    warning([num2str(nzero_final_tmp),' shots had zero-counts after processing. Removing from analysis...']);
end
halo.k=halo.k(~zero_final_halo,:);
clear nzero_final_tmp;

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
if do_corr_analysis
    % TODO *_tmp variables should be saved to analysis.corr. at the end
    %       code should be generic

    %% Do tasked correlation analysis 
    for iCorr=1:length(analysis.corr.type)        
        % Set up bins
        bin_dim=length(analysis.corr.lim{iCorr});  % get binning dims
        bin_edge_tmp=cell(1,bin_dim);
        bin_cent_tmp=cell(1,bin_dim);
        for i=1:bin_dim
            % make bin edge and centre vectors
            bin_edge_tmp{i}=linspace(analysis.corr.lim{iCorr}{i}(1),...
                analysis.corr.lim{iCorr}{i}(2),analysis.corr.nBin{iCorr}(i)+1);
            bin_cent_tmp{i}=0.5*(bin_edge_tmp{i}(1:end-1)+bin_edge_tmp{i}(2:end));
        end
        
        % Evaluate G2 correlation
        % TODO - insert debug code for G2 and halo.k manipulator
        [G2_shot_tmp,G2_all_tmp]=G2_caller(halo.k(:,analysis.corr.type{iCorr}.comp),...
            bin_edge_tmp,analysis.corr.type{iCorr}.coord,analysis.corr.type{iCorr}.opt,verbose);
        g2_tmp=size(halo.k,1)*G2_shot_tmp./G2_all_tmp;      % normalised g2
        
        % Get results
        analysis.corr.bEdge{iCorr}=bin_edge_tmp;
        analysis.corr.bCent{iCorr}=bin_cent_tmp;
        analysis.corr.G2shot{iCorr}=G2_shot_tmp;
        analysis.corr.G2all{iCorr}=G2_all_tmp;
        analysis.corr.g2{iCorr}=g2_tmp;
    end
    % clear workspace
    clear bin_dim bin_edge_tmp bin_cent_tmp G2_shot_tmp G2_all_tmp g2_tmp;
    
    %% Plot the original g2 correlation function
    for iCorr=1:length(analysis.corr.type)
        nfig_tmp=10+iCorr;  % g2 figures start from figure 11
        hfig=plotCorr(nfig_tmp,analysis.corr,iCorr);
        
        % save figs
        fname_str=['corr_',num2str(iCorr)];
        saveas(hfig,[dir_output,fname_str,'.fig']);
        saveas(hfig,[dir_output,fname_str,'.png']);
    end
    
    %% 1D correlation profile and Gaussian fit
    for iCorr=1:length(analysis.corr.type)
        nfig_tmp=20+iCorr;  % 1D starts from fig 21
        hfig=figure(nfig_tmp);
        ax=gca;
        
        this_corr_type=analysis.corr.type{iCorr}.coord;
        % Get integrated or sliced 1D correlation profile
        if isequal(this_corr_type,'angular')
            % integrate dk
            g2_1d_tmp=size(halo.k,1)*sum(analysis.corr.G2shot{iCorr},1)./sum(analysis.corr.G2all{iCorr},1);
            
            plot(analysis.corr.bCent{iCorr}{2},g2_1d_tmp,'*');
            
            hold(ax,'on');
            title_str=['$\Delta k$-integrated, ',...
                '(',num2str(analysis.corr.type{iCorr}.comp),') halos'];
            title(title_str);
            xlabel('$\Delta\theta$'); ylabel('$\bar{g}^{(2)}$');
            xlim([0,pi]); ylim auto;
            
            % Gaussian fit
            param0=[4,pi,0.1,1];     % fit estimate [amp,mu,sigma,offset]
            [fitparam_tmp,fit_g2_tmp]=gaussfit(analysis.corr.bCent{iCorr}{2},g2_1d_tmp,param0,verbose);
            plot(ax,fit_g2_tmp.x,fit_g2_tmp.y,'r');     % plot the fit
            hold(ax,'off');
            
        elseif isequal(this_corr_type,'cart')
            % TODO - do in X/Y?
            % Get line through Z-axis
            ind_zero_tmp=round((analysis.corr.nBin{iCorr}+1)/2);    % zero-cent'd bin index for sampling 3D-g2
            plot(analysis.corr.bCent{iCorr}{1},analysis.corr.g2{iCorr}(:,ind_zero_tmp(2),ind_zero_tmp(3)),'*');
            
            hold(ax,'on');
            title_str=['(',num2str(analysis.corr.type{iCorr}.comp),') halos, ',...
                analysis.corr.type{iCorr}.opt,', ','$Z$-axis'];
            title(title_str);
            xlabel('$\Delta K_z$'); ylabel(['$g^{(2)}_{',analysis.corr.type{iCorr}.opt,'}$']);
            
            % Gaussian fit
            param0=[4,0,0.1,1];     % fit estimate [amp,mu,sigma,offset]
            [fitparam_tmp,fit_g2_tmp]=gaussfit(analysis.corr.bCent{iCorr}{1},analysis.corr.g2{iCorr}(:,ind_zero_tmp(2),ind_zero_tmp(3)),param0,verbose);
            plot(ax,fit_g2_tmp.x,fit_g2_tmp.y,'r');     % plot the fit
            hold(ax,'off');
            
        else
            warning('SOMETHING IS WRONG!');
        end
        
        % Get fit params
        analysis.corr.fit{iCorr}=fitparam_tmp;
        
        % Save figs
        fname_str=['corr1d_',num2str(iCorr)];
        saveas(hfig,[dir_output,fname_str,'.fig']);
        saveas(hfig,[dir_output,fname_str,'.png']);
    end
    clear g2_1d_tmp param0 fitparam_tmp fit_g2_tmp ax this_corr_type;
    
    %% Save data
    save([configs.files.path,'data.mat'],'analysis','-append');
end
clear iCorr nfig_tmp fname_str;

%% end of code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_main_end=toc(t_main_start);   % end of code
disp('-----------------------------------------------');
fprintf('Total elapsed time (s): %7.1f\n',t_main_end);
disp('===================ALL TASKS COMPLETED===================');