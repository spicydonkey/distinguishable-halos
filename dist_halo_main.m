% Main script for the analysis of distinguishable s-wave scattered halos
% DKS

clear all; close all; clc;

%% USER CONFIG
% GENERAL
use_saved_data=1;   %if false will remake the fully processed data files used in analysis
use_txy=1;          %if false will remake the txy_forc files

verbose=2;

% ANALYSIS
% TODO - usr configurable switches for analysis
corr.rad_theta.nBin_pol=[10,10];     % number of bins for rad-angular (dr,dtheta)
corr.rad_theta.nBin_cart=50*[1,1,1];    % bins to use for g2 in cartesian (Z,X,Y)

% mode
% TODO
% viewall: overlay all shots for first visualisation
%usrconfigs.mode='viewall';


% IN/OUTPUTS
% files -  data file
usrconfigs.files.path='..\data\test\d';    % path to unindexed data file (e.g. 'a\b\datadir\datafile')
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
usrconfigs.halo.dR{1}=0.25;      % halo fractional thickness each dir (in/out)
usrconfigs.halo.R{2}=10e-3;
usrconfigs.halo.dR{2}=0.25;

%% PLOTS
% 3D real space
doplot.real.all=1;      % real space
doplot.real.ind=1:3;    % plots the selection of shots

% 3D k-space (normed)   TODO
doplot.kspace.all=1;    % k-space
doplot.kspace.ind=1:3;  % plots the selection of shots

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
    figure(figN);
    scatter_zxy(figN,vertcat(halo.zxy{:,1}),dotSize,'r');
    scatter_zxy(figN,vertcat(halo.zxy{:,2}),dotSize,'b');
    scatter_zxy(figN,vertcat(bec.zxy{:,1}),dotSize,'m');
    scatter_zxy(figN,vertcat(bec.zxy{:,2}),dotSize,'c');
    scatter_zxy(figN,vertcat(culled.tail.zxy{:,1}),dotSize,'k');
    scatter_zxy(figN,vertcat(culled.tail.zxy{:,2}),dotSize,'k');
    scatter_zxy(figN,vertcat(culled.fuzz.zxy{:}),dotSize,'k');
    
    title('All shots (real space)');
    xlabel('X'); ylabel('Y'); zlabel('Z');
end
if ~isempty(doplot.real.ind)
    figN=102; dotSize=100;
    figure(figN);
    scatter_zxy(figN,vertcat(halo.zxy{doplot.real.ind,1}),dotSize,'r');
    scatter_zxy(figN,vertcat(halo.zxy{doplot.real.ind,2}),dotSize,'b');
    scatter_zxy(figN,vertcat(bec.zxy{doplot.real.ind,1}),dotSize,'m');
    scatter_zxy(figN,vertcat(bec.zxy{doplot.real.ind,2}),dotSize,'c');
    scatter_zxy(figN,vertcat(culled.tail.zxy{doplot.real.ind,1}),dotSize,'k');
    scatter_zxy(figN,vertcat(culled.tail.zxy{doplot.real.ind,2}),dotSize,'k');
    scatter_zxy(figN,vertcat(culled.fuzz.zxy{doplot.real.ind}),dotSize,'k');
    
    title(['Selected ',num2str(length(doplot.real.ind)),' shots (real space)']);
    xlabel('X'); ylabel('Y'); zlabel('Z');
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
    figure(figN); title('All shots (k"-space)');
    scatter_zxy(figN,vertcat(halo.k{:,1}),dotSize,'r');
    scatter_zxy(figN,vertcat(halo.k{:,2}),dotSize,'b');
    scatter_zxy(figN,vertcat(bec.k{:,1}),dotSize,'m');
    scatter_zxy(figN,vertcat(bec.k{:,2}),dotSize,'c');
    
    title('All halos and BEC (k-space)');
    xlabel('$K_{X}$'); ylabel('$K_{Y}$'); zlabel('$K_{Z}$');
end
if ~isempty(doplot.kspace.ind)
    figN=112; dotSize=100;
    figure(figN); title('All shots (k"-space)');
    scatter_zxy(figN,vertcat(halo.k{doplot.kspace.ind,1}),dotSize,'r');
    scatter_zxy(figN,vertcat(halo.k{doplot.kspace.ind,2}),dotSize,'b');
    scatter_zxy(figN,vertcat(bec.k{doplot.kspace.ind,1}),dotSize,'m');
    scatter_zxy(figN,vertcat(bec.k{doplot.kspace.ind,2}),dotSize,'c');
    
    title(['Selected ',num2str(length(doplot.kspace.ind)),' shots (k-space)']);
    xlabel('$K_{X}$'); ylabel('$K_{Y}$'); zlabel('$K_{Z}$');
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
% define bin edges for correlation binning
%   {1}: diff radius, {2}: diff angle
bin_edge_pol{1}=linspace(-2*configs.halo.dR{1},2*configs.halo.dR{1},corr.rad_theta.nBin_pol(1)+1);      % TODO: max set to 1 since halo is relatively thin and normalised in k-space
bin_edge_pol{2}=linspace(0,pi,corr.rad_theta.nBin_pol(2)+1);

for i=1:3
    bin_edge_cart{i}=linspace(-2*(1+configs.halo.dR{1}),2*(1+configs.halo.dR{1}),corr.rad_theta.nBin_cart(i)+1);     % TODO: 
end

%% Cross-halo back-to-back: in (dk,dtheta)
nHalo=size(halo.k_pol,1);

G2{1}=zeros(corr.rad_theta.nBin_pol);  % initialise G2 (unnormalised)
for i=1:nHalo   % iterate shot-to-shot
    nAtom=size(halo.k_pol{i,1},1);      % number of counts in this halo
    dk_BB_tmp=[];
    for j=1:nAtom     % iterate through each atom in halo#1
        % back-to-back condition
        atom_ref=halo.k_pol{i,1}(j,:);  % this reference atom in k-pol
        
        dk_BB_tmp(:,1)=halo.k_pol{i,2}(:,1)-atom_ref(1);    % dk
        dk_BB_tmp(:,2)=acos(cos(-atom_ref(3)).*cos(halo.k_pol{i,2}(:,3)).*cos(halo.k_pol{i,2}(:,2)-(atom_ref(2)+pi)) ...
            +sin(-atom_ref(3)).*sin(halo.k_pol{i,2}(:,3)));     % dtheta
        
        count_tmp=histcn(dk_BB_tmp,bin_edge_pol{1},bin_edge_pol{2});
        G2{1}=G2{1}+count_tmp;
    end
end

% Normalisation
G2_all{1}=zeros(corr.rad_theta.nBin_pol);  % normalisation
halo_all{1}=vertcat(halo.k_pol{:,1});
halo_all{2}=vertcat(halo.k_pol{:,2});

nAtom=size(halo_all{1},1);
dk_BB_tmp=[];
for j=1:nAtom     % iterate through each atom in halo#1
    % back-to-back condition
    atom_ref=halo_all{1}(j,:);  % this reference atom in k-pol
    
    dk_BB_tmp(:,1)=halo_all{2}(:,1)-atom_ref(1);    % dk
    dk_BB_tmp(:,2)=acos(cos(-atom_ref(3)).*cos(halo_all{2}(:,3)).*cos(halo_all{2}(:,2)-(atom_ref(2)+pi)) ...
        +sin(-atom_ref(3)).*sin(halo_all{2}(:,3)));     % dtheta
    
    count_tmp=histcn(dk_BB_tmp,bin_edge_pol{1},bin_edge_pol{2});
    G2_all{1}=G2_all{1}+count_tmp;
end

g2{1}=G2{1}./G2_all{1}*nHalo;

% DEBUG PLOT
figure(11);
X=0.5*(bin_edge_pol{1}(1:end-1)+bin_edge_pol{1}(2:end));    % bin centers
Y=0.5*(bin_edge_pol{2}(1:end-1)+bin_edge_pol{2}(2:end));
[X,Y]=meshgrid(X,Y);
surf(X',Y',G2{1},'edgecolor','none');
title('');
xlabel('$\delta k$'); ylabel('$\delta\theta$'); zlabel('$G^{(2)}_{BB(0,1)}$');

figure(12);
surf(X',Y',g2{1},'edgecolor','none');
title('');
xlabel('$\delta k$'); ylabel('$\delta\theta$'); zlabel('$g^{(2)}_{BB(0,1)}$');

%% Cross-halo back-to-back: in Cartesian delta_k
nHalo=size(halo.k_pol,1);

G2{3}=zeros(corr.rad_theta.nBin_cart);  % initialise G2 (unnormalised)
for i=1:nHalo   % iterate shot-to-shot
    nAtom=size(halo.k{i,1},1);      % number of counts in this halo
    dk_BB_tmp=[];   
        atom_ref=halo.k{i,1}(j,:);  % this reference atom in k(Z,X,Y)
        
        % set-up dk_BB_tmp in cart coord
        dk_BB_tmp=halo.k{i,2}-repmat(atom_ref,[size(halo.k{i,2},1),1]);    % dk in cart coord
        
        count_tmp=histcn(dk_BB_tmp,bin_edge_cart{1},bin_edge_cart{2},bin_edge_cart{3});
        G2{3}=G2{3}+count_tmp;
    end
end

% Normalisation
G2_all{3}=zeros(corr.rad_theta.nBin_cart);  % normalisation
halo_all_cart{1}=vertcat(halo.k{:,1});
halo_all_cart{2}=vertcat(halo.k{:,2});

nAtom=size(halo_all_cart{1},1);
dk_BB_tmp=[];
for j=1:nAtom     % iterate through each atom in halo#1
    % back-to-back condition
    atom_ref=halo_all_cart{1}(j,:);  % this reference atom in k(Z,X,Y)
    
    % set-up dk_BB_tmp in cart coord
    dk_BB_tmp=halo_all_cart{2}-repmat(atom_ref,[size(halo_all_cart{2},1),1]);    % dk in cart coord
    
    count_tmp=histcn(dk_BB_tmp,bin_edge_cart{1},bin_edge_cart{2},bin_edge_cart{3});
    G2_all{3}=G2_all{3}+count_tmp;
end

g2{3}=G2{3}./G2_all{3}*nHalo;

% DEBUG PLOT
Xcart=0.5*(bin_edge_cart{1}(1:end-1)+bin_edge_cart{1}(2:end));    % bin centers
Ycart=0.5*(bin_edge_cart{2}(1:end-1)+bin_edge_cart{2}(2:end));
[Xcart,Ycart]=meshgrid(Xcart,Ycart);

mid_slice=round(corr.rad_theta.nBin_cart(1)/2);

figure(31);
surf(Xcart',Ycart',squeeze(G2{3}(mid_slice,:,:)),'edgecolor','none');
title('');
xlabel('$\delta ki$'); ylabel('$\delta kj$'); zlabel('$G^{(2)}_{BB(0,1)}$');

figure(32);
surf(Xcart',Ycart',squeeze(g2{3}(mid_slice,:,:)),'edgecolor','none');
title('');
xlabel('$\delta ki$'); ylabel('$\delta kj$'); zlabel('$g^{(2)}_{BB(0,1)}$');

% %% Cross-halo colinear
% G2{2}=zeros(corr.rad_theta.nBin_pol);  % initialise G2 (unnormalised)
% for i=1:nHalo   % iterate shot-to-shot
%     nAtom=size(halo.k_pol{i,1},1);      % number of counts in this halo
%     dk_CL_tmp=[];
%     for j=1:nAtom     % iterate through each atom in halo#1
%         % back-to-back condition
%         atom_ref=halo.k_pol{i,1}(j,:);  % this reference atom in k-pol
%         
%         dk_CL_tmp(:,1)=halo.k_pol{i,2}(:,1)-atom_ref(1);    % dk
%         dk_CL_tmp(:,2)=acos(cos(atom_ref(3)).*cos(halo.k_pol{i,2}(:,3)).*cos(halo.k_pol{i,2}(:,2)-atom_ref(2)) ...
%             +sin(atom_ref(3)).*sin(halo.k_pol{i,2}(:,3)));  % dtheta
%         
%         count_tmp=histcn(dk_CL_tmp,bin_edge{1},bin_edge{2});
%         G2{2}=G2{2}+count_tmp;
%     end
% end
% % DEBUG PLOT
% figure(21);
% surf(X',Y',G2{2},'edgecolor','none');
% title('');
% xlabel('$\delta k$'); ylabel('$\delta\theta$'); zlabel('$G^{(2)}_{CL(0,1)}$');

%% end of code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc;