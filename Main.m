%======================================90char=============================================
%+++++++++++++++++++++++++++++++++++++++ABOUT+++++++++++++++++++++++++++++++++++++++++++++
%Halo analysis for investigating the spontanious/ stimulated transistion

%This program will import the halo data using HaloImportData then
%plot
%for everything inc the bec's (turn on with save_all_points)
%   the 3d dist
%   tof dist
%the Time integrated histogram
%the radial distributionsin the halo
%the halo 3d sidtribution (plot3d_halo)

%then bin will bin  the halos by the number of counts for each bin we calc
%then correlation function (turn on with find_correlation)
%correlation
%squezing
%radial width

%then plot the output of these as a function of number bin

%++++++++++++++++++++++++++++++++++++To Be Improved+++++++++++++++++++++++++++++++++++++++

%look at pos correction as it seems to be dancing arround a fair bit

%take radial fft

%radial correlations

%calculate center of mass of each halo, look for trend with number in halo

%should find the halo radius from the distance between bec and condensate
%and compare with fit to get the velocity right on


%sqz
    %use proper expresssion for the uncert
    
%move to arrays for xy limits

%move to array for spherical cordinates instead of sep variables

%% -----------------------START user var---------------------------
use_txy=1;                  %if false will remake the txy_forc files
use_saved_halos=1;          %if false will remake the halo files (useful for any change to the halo cut)

clean_radialy=1;            %only keep the halo values within a certian radius

isverbose=1;                %print progress ect
    progress_scaling=50;    %what fraction of the time the progress bar should update

plot_sph_dist=1;            %plot the spherical distibution hitograms( radial, azm ,elev)
    fit_rad_dist=0;         %fit the radial distribution with a gaussian +offset & linear or just est from avg and sd
    azm_bins=300;            %number of bins for histogram and azm bin
    radial_width_azm_dep=0; %plot the radial width as a function of azm angle
        radial_width_plots=0;%plot each azm bin and delay by 1s to show
    fit_sine_azm=0;              %fit a sine wave to the counts and width
plot3d_hist=0;              %plot a 3d histogram of the halo collapsed in time
plot2d_hist=0;              %plot a 2d histogram(image) of the halo colapsed in time
    plot2d_hist_binw=0.0001; %bin width for hist in meters
    plot2d_hist_gauss=0.000;%apply gaussian blur in meters, for off set to zero otherwise specify blur radius
    plot2d_hist_each_shot=0;%plot for each shot indiv. to look for phase grains, will save
plot3d_halo=1;              %display all the combined halos as a 3d plot
plot_counts_dist=0;         %plot a histogram of the counts in the halo
movies3d=0;                 %make movies of the 3d plots %TODO don't choke analysis & do movie of final data
    
find_sqz=0;
find_correlation=1;         %find the correlation in x,y,z
    corr.norm=1;            %normalize the correlations
    corr.fit=1;             %fit the correlaton with a gaussian
    %params for correlations
    corr.yy=linspace(-0.005,0.005,100);    %define the bin size
    corr.dx=0.0005;         %define the condition in other axes
    corr.dt=corr.dx;
    
split_by_halocounts=0;
    halocounts_min=3;
    halocounts_max=1000;
    halocounts_bins=5;
    plots_for_each_bin=1;
    
% files
files.count_min=100;        %min mumber of counts to read in 
files.rot_angle=0.64;       %rotation anle from the txy data 
files.do_pos_correction=1;	%find the halo pos from the ceneter of bragg orders (position correction)
                            %if false will not radialy mask
files.save_all_points=0;    %create array that off all the points usefull for intial investigation and plot TOF and 3d

files.path='..\data\multihalos\multihalos_';    % path to unindexed data file (e.g. 'a\b\datadir\datafile')
files.numstart=1;           %start file num
files.numtoimport=3700;       %number of files to import
files.velocity=9.8*0.430;   %should be 2*9.8*0.6
%files.velocity=sqrt(2*9.8*0.7);%*50;

% windows
windows.bragg.tmin=0.204;
windows.bragg.tmax=0.207;

windows.bec.tmin=0.217;
windows.bec.tmax=0.22;

windows.halo.tmin=0.208;
windows.halo.tmax=0.216;

windows.all.xmin=-0.04;
windows.all.xmax=0.04;
windows.all.ymin=-0.04;
windows.all.ymax=0.04;

windows.halo.rmin=0.025;
windows.halo.rmax=0.03;

windows.reflections=0; %allows reflections throguh the mask
windows.halo.reflecrmin=0.010;

%for the all points the time windows can be different, this is usefull for
%initialy finding where the bec is in time
tmin_allpoints=windows.bragg.tmin;
tmax_allpoints=windows.bec.tmax;

%------------------------END user var------------------------------


%% Path setting
dir_repo = fileparts(which(mfilename));     % directory of main code repo
addpath(genpath(dir_repo));     % add all subdirectories to search path


%% Main code
tic
close all; clc;

%import the halo data
[halo_centered_cells,bec_bragg,all_points]=HaloReflecImportData(files,windows,use_txy,use_saved_halos,isverbose,progress_scaling);
halo_centered=vertcat(halo_centered_cells{:});  % collated halos

%if all the points were saved then a TOF and 3d plot will be done
if files.save_all_points
    all_points=vertcat(all_points{:}); 
    
    %mask the dat 
    figure(4);
    mask=all_points(:,1)>tmin_allpoints & all_points(:,1)<tmax_allpoints;
    all_points=all_points(mask,:);
    %plot TOF
    hist(all_points(:,1),10000)
    %set(gcf,'Position',[400 100 600 600])
    set(gcf,'Color',[1 1 1]);
    
    %here i multiply by velocity for time plot comment this out
    all_points=all_points.*repmat([files.velocity 1 1],size(all_points,1),1);
    
    %plot the halo and bec
    figure(5)
    scatter3(all_points(:,1),all_points(:,2),all_points(:,3),1,'k.','SizeData',1)
    axis equal;
    axis vis3d;
    set(gcf,'Position',[200 50 800 600])
    set(gcf,'Color',[1 1 1]);
    xlabel('time*vel')
    ylabel('X(m)')
    zlabel('y(m)')
    if movies3d
        OptionZ.FrameRate=15;OptionZ.Duration=25;OptionZ.Periodic=true;
        CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'halo_and_bec',OptionZ)
    end
end

%plots for all the count in the halo (all shots combined)
%plot the halo
if plot3d_halo
    figure(1)
    scatter3(halo_centered(:,1),halo_centered(:,2),halo_centered(:,3),1,'k.')
    %set(gcf,'Position',[200 50 800 600])
    set(gcf,'Color',[1 1 1]);
    axis vis3d;
    axis equal;
    xlabel('time*vel (m)')
    ylabel('X(m)')
    zlabel('y(m)')
    if movies3d
        OptionZ.FrameRate=15;OptionZ.Duration=25;OptionZ.Periodic=true;
        CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'halo',OptionZ)
    end
end

% halo histogram in (x,y)
if plot3d_hist
    figure(2)
    %hist3(halo_centered(:,2:3),[100 100])
    hist3(halo_centered(:,2:3),[100 100])
    set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
    %set(gcf,'Position',[400 100 600 600])
    set(gcf,'Color',[1 1 1]);
    xlabel('X(m)')
    ylabel('Y(m)')
    zlabel('Counts/Files')
end

% halo histogram (2D figure)
if plot2d_hist
    if ~plot2d_hist_each_shot
        figure(11)
        
        XEdges=windows.all.xmin:plot2d_hist_binw:windows.all.xmax;
        YEdges=windows.all.ymin:plot2d_hist_binw:windows.all.ymax;
        
        [counts,centers]=hist3(halo_centered(:,2:3),'edges',{XEdges,YEdges});
        
        % gaussian blur
        if  plot2d_hist_gauss~=0
            filter=fspecial('gaussian',round(10*plot2d_hist_gauss/plot2d_hist_binw),plot2d_hist_gauss/plot2d_hist_binw);
            counts=imfilter(counts, filter, 'replicate');
        end
        
        %imagesc seems to plot the wrong way round so we transpose here
        imagesc(centers{1},centers{2},transpose(counts));
        set(gcf,'Color',[1 1 1]);
        title(['Halo Count Dist Combined Blur=',num2str(plot2d_hist_gauss)])
        
        % plot halo histogram for individual shot and save
    else
        %calc the filter kernel out of the loop
        if  plot2d_hist_gauss~=0
            filter=fspecial('gaussian',round(10*plot2d_hist_gauss/plot2d_hist_binw),plot2d_hist_gauss/plot2d_hist_binw);
        end
        XEdges=windows.all.xmin:plot2d_hist_binw:windows.all.xmax;
        YEdges=windows.all.ymin:plot2d_hist_binw:windows.all.ymax;
        figure(11)
        for n=1:length(halo_centered_cells)
            current_file_str = num2str(files.numstart+n-1);
            temp_halo=halo_centered_cells{n};
            
            if length(temp_halo)>300
                [counts,centers]=hist3(temp_halo(:,2:3),'edges',{XEdges,YEdges});
                if  plot2d_hist_gauss~=0
                    counts=imfilter(counts, filter, 'replicate');
                end
                imagesc(centers{1},centers{2},transpose(counts))
                set(gcf,'Color',[1 1 1]);
                title(['Halo Count Dist ',int2str(n),' Blur=',num2str(plot2d_hist_gauss),'m'])
                saveas(gcf,[files.path,current_file_str,'_2d_hist.jpg'])
            else
                disp(['file number ',current_file_str,' has too few counts to plot 2dhist'])
            end
        end
        clear temp_halo filter;
    end
end
clear XEdges YEdges counts centers;

% get atom counts per halo
halo_counts=cellfun(@(x) size(x,1),halo_centered_cells);

% plot histogram of counts in halo
if plot_counts_dist
    figure(12);
    set(gcf,'Color',[1 1 1]);
    hist(halo_counts,100);
    xlabel('Halo Counts');
    ylabel('Number of shots');
    title('Halo Count Dist');
end

% TODO - generalise this process
% categorise halos into different atom numbers (for squeezing?)
if split_by_halocounts
	if isverbose
        fprintf('Categorising halos into atom number bins... ');
    end
    
    halo_numcat={};
    
    %define the bins
    halo_bin_edge=linspace(halocounts_min,halocounts_max,halocounts_bins+1);  % bin edges
    halo_bin_cent=halo_bin_edge(1:end-1)+diff(halo_bin_edge)/2;     % bin centers
    
    % TODO - halo_centered_cells_count_bined really necessary for analysis?
    % bin halo_counts
    halo_numcat={};     %this will contan the halo_bins for each count range
    for n=1:length(halo_bin_cent)
        halo_numcat{n}=halo_centered_cells(halo_counts<halo_bin_edge(n+1) & halo_counts>halo_bin_edge(n));
    end
%     clear mask counts_min counts_max;
    verbose_corr=0;
    
    if isverbose
        fprintf('Done!\n--------------------------------------------------\n');
    end
else
    halo_numcat={halo_centered_cells};  % TODO - if no filtering req'd then shouldn't create this var - it's not categorised
    halo_bin_cent=[NaN];   % TODO - this var is dangerous to refer to, if no binning was done at all
    verbose_corr=1;
end

% Correlation analysis
if find_correlation   
    if isverbose
        fprintf('Starting correlation analysis\n');
    end
    
    corr_params=[];
    
    if split_by_halocounts
        % correlation analysis is conducted separately for halos of similar counts
        if ~verbose_corr
            parfor_progress(length(halo_numcat));
        end
        
        for n=1:length(halo_numcat)
            [~,corr_params(n,:,:,:)]=CalcCorr(halo_numcat{n}, corr, verbose_corr);
            saveas(gcf,[files.path,'_Correlations_Bin',num2str(halo_bin_cent(n),'%5.0f'),'_counts.jpg'])
            saveas(gcf,[files.path,'_Correlations_Bin',num2str(halo_bin_cent(n),'%5.0f'),'_counts.fig'])
            if ~verbose_corr
                parfor_progress;
            end
        end
        
        % TODO - make user control for generation of below figures
        if length(halo_numcat)>2 && corr.fit
            figure(100);
            set(gcf,'Color',[1 1 1]);
            subplot(2,1,1)
            errorbar(halo_bin_cent,corr_params(:,1,1,1),corr_params(:,1,1,2))
            hold on;
            errorbar(halo_bin_cent,corr_params(:,2,1,1),corr_params(:,2,1,2),'r')
            errorbar(halo_bin_cent,corr_params(:,3,1,1),corr_params(:,3,1,2),'g')
            hold off;
            xlabel('Counts')
            ylabel('Corr Amp')
            legend('Y','X','Z')     % TODO - figure out this var ordering
            
            subplot(2,1,2)
            errorbar(halo_bin_cent,abs(corr_params(:,1,3,1)),corr_params(:,1,3,2))
            hold on;
            errorbar(halo_bin_cent,abs(corr_params(:,2,3,1)),corr_params(:,2,3,2),'r')
            errorbar(halo_bin_cent,abs(corr_params(:,3,3,1)),corr_params(:,3,3,2),'g')
            hold off;
            xlabel('Counts')
            ylabel('Corr Width (m)')
            legend('Y','X','Z')
        end
        
        if ~verbose_corr
            parfor_progress(0);
        end
    
    % DKS
    else
        % correlation analysis done across all halos
        [~,corr_params(1,:,:,:)]=CalcCorr(halo_numcat{1}, corr, verbose_corr);
        saveas(gcf,[files.path,'_corr.jpg'])
        saveas(gcf,[files.path,'_corr.fig'])
    end

    if isverbose
        fprintf('Done!\n--------------------------------------------------\n');
    end
end

% Squeezing
if find_sqz
    disp('Starting analysis for squeezing...');
    
    sqz_params={};
    if split_by_halocounts
        parfor_progress(length(halo_numcat));
        for n=1:length(halo_numcat)
            sqz_params{n}=squezing(halo_numcat{n},0);
            saveas(gcf,[files.path,'_Sqz_Bin_',num2str(halo_bin_cent(n),'%5.0f'),'_counts.jpg'])
            saveas(gcf,[files.path,'_Sqz _Bin',num2str(halo_bin_cent(n),'%5.0f'),'_counts.fig'])
            parfor_progress;
        end
        parfor_progress(0);
        
        % TODO - make user control for generation of below figures
        if length(halo_numcat)>2
            figure(101)
            set(gcf,'Color',[1 1 1]);
            errorbar(halo_bin_cent,cellfun(@(x) x{2}(1,1),sqz_params),cellfun(@(x) x{2}(1,2),sqz_params))
            xlabel('Counts')
            ylabel('Mean Norm Var (opst Bins)')
            title('Sqz')
            hold on
            plot(halo_bin_cent,cellfun(@(x) x{2}(1,3),sqz_params),'r')
        end
        
    % DKS
    else
        sqz_params{1}=squezing(halo_numcat{1},isverbose);
        saveas(gcf,[files.path,'_sqz.jpg']);
        saveas(gcf,[files.path,'_sqz.fig']);
    end
    
    if isverbose
        fprintf('Done!\n--------------------------------------------------\n');
    end
end

% Spherical plot of number distribution
if plot_sph_dist 
    rad_fit_params=[];
    for n=1:length(halo_numcat)
%         halo_centered_cells_count_bined_single=halo_numcat{n};    $ DKS -
%         rename variable appropriately
%         halo_single_tmp=halo_numcat{n};   % waste of tmp variable
        
        %here we plot the radial distribution of all halos combined
%         halo_single_comb=vertcat(halo_single_tmp{:});     % this should
%         be the temp var to use
        halo_numcat_tmp = vertcat(halo_numcat{n}{:});
%         bin_counts=length(halo_numcat_tmp);   % usused variable
        [halo_radial,halo_azm,halo_elev]=ConvToSph(halo_numcat_tmp);    % could use cart2sph and solve wrapping problem
        
        figure(6);
        set(gcf,'Color',[1 1 1]);
        %set(gcf,'Position',[400 100 600 600])
        
        subplot(2,2,1);
        hist(halo_elev/pi,100);
        title('Elevation');
        set(gcf,'Color',[1 1 1]);
        xlabel('Angle (Rad)/pi');
        ylabel('Counts');
        
        subplot(2,2,2)
        if azm_bins>1
            [azmcounts,azmedges]=histcounts(halo_azm/pi,azm_bins);
            azmcenters=0.5*(azmedges(1:end-1)+azmedges(2:end));
            plot(azmcenters,azmcounts);
            title('Azm Bin Counts');
            set(gcf,'Color',[1 1 1]);
            xlabel('Angle (Rad)/pi');
            ylabel('Counts');
            if fit_sine_azm
                disp('fit sine');
                %fit_params=fit_sine(xdata,ydata,amp_guess,phase_guess,freq_guess,isverbose)
                azm_counts_fit_params(n,:,:)=fit_azm_sine(azmcenters,azmcounts,range(azmcounts),mean(azmcounts),0,1,1);
            end
        else
            warning('azm_bins should be set to greater than 1 for histogram.');
        end
        
        subplot(2,2,3)
        [radcounts,radcenters]=hist(halo_radial,100);
        plot(radcenters,radcounts);
        title('Com. Radial Width')
        set(gcf,'Color',[1 1 1]);
        xlabel('Distance From Cen. (m)')
        ylabel('Counts')
        if fit_rad_dist
            %coef*SE of {'amp','mu','sig','off','grad'}=rad_fit(xdata,ydata,mu_guess,width_guess,amp_guess)
            %estimate the fit params by assuming that it is gaussian like
            rad_fit_params(n,:,:)=rad_fit(radcenters,radcounts,mean(halo_radial),std(halo_radial),max(radcounts),1);
        end
        
        %want to calc the radial width as a function of the azm width
        %to start with will just calc the std radialy then maybe do a fit
        %if there are enough counts
        if radial_width_azm_dep
            halo_radial_std_azmbin=zeros(1,length(azmedges)-1);
            for m=1:(length(azmedges)-1)
                mask_tmp=halo_azm/pi>azmedges(m) & halo_azm/pi<azmedges(m+1);
                halo_radial_azmbin=halo_radial(mask_tmp);
                %can either find the sd of the radial for this azm bin or
                %better yet do a fit
                if fit_rad_dist  
                    [radcounts,radcenters]=hist(halo_radial_azmbin,100);
                    if radial_width_plots
                        figure(7)
                        plot(radcenters,radcounts,'+')   
                    end
                    %coef*SE of {'amp','mu','sig','off','grad'}=rad_fit(xdata,ydata,mu_guess,width_guess,amp_guess,isverbose)
                    %estimate the fit params by assuming that it is gaussian like
                    rad_fit_azmbin_params(n,m,:,:)=rad_fit(radcenters,radcounts,mean(halo_radial_azmbin),std(halo_radial_azmbin),max(radcounts),radial_width_plots);
                    if radial_width_plots
                        pause(0.1) 
                    end
                else
                   rad_fit_azmbin_params(n,m,:,1)=[max(radcounts),mean(halo_radial_azmbin),std(halo_radial_azmbin),0,0];
                    %should set params to best gueses from 
                    %rad_fit_azmbin_params(n,m,:,:)
                    %halo_radial_std_azmbin(m)=std(halo_radial_azmbin);
                end               

            end
            figure(6)
            subplot(2,2,4)
            errorbar(azmcenters,rad_fit_azmbin_params(n,:,3,1),rad_fit_azmbin_params(n,:,3,2))
            title('Azm. Bin Rad. Width')
            xlabel('Angle(rad)/pi')
            ylabel('Width')
            if fit_sine_azm
                disp('fit sine to rad width of azm bins')
                %fit_params=fit_sine(xdata,ydata,amp_guess,phase_guess,freq_guess,isverbose)
                azm_rad_fit_params(n,:,:)=fit_azm_sine(azmcenters,rad_fit_azmbin_params(n,:,3,1),range(rad_fit_azmbin_params(n,:,3,1)),mean(rad_fit_azmbin_params(n,:,3,1)),0,1,1);
            end
        end
        
        saveas(gcf,[files.path,'_RadDist_Bin_',num2str(halo_bin_cent(n),'%5.0f'),'_counts.jpg'])
        saveas(gcf,[files.path,'_RadDist _Bin',num2str(halo_bin_cent(n),'%5.0f'),'_counts.fig']) 
    end%loop over halo count bins
    
%     clear halo_centered_cells_count_bined_single
%     clear halo_centered_cells_count_bined_single_comb
%     clear bin_counts
    clear mask_tmp
    
    if size(halo_numcat,2)>2 && fit_rad_dist
        figure(8)
        set(gcf,'Color',[1 1 1]);
        subplot(2,1,1)
        
        errorbar(halo_bin_cent,rad_fit_params(:,3,1),rad_fit_params(:,3,2),'x');
        xlabel('Counts')
        ylabel('Halo Radial Width (m)')
        title('Halo Width')
        subplot(2,1,2)
        fracwidth=rad_fit_params(:,3,1)./rad_fit_params(:,2,1);
        fracunc=fracwidth.*sqrt( (rad_fit_params(:,3,2)./rad_fit_params(:,3,1)).^2 +...
            (rad_fit_params(:,2,2)./rad_fit_params(:,2,1)).^2);
        errorbar(halo_bin_cent,fracwidth,fracunc,'x');
        xlabel('Counts')
        ylabel('Halo Radial Width Fraction')
        
        clear fracunc
    end
    
end

toc
if isverbose
   disp('ALL TASKS COMPLETED');
end

%the halo width as frac of the halo radius
%rad_fit_params(:,3,1)./rad_fit_params(:,2,1)

%the bec width is given by
%bec_bragg(1,1,5) in x
%bec_bragg(1,1,6) in y
%mean(bec_bragg(:,1,5))
%mean(bec_bragg(:,1,6))

%as fraction of halo radius
%mean(bec_bragg(:,1,5))/rad_fit_params(:,2,1)
%mean(bec_bragg(:,1,6))/rad_fit_params(:,2,1)