function [halo_centered_cells,halo_centered,bec_bragg,all_points]=HaloReflecImportData(files,windows,use_txy,use_saved_halos,isverbose,progress_scaling)
%======================================90char=============================================
%+++++++++++++++++++++++++++++++++++++++ABOUT+++++++++++++++++++++++++++++++++++++++++++++
%This function reads in txy data for listed in the LOG_parameters.txt file
%(will convert DLD data if TXY does not exist) checks if there are at least files.count_min
%counts in the file (returns empty otherwise). The bec and bragg is
%windowed and the average is used to find the center of the halo. The halo
%counts have this center subtracted and multiplied by the velocity in the z-axis to
%convert to position. The halo counts are then radially selected to reduce
%the background. These counts and the window/file parameters are then saved
%to $files.path$_all_halos_saved.mat which can just be imported the next time.
%to allow the reflections in I add a second criterion that is a halo
%shifted by ± one radius

% INPUTS
% $file$ structure:
% file.path             % path to unindexed data file
%
% files.numstart        % start file num
% files.numtoimport     % number of files to import
%
% windows.bragg.tmin	% init t of Bragg scattered condensate
% windows.bragg.tmax
%
% windows.bec.tmin      % init t of unscattered BEC
% windows.bec.tmax
%
% windows.halo.tmin     % init t of halo window
% windows.halo.tmax
% windows.halo.rmin     % halo radial mask
% windows.halo.rmax
%
% windows.all.xmin      % common window
% windows.all.xmax
% windows.all.ymin
% windows.all.ymax
% windows.reflections   % allows reflections throguh the mask


% OUTPUTS
% bec_bragg files*2*7 matrix which contains each file contains bragg avg (t,x,y)
% sd(t x y) numcounts as flat vector and then the same for bec


% TO IMPROVE:
% if files =0 import all
% should pass a vector of paths for multi folder import
% should import params file for scans
% should have a xy radial mask

% measure widths and number in BEC, bragg

%================================ START USER INPUT=====================================
min_counts_in_bragg_or_BEC=5;
%===================================END USER INPUT========================================

%% Variable manager
% variables to save in $files.path$_all_halos_saved.mat file
vars_saved = {'halo_centered_cells','windows','files','bec_bragg',...
    'halo_radius','counts_plus','counts_minus','counts_zero',...
    'lowcountfiles','missingfiles','filestotxy','bec_or_bragg_zero','importokfiles'};

%% Use saved settings
if use_saved_halos
    if ~fileExists([files.path,'_all_halos_saved.mat'])
        disp(['WARNING: Halo data "', [files.path,'_all_halos_saved.mat'], '" does not exist - setting use_saved_halos=0 and continue...']);
        use_saved_halos=0;
    else
        % Error checks on prev saved file
        S = load([files.path,'_all_halos_saved.mat']);
        
        is_complete=1;
        for i=1:length(vars_saved)
            is_complete=is_complete&isfield(S,vars_saved{i});
        end
        
        % Completeness of information
        if ~is_complete
            disp('WARNING: Previously saved file is incomplete - setting use_saved_halos=0 and continue...');
            use_saved_halos=0;
            
            % Check consistency of requested analysis settings and saved data
        elseif ~isequal(windows,S.windows) || ~isequal(files,S.files)
            disp('WARNING: Previously saved data contains settings (windows, files) different to requested - setting use_saved_halos=0 and continue...');
            use_saved_halos=0;
            
            % OK load all vars in file
        else
            load([files.path,'_all_halos_saved.mat']);
        end
        
        clear S;
    end
end

%% Import data files and do pre-processing
if ~use_saved_halos
    %initalize variables
    halo_centered_cells={};
    
    %initialize the parfor compatible counters
    lowcountfiles=zeros(files.numtoimport,1);
    missingfiles=zeros(files.numtoimport,1);
    filestotxy=zeros(files.numtoimport,1);
    bec_or_bragg_zero=zeros(files.numtoimport,1);
    importokfiles=zeros(files.numtoimport,1);
    counts_zero=zeros(files.numtoimport,1);
    counts_plus=zeros(files.numtoimport,1);
    counts_minus=zeros(files.numtoimport,1);
    halo_radius=zeros(files.numtoimport,1);
    
    if isverbose
        disp('Importing data...');
        parfor_progress(round(files.numtoimport/progress_scaling));
    end
    
    bec_bragg=zeros(files.numtoimport,2,7);
    
    % use parfor below? - DKS
    for n=1:files.numtoimport
        current_file_str = num2str(files.numstart+n-1);
        
        % BELOW IS BUGGY for use_txy switch
        %TODO
        %check if any (raw or txy) file exists
        if ~fileExists([files.path,current_file_str,'.txt'])&&~fileExists([files.path,'_txy_forc',current_file_str,'.txt'])
            missingfiles(n)=1;  % a missing file
            if isverbose
                disp(['WARNING: missing file - file #', curr_file_str]);
            end
            
        % user reqs to use_txy but txy data doesn't exist
        elseif use_txy && ~fileExists([files.path,'_txy_forc',current_file_str,'.txt'])
            missingfiles(n)=1;  % treat as missing file
            if verbose
                disp(['WARNING: use_txy=1 but TXY-data does not exist - file #', curr_file_str]);
            end
            
        % user reqs to build txy data but raw doesn't exist
        elseif use_txy && ~fileExists([files.path,current_file_str,'.txt'])
            missingfiles(n)=1;  %treat as missing file
            if verbose
                disp(['WARNING: use_txy=0 but raw DLD data does not exist - file #', curr_file_str]);
            end
            
        % Import from file
        else                         
            %three_channel_output=importdata([files.path,'_txy_forc',current_file_str,'.txt'],',');
            %three_channel_output=load('-ascii',[files.path,'_txy_forc',current_file_str,'.txt']);
            three_channel_output=txy_importer(files.path,current_file_str);
            
            if size(three_channel_output,1)<files.count_min
                %disp('selected shot has too few counts')
                lowcountfiles(n)=1;
                
            else
                %                 if sum(sum(isnan(three_channel_output)))>0
                %                     disp('NAN values found')
                %                     return
                %                 end
                
                three_channel_output_rot=zeros(size(three_channel_output));
                %rotate the spatial cord
                sin_theta = sin(files.rot_angle);
                cos_theta = cos(files.rot_angle);
                three_channel_output_rot(:,1) = three_channel_output(:,1);
                three_channel_output_rot(:,2) = three_channel_output(:,2)*cos_theta...
                    - three_channel_output(:,3)*sin_theta;
                three_channel_output_rot(:,3) = three_channel_output(:,2)*sin_theta...
                    + three_channel_output(:,3)*cos_theta;
                
                %which counts in the file are within the xy limits
                mask_XY=three_channel_output_rot(:,2)>windows.all.xmin &...
                    three_channel_output_rot(:,2)<windows.all.xmax & ...
                    three_channel_output_rot(:,3)>windows.all.ymin & ...
                    three_channel_output_rot(:,3)<windows.all.ymax;
                
                %save all points is mainly just usefull for a diagnostic
                if files.save_all_points
                    all_points{n}=three_channel_output_rot(mask_XY,:);
                end
                
                %windows in time for bragg and BEC
                %TO IMPROVE: little point finding whats mask and bragg for
                %pos correction off,but this is a rare use case..
                mask_BEC=three_channel_output_rot(:,1)>windows.bec.tmin &...
                    three_channel_output_rot(:,1)<windows.bec.tmax &...
                    mask_XY;
                
                mask_bragg=three_channel_output_rot(:,1)>windows.bragg.tmin &...
                    three_channel_output_rot(:,1)<windows.bragg.tmax &...
                    mask_XY;
                %if the bragg or bec is empty give up and set
                %halo_centered_cells_temp to empty
                if  (sum(mask_bragg)<min_counts_in_bragg_or_BEC || sum(mask_BEC)<min_counts_in_bragg_or_BEC) ...
                        && files.do_pos_correction
                    bec_or_bragg_zero(n)=1;
                    %write a halo file to prevent it getting here again
                    halo_centered_temp=[];
                else
                    three_channel_output_rot_bragg=three_channel_output_rot(mask_bragg,:);
                    three_channel_output_rot_BEC=three_channel_output_rot(mask_BEC,:);
                    %find the center of the bragg and bec and then subtract
                    %this position off the couns in the halo
                    if files.do_pos_correction
                        mean_bec=mean(three_channel_output_rot_BEC);
                        mean_bragg=mean(three_channel_output_rot_bragg);
                        bec_bragg(n,:,:)=[[mean_bec,std(three_channel_output_rot_BEC),sum(mask_BEC)];[mean_bragg,std(three_channel_output_rot_bragg),sum(mask_bragg)]];
                        halo_center=(mean_bragg+mean_bec)/2;
                        halo_radius(n)=norm((mean_bragg-mean_bec).*[files.velocity 1 1])/2;
                        
                        bec_centered=mean_bec-halo_center;
                        bragg_centered=mean_bragg-halo_center;
                    else
                        %same as not subtracting pos
                        halo_center=[0 0 0];
                    end
                    
                    %mask the halo
                    mask_halo=three_channel_output_rot(:,1)>windows.halo.tmin &...
                        three_channel_output_rot(:,1)<windows.halo.tmax &...
                        mask_XY;
                    three_channel_output_rot_halo=three_channel_output_rot(mask_halo,:);
                    
                    %subtract the t,x,y center from every row
                    %repmat is apparently the fastest http://au.mathworks.com/matlabcentral/newsreader/view_thread/267082
                    %then convert the time axis into z pos by multiplying time by
                    %files.velocity
                    halo_centered_temp_zero=(three_channel_output_rot_halo-repmat(halo_center,size(three_channel_output_rot_halo,1),1)).*...
                        repmat([files.velocity 1 1],size(three_channel_output_rot_halo,1),1);
                    
                    %here i mask radialy by computing the radius and seeing if it
                    %is within rmax,windows.halo.rmin
                    %will not do if files.do_pos_correction=0 as the halo is not zero
                    %centered
                    importokfiles(n)=1;
                    if files.do_pos_correction
                        %radius for the ceneter halo
                        radial_temp_zero=sqrt(sum(abs(halo_centered_temp_zero).^2,2));
                        %radial_temp_zero=sqrt((:,1).^2+halo_centered_temp_zero(:,2).^2+halo_centered_temp_zero(:,3).^2);
                        %radius for the halo reflections
                        radial_mask=radial_temp_zero<windows.halo.rmax & radial_temp_zero>windows.halo.rmin;
                        counts_zero(n)=sum(radial_mask);
                        
                        %if reflections are allowed then move all the
                        %points, cacluate the radius and then add points
                        %within the radius to the radial mask
                        if windows.reflections
                            %this creates two shifted version of the counts
                            %does it by scratch (not using above computer centered)
                            %which is easy but perhaps a little wastefull
                            
                            halo_centered_temp_plus=(three_channel_output_rot_halo-repmat(mean_bragg,size(three_channel_output_rot_halo,1),1)).*...
                                repmat([files.velocity 1 1],size(three_channel_output_rot_halo,1),1);
                            halo_centered_temp_minus=(three_channel_output_rot_halo-repmat(mean_bec,size(three_channel_output_rot_halo,1),1)).*...
                                repmat([files.velocity 1 1],size(three_channel_output_rot_halo,1),1);
                            
                            radial_temp_plus=sqrt(sum(abs(halo_centered_temp_plus).^2,2));
                            radial_temp_minus=sqrt(sum(abs(halo_centered_temp_minus).^2,2));
                            radial_mask_plus=radial_temp_plus<windows.halo.rmax & radial_temp_plus>windows.halo.rmin;
                            counts_plus(n)=sum(radial_mask_plus);
                            radial_mask_minus=radial_temp_minus<windows.halo.rmax & radial_temp_minus>windows.halo.rmin;
                            counts_minus(n)=sum(radial_mask_minus);
                            reflec_rmin_mask=radial_temp_zero>windows.halo.reflecrmin;
                            radial_mask=(radial_mask | radial_mask_plus | radial_mask_minus)& reflec_rmin_mask;
                            %radial_mask=radial_mask_plus | radial_mask_minus;
                        end
                        
                        halo_centered_temp=halo_centered_temp_zero(radial_mask,:);
                    end
                end
                
                %save the halo counts
                %use cells so that the output can be ragged
                halo_centered_cells{n}=halo_centered_temp;
                
                %clean up variables
                
            end %file size condition
        end %file exists condition
        
        if isverbose
            % update parfor only occasionally
            if rand(1)<(1/progress_scaling)
                parfor_progress;
            end
        end
    end %file loop
    
    %save the imported data
    for i = 1:length(vars_saved)
        save([files.path,'_all_halos_saved.mat'],vars_saved{i},'-append');
    end
    
    %clean up, this is only important if the loop is serial not parfor
    clear('mask_BEC','mask_bragg','mask_XY','mask_halo','mean_bec','mean_bragg')
    clear('three_channel_output','three_channel_output_rot_BEC','three_channel_output_rot_bragg')
    clear('three_channel_output_rot_halo','halo_centered_temp','radial_temp','three_channel_output_rot')
    
    if isverbose
        parfor_progress(0);
        disp('Import Done')
    end
end

if ~files.save_all_points
    all_points=[];
end

%% Summary
if isverbose
    if use_saved_halos
        disp(['Data loaded from: ',[files.path,'_all_halos_saved.mat']]);
    else
        % New file created
        disp(['Data imported and saved: ',[files.path,'_all_halos_saved.mat']]);
    end
    
    disp('=======HaloReflecImportData========');
    disp([int2str(sum(lowcountfiles)),' files with low counts'])
    disp([int2str(sum(missingfiles)),' files missing'])
    disp([int2str(sum(filestotxy)),' files converted to txy'])
    disp([int2str(sum(bec_or_bragg_zero)),' files with no counts in bec or bragg'])
    disp([int2str(sum(importokfiles)),' imported ok'])
    
    halo_centered=vertcat(halo_centered_cells{:});  % var for summary
    
    disp(['avergage halo radius: ',num2str(mean(halo_radius))]);
    disp([int2str(sum(counts_zero)),',',int2str(sum(counts_plus)),',',int2str(sum(counts_minus)),' counts in zero,plus,minus'])
    disp([int2str(size(halo_centered,1)),' total counts in halo'])
    disp([num2str(size(halo_centered,1)/(files.numtoimport-(sum(missingfiles)+sum(bec_or_bragg_zero)))),' counts per file in halo'])
    disp(['mean x,y,t pos halo ',num2str(mean(halo_centered(:,2))),' , ', ...
        num2str(mean(halo_centered(:,3))),' , ', num2str(mean(halo_centered(:,1)))]);
    disp('===================================');
end

end