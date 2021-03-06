function [G2_corr,G2_uncorr]=G2_cart(DATA,BIN_EDGE,CORR_INFO,VERBOSE)
% Unnormalised G2 calculator for more general applicaton. e.g. cross-halo
% DKS 06/11/2016
%
% INPUT
% DATA - N(shots)xM(1 or 2) cell of counts in ZXY format
% CORR_INFO - 'CL': collinear
%             'BB': back-to-back
%             (default is BB)
% BIN_EDGE - 1x3 cell, bin-edges deltaK in ZXY axis
%
% OUTPUT
% G2_corr - G2 for correlated pairs
% G2_uncorr - G2 for all uncorrelated pairs
%

if VERBOSE>0
    t_fun_start=tic;
end

% Input check
% DATA
if ~iscell(DATA)
    error('DATA should be a cell of counts.');
end
% BIN_EDGE
if ~iscell(BIN_EDGE)
    error('BIN_EDGE should be a cell of bin edges');
end
if ~isequal(size(BIN_EDGE),[1,3])
    error('BIN_EDGE should be a 1x3 cell of bin edges in Z,X,Y axis.');
end

% CORR_INFO
if ~exist('CORR_INFO','var')  % missing arg
    error('CORR_INFO must be "BB" or "CL". Setting to default "BB"');
elseif ~(isequal(CORR_INFO,'BB')||isequal(CORR_INFO,'CL'))  % invalid input
    error('CORR_INFO must be "BB" or "CL". Setting to default "BB"');
end

% VERBOSE
if ~exist('VERBOSE','var')
    VERBOSE=0;  % default is quiet
end

% Parse inputs
if size(DATA,2)==2
    % cross-data G2
    data1=DATA(:,1);
    data2=DATA(:,2);
    ncomp=2;    % number of components
elseif size(DATA,2)==1
    % single-component G2
    data1=DATA;
    ncomp=1;
else
    error('DATA must be a 1 or 2 column cell');
end
clear DATA;     % clean workspace

nBin=zeros(1,length(BIN_EDGE));
for i=1:length(BIN_EDGE)
    nBin(i)=length(BIN_EDGE{i})-1;  % number of bins from bin edges
end

% Initialise variables
nShot=size(data1,1);     % number of shots in data
if VERBOSE>0
    disp('----------------------------------------------');
    disp([num2str(nShot),' shots to analyse for G2_cart...']);
end
G2_corr=zeros(nBin);
G2_uncorr=zeros(nBin);

% Branch G2 analysis so that condition is out of loop: just add other
%   conditions as a conditional branch following one as template
erase_str='';   % initialise for progress reporting
if ncomp==2
    %% 2 component G(2) analysis
    if isequal(CORR_INFO,'BB')
        %% Back-to-back G2 analysis - x-state
        % correlated pairs G2
%         for i=1:nShot
        parfor i=1:nShot
            nAtom=size(data1{i},1);         % number of counts in shot - comp #1
            Npairs=size(data2{i},1);        % number of counts in shot - comp #2
            
            for j=1:nAtom
                % back-to-back condition
                this_atom=data1{i}(j,:);    % ZXY-vector for this atom (to find pairs)
                diff_tmp=data2{i}+repmat(this_atom,[Npairs,1]);   % sum for diff_BB in k-space
                
                G2_corr=G2_corr+nhist(diff_tmp,BIN_EDGE);	% update G2
            end
            
%             % Progress report
%             if VERBOSE>1
%                 prog_msg=sprintf('[1/2]: %4.1f',100*i/nShot);
%                 fprintf([erase_str,prog_msg]);
%                 erase_str=repmat(sprintf('\b'),1,length(prog_msg));
%             end
        end
        
%         if VERBOSE>1
%             fprintf('\n');
%             erase_str='';   % reset to null
%         end
        
        % UNcorrelated pairs G2
%         for i=1:nShot
        parfor i=1:nShot
            data_collated=vertcat(data2{[1:i-1,i+1:end]});  % all shots except self
            Ntotpair=size(data_collated,1);  % total number of counts in the cross-species
            nAtom=size(data1{i},1);
            
            for j=1:nAtom
                % back-to-back condition
                this_atom=data1{i}(j,:);
                diff_tmp=data_collated+repmat(this_atom,[Ntotpair,1]);   % sum for BB
                
                G2_uncorr=G2_uncorr+nhist(diff_tmp,BIN_EDGE);     % update G2
            end
            
%             % Progress report
%             if VERBOSE>1
%                 prog_msg=sprintf('[2/2]: %4.1f',100*i/nShot);
%                 fprintf([erase_str,prog_msg]);
%                 erase_str=repmat(sprintf('\b'),1,length(prog_msg));
%             end
        end
        
    elseif isequal(CORR_INFO,'CL')
        %% Colinear G2 analysis - x-state
        % correlated pairs G2
%         for i=1:nShot
        parfor i=1:nShot
            nAtom=size(data1{i},1); % number of counts in DATA1
            Npairs=size(data2{i},1);
            
            for j=1:nAtom
                % colinear condition
                this_atom=data1{i}(j,:);    % ZXY-vector for this atom (to find pairs)
                diff_tmp=abs(data2{i}-repmat(this_atom,[Npairs,1]));   % "absolute valued-DIFF" for diff_CL in k-space
                
                G2_corr=G2_corr+nhist(diff_tmp,BIN_EDGE);	% update G2
            end
            
%             % Progress report
%             if VERBOSE>1
%                 prog_msg=sprintf('[1/2]: %4.1f',100*i/nShot);
%                 fprintf([erase_str,prog_msg]);
%                 erase_str=repmat(sprintf('\b'),1,length(prog_msg));
%             end
        end
        
%         if VERBOSE>1
%             fprintf('\n');  % need to return carriage in output
%             erase_str='';   % reset to null
%         end
        
        % UNcorrelated pairs G2
%         for i=1:nShot
        parfor i=1:nShot
            data_collated=vertcat(data2{[1:i-1,i+1:end]});  % all shots except self
            Ntotpair=size(data_collated,1);  % total number of counts in the cross-species
            nAtom=size(data1{i},1);
            
            for j=1:nAtom
                % colinear condition
                this_atom=data1{i}(j,:);
                diff_tmp=abs(data_collated-repmat(this_atom,[Ntotpair,1]));   % abs-diff for CL
                
                G2_uncorr=G2_uncorr+nhist(diff_tmp,BIN_EDGE);     % update G2
            end
            
%             % Progress report
%             if VERBOSE>1
%                 prog_msg=sprintf('[2/2]: %4.1f',100*i/nShot);
%                 fprintf([erase_str,prog_msg]);
%                 erase_str=repmat(sprintf('\b'),1,length(prog_msg));
%             end
        end
    else
        error('BUG: CORR_INFO must be BB or CL at this point: this line should never be called.');
    end
    
    %% normalise G2 function
    % get number of detected atoms in each shot/species
    nn1=cell2mat(cellfun(@(C) size(C,1),data1,'UniformOutput',false));      % counts in shot (species 1)
    nn2=cell2mat(cellfun(@(C) size(C,1),data2,'UniformOutput',false));      % counts in shot (species 2)
    
    % X-state unique/unordered pair counting
    N_pair_corr=sum(nn1.*nn2);
    N_pair_uncorr=sum(nn1)*sum(nn2)-sum(nn1.*nn2);
    
    % normalise with respect to all possible pairs
    G2_corr=G2_corr/N_pair_corr;
    G2_uncorr=G2_uncorr/N_pair_uncorr;
    
%     if VERBOSE>1
%         % new-line at the end of report
%         fprintf('\n');
%     end
    
elseif ncomp==1
    % TODO - BUG with CL because counts are ordered in Z (in a shot, CL-diff is always positive)
    %% Single component G(2) analysis
    if isequal(CORR_INFO,'BB')
        %% Back-to-back G2 analysis - in-halo
        % correlated pairs G2
%         for i=1:nShot
        parfor i=1:nShot
            shot_tmp=data1{i};
            nAtom=size(shot_tmp,1); % number of counts in this shot
            
            for j=1:(nAtom-1)   % skip last atom - no unique pairs left
                % back-to-back condition
                this_atom=shot_tmp(j,:);    % ZXY-vector for this atom (to find pairs)
                % Get BB-diff vectors for all pairs except self
                diff_tmp=shot_tmp((j+1):end,:)+repmat(this_atom,[nAtom-j,1]);  % BB-diff for unique pairs
                
                G2_corr=G2_corr+nhist(diff_tmp,BIN_EDGE);   % update G2
            end
            
%             % Progress report
%             if VERBOSE>1
%                 prog_msg=sprintf('[1/2]: %4.1f',100*i/nShot);
%                 fprintf([erase_str,prog_msg]);
%                 erase_str=repmat(sprintf('\b'),1,length(prog_msg));
%             end
        end
%         
%         if VERBOSE>1
%             fprintf('\n');
%             erase_str='';   % reset to null
%         end
        
        % UNcorrelated pairs G2
%         for i=1:(nShot-1)     % skip last shot - no unique pairs left
        parfor i=1:(nShot-1)     % skip last shot - no unique pairs left
            data_collated=vertcat(data1{(i+1):end});  % all shots except self for unique unordered pairing
            Ntotpair=size(data_collated,1);     % total number of counts to search pairs
            nAtom=size(data1{i},1);     % counts in this shot
            
            for j=1:nAtom
                % back-to-back condition
                this_atom=data1{i}(j,:);
                diff_tmp=data_collated+repmat(this_atom,[Ntotpair,1]);   % BB-diff
                
                G2_uncorr=G2_uncorr+nhist(diff_tmp,BIN_EDGE); % update G2
            end
            
%             % Progress report
%             if VERBOSE>1
%                 prog_msg=sprintf('[2/2]: %4.1f',100*i/nShot);
%                 fprintf([erase_str,prog_msg]);
%                 erase_str=repmat(sprintf('\b'),1,length(prog_msg));
%             end
        end
    elseif isequal(CORR_INFO,'CL')
        %% Colinear G2 analysis - in-halo
        % correlated pairs G2
        parfor i=1:nShot
%             for i=1:nShot
            shot_tmp=data1{i};
            nAtom=size(shot_tmp,1); % number of counts in this shot
            
            for j=1:(nAtom-1) % skip last atom - no unique pairs left
                % colinear condition
                this_atom=shot_tmp(j,:);    % ZXY-vector for this atom (to find pairs)
                % Get "CL"-diff vectors for all pairs except self
                diff_tmp=abs(shot_tmp((j+1):end,:)-repmat(this_atom,[nAtom-j,1]));  % remove ordering by subtraction with abs
                G2_corr=G2_corr+nhist(diff_tmp,BIN_EDGE);   % update G2
            end
            
%             % Progress report
%             if VERBOSE>1
%                 prog_msg=sprintf('[1/2]: %4.1f',100*i/nShot);
%                 fprintf([erase_str,prog_msg]);
%                 erase_str=repmat(sprintf('\b'),1,length(prog_msg));
%             end
        end
        
%         if VERBOSE>1
%             fprintf('\n');
%             erase_str='';   % reset to null
%         end
        
        % UNcorrelated pairs G2
        parfor i=1:(nShot-1)  % skip last shot - no unique pairs left
%             for i=1:(nShot-1)  % skip last shot - no unique pairs left
            data_collated=vertcat(data1{(i+1):end});  % all shots except self for unique unordered pairing
            Ntotpair=size(data_collated,1);     % total number of counts to search pairs
            nAtom=size(data1{i},1);     % counts in this shot
            
            for j=1:nAtom
                % colinear condition
                this_atom=data1{i}(j,:);
                diff_tmp=abs(data_collated-repmat(this_atom,[Ntotpair,1]));
                
                G2_uncorr=G2_uncorr+nhist(diff_tmp,BIN_EDGE); % update G2
            end
            
%             % Progress report
%             if VERBOSE>1
%                 prog_msg=sprintf('[2/2]: %4.1f',100*i/nShot);
%                 fprintf([erase_str,prog_msg]);
%                 erase_str=repmat(sprintf('\b'),1,length(prog_msg));
%             end
        end
    else
        error('BUG: CORR_INFO must be BB or CL at this point: this line should never be called.');
    end
    
	%% normalise G2 function
    % get number of detected atoms in each shot/species
    nn1=cell2mat(cellfun(@(C) size(C,1),data1,'UniformOutput',false));      % counts in shot
    
    % single-species unique/unordered pair counting
    N_pair_corr=sum(nn1.*(nn1-1)/2);
    N_pair_uncorr=(sum(nn1)^2-sum(nn1.^2))/2;
    
    % normalise with respect to all possible pairs
    G2_corr=G2_corr/N_pair_corr;
    G2_uncorr=G2_uncorr/N_pair_uncorr;
    
%     if VERBOSE>1
%         % new-line at the end of report
%         fprintf('\n');
%     end
end

if VERBOSE>0
    t_fun_end=toc(t_fun_start);
    fprintf('Total elapsed time (s): %7.2f\n',t_fun_end);
    disp('----------------------------------------------');
end