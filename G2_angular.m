function [G2_corr,G2_uncorr]=G2_angular(DATA,BIN_EDGE,VERBOSE)
% Unnormalised G2 calculator in 2D radial-angular separation with cartesian
%   count format
% DKS 14/11/16
%
% INPUT
% DATA - N (shot) x M(1 or 2) cell of counts in ZXY format
% BIN_EDGE - 1x2 cell of bin edges in radial and angular separation
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
if ~isequal(size(BIN_EDGE),[1,2])
    error('BIN_EDGE should be a 1x2 cell of bin edges in radial and angular separation.');
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
clear DATA;

nBin=zeros(1,length(BIN_EDGE));
for i=1:length(BIN_EDGE)
    nBin(i)=length(BIN_EDGE{i})-1;  % number of bins from bin edges
end

% Initialise variables
nShot=size(data1,1);     % number of shots in data
if VERBOSE>0
    disp('----------------------------------------------');
    disp([num2str(nShot),' shots to analyse in G2_angular...']);
end
G2_corr=zeros(nBin);
G2_uncorr=zeros(nBin);


% Branch G2 analysis so that condition is out of loop: just add other
%   conditions as a conditional branch following one as template
erase_str='';   % initialise for progress reporting
if ncomp==2
    %% 2 component G(2) analysis
    % correlated pairs G2   
    for i=1:nShot
        nAtom=size(data1{i},1);     % number of counts in DATA1
        Npairs=size(data2{i},1);    % number of 'pairable' counts in counterpart
        diff_tmp=[];   % diff vectors for pair search (dradius,dangle)
        
        r_shot=sqrt(sum(data1{i}.^2,2));    % radial dist for counts in halo1
        r_pairs=sqrt(sum(data2{i}.^2,2));   % radial dist for counts in halo2 (pair)
        
        for j=1:nAtom
            zxy_this=data1{i}(j,:);    % ZXY-vector for this atom (to find pairs)
            r_this=r_shot(j);   % norm of this atom
%             r_pairs=sqrt(sum(data2{i}.^2,2));  % norm of counterpart atoms (same shot)
            
            diff_tmp(:,1)=abs(r_pairs-r_this);   % radial diff - remove ordering with abs
            dotp_tmp=sum(data2{i}.*repmat(zxy_this,[Npairs,1]),2); % dot product
            diff_tmp(:,2)=real(acos(dotp_tmp./(r_this*r_pairs)));    % angular diff
            
            G2_corr=G2_corr+nhist(diff_tmp,BIN_EDGE);	% update G2
        end
        
        % Progress report
        if VERBOSE>1
            prog_msg=sprintf('[1/2]: %4.1f',100*i/nShot);
            fprintf([erase_str,prog_msg]);
            erase_str=repmat(sprintf('\b'),1,length(prog_msg));
        end
    end
    
    if VERBOSE>1
        fprintf('\n');
        erase_str='';   % reset to null
    end
    
    % UNcorrelated pairs G2
    for i=1:nShot
        data_collated=vertcat(data2{[1:i-1,i+1:end]});  % except self
        Ntotpair=size(data_collated,1);  % total number of counts in the cross-species
        nAtom=size(data1{i},1);
        diff_tmp=[];
        
        r_shot=sqrt(sum(data1{i}.^2,2));        % radial dist for counts in halo1
        r_pairs=sqrt(sum(data_collated.^2,2));  % radial dist for all counts - uncorr pairs
        
        for j=1:nAtom
            zxy_this=data1{i}(j,:);
            r_this=r_shot(j);   % norm of this atom
%             r_pairs=sqrt(sum(data_collated.^2,2));  % norm of counterpart atoms (same shot)
            
            diff_tmp(:,1)=abs(r_pairs-r_this);   % radial diff - remove ordering with abs
            dotp_tmp=sum(data_collated.*repmat(zxy_this,[Ntotpair,1]),2); % dot product           
            diff_tmp(:,2)=real(acos(dotp_tmp./(r_this*r_pairs)));    % angular diff
            
            G2_uncorr=G2_uncorr+nhist(diff_tmp,BIN_EDGE);     % update G2
        end
        
        % Progress report
        if VERBOSE>1
            prog_msg=sprintf('[2/2]: %4.1f',100*i/nShot);
            fprintf([erase_str,prog_msg]);
            erase_str=repmat(sprintf('\b'),1,length(prog_msg));
        end
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
    
    if VERBOSE>1
        % new-line at the end of report
        fprintf('\n');
    end
    
elseif ncomp==1
    %% Single component G(2) analysis
    % correlated pairs G2
    for i=1:nShot
        shot_tmp=data1{i};
        nAtom=size(shot_tmp,1); % number of counts in this shot
        
        r_all=sqrt(sum(shot_tmp.^2,2));    % radial dist for all counts in this shot
        
        for j=1:(nAtom-1)   % skip last atom - no unique pairs
            diff_tmp=[];   % initialise diff vectors for pair search 
            
            zxy_this=shot_tmp(j,:);    % ZXY-vector for this atom (to find pairs)
            zxy_pairs=shot_tmp((j+1):end,:);    % rem count available for unique pairing
            
            % Get rad/angular diff vectors for all pairs except self
            % norms of zxy vector
%             r_this=sqrt(sum(zxy_this.^2,2));
%             r_pairs=sqrt(sum(zxy_pairs.^2,2));
            r_this=r_all(j);
            r_pairs=r_all(j+1:end);
            
            % radial diff - remove ordering with abs
            diff_tmp(:,1)=abs(r_pairs-r_this);      
            
            % calculate angular difference
            dotp_tmp=sum(zxy_pairs.*repmat(zxy_this,[nAtom-j,1]),2);    % vector inner-product
            diff_tmp(:,2)=real(acos(dotp_tmp./(r_this*r_pairs)));

            G2_corr=G2_corr+nhist(diff_tmp,BIN_EDGE);   % update G2
        end
        
        % Progress report
        if VERBOSE>1
            prog_msg=sprintf('[1/2]: %4.1f',100*i/nShot);
            fprintf([erase_str,prog_msg]);
            erase_str=repmat(sprintf('\b'),1,length(prog_msg));
        end
    end
    
    if VERBOSE>1
        fprintf('\n');
        erase_str='';   % reset to null
    end
    
    % UNcorrelated pairs G2
    for i=1:(nShot-1)   % skip last shot - no unique pairs
        data_collated=vertcat(data1{(i+1):end});  % all shots except self
        Ntotpair=size(data_collated,1);     % total number of counts to search pairs
        nAtom=size(data1{i},1);     % counts in this shot
        
        r_shot=sqrt(sum(data1{i}.^2,2));            % radial dist for counts in this shot
        r_pairs=sqrt(sum(data_collated.^2,2));      % radial dist for all counts - uncorr pairs
        
        for j=1:nAtom
            diff_tmp=[];   % initialise diff vectors
            
            zxy_this=data1{i}(j,:);
            
            % Get rad/angular diff vectors for all pairs
%             r_this=sqrt(sum(zxy_this.^2,2));
%             r_pairs=sqrt(sum(data_collated.^2,2));
            r_this=r_shot(j);
            
            diff_tmp(:,1)=abs(r_pairs-r_this);
            dotp_tmp=sum(data_collated.*repmat(zxy_this,[Ntotpair,1]),2);
            diff_tmp(:,2)=real(acos(dotp_tmp./(r_this*r_pairs)));
            
            G2_uncorr=G2_uncorr+nhist(diff_tmp,BIN_EDGE); % update G2
        end
        
        % Progress report
        if VERBOSE>1
            prog_msg=sprintf('[2/2]: %4.1f',100*i/nShot);
            fprintf([erase_str,prog_msg]);
            erase_str=repmat(sprintf('\b'),1,length(prog_msg));
        end
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
    
    if VERBOSE>1
        % new-line at the end of report
        fprintf('\n');
    end
    
end

if VERBOSE>0
    t_fun_end=toc(t_fun_start);
    fprintf('Total elapsed time (s): %7.2f\n',t_fun_end);
    disp('----------------------------------------------');
end