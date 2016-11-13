function [G2_SINGLE,G2_ALL]=G2_angular(DATA,BIN_EDGE,VERBOSE)
% Unnormalised G2 calculator in 2D radial-angular separation
% DKS 14/11/16
%
% INPUT
% DATA - N (shot) x M(1 or 2) cell of counts in ZXY format
% BIN_EDGE - 1x2 cell of bin edges in radial and angular separation
%
% OUTPUT
% G2_SINGLE - G2 summed for each shot
% G2_ALL - G2 for all data collated (for normalisation)

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
    error('BIN_EDGE should be a 1x2 cell of radial and angular separation bin edges.');
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
if VERBOSE>1,disp([num2str(nShot),' shots to analyse in G2_angular...']);,end;
G2_SINGLE=zeros(nBin);
G2_ALL=zeros(nBin);


% Branch G2 analysis so that condition is out of loop: just add other
%   conditions as a conditional branch following one as template
if ncomp==2
    %% 2 component G(2) analysis
    % Single-shot correlations
    for i=1:nShot
        nAtom=size(data1{i},1);     % number of counts in DATA1
        Npairs=size(data2{i},1);    % number of 'pairable' counts in counterpart
        diff_tmp=[];   % diff vectors for pair search (dradius,dangle)
        
        for j=1:nAtom
            this_atom=data1{i}(j,:);    % ZXY-vector for this atom (to find pairs)
            norm_this_atom=sqrt(sum(this_atom.^2,2));   % norm of this atom
            norm_tmp=sqrt(sum(data2{i}.^2,2));  % norm of counterpart atoms (same shot)
            
            diff_tmp(:,1)=norm_tmp-norm_this_atom;   % radial diff
            dotp_tmp=sum(norm_tmp.*repmat(this_atom,[Npairs,1]),2); % dot product
            diff_tmp(:,2)=acos(dotp_tmp./(norm_this_atom*norm_tmp));    % angular diff
            
            G2_SINGLE=G2_SINGLE+nhist(diff_tmp,BIN_EDGE);	% update G2
        end
    end
    
    % all shots - except self
    for i=1:nShot
        data_collated=vertcat(data2{[1:i-1,i+1:end]});  % except self
        %data_collated=vertcat(data2{:}); % collate all shots inc. self
        Ntotpair=size(data_collated,1);  % total number of counts in the cross-species
        nAtom=size(data1{i},1);
        diff_tmp=[];
        
        for j=1:nAtom
            this_atom=data1{i}(j,:);
            norm_this_atom=sqrt(sum(this_atom.^2,2));   % norm of this atom
            norm_tmp=sqrt(sum(data_collated.^2,2));  % norm of counterpart atoms (same shot)
            
            diff_tmp(:,1)=norm_tmp-norm_this_atom;   % radial diff
            dotp_tmp=sum(norm_tmp.*repmat(this_atom,[Ntotpair,1]),2); % dot product
            diff_tmp(:,2)=acos(dotp_tmp./(norm_this_atom*norm_tmp));    % angular diff
            
            G2_ALL=G2_ALL+nhist(diff_tmp,BIN_EDGE);     % update G2
        end
    end
elseif ncomp==1
    %TODO
end