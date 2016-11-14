function [G2_SINGLE,G2_ALL]=G2_polar(DATA,BIN_EDGE,CORR_INFO,VERBOSE)
% Unnormalised G2 calculator for more general applicaton. e.g. cross-halo
% Retrieves G2 in polar system: diff_norm, diff_angle
% DKS 06/11/2016
%
% INPUT
% DATA - N(shots)xM(1 or 2) cell of counts in R,azim,elev (polar coord)
% CORR_INFO - 'CL': collinear
%             'BB': back-to-back
%             (default is BB)
% BIN_EDGE - 1x2 cell, bin-edges deltaK in (dr,dtheta) axis
%
% OUTPUT
% G2_SINGLE - G2 summed for each shot
% G2_ALL - G2 for all data collated (for normalisation)

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

% CORR_INFO
if ~exist('CORR_INFO','var')  % missing arg - default
    warning('CORR_INFO is not set by user: default to back-to-back (BB) G(2) analysis');
    CORR_INFO='BB';
elseif ~(isequal(CORR_INFO,'BB')||isequal(CORR_INFO,'CL'))  % invalid input
    warning('CORR_INFO must be "BB" or "CL". Setting to default "BB"');
    CORR_INFO='BB';
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
if VERBOSE>1,disp([num2str(nShot),' shots to analyse for G2_polar...']);,end;
G2_SINGLE=zeros(nBin);
G2_ALL=zeros(nBin);

% Branch G2 analysis so that condition is out of loop: just add other
%   conditions as a conditional branch following one as template
if ncomp==2
    %% 2 component G(2) analysis
    if isequal(CORR_INFO,'BB')
        % Back-to-back G2 analysis
        for i=1:nShot
            nAtom=size(data1{i},1); % number of counts in DATA1
            diff_tmp=[];   % diff vectors for pair search
            for j=1:nAtom
                % back-to-back condition in k-theta is dk=0,dtheta=pi
                this_atom=data1{i}(j,:);    % polar-vector for this atom (to find pairs)
                
                diff_tmp(:,1)=data2{i}(:,1)-this_atom(1);   % diff in norm=k2-k1
                diff_tmp(:,2)=real(acos(cos(this_atom(3)).*cos(data2{i}(:,3)).*cos(data2{i}(:,2)-this_atom(2)) ...
                    + sin(this_atom(3)).*sin(data2{i}(:,3))));  % diff angle (see dot-product in polar coord)
                
                G2_SINGLE=G2_SINGLE+nhist(diff_tmp,BIN_EDGE);	% update G2
            end
        end
        
        % all shots - except self
        for i=1:nShot
            data_collated=vertcat(data2{[1:i-1,i+1:end]}); % collate all except self
            %data_collated=vertcat(data2{:}); % collate all shots inc. self
            nAtom=size(data1{i},1);
            diff_tmp=[];
            
            for j=1:nAtom
                this_atom=data1{i}(j,:);
                
                diff_tmp(:,1)=data_collated(:,1)-this_atom(1);
                diff_tmp(:,2)=real(acos(cos(this_atom(3)).*cos(data_collated(:,3)).*cos(data_collated(:,2)-this_atom(2)) ...
                    + sin(this_atom(3)).*sin(data_collated(:,3))));  % diff angle (see dot-product in polar coord)
                
                G2_ALL=G2_ALL+nhist(diff_tmp,BIN_EDGE);        % update G2
            end
        end
    elseif isequal(CORR_INFO,'CL')
        error('CL is not set up yet');
    else
        error('BUG: CORR_INFO must be BB or CL at this point: this line should never be called.');
    end
    
elseif ncomp==1
    %% Single component G(2) analysis
    if isequal(CORR_INFO,'BB')
        % Back-to-back G2 analysis
        for i=1:nShot
            shot_tmp=data1{i};
            nAtom=size(shot_tmp,1); % number of counts in this shot
            diff_tmp=[];   % diff vectors for pair search
            
            for j=1:nAtom
                % back-to-back condition
                this_atom=shot_tmp(j,:);    % polar-vector for this atom (to find pairs)
                other_tmp=shot_tmp([1:j-1,j+1:end],:);  % all other counts
                
                % diff for BB in polar
                diff_tmp(:,1)=other_tmp(:,1)-this_atom(1);   % diff in norm
                diff_tmp(:,2)=real(acos(cos(this_atom(3)).*cos(other_tmp(:,3)).*cos(other_tmp(:,2)-this_atom(2)) ...
                    + sin(this_atom(3)).*sin(other_tmp(:,3))));  % diff angle (see dot-product in polar coord)
                
                G2_SINGLE=G2_SINGLE+nhist(diff_tmp,BIN_EDGE);	% update G2
            end
        end
        
        % all shots - except self
        for i=1:nShot
            data_collated=vertcat(data1{[1:i-1,i+1:end]}); % collate all except self
            %data_collated=vertcat(data2{:}); % collate all shots inc. self
            nAtom=size(data1{i},1);
            diff_tmp=[];
            
            for j=1:nAtom
                % back-to-back condition
                this_atom=data1{i}(j,:);
                
                % diff for BB in polar
                diff_tmp(:,1)=data_collated(:,1)-this_atom(1);   % diff in norm
                diff_tmp(:,2)=real(acos(cos(this_atom(3)).*cos(data_collated(:,3)).*cos(data_collated(:,2)-this_atom(2)) ...
                    + sin(this_atom(3)).*sin(data_collated(:,3))));  % diff angle (see dot-product in polar coord)
                
                G2_ALL=G2_ALL+nhist(diff_tmp,BIN_EDGE);     % update G2
            end
        end
        
    elseif isequal(CORR_INFO,'CL')
        error('CL is not set up yet');
    else
        error('BUG: CORR_INFO must be BB or CL at this point: this line should never be called.');
    end
end

if VERBOSE>0
    t_fun_end=toc(t_fun_start);
    fprintf('Total elapsed time (s): %7.2f\n',t_fun_end);
end