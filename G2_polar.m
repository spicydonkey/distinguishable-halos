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
%
% TODO: CRITICAL - self-G2 includes own count

% Input check
% DATA
if ~iscell(DATA)
    error('DATA should be a cell of counts.');
end
% BIN_EDGE
if ~iscell(BIN_EDGE)
    error('BIN_EDGE should be a cell of bin edges');
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
data1=DATA(:,1);
if size(DATA,2)==2
     % cross-data G2
    data2=DATA(:,2);
else
    % ordinary G2 for single component
    warning('CRITICAL BUG: Currently includes self to do G2!');
    data2=data1;
end
clear DATA;

nBin=[];
for i=1:length(BIN_EDGE)
    nBin(i)=length(BIN_EDGE{i})-1;  % number of bins from bin edges
end

% Initialise variables
nShot=size(data1,1);     % number of shots in data
if VERBOSE>1,disp([num2str(nShot),' shots to analyse for G2 (polar)...']);,end
G2_SINGLE=zeros(nBin);
G2_ALL=zeros(nBin);

% Branch G2 analysis so that condition is out of loop: just add other
%   conditions as a conditional branch following one as template
if isequal(CORR_INFO,'BB')
    % Back-to-back G2 analysis
    for i=1:nShot
        nAtom=size(data1{i},1); % number of counts in DATA1
        diff_tmp=[];   % diff vectors for pair search
        
        for j=1:nAtom
            % back-to-back condition
            this_atom=data1{i}(j,:);
            
            % diff for BB in polar
            diff_tmp(:,1)=data2{i}(:,1)-this_atom(1);   % diff in norm
            diff_tmp(:,2)=acos(cos(this_atom(3)).*cos(data2{i}(:,3)).*cos(data2{i}(:,2)-(this_atom(2)+pi)) ...
                + sin(-this_atom(3)).*sin(data2{i}(:,3)));  % diff angle
            
            count_tmp=histcn(diff_tmp,BIN_EDGE{1},BIN_EDGE{2});
            G2_SINGLE=G2_SINGLE+count_tmp;
        end
    end
    
    % all shots - (includes correlated shots)
    data_ncorr=vertcat(data2{:}); % collate all shots - including corr
    for i=1:nShot
        nAtom=size(data1{i},1);
        
        diff_tmp=[];
        for j=1:nAtom
            this_atom=data1{i}(j,:);
            
            % diff for BB in polar
            diff_tmp(:,1)=data_ncorr(:,1)-this_atom(1);   % diff in norm
            diff_tmp(:,2)=acos(cos(this_atom(3)).*cos(data_ncorr(:,3)).*cos(data_ncorr(:,2)-(this_atom(2)+pi)) ...
                + sin(-this_atom(3)).*sin(data_ncorr(:,3)));  % diff angle
            
            count_tmp=histcn(diff_tmp,BIN_EDGE{1},BIN_EDGE{2});
            G2_ALL=G2_ALL+count_tmp;
        end
    end
    
elseif isequal(CORR_INFO,'CL')
    error('CL is not set up yet');
else
    error('BUG: CORR_INFO must be BB or CL at this point: this line should never be called.');
end