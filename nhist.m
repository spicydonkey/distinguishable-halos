function COUNT = nhist(X, EDGE)
% Edited by DKS 07/11/2016 (original: histcn.m by Bruno Luong)
% function count = histcn(X, EDGE)
%  edge1, edge2, ..., edgeN
% Purpose: compute n-dimensional histogram
%
% INPUT
%   - X: (M x N) array, represents M-many, N-dimensional vectors
%   - EDGE: (1 x N) cell of edge bin vectors on dimension k, k=1...N.
%
% OUTPUT
%   - COUNT: N-dimensional array count of X on the bins, i.e.,
%         COUNT(i1,i2,...,iN) = cardinal of X such that
%                  EDGE{1}(i1) <= X(:,i1) < EDGE{1}(i1)+1
%                       ...
%                  EDGE{N}(iN) <= X(:,iN) < EDGE{N}(iN)+1
%
% TODO: Usage examples:
%
% 

% Error checks/Parse inputs
if ndims(X)>2
    error('nhist: X requires to be an (M x N) array of M points in R^N');
end

ndim = size(X,2);     % get the dimension of vectors

if ~iscell(EDGE)
    error('nhist: EDGE must be a "cell" of bin edges');
else
    if ndim~=length(EDGE)
        error('nhist: EDGE must be a (1 x Ndim) cell of bin edges.');
    end
end

% Do histogram
loc = zeros(size(X));   % Allocation of array bin location: index location of X in the bins
sz = zeros(1,ndim);     % size of COUNT array: number of bins in each dim
% Loop in the dimension
for d=1:ndim
    ed = EDGE{d};   % this dim's bin-edge vector
    sz(d) = length(ed)-1;   % num bin is length of edge vector - 1
    Xd = X(:,d);    % get vectors' component in this dim
    
    [~,~,loc(:,d)]=histcounts(Xd,ed);   % get bin location for this dimension
end % for-loop

if ndim==1
    sz = [sz 1];    % special case for a 1D-vectors
end

hasdata = all(loc>0,2);     % successfully binned vectors
COUNT = accumarray(loc(hasdata,:),1,sz);    % build COUNT

return

end % histcn