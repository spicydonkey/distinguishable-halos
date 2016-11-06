function [G2_SINGLE,G2_ALL]=G2_cart(DATA,BIN_EDGE,CORR_INFO)
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
% CORR_INFO
if ~isexist('CORR_INFO','var')  % missing arg - default
    warning('CORR_INFO is not set by user: default to back-to-back (BB) G(2) analysis');
    CORR_INFO='BB';
elseif ~(isequal(CORR_INFO,'BB')||isequal(CORR_INFO,'CL'))  % invalid input
    warning('CORR_INFO must be "BB" or "CL". Setting to default "BB"');
    CORR_INFO='BB';
end

nShot=size(DATA,1);     % number of shots in data
