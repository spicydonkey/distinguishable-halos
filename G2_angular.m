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

