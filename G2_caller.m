function [G2_SINGLE,G2_ALL]=G2_caller(DATA,BIN_EDGE,TYPE,OPTIONAL,VERBOSE)
% Wrapper to call different G2 routines
%
% INPUT
% TYPE: must be a string 'cart' or 'polar'
%
% OUTPUT: see G2_cart.m or G2_angular.m
%

if ~exist('TYPE','var')
    error('Specify TYPE for G2_caller: must be either "cart" of "angular".');
end
if ~(isequal(TYPE,'cart')||isequal(TYPE,'angular'))
    error('TYPE must be either "cart" or "angular".');
end

if isequal(TYPE,'cart')
    [G2_SINGLE,G2_ALL]=G2_cart(DATA,BIN_EDGE,OPTIONAL,VERBOSE);
elseif isequal(TYPE,'angular')
    [G2_SINGLE,G2_ALL]=G2_angular(DATA,BIN_EDGE,VERBOSE);
else
    % some other routine would go here
    return
end