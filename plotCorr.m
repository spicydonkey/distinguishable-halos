function [HFIG]=plotCorr(FIGNUM,CORRDATA,INDEX)
% wrapper to plot correlation functions
% Plots correlation function data in CORRDATA - INDEX locates which config
%
% INPUT
% CORRDATA: correlation data with configuration and results
%   contains fields including: type, nBin, bEdge, bCent, G2shot, G2all, g2,
%   etc.
% INDEX: which one?
%
% OUTPUT
% HFIG: generated figure handle
%

if isequal(CORRDATA.type{INDEX}.coord,'angular')
    HFIG=plotCorrAngle(FIGNUM,CORRDATA,INDEX);
    return
elseif isequal(CORRDATA.type{INDEX}.coord,'cart')
    HFIG=plotCorrCart(FIGNUM,CORRDATA,INDEX); 
    return
else
    error('Correlation must in evaluated in "angular" or "cart".');
end