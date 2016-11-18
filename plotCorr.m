function [HFIG]=plotCorr(FIGNUM,CORRDATA,CORRCONFIG,INDEX)
% wrapper to plot correlation functions
% Plots correlation function data in CORRDATA - INDEX locates which config
%
% INPUT
% CORRDATA: correlation results
%   contains fields including:  bEdge, bCent, G2shot, G2all, g2,
%
% CORRCONFIG: type, nBin
% INDEX: which one?
%
% OUTPUT
% HFIG: generated figure handle
%

if isequal(CORRCONFIG.type{INDEX}.coord,'angular')
    HFIG=plotCorrAngle(FIGNUM,CORRDATA,CORRCONFIG,INDEX);
    return
elseif isequal(CORRCONFIG.type{INDEX}.coord,'cart')
    HFIG=plotCorrCart(FIGNUM,CORRDATA,CORRCONFIG,INDEX); 
    return
else
    error('Correlation must in evaluated in "angular" or "cart".');
end