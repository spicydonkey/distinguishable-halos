function plotCorr(CORRDATA,CORRCONFIG)
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

if isequal(CORRCONFIG.type.coord,'angular')
    plotCorrAngle(CORRDATA,CORRCONFIG);
    return
elseif isequal(CORRCONFIG.type.coord,'cart')
    plotCorrCart(CORRDATA,CORRCONFIG); 
    return
else
    error('Correlation must in evaluated in "angular" or "cart".');
end