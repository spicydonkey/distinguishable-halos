%% <template> CONFIG for capturing halo from TXY

%-----------------------------------------------------------------
% PARAMETERS
%-----------------------------------------------------------------
configs.bec.pos{1}=[20.7024,4.74e-3,2.72e-3];   % approx condensate locations (z,x,y)
configs.bec.Rmax{1}=7e-3;       % max condensate sph radius
configs.bec.dR_tail{1}=0.3;     % BEC tail radial frac diff

configs.bec.pos{2}=[20.7005,-7.38e-3,6.55e-3];
configs.bec.Rmax{2}=7e-3;
configs.bec.dR_tail{2}=0.3;

configs.halo.R{1}=11e-3;        % estimated radius of halo
configs.halo.dR{1}=0.15;        % radial mask fractional width (bi-dir)
configs.halo.R{2}=10e-3;
configs.halo.dR{2}=0.15;

configs.halo.zcap=0.7;   % z-cutoff (fractional wrt radius)