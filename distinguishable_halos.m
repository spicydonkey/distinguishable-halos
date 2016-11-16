% Distinguishable halo
% DKS 16/11/2016

%% CONFIGS
%   <<Group into each stages>>
% raw - raw-data handling
% preproc - broad capture
% postproc - manipulations
% analysis

%% Load raw-data DLD/TXY
%   <<Skip if already done and configs are same>>
%   Process rawDLD/TXY archive to:
%     - get counts separated by shot# (ditch low-counts)
%     - save (ZXY)-counts to .mat file



%% Preprocess: Get halo counts
%   <<Skip if already done AND configs are same>>
%   Process loaded counts to:
%     - locate condensate accurately in each shot
%     - capture all counts in the halo with approx centre location and
%       a broad radial mask
%     - collate all shots with BEC oscillations cancelled (esp. halo)
%     - save to file (think about not saving BEC's - too many counts)

%% Fit halos
%   Process the broadly captured halos (ZXY):
%     - fit the point cloud to some shape (sphere/oblate spheroid/tri-axial
%       ellipsoid, etc.)
%     - get "centre", "radii", etc. for each halo
%     - reference to centre of mass (halo centres)
%     - extract well-fitted counts (tight capture)
%     - save to file (fitted halo with shape and fit params)

%% Postprocess: Manipulate halos
%   <<Skip if already done and configs are same>>
%   Manipulate the halos to account for distortions:
%     - e.g. spin in Z-axis, isotropic/anisotropic scaling, etc.
%     - save to file

%% Analysis
%   <<Skip if already done and configs are same>>
%   Correlation analysis (in angular, cartesian):
%     - cross-halo
%     - single-halo