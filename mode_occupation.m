% % RUN 1
% kr=1;       % halo radius
% dk=0.027;   % halo width rms
% Nsc=45;  % number scattered
% wbb=0.03;   % bb correlation length

% % RUN2
% kr=1;       % halo radius
% dk=0.056;   % halo width rms
% Nsc=140;  % number scattered
% wbb=0.025;   % bb correlation length

% % 20170622 --> VERY HOT!
% kr=1;
% dk=0.044;
% Nsc=710;
% wbb=0.023;

% % 20170624 --> BAD
% kr=1;
% dk=0.1;
% Nsc=110;
% wbb=0.03;   % TODO

% % 20170627 --> BAD
% kr=1;
% dk=0.12;
% Nsc=38;
% wbb=0.03;       % TODO

% % RUN4
% kr=1;
% dk=0.034;
% Nsc=170;
% wbb=0.028;

% % RUN5 - high n!
% kr=1;
% dk=0.03;
% Nsc=1900;
% wbb=0.033;

% % RUN6
% kr=1;
% dk=0.033;
% Nsc=320;
% wbb=0.029;

% % RUN7 - very low number
% kr=1;
% dk=0.034;
% Nsc=60;
% wbb=0.02;


Vm=(2*pi)^(3/2)*(wbb/1.1)^3        % mode/phase grain volume
Vs=4*pi*sqrt(2*pi)*kr^2*dk         % total scattering volume

M=Vs/Vm    % number of modes in halo

n=Nsc/M      % halo mode occupancy

g2_BB_max=1+1/n    % predicted g2 peak
