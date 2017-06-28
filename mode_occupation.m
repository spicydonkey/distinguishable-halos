% % RUN 1
% kr=1;       % halo radius
% wbb=0.03;   % bb correlation length
% dk=0.027;   % halo width rms
% Nsc=45;  % number scattered

% % RUN2
% kr=1;       % halo radius
% wbb=0.028;   % bb correlation length
% dk=0.04;   % halo width rms
% Nsc=120;  % number scattered

% % 20170622
% kr=1;
% wbb=0.045;
% dk=0.06;
% Nsc=315;

% % 20170624 --> BAD
% kr=1;
% wbb=0.03;   % TODO
% dk=0.1;
% Nsc=110;

% % 20170627
% kr=1;
% wbb=0.03;       % TODO
% dk=0.12;
% Nsc=38;

% RUN4
kr=1;
wbb=0.03;
dk=0.033;
Nsc=160;


Vm=(2*pi)^(3/2)*(wbb/1.1)^3        % mode/phase grain volume
Vs=4*pi*sqrt(2*pi)*kr^2*dk         % total scattering volume

M=Vs/Vm    % number of modes in halo

n=Nsc/M      % halo mode occupancy

g2_BB_max=1+1/n    % predicted g2 peak
