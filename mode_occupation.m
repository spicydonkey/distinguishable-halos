kr=1;       % halo radius
wbb=0.04;   % bb correlation length
dk=0.053;   % halo width rms
Nsc=950*2;  % number scattered

Vm=(2*pi)^(3/2)*(wbb/1.1)^3        % mode/phase grain volume
Vs=4*pi*sqrt(2*pi)*kr^2*dk         % total scattering volume

M=Vs/Vm    % number of modes in halo

n=Nsc/M      % halo mode occupancy

g2_BB_max=1+1/n    % predicted g2 peak
