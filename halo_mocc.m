% s-wave scattering halo mode occupancy
% From SOMs in https://arxiv.org/abs/1702.03617

% n = halo_mocc(kr,w,N,sigk)
%   n : halo mode occupancy
%   kr : shell radius
%   w : rms width of shell
%   N : total number of atoms in the halo
%   sigk : rms width of the momentum distribution of the source BEC
%   (assuming isotropic Gaussian distribution) [sigk ~ wbb/1.1]

function n = halo_mocc(kr,w,N,sigk)
Vs=4*pi*sqrt(2*pi)*(kr^2)*w;        % total scattering volume
Vm=((2*pi)^(3/2))*(sigk^3);         % mode/phase grain volume

M=Vs/Vm;        % number of modes in halo

n=N/M;      % mode occupancy
end