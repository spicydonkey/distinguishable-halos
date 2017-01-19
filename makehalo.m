function p_out = makehalo(n_collide,det_qe,p_dist,bgnd)
% function [p_out,is_detected,p_in,p_0_abs,scat_angle,p_0_scat] = makehalo(n_collide,det_qe,p_dist,bgnd)
%% Function to generate a pair of s-wave scattered halos
% DK SHIN
% 18/01/2017

%% I/O
% Inputs:
%   n_collide   {num atoms in collision; per colliding source}
%   det_qe      {detection efficiency}
%   p_dist      {momentum width/uncert of source; Gaussian, symmetric}
%   bgnd        {background counts}

% Outputs:
%   scat_prod   {scattering proudcts}   %data structure?

is_detected=(rand(n_collide,2)<det_qe);   % boolean detectable table
% TODO: cull all non-detected pairs? - would save much memory
n_pair=n_collide;   % num of collision pairs to simulate

% Assign collision momenta from given distribution
p_in=cell(2,1);
for ii=1:2  % collision source
    for jj=1:3  % momentum xyz components
        p_in{ii}(:,jj)=random('norm',(-1)^ii*p_dist(jj,1),p_dist(jj,2),[n_pair,1]);
    end
end

%% Two-body collision: s-wave
P_com=(p_in{1}+p_in{2})/2;      % collision com frame
p_0=p_in{1}-P_com;              % source-1 momentum in com
p_0_abs=sqrt(sum(p_0.^2,2));    % abs momentum mag of colliding atoms in com

scat_angle=zeros(n_pair,2);
scat_angle(:,1)=acos(2*rand([n_pair,1])-1);     % s-wave polar scat angle
scat_angle(:,2)=2*pi*rand([n_pair,1]);      % azimuthal scattering angle

% scattered momenta in com
p_0_scat=zeros(n_pair,3);   
p_0_scat(:,1)=p_0_abs.*sin(scat_angle(:,1)).*cos(scat_angle(:,2));
p_0_scat(:,2)=p_0_abs.*sin(scat_angle(:,1)).*sin(scat_angle(:,2));
p_0_scat(:,3)=p_0_abs.*cos(scat_angle(:,1));

% collision products
p_out=cell(2,1);
p_out{1}=p_0_scat+P_com;    % transform scattered momentum back to orig ref frame
p_out{2}=-p_0_scat+P_com;

% build detected counts
for ii=1:2
    p_out{ii}=p_out{ii}(is_detected(:,ii),:);
end