% Example script to test makehalo.m
% DK Shin
% 19/1/2017

clear all; clc; close all;

n_collide=10000;
det_qe=0.1;

% source momentum distribution (COM): [Pcollision,Pwidth]
p_x=[0,0.1];
p_y=[0,0.1];
p_z=[1,0.1];
p_dist=[p_x;p_y;p_z];

bgnd=0; % TODO no functionality

% Simulate collision
p_out=makehalo(n_collide,det_qe,p_dist,bgnd);

% Plot result
figure(11);
scatter3(p_out{1}(:,1),p_out{1}(:,2),p_out{1}(:,3),'b.');
hold on;
scatter3(p_out{2}(:,1),p_out{2}(:,2),p_out{2}(:,3),'r.');
axis equal;
title('S-wave scattering simulation in COM frame');
xlabel('p_x');
ylabel('p_y');
zlabel('p_z');