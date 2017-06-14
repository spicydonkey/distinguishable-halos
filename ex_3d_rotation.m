% Rotational matrix
npaxis=100;     % number of points along basis
ubasis=linspace(-1,1,npaxis)';      % 1D vector - linearly spaced coords

basis3d=cell(3,1);  % lin-spaced points along each 3D basis
for ii=1:3
    % create lin-spaced points along each basis vectors
    basis3d{ii}=zeros(npaxis,3);
    basis3d{ii}(:,ii)=ubasis;
end

% plotting
colors=distinguishable_colors(3);
hfig=figure();
subplot(1,2,1);
hold on;
for ii=1:3
    scatter3(basis3d{ii}(:,1),basis3d{ii}(:,2),basis3d{ii}(:,3),'.','MarkerFaceColor',colors(ii,:));
end
axis equal;
view(3);

% Rotate original Cart coord by Euler angle (alpha,beta,gamma)
    % 3D Euler angles --> 3D rotation matrix
eul_angle=[0,0,pi/4];
R=euler2rotm(eul_angle);
rot_basis3d=basis3d;
for ii=1:3
    rot_basis3d{ii}=(R*rot_basis3d{ii}')';
end    

% plot rotated basis
figure(hfig);
subplot(1,2,2);
hold on;
for ii=1:3
    scatter3(rot_basis3d{ii}(:,1),rot_basis3d{ii}(:,2),rot_basis3d{ii}(:,3),'.','MarkerFaceColor',colors(ii,:));
end
axis equal;
view(3);