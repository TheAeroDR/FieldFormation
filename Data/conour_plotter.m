% files = dir(['./10cm_res/contour_part_*.csv']);
% 
% x_all = [];
% y_all = [];
% z_all = [];
% Bz_all = [];
% Ez_all = [];
% 
% % Read data from each file and store all x, z, Bz values
% for i = 1:length(files)
%     % Read data from file
%     data = readtable(fullfile(files(i).folder, files(i).name));
% 
%     % Append to arrays
%     x_all = [x_all; data.x]; 
%     y_all = [y_all; data.y]; 
%     z_all = [z_all; data.z];
%     Bz_all = [Bz_all; data.Bz];
%     Ez_all = [Ez_all; -1*data.Ez];
% end
% 
% [x_all,y_all,z_all,Bz_all,Ez_all];
% array2table(ans);
% ans.Properties.VariableNames{1} = 'x';
% ans.Properties.VariableNames{2} = 'y';
% ans.Properties.VariableNames{3} = 'z';
% ans.Properties.VariableNames{4} = 'Bz';
% ans.Properties.VariableNames{5} = 'Ez';
% writetable(ans,"contour_3ed_10cm.csv");

%% vertical plane
load('2ed5_10cm.mat')

x_unique = unique(x_all, 'sorted');
y_unique = unique(y_all, 'sorted');
z_unique = unique(z_all, 'sorted');

[X,Z] = meshgrid(x_unique,z_unique);
Bz = griddata(x_all,z_all,Bz_all,X,Z);
Ez = griddata(x_all,z_all,Ez_all,X,Z);

Bz(Bz>1e-9)=1e-9;
Bz(Bz<-1e-9)=-1e-9;
Ez(Ez>5e6)=5e6;
Ez(Ez<-5e6)=-5e6;

figure(1)
tiledlayout(1,2)
nexttile(1)
surf(x_unique, z_unique, Bz, 'EdgeColor', 'none');
view(2);
hold on

diam = 7;
x_rect = [diam/4, diam/4 + diam/4, diam/4 + diam/4, diam/4, diam/4]; % x-coordinates
z_rect = [0, 0, 14.5, 14.5, 0]; % z-coordinates
z_top = max(Bz(:)) + 1e-10; % Place rectangle slightly above max Bz
plot3(x_rect, z_rect, repmat(z_top, size(x_rect)), 'k--', 'LineWidth', 2);
x_rect2 = [-diam/2, -diam/2 + diam/4, -diam/2 + diam/4, -diam/2, -diam/2];
plot3(x_rect2, z_rect, repmat(z_top, size(x_rect2)), 'k--', 'LineWidth', 2);

axis([-25,25,-inf,25])

axis square
xlabel('x [m]');
ylabel('z [m]');
c = colorbar();
c.Label.String = '$B_z$ [T]';
c.Label.Interpreter = 'latex';

figure(2)
tiledlayout(1,2)
nexttile(1)
surf(x_unique, z_unique, Ez, 'EdgeColor', 'none');
view(2);
hold on

diam = 7;
x_rect = [diam/4, diam/4 + diam/4, diam/4 + diam/4, diam/4, diam/4]; % x-coordinates
z_rect = [0, 0, 14.5, 14.5, 0]; % z-coordinates
z_top = max(Ez(:)) + 1; % Place rectangle slightly above max Bz
plot3(x_rect, z_rect, repmat(z_top, size(x_rect)), 'k--', 'LineWidth', 2);
x_rect2 = [-diam/2, -diam/2 + diam/4, -diam/2 + diam/4, -diam/2, -diam/2];
plot3(x_rect2, z_rect, repmat(z_top, size(x_rect2)), 'k--', 'LineWidth', 2);

axis([-25,25,-inf,25])

axis square
xlabel('x [m]');
ylabel('z [m]');
c = colorbar();
c.Label.String = '$E_z$ [V/m]';
c.Label.Interpreter = 'latex';

%%
load('2ed5_10cmh.mat')

x_unique = unique(x_all, 'sorted');
y_unique = unique(y_all, 'sorted');
z_unique = unique(z_all, 'sorted');

% Reshape Bz into a 2D matrix
[X,Y] = meshgrid(x_unique,y_unique);
Bzh = griddata(x_all,y_all,Bz_all,X,Y);
Ezh = griddata(x_all,y_all,Ez_all,X,Y);

Bzh(Bzh>1e-9)=1e-9;
Bzh(Bzh<-1e-9)=-1e-9;
Ezh(Ezh>5e6)=5e6;
Ezh(Ezh<-5e6)=-5e6;

figure(1)
nexttile(2)
surf(x_unique, y_unique, Bzh, 'EdgeColor', 'none');
view(2);
hold on

contour(x_unique, y_unique, Bzh, 200);

theta = linspace(0, 2*pi, 500); % Circle angles

diam = 7;

x_circle1 = 0 + diam/2 * cos(theta);
y_circle1 = 0 + diam/2 * sin(theta);

x_circle2 = 0 + diam/4 * cos(theta);
y_circle2 = 0 + diam/4 * sin(theta);

% Plot circles on top of the surface
z_top = max(Bzh(:)) + 1e-10; % Place rectangle slightly above max Bz

plot3(x_circle1, y_circle1, repmat(z_top, size(x_circle1)), 'k--', 'LineWidth', 2);
plot3(x_circle2, y_circle2, repmat(z_top, size(x_circle2)), 'k--', 'LineWidth', 2);

axis([-25,25,-25,25])

axis square
xlabel('x [m]');
ylabel('y [m]');
c = colorbar();
c.Label.String = '$B_z$ [T]';
c.Label.Interpreter = 'latex';

figure(2)
nexttile(2)
surf(x_unique, y_unique, Ezh, 'EdgeColor', 'none');
view(2);
hold on


theta = linspace(0, 2*pi, 500); % Circle angles

x_circle1 = 0 + diam/2 * cos(theta);
y_circle1 = 0 + diam/2 * sin(theta);

x_circle2 = 0 + diam/4 * cos(theta);
y_circle2 = 0 + diam/4 * sin(theta);

% Plot circles on top of the surface
z_top = max(Ezh(:)) + 1; % Place rectangle slightly above max Bz

plot3(x_circle1, y_circle1, repmat(z_top, size(x_circle1)), 'k--', 'LineWidth', 2);
plot3(x_circle2, y_circle2, repmat(z_top, size(x_circle2)), 'k--', 'LineWidth', 2);

axis([-25,25,-25,25])

axis square
xlabel('x [m]');
ylabel('y [m]');
c = colorbar();
c.Label.String = '$E_z$ [V/m]';
c.Label.Interpreter = 'latex';

%%
load("farrell_digitisation.mat");

temp = sortrows(farrell_m);
temp_fit = smooth(farrell_m(:,1), farrell_m(:,2), 0.01, 'lowess');
transect = 0;

% Find the closest z index
[~, t_idx] = min(abs(y_unique - transect));

% Extract Bz values along this z index
Bz_transect = Bzh(t_idx, :);
Ez_transect = Ezh(t_idx, :);

dd_vel = 3;
dd_centre = 53030.5;
% Plot the transect
figure(3)
tiledlayout(2,1)
nexttile(1)
plot(x_unique, Bz_transect*1e9, '-x', 'LineWidth', 1.5);
xlabel('x [m]');
ylabel('$B_z$ [nT]');
hold on
plot(dd_vel * (farrell_m(:,1)-dd_centre),farrell_m(:,2),'.','MarkerSize',20)

nexttile(2)
Ez_transect(Ez_transect>4300)=4300;
Ez_transect(Ez_transect<-4300)=-4300;
plot(x_unique, Ez_transect*1e-3, '-x', 'LineWidth', 1.5);
xlabel('x [m]');
ylabel('$E_z$ [kV/m]');
hold on
plot(dd_vel * (farrell_e(:,1)-dd_centre),farrell_e(:,2),'.','MarkerSize',20)

