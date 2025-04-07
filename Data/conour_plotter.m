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
load('2ed5_10cmv.mat')

x_unique = unique(x_all, 'sorted');
y_unique = unique(y_all, 'sorted');
z_unique = unique(z_all, 'sorted');

[X,Z] = meshgrid(x_unique,z_unique);
Bz = griddata(x_all,z_all,Bz_all,X,Z);
Ez = griddata(x_all,z_all,Ez_all,X,Z);

Bz_plim = prctile(Bz_all, 99.999);
Bz_mlim = prctile(Bz_all, 0.001);
Ez_plim = prctile(Ez_all, 99.999);
Ez_mlim = prctile(Ez_all, 0.001);

Bz(Bz>Bz_plim)=Bz_plim;
Bz(Bz<Bz_mlim)=Bz_mlim;
Ez(Ez>Ez_plim)=Ez_plim;
Ez(Ez<Ez_mlim)=Ez_mlim;

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

Bzh_plim = prctile(Bz_all, 99.999);
Bzh_mlim = prctile(Bz_all, 0.001);
Ezh_plim = prctile(Ez_all, 99.999);
Ezh_mlim = prctile(Ez_all, 0.001);

Bzh(Bzh>Bzh_plim)=Bzh_plim;
Bzh(Bzh<Bzh_mlim)=Bzh_mlim;
Ezh(Ezh>Ezh_plim)=Ezh_plim;
Ezh(Ezh<Ezh_mlim)=Ezh_mlim;

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
temp2 = sortrows(farrell_e);
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
plot(dd_vel * (temp(:,1)-dd_centre),temp(:,2),'.-','MarkerSize',20)
plot(5.7 * (temp(:,1)-dd_centre),temp(:,2),'s-','MarkerSize',8,'MarkerFaceColor',[0.9290 0.6940 0.1250]	)
plot(1.58 * (temp(:,1)-dd_centre),temp(:,2),'^-','MarkerSize',5,'MarkerFaceColor',[0.4940 0.1840 0.5560])
axis([-100,100,-inf,inf])

nexttile(2)
Ez_transect(Ez_transect>4300)=4300;
Ez_transect(Ez_transect<-4300)=-4300;
plot(x_unique, Ez_transect*1e-3, '-x', 'LineWidth', 1.5);
xlabel('x [m]');
ylabel('$E_z$ [kV/m]');
hold on
plot(dd_vel * (temp2(:,1)-dd_centre),temp2(:,2),'.-','MarkerSize',20)
plot(5.7 * (temp2(:,1)-dd_centre),temp2(:,2),'s-','MarkerSize',8,'MarkerFaceColor',[0.9290 0.6940 0.1250])
plot(1.58 * (temp2(:,1)-dd_centre),temp2(:,2),'^-','MarkerSize',5,'MarkerFaceColor',[0.4940 0.1840 0.5560])
axis([-100,100,-inf,inf])

legend('Simulated Field','Farrell Field, v=3m/s','Farrell Field, v=35.7/s','Farrell Field, v=1.58m/s','location','none')