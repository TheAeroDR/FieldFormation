x_all = [];
y_all = [];
z_all = [];
Bz_all = [];
Ez_all = [];

filestring = sprintf('./Farrell_2ed_CUDA_h1/field_t*.csv');

files = dir(filestring);

% Sort files by name to ensure correct temporal order
[~, idx] = sort({files.name});
files = files(idx);

for j = 1:length(files)
    data = readtable(fullfile(files(j).folder, files(j).name),'ReadVariableNames',true);

    x_all(:,j) = data.x;
    y_all(:,j) = data.y;
    z_all(:,j) = data.z;
    Ex_all(:,j) = data.Ex;
    Ey_all(:,j) = data.Ey;
    Ez_all(:,j) = data.Ez;
    Bx_all(:,j) = data.Bx;
    By_all(:,j) = data.By;
    Bz_all(:,j) = data.Bz;
end
%%
grid_size_x = length(unique(x_all(:,1)));
grid_size_y = length(unique(y_all(:,1)));

%%
x_reshaped = reshape(x_all(:,1), grid_size_x, grid_size_y)';
x_unique = x_reshaped(1,:);
y_reshaped = reshape(y_all(:,1), grid_size_x, grid_size_y)'; 
y_unique = y_reshaped(:,1);


[X,Y] = meshgrid(x_unique, y_unique);
for i = 1:size(Bz_all,2)
    Bx(:,:,i) = reshape(Bx_all(:,i), grid_size_x, grid_size_y)';
    By(:,:,i) = reshape(By_all(:,i), grid_size_x, grid_size_y)';
    Bz(:,:,i) = reshape(Bz_all(:,i), grid_size_x, grid_size_y)';
    Ex(:,:,i) = reshape(Ex_all(:,i), grid_size_x, grid_size_y)';
    Ey(:,:,i) = reshape(Ey_all(:,i), grid_size_x, grid_size_y)';
    Ez(:,:,i) = reshape(Ez_all(:,i), grid_size_x, grid_size_y)';
end

%%
figure;

bz_min = -5e-10;
bz_max = 5e-10;
ez_min = -6e6;
ez_max = 6e6;

t = 0:0.5:37.5;
x_centre = -102+5.38*t;
while true

    for i = 1:size(Bz,3)
        surf(X, Y, Bz(:,:,i), 'EdgeColor', 'none');
        colormap jet;
        colorbar;
        title(['Bz at frame ', num2str(i)]);
        xlabel('X'); ylabel('Y'); zlabel('Bz');
        axis tight;
        zlim([bz_min, bz_max]);
        xlim([x_centre(i)-10,x_centre(i)+10]);
        ylim([-10,10]);
        caxis([bz_min, bz_max]);
        view(2);
        pause(0.1);
    end
end

%%
field_names = {'$B_x$ [T]', '$B_y$ [T]', '$B_z$ [nT]', '$E_x$ [V/m]', '$E_y$ [V/m]', '$E_z$ [V/m]'};

figure
tiledlayout(2,3)
t_index = 39
for i = 1:length(field_names)
    if i == 1
        data = Bx(:,:,t_index);
        c_min = -5e-11;
        c_max = 5e-11;
    elseif i == 2
        data = By(:,:,t_index);
        c_min = -5e-10;
        c_max = 5e-10;
    elseif i == 3
        data = Bz(:,:,t_index);
        c_min = -5e-10;
        c_max = 5e-10;
    elseif i == 4
        data = Ex(:,:,t_index);
        c_min = -5e6;
        c_max = 5e6;
    elseif i == 5
        data = Ey(:,:,t_index);
        c_min = -5e6;
        c_max = 5e6;
    elseif i == 6
        data = Ez(:,:,t_index);
        c_min = -6e6;
        c_max = 6e6;
    end
    nexttile(i)
    surf(X, Y, data, 'EdgeColor', 'none');
    colormap jet;
    c = colorbar;
    c.Label.String = field_names(i);
    c.Label.Interpreter = "latex";
    c.TickLabelInterpreter = "latex";
    if i == 4 || i == 5 || i == 6
        xlabel('X');
    else
        xlabel('');
    end
    if i == 1 || i == 4
        ylabel('Y');
    else
        ylabel('');
    end
    zlabel(field_names);
    axis tight;
    xlim([-10,10]);
    ylim([-10,10]);
    caxis([c_min, c_max]);
    view(2);
    hold on

    diam = 7;
    theta = linspace(0,2*pi,100);
    x_circle1 = (-102+5.38*(19.0)) + diam/2 * cos(theta);
    y_circle1 = 0 + diam/2 * sin(theta);

    x_circle2 = (-102+5.38*(19.0)) + diam/4 * cos(theta);
    y_circle2 = 0 + diam/4 * sin(theta);

    % Plot circles on top of the surface
    z_top = max(max(data(:,:))) + 1e-10; % Place rectangle slightly above max

    plot3(x_circle1, y_circle1, repmat(z_top, size(x_circle1)), 'k--', 'LineWidth', 2);
    plot3(x_circle2, y_circle2, repmat(z_top, size(x_circle2)), 'k--', 'LineWidth', 2);

end

%%
x_all = [];
y_all = [];
z_all = [];
Bz_all = [];
Ez_all = [];

filestring = sprintf('./Farrell_2ed_CUDA_v0/field_t*.csv');

files = dir(filestring);

% Sort files by name to ensure correct temporal order
[~, idx] = sort({files.name});
files = files(idx);

for j = 1:length(files)
    data = readtable(fullfile(files(j).folder, files(j).name),'ReadVariableNames',true);

    x_all(:,j) = data.x;
    y_all(:,j) = data.y;
    z_all(:,j) = data.z;
    Ex_all(:,j) = data.Ex;
    Ey_all(:,j) = data.Ey;
    Ez_all(:,j) = data.Ez;
    Bx_all(:,j) = data.Bx;
    By_all(:,j) = data.By;
    Bz_all(:,j) = data.Bz;
end
%%
grid_size_x = length(unique(x_all(:,1)));
grid_size_y = length(unique(z_all(:,1)));

%%
x_reshaped = reshape(x_all(:,1), grid_size_x, grid_size_y)';
x_unique = x_reshaped(1,:);
y_reshaped = reshape(z_all(:,1), grid_size_x, grid_size_y)'; 
y_unique = y_reshaped(:,1);


[X,Y] = meshgrid(x_unique, y_unique);
for i = 1:size(Bz_all,2)
    Bx(:,:,i) = reshape(Bx_all(:,i), grid_size_x, grid_size_y)';
    By(:,:,i) = reshape(By_all(:,i), grid_size_x, grid_size_y)';
    Bz(:,:,i) = reshape(Bz_all(:,i), grid_size_x, grid_size_y)';
    Ex(:,:,i) = reshape(Ex_all(:,i), grid_size_x, grid_size_y)';
    Ey(:,:,i) = reshape(Ey_all(:,i), grid_size_x, grid_size_y)';
    Ez(:,:,i) = reshape(Ez_all(:,i), grid_size_x, grid_size_y)';
end
%%
%%
figure;

bz_min = -5e-10;
bz_max = 5e-10;
ez_min = -6e6;
ez_max = 6e6;

t = 0:0.5:37.5;
x_centre = -102+5.38*t;
while true

    for i = 1:size(Ez,3)
        surf(X, Y, Ez(:,:,i), 'EdgeColor', 'none');
        colormap jet;
        colorbar;
        title(['Ez at frame ', num2str(i)]);
        xlabel('X'); ylabel('Y'); zlabel('Ez');
        axis tight;
        zlim([ez_min, ez_max]);
        xlim([x_centre(i)-10,x_centre(i)+10]);
        ylim([0,20]);
        %caxis([ez_min, ez_max]);
        caxis([bz_min, bz_max]);
        view(2);
        pause(0.1);
    end
end

%%
field_names = {'$B_x$ [T]', '$B_y$ [T]', '$B_z$ [nT]', '$E_x$ [V/m]', '$E_y$ [V/m]', '$E_z$ [V/m]'};

figure
tiledlayout(2,3)
t_index = 39
for i = 1:length(field_names)
    if i == 1
        data = Bx(:,:,t_index);
        c_min = -5e-11;
        c_max = 5e-11;
    elseif i == 2
        data = By(:,:,t_index);
        c_min = -5e-10;
        c_max = 5e-10;
    elseif i == 3
        data = Bz(:,:,t_index);
        c_min = -5e-10;
        c_max = 5e-10;
    elseif i == 4
        data = Ex(:,:,t_index);
        c_min = -5e6;
        c_max = 5e6;
    elseif i == 5
        data = Ey(:,:,t_index);
        c_min = -5e6;
        c_max = 5e6;
    elseif i == 6
        data = Ez(:,:,t_index);
        c_min = -6e6;
        c_max = 6e6;
    end
    nexttile(i)
    surf(X, Y, data, 'EdgeColor', 'none');
    colormap jet;
    c = colorbar;
    c.Label.String = field_names(i);
    c.Label.Interpreter = "latex";
    c.TickLabelInterpreter = "latex";
    if i == 4 || i == 5 || i == 6
        xlabel('X');
    else
        xlabel('');
    end
    if i == 1 || i == 4
        ylabel('Y');
    else
        ylabel('');
    end
    zlabel(field_names);
    axis tight;
    xlim([-10,10]);
    ylim([0,20]);
    caxis([c_min, c_max]);
    view(2);
    hold on

    diam = 7;
    x_rect = [diam/4, diam/4 + diam/4, diam/4 + diam/4, diam/4, diam/4]; % x-coordinates
    z_rect = [0, 0, 14.5, 14.5, 0]; % z-coordinates
    z_top = max(Ez(:)) + 1; % Place rectangle slightly above max Bz
    plot3(x_rect, z_rect, repmat(z_top, size(x_rect)), 'k--', 'LineWidth', 2);
    x_rect2 = [-diam/2, -diam/2 + diam/4, -diam/2 + diam/4, -diam/2, -diam/2];
    plot3(x_rect2, z_rect, repmat(z_top, size(x_rect2)), 'k--', 'LineWidth', 2);

end
%%

B = sqrt(Bx(:,:,39).^2 + By(:,:,39).^2 + Bz(:,:,39).^2);
E = sqrt(Ex(:,:,39).^2 + Ey(:,:,39).^2 + Ez(:,:,39).^2);