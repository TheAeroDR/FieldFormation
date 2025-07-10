% MATLAB script to load convergence results
% Run this after convergence_cuda.exe completes

clear; clc;

fprintf('Loading convergence study data...\n');

% Load individual particle count data
particle_counts = [1000, 10000, 100000, 1000000, 10000000];
loaded_data = struct();

for i = 1:length(particle_counts)
    N = particle_counts(i);
    
    % Load time series data
    timeseries_file = sprintf('convergence_%d_particles_timeseries.csv', N);
    
    if exist(timeseries_file, 'file')
        loaded_data.(sprintf('N%d_timeseries', N)) = readtable(timeseries_file);
        fprintf('Loaded: %s (%d data points)\n', timeseries_file, height(loaded_data.(sprintf('N%d_timeseries', N))));
    else
        fprintf('Missing: %s\n', timeseries_file);
    end
    
end


%% Create IOSR boxplots with tiled layout
% First, create figures with tiledlayouts for each field component
field_names = {'Bx', 'By', 'Bz', 'Ex', 'Ey', 'Ez'};
field_units = {'$B_x$ [T]', '$B_y$ [T]', '$B_z$ [T]', '$E_x$ [V/m]', '$E_y$ [V/m]', '$E_z$ [V/m]'};

% Create figures and tiledlayouts first
fig_handles = [];
for i = 1:length(field_names)
    fig_handles(i) = figure(i);
    tiledlayout(length(particle_counts), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

end

n_runs = unique(loaded_data.N1000_timeseries.run);

boxy_bx = NaN(n_runs, length(particle_counts));
boxy_by = NaN(n_runs, length(particle_counts));
boxy_bz = NaN(n_runs, length(particle_counts));
boxy_ex = NaN(n_runs, length(particle_counts));
boxy_ey = NaN(n_runs, length(particle_counts));
boxy_ez = NaN(n_runs, length(particle_counts));

for j = 1:length(particle_counts)
    N = particle_counts(j);
    field_name = sprintf('N%d_timeseries', N);
    
    % Check if data exists for this particle count
    if ~isfield(loaded_data, field_name)
        fprintf('Skipping N=%d - no timeseries data found\n', N);
        continue;
    end
    
    % Get the data for current count
    data = loaded_data.(field_name);
    
    % Get unique timesteps
    unique_timesteps = unique(data.timestep);
    unique_times = unique(data.time);
      
    % Prepare data matrices for boxplots (timesteps x runs)
    n_timesteps = length(unique_timesteps);
    n_runs = length(unique(data.run));
    
    % Check data distribution across timesteps
    for i = 1:min(5, n_timesteps)
        timestep = unique_timesteps(i);
        mask = data.timestep == timestep;
        n_points = sum(mask);
    end
    
    % Initialize matrices
    Bx_matrix = NaN(n_runs, n_timesteps);
    By_matrix = NaN(n_runs, n_timesteps);
    Bz_matrix = NaN(n_runs, n_timesteps);
    Ex_matrix = NaN(n_runs, n_timesteps);
    Ey_matrix = NaN(n_runs, n_timesteps);
    Ez_matrix = NaN(n_runs, n_timesteps);
    
    % Fill matrices - handle variable number of data points per timestep
    for t_idx = 1:n_timesteps
        timestep = unique_timesteps(t_idx);
        mask = data.timestep == timestep;
        data_subset = data(mask, :);
        
        n_points = height(data_subset);
        if n_points > 0
            % Fill available data points, leave rest as NaN
            Bx_matrix(1:n_points, t_idx) = data_subset.Bx;
            By_matrix(1:n_points, t_idx) = data_subset.By;
            Bz_matrix(1:n_points, t_idx) = data_subset.Bz;
            Ex_matrix(1:n_points, t_idx) = data_subset.Ex;
            Ey_matrix(1:n_points, t_idx) = data_subset.Ey;
            Ez_matrix(1:n_points, t_idx) = data_subset.Ez;
        end
    end
    
    % Remove rows that are all NaN (if some timesteps have fewer runs)
    all_nan_rows = all(isnan(Bx_matrix), 2);
    if any(all_nan_rows)
        Bx_matrix(all_nan_rows, :) = [];
        By_matrix(all_nan_rows, :) = [];
        Bz_matrix(all_nan_rows, :) = [];
        Ex_matrix(all_nan_rows, :) = [];
        Ey_matrix(all_nan_rows, :) = [];
        Ez_matrix(all_nan_rows, :) = [];
    end
    
    % Create plots for each field component
    field_matrices = {Bx_matrix, By_matrix, Bz_matrix, Ex_matrix, Ey_matrix, Ez_matrix};
    
    for i = 1:length(field_names)
        figure(fig_handles(i));
        nexttile;
        
        try
            iosr.statistics.boxPlot(unique_times, field_matrices{i}, 'scaleWidth', true, 'notch', true, 'showOutliers', false);
        catch
            boxplot(field_matrices{i}, unique_times);
        end

        if j == 5
            xlabel('Time (s)', 'FontSize', 10);
        elseif j == 3
            ylabel(sprintf('%s',field_units{i}))
            xticklabels([])
        else
            xticklabels([])
        end      
       
    end

        boxy_bx(:,j) = Bx_matrix(:,9);
        boxy_by(:,j) = By_matrix(:,9);
        boxy_bz(:,j) = Bz_matrix(:,9);
        boxy_ex(:,j) = Ex_matrix(:,9);
        boxy_ey(:,j) = Ey_matrix(:,9);
        boxy_ez(:,j) = Ez_matrix(:,9); 
end

%%
particles_10_form = {'$10^3$','$10^4$','$10^5$','$10^6$','$10^7$'};
figure
iosr.statistics.boxPlot(boxy_bx, 'scaleWidth', true, 'notch', true, 'showOutliers', false);
xticklabels(particles_10_form)
xlabel("N")
ylabel(field_units{1})
figure
iosr.statistics.boxPlot(boxy_by, 'scaleWidth', true, 'notch', true, 'showOutliers', false);
xticklabels(particles_10_form)
xlabel("N")
ylabel(field_units{2})
figure
iosr.statistics.boxPlot(boxy_bz, 'scaleWidth', true, 'notch', true, 'showOutliers', false);
xticklabels(particles_10_form)
xlabel("N")
ylabel(field_units{3})
figure
iosr.statistics.boxPlot(boxy_ex, 'scaleWidth', true, 'notch', true, 'showOutliers', false);
xticklabels(particles_10_form)
xlabel("N")
ylabel(field_units{4})
figure
iosr.statistics.boxPlot(boxy_ey, 'scaleWidth', true, 'notch', true, 'showOutliers', false);
xticklabels(particles_10_form)
xlabel("N")
ylabel(field_units{5})
figure
iosr.statistics.boxPlot(boxy_ez, 'scaleWidth', true, 'notch', true, 'showOutliers', false);
xticklabels(particles_10_form)
xlabel("N")
ylabel(field_units{6})