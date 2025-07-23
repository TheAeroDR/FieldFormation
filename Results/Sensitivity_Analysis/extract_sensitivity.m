filestring_base = sprintf('sensitivity100000_timeseries_');
filestring_search = sprintf('./results/%s*.csv',filestring_base);


% get file names
files = dir(filestring_search);

%extract unique string identifiers
for i = 1:length(files)
    data_name{i} = erase(files(i).name, filestring_base);
    data_name{i} = erase(data_name{i}, '.csv');
    data_name{i} = strrep(data_name{i},'_',' ');
    all_data{i} = readtable([files(i).folder,'/',files(i).name]);
end

%%
timesteps = unique(all_data{1}.timestep);

field_names = {'$B_x$ [T]', '$B_y$ [T]', '$B_z$ [V/m]', '$E_x$ [V/m]', '$E_y$ [V/m]', '$E_z$ [V/m]'};

figure(1)
tiledlayout(2,3);
figure(2)
tiledlayout(2,3);


n_runs = 1000;
n_timesteps = length(timesteps);

for i = 1:length(all_data)
    data = all_data{i};

    Bx_matrix = NaN(n_runs, n_timesteps);
    By_matrix = NaN(n_runs, n_timesteps);
    Bz_matrix = NaN(n_runs, n_timesteps);
    Ex_matrix = NaN(n_runs, n_timesteps);
    Ey_matrix = NaN(n_runs, n_timesteps);
    Ez_matrix = NaN(n_runs, n_timesteps);

    for j = 1:length(timesteps)
        timestep = timesteps(j);
        mask = data.timestep == timestep;
        data_subset = data(mask, :);

        n_points = height(data_subset);

        Bx_matrix(:,j) = data_subset.Bx;
        By_matrix(:,j) = data_subset.By;
        Bz_matrix(:,j) = data_subset.Bz;
        Ex_matrix(:,j) = data_subset.Ex;
        Ey_matrix(:,j) = data_subset.Ey;
        Ez_matrix(:,j) = data_subset.Ez;
    end
    field_matrices = {Bx_matrix, By_matrix, Bz_matrix, Ex_matrix, Ey_matrix, Ez_matrix};
    
    for j = 1:length(field_names)
        figure(1)
        nexttile(j);
        
        temp = iosr.statistics.boxPlot(timesteps, field_matrices{j}, 'scaleWidth', true, 'notch', true, 'showOutliers', false);
        if j<0
            xticklabels([])
        end
        ylabel(field_names(j))

        figure(2)
        nexttile(j)
        hold on
        plot(timesteps,temp.statistics.median)
        ylabel(field_names(j))

        
   
    end
    display(data_name(i));
end
%%
[~,t_index] = min(abs(timesteps-20));
[~,baseline_index] = max(matches(data_name,'baseline'));
medians = NaN([length(all_data),length(field_names)]);
for i = 1:length(data_name)
    data = all_data{i};
    mask = data.timestep == timesteps(t_index);
    data_sub = data(mask, :);
    medians(i,:) = [median(data_sub.Bx) median(data_sub.By) median(data_sub.Bz) median(data_sub.Ex) median(data_sub.Ey) median(data_sub.Ez)];
end

for i =1:length(data_name)
    perc_change(i,:) = (medians(i,:) - medians(baseline_index,:))./medians(baseline_index,:);
end

minus_set = contains(data_name,'minus');
minus_set(22) = 1;
plus_set = contains(data_name,'plus');
plus_set(23) = 1;

grouped = [];
grouped(1,:,:) = perc_change(minus_set,:);
grouped(2,:,:) = perc_change(plus_set,:);


field_names_stringy = {'B_x$', 'B_y$', 'B_z$', 'E_x$', 'E_y$', 'E_z$'};

figure
tiledlayout(6,1)
for i = 1:length(field_names)
    nexttile(i)
    temp = squeeze(grouped(:,:,i));
    bar(temp');
    ylabel(['$\Delta{}' field_names_stringy{i} '[\%]'])

end
legend('negative change in parameter','positive change in parameter','Location','none','Orientation','horizontal')