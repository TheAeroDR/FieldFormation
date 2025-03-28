close all
clear all
clc

% Get a list of all files in the directory
files = dir(['lack_model_output_1ed*.txt']);

% Create a cell array to store tables
dataTables = cell(1, length(files));
% %%  
% Loop through each file
for i = 1:length(files)
    % Get the current file name
    currentFileName = files(i).name;
    
    % Read the file into a table
    dataTable = readtable(currentFileName, 'Delimiter', ',', 'ReadVariableNames', true, 'HeaderLines', 9); %old files are 17
    
    % Store the table in the cell array with a sequential name
    dataTables{i} = dataTable;
end
%%
% Combine all data points for the fit
allRadii = [];
allCharges = [];
for i = 1:length(dataTables)
    allRadii = [allRadii; dataTables{i}.radii];
    allCharges = [allCharges; dataTables{i}.charge];
end

% Sort the data based on radii
[xData, yData] = prepareCurveData(allRadii, allCharges);

% Set up fittype and options.
ft = fittype( 'a*x^2 +c*x^-0.8', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 -Inf];
opts.MaxIter = 4000;
opts.Robust = 'Bisquare';
opts.StartPoint = [0.711594920483318 0.313254284887248];
opts.Upper = [Inf 0];

% Fit model to data.
[fitresult, gof] = fit(xData, yData, ft, opts);

figure
hold on
for i = 1:length(dataTables)
    % Plot individual data points
    plot(dataTables{i}.radii, dataTables{i}.charge, '.', 'MarkerSize', 30);
end
% Plot the fitted curve
plot(sort(xData),fitresult(sort(xData)));

% Add labels and legend
xlabel('Radii [$\mu$m]');
ylabel('Charge [e]');
legend('n = 150','n = 250','n = 500', 'Fit');
hold off

%%
figure
[xS, yS, zS] = sphere;
cmap=parula(length(dataTable.radii));
% Plot spheres:
for pt = 1:numel(dataTable.radii)
    hs=surf(dataTable.position_x_(pt)+xS.*dataTable.radii(pt), ...  % Shift and scale x data
       dataTable.position_y_(pt)+yS.*dataTable.radii(pt), ...  % Shift and scale y data
       dataTable.position_z_(pt)+zS.*dataTable.radii(pt), ...  % Shift and scale z data
       'EdgeColor', 'none');
       set(hs,'FaceColor',cmap(pt,:))
       set(hs,'FaceAlpha',1)
  hold on;
end
material shiny;            % Set material properties to shiny
set(gca, 'DataAspectRatio', [1,1,1]);  % Ensure equal aspect ratio for 3D plot
camlight('headlight');     % Add headlight effect for better visibility
axis square



%%

x = log(allRadii); 
y = sign(allCharges) .* log(abs(allCharges)); 
y = log(abs(allCharges)); 


[x, idx] = sort(x);

y = y(idx);

% Define the piecewise function
customModel = @(p, x) (x >= 2.7).*(p(1) + p(2)*x) + (x < 2.4).*(p(3) + p(4)*x);

% Initial parameter guesses [a, b, c, d]
p0 = [1, 1, -1, -1];

% Objective function for lsqcurvefit
objectiveFcn = @(p, x) customModel(p, x);

% Fit the model using lsqcurvefit
options = optimset('Display', 'iter');
params = lsqcurvefit(objectiveFcn, p0, x, y, [], [], options);

% Plot the results
figure;
scatter(x, y, 'filled'); hold on;
fittedY = customModel(params, x);
fittedY(fittedY==0)=NaN;
plot(x, fittedY, 'LineWidth', 2);
legend('Data', 'Fit');
xlabel('$\log_{10}(r)$'); ylabel('$\log_{10}(|Q|)$');
grid on;
