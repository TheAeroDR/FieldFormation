close all

load('Stanzel_data.mat')
D = table2array(data(:,1));
D_err = table2array(data(:,2)); 
H = table2array(data(:,3));
H_err = table2array(data(:,4));

f1 = figure(1);

%l1 = errorbar(D,H,H_err,H_err,D_err,D_err,'.');
%l1.MarkerSize = 16;
%l1.Color = [0.80,0.80,0.80];

l1 = plot(H,D/2,'.','DisplayName','Stanzel MX Data');
l1.MarkerSize = 16;

xlabel('Height [m]')
ylabel('Radius [m]')
f1.Children.XScale = 'log';
f1.Children.YScale = 'log';
hold on

model = fittype('a*x + b', 'coeff', {'a', 'b'});
model = fittype('0.5*x + b', 'coeff', {'b'});

opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.MaxIter = 4000;
opts.StartPoint = [0.7];

% Fit model to data.
[fitresult, gof] = fit(log10(H), log10(D/2), model, opts);

ci = confint(fitresult, 0.95);

fitH = linspace(70,5000);
fitR = (10^fitresult(0)) * fitH.^0.5;
fitRl = (10^ci(1)) * fitH.^0.5;
fitRh = (10^ci(2)) * fitH.^0.5;

plot(fitH,fitR,'DisplayName','Fit');

patch([fitH flip(fitH)],[fitRh flip(fitRl)],[0.5 0.5 0.5],'FaceAlpha',0.3,'EdgeColor', 'none','DisplayName', "95\% Confidence")

legend('Location','best')

xlim([70,5000])


%%
fit_data = table(log10(D), log10(H), log10(D_err), log10(H_err), 'VariableNames', {'x', 'y', 'x_error', 'y_error'});

ymxc_model = fitlm(fit_data, 'y ~ x', 'Weights', 1./(fit_data.y_error.^2));

xnyd_model = fitlm(fit_data, 'x ~ y', 'Weights', 1./(fit_data.x_error.^2));

x1_fit = linspace(min(fit_data.x), max(fit_data.x), 100);
y1_fit= (ymxc_model.Coefficients.Estimate(2) .* x1_fit) + ymxc_model.Coefficients.Estimate(1);

y2_fit = linspace(min(fit_data.y), max(fit_data.y),100);
x2_fit = (xnyd_model.Coefficients.Estimate(2) .* y2_fit) + xnyd_model.Coefficients.Estimate(1);

y3_fit = (0.6754 .* x1_fit) + 1.154;


plot(10.^x1_fit, 10.^y1_fit, 'Color', [0.00,0.45,0.74]);
plot(10.^x2_fit, 10.^y2_fit, 'Color', [0.47,0.67,0.19]);
plot(10.^x1_fit, 10.^y3_fit, 'k--');

legend('Data', 'Fit $y = 0.7409x + 1.0054$', 'Fit $x = 0.3806y + 1.2372$','Jackson Fit');