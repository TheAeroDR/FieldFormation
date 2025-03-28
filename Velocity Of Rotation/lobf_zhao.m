load("zhao_data.mat")

r = [100 90 80 70 60 50 40 30];
vel = [data(:,1) data(:,3) data(:,5) data(:,7) data(:,9) data(:,11) data(:,13) data(:,15)];
height = [data(:,2) data(:,4) data(:,6) data(:,8) data(:,10) data(:,12) data(:,14) data(:,16)];
av_vel = mean(vel(8:20,:));

plot(r,av_vel,'.','Color', [0.85,0.33,0.10],'MarkerSize',16)
hold on

model = fittype('a*x.^2 + b*x +c','coeff', {'a','b','c'});

opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.MaxIter = 4000;
opts.StartPoint = [1,1,1];


[fitresult, gof] = fit(r', av_vel', model, opts);

model2 = fittype('d*x','coeff', {'d'});

opts2 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts2.Display = 'Off';
opts2.MaxIter = 4000;
opts2.StartPoint = [1];

[fitresult2, gof2] = fit(r', av_vel', model2, opts2);

x = linspace(0,100,101);

fit = fitresult(x);

plot(x,fit,'Color', [0.00,0.45,0.74])
plot(x,fit2)


xlabel('Radius $[m]$')
ylabel('Tangential Velocity $[ms^{-1}]$')
legend('Data', 'Fit $y = 0.0007982x^2 -0.1637x +10.68$');