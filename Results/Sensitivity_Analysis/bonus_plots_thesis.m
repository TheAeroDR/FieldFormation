a = rand(100000,1);
R = 3.5;

baseline = R*(1-0.5*(a));

increased = R*(1-0.5*(a.^1.1));

reduced = R*(1-0.5*(a.^0.9));

edges = 1.75:0.05:3.5;
mids = 1.775:0.05:3.475;
bin_n = histcounts(baseline,edges);
bin_r = histcounts(reduced,edges);
bin_i = histcounts(increased,edges);

figure
plot(mids,bin_n)
hold on
plot(mids,bin_r)
plot(mids,bin_i)
legend('Baseline', 'Reduced Radius, 0.9', 'Increased Radius, 1.1')

%%
a = rand(10000,1);
R = 3.5;
baseline = R*(1-0.5*(a));
b = rand(10000,1);

colours = lines(3);
figure
tiledlayout(1,3)
nexttile(1)
h1 = plot(baseline.*sin(2*pi*b),baseline.*cos(2*pi*b),'.','Color',colours(1,:));
xlabel('$x$ [m]')
ylabel('$y$ [m]')
title('Baseline, $2\pi$')
axis square
nexttile(2)
h2 = plot(baseline.*sin(1.8*pi*b),baseline.*cos(1.8*pi*b),'.','Color',colours(2,:));
xlabel('$x$ [m]')
title('Void, $1.8\pi$')
axis square
nexttile(3)
h3 = plot(baseline.*sin(2.2*pi*b),baseline.*cos(2.2*pi*b),'.','Color',colours(3,:));
xlabel('$x$ [m]')
title('Overlap, $2.2\pi$')
axis square