clear all
load('Raack-deep.mat')
%%
figure(1)
tiledlayout(2,3)
nexttile([1,3])
plot(mean_data,data.height,'Color',[0 0.4470 0.7410],'Marker','square','MarkerFaceColor',[0 0.4470 0.7410])
hold on
plot(median_data,data.height,'Color',[0.8500 0.3250 0.0980],'Marker','^','MarkerFaceColor',[0.8500 0.3250 0.0980])
y_full = linspace(0,0.2,101);
x_full = 81*exp(-((y_full)/0.03304).^2)+4;
plot(x_full,38 * y_full)
legend('Mean Value $[\mu m]$','Median Value $[\mu m]$','Fit $D =  81 \exp\left[-\left(\frac{H}{0.03304}\right)^2\right] + 4$')
ylabel("Height $[m]$")
xlabel("Diameter $[\mu m]$")

nexttile(4)
plot(100 * data.("<2um")/sum(data.("<2um")),data.height,'Color',[0 0.4470 0.7410],'Marker','square','MarkerFaceColor',[0 0.4470 0.7410])
hold on
plot(100 * data.("2-4um")/sum(data.("2-4um")),data.height,'Color',[0.8500 0.3250 0.0980],'Marker','^','MarkerFaceColor',[0.8500 0.3250 0.0980])
legend('$<2 \mu m$', '$2\leq x<4 \mu m$')
ylabel("Height $[m]$")
xlabel("Percentage $[\%]$")

nexttile(5)
plot(100 * data.("4-8um")/sum(data.("4-8um")),data.height,'Color',[0 0.4470 0.7410],'Marker','square','MarkerFaceColor',[0 0.4470 0.7410])
hold on
plot(100 * data.("8-16um")/sum(data.("8-16um")),data.height,'Color',[0.8500 0.3250 0.0980],'Marker','^','MarkerFaceColor',[0.8500 0.3250 0.0980])
plot(100 * data.("16-31um")/sum(data.("16-31um")),data.height,'Color',[0.9290 0.6940 0.1250],'Marker','o','MarkerFaceColor',[0.9290 0.6940 0.1250])
plot(100 * data.("31-63um")/sum(data.("31-63um")),data.height,'Color',[0.4940 0.1840 0.5560],'Marker','diamond','MarkerFaceColor',[0.4940 0.1840 0.5560])
legend('$4\leq x < 8 \mu m$', '$8\leq x<16 \mu m$', '$16\leq x<31 \mu m$','$31\leq x<63 \mu m$')
ylabel("Height $[m]$")
xlabel("Percentage $[\%]$")

nexttile(6)
plot(100 * data.("63-125um")/sum(data.("63-125um")),data.height,'Color',[0 0.4470 0.7410],'Marker','square','MarkerFaceColor',[0 0.4470 0.7410])
hold on
plot(100 * data.("125-250um")/sum(data.("125-250um")),data.height,'Color',[0.8500 0.3250 0.0980],'Marker','^','MarkerFaceColor',[0.8500 0.3250 0.0980])
plot(100 * data.("250-500um")/sum(data.("250-500um")),data.height,'Color',[0.9290 0.6940 0.1250],'Marker','o','MarkerFaceColor',[0.9290 0.6940 0.1250])
legend('$63\leq x < 125 \mu m$', '$125\leq x<250 \mu m$', '$250\leq x<500 \mu m$')
ylabel("Height $[m]$")
xlabel("Percentage $[\%]$")

%%
bins = {'<2', '2 \leq x < 4', '4 \leq x < 8', '8 \leq x < 16', '16 \leq x < 31', '31 \leq x < 63', '63 \leq x < 125', '125 \leq x < 250', '250 \leq x < 500'};
bin_widths = [2,2,4,8,15,32,62,125,250];
f_data = table2array(data);
f_data = f_data(:,2:end);
f_data = f_data./sum(f_data,2);
%f_data(f_data<0.09)=0;
figure(2)
tiledlayout(2,1)
nexttile(1)
h = heatmap(bins, flip(y), flip(f_data));
h.CellLabelFormat = '%.2f';
xlabel("Binned Diameter [\mu{}m]");
ylabel("Height [m]");
colormap('jet');
ylabel(struct(h).Colorbar,"Percentage of particulate at height")
%% 
nexttile(2)
markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
semilogx(x,ext(1,:),'LineWidth',0.5,'LineStyle','--','Marker',markers{1},'MarkerSize',20)
hold on
for i = 2:length(y)
    semilogx(x,ext(i,:),'LineWidth',0.5,'LineStyle','--','Marker',markers{i},'MarkerSize',10)
end
xlabel('Diameter $[\mu m]$')
ylabel('Number of Particles $[-]$')
l = legend(num2str(y),'Location','eastoutside');
 %%
for i = 1:2
    if i == 2
        i = 9;
    end
    freq = zeros(9);
    freq(1,:) = ([0,0,0,0,ext(1,5:end)]);
    freq(9,:) = ([ext(9,1:5),0,0,0,0]);
    x_full = repelem(x,ext(i,:));
    
    a{i} = fitdist(log(x)','Normal','Frequency',freq(i,:));
    b{i} = fitdist(x','Lognormal','Frequency',freq(i,:));

end
normal_func = @(params, x) ((1/(params(1) .* sqrt(2 .* pi))) .* exp(-0.5 .* ((x - params(2))/params(1)).^2));
    
log_normal_func = @(params, x) ((1./(params(1) .* x .* sqrt(2 .* pi))) .* exp(-1 .* (((log(x) - params(2)).^2) ./(2 .* params(1).^2))));
%%
figure(3)
x_fit = logspace(-4,3,100001);
height = linspace(0,5,51);
low_params=a{1};
high_params = a{9};
pDisFunc = @(params, x) (1./x).*((normal_func([low_params.sigma,low_params.mu],log(x))*(1-Heavi_alt(params,2))) + (normal_func([high_params.sigma,high_params.mu],log(x))*(Heavi_alt(params,2))));

for i = 1:length(height)
    h = height(i);
    pDis(i,:) = pDisFunc(h,x_fit);
end
%plot3(x_fit,ones(length(x_fit),1)*height,pDisFunc)

[X, H] = meshgrid(x_fit, height);
surf(X, H, pDis,'LineStyle','none','FaceColor','interp');

set(gca,'XScale','log')
xlabel('Diameter $[\mu m]$')
ylabel('Height $[m]$')
zlabel('PDF $[-]$')

%%
nexttile(2)
%%
markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
% plot3(x,ones(length(x),1)*data.height(1),ext(1,:)./sum(ext(1,:)),'LineWidth',0.5,'LineStyle','--','Marker',markers{1},'MarkerSize',20)
% hold on
% for i = 2:length(y)
% end
heights = data.height';
bins_low = [0,2,4,8,16,31,63,125,250];
bins_high = [2,4,8,16,31,63,125,250,500];
for i = 1:length(heights)
    for j = 1:length(bins_low)
        x1 = bins_low(j);
        x2 = bins_high(j);
        idxl = find(x_fit>=x1,1);
        idxh = find(x_fit>=x2,1);
        p(i,j) = trapz((x_fit(idxl:idxh)),pDisFunc(heights(i),x_fit(idxl:idxh)));
    end
    plot3(x,ones(length(x),1)*heights(i),ext(i,:)./sum(ext(i,:)),'LineWidth',0.5,'LineStyle','--','Marker',markers{i},'MarkerSize',10)
    hold on
    %plot3(x,ones(length(x),1)*heights(i),p(i,:))
end
[X, H] = meshgrid(x, heights);
surf(X, H, p,'LineStyle','none','FaceColor','interp','FaceAlpha',0.5);

set(gca,'XScale','log')
xlabel('Diameter $[\mu m]$')
ylabel('Height $[m]$')
zlabel('Probability $[-]$')
%l = legend(num2str(y),'Location','eastoutside');

%%
function y = Heavi_alt(x,a)
    y = 0.5*tanh(x-a) +1/2;
end
function y = Heavi(x,a)
    if x - a < 0
        y = 0;
    elseif x - a > 0
        y = 1;
    else
        y = 0.5;
    end
end