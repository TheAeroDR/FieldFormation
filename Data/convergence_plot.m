dat = readtable('all_convergence.csv');

N_values = unique(dat.N);

% Group data by N values
for i = 1:length(N_values)
    idx = dat.N == N_values(i); % Get indices f456104531764.0or each unique N
    Bx_groups(:,i) = padarray(dat.Bx(idx),1000-length(dat.Bx(idx)),NaN,'post');
    By_groups(:,i) = padarray(dat.By(idx),1000-length(dat.Bx(idx)),NaN,'post');
    Bz_groups(:,i) = padarray(dat.Bz(idx),1000-length(dat.Bx(idx)),NaN,'post');
    Ex_groups(:,i) = padarray(dat.Ex(idx),1000-length(dat.Bx(idx)),NaN,'post');
    Ey_groups(:,i) = padarray(dat.Ey(idx),1000-length(dat.Bx(idx)),NaN,'post');
    Ez_groups(:,i) = padarray(dat.Ez(idx),1000-length(dat.Bx(idx)),NaN,'post');
end

% Create box plots
figure;
%subplot(2,3,1);
iosr.statistics.boxPlot(N_values, Bx_groups,'scaleWidth',true,'notch',true,'xspacing','equal');
xticklabels({'$10^3$','$10^4$','$10^5$','$10^6$','$10^7$'})
ylabel('$B_x$ [T]')
xlabel('N')
%subplot(2,3,2);
figure
iosr.statistics.boxPlot(N_values, By_groups,'scaleWidth',true,'notch',true,'xspacing','equal');
xticklabels({'$10^3$','$10^4$','$10^5$','$10^6$','$10^7$'})
ylabel('$B_y$ [T]')
xlabel('N')
%subplot(2,3,3);
figure
iosr.statistics.boxPlot(N_values, Bz_groups,'scaleWidth',true,'notch',true,'xspacing','equal');
xticklabels({'$10^3$','$10^4$','$10^5$','$10^6$','$10^7$'})
ylabel('$B_z$ [T]')
xlabel('N')
%subplot(2,3,4);
figure
iosr.statistics.boxPlot(N_values, Ex_groups,'scaleWidth',true,'notch',true,'xspacing','equal');
xticklabels({'$10^3$','$10^4$','$10^5$','$10^6$','$10^7$'})
ylabel('$E_x$ [V/m]')
xlabel('N')
%subplot(2,3,5);
figure
iosr.statistics.boxPlot(N_values, Ey_groups,'scaleWidth',true,'notch',true,'xspacing','equal');
xticklabels({'$10^3$','$10^4$','$10^5$','$10^6$','$10^7$'})
ylabel('$E_y$ [V/m]')
xlabel('N')
%subplot(2,3,6);
figure
iosr.statistics.boxPlot(N_values, Ez_groups,'scaleWidth',true,'notch',true,'xspacing','equal');
xticklabels({'$10^3$','$10^4$','$10^5$','$10^6$','$10^7$'})
ylabel('$E_z$ [V/m]')
xlabel('N')