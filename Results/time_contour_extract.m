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

load('timeset_10ed_200Hz.mat');

%%
grid_size_x = length(unique(x_all(:,1)));
grid_size_y = length(unique(y_all(:,1)));
grid_size_y = length(unique(z_all(:,1)));

%%
x_reshaped = reshape(x_all(:,1), grid_size_x, grid_size_y)';
x_unique = x_reshaped(1,:);
y_reshaped = reshape(y_all(:,1), grid_size_x, grid_size_y)'; 
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
% %%
% for i = 1:size(Bz,3)
%     Bz_plim = prctile(Bz_all(:,i), 99.999);
%     Bz_mlim = prctile(Bz_all(:,i), 0.001);
%     Ez_plim = prctile(Ez_all(:,i), 99.999);
%     Ez_mlim = prctile(Ez_all(:,i), 0.001);
% 
% 
%     Bz_plim = 1e-8;
%     Bz_mlim = -1e-9;
%     Ez_plim = prctile(Ez_all(:,i), 99.999);
%     Ez_mlim = prctile(Ez_all(:,i), 0.001);
% 
%     Bz(:,:,i) = min(max(Bz(:,:,i), Bz_mlim), Bz_plim);
%     Ez(:,:,i) = min(max(Ez(:,:,i), Ez_mlim), Ez_plim);
% end
%%
figure;
bz_min = min(Bz(:));
bz_max = max(Bz(:));

bz_min = -5e-10;
bz_max = 5e-10;
while true

    for i = 1:size(Bz,3);
        surf(X, Y, Bz(:,:,i), 'EdgeColor', 'none');
        colormap jet;
        colorbar;
        title(['Bz at frame ', num2str(i)]);
        xlabel('X'); ylabel('Y'); zlabel('Bz');
        axis tight;
        zlim([bz_min, bz_max]);
        %xlim([-101,-80]);
        %ylim([-10,10]);
        caxis([bz_min, bz_max]);
        view(2);
        pause(0.1);
    end
end

%%

S = squeeze(Bz(511,91,:));

S = S - mean(S);

Fs = 20;                
T = 1/Fs;      
L = length(S);
f = Fs/L*(0:(L/2));
t = (0:L-1)*T;

figure
plot(t,S);

Y = fft(S);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure
plot(f,P1,"LineWidth",3) 
title("Single-Sided Amplitude Spectrum of S(t)")
xlabel("f (Hz)")
ylabel("|P1(f)|")
% %%
% clear data
% sin_fs = [3.2,3.4,5,6,2,3.2,3.6,3.25,3.679];
% sin_amps = [0.8,1,1,5,1,1,3,1,1];
% t  = 0:5e-3:10;
% for i=1:length(sin_fs)
%     data(i,:) = sin_amps(i)*sin(2*pi*sin_fs(i)*t);
% end
% M = sum(data,1);
% figure
% plot(t,M);
% 
% Fs = 200;                
% T = 1/Fs;      
% L = length(M);
% f = Fs/L*(0:(L/2));
% 
% Y = fft(M);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% 
% figure
% plot(f,P1,"LineWidth",3) 
% title("Single-Sided Amplitude Spectrum of S(t)")
% xlabel("f (Hz)")
% ylabel("|P1(f)|")

%%
load("timeset_10ed_200Hz_processed_L1.mat")
%%
for i = 1:size(Bz,1)
    for j = 1:size(Bz,2)
        S = squeeze(Bz(i,j,:));
        DC_comp = mean(S);
        S = S - DC_comp;
        Bz_DC(i,j) = DC_comp;
        Bz_AC(i,j,:) = S;

    end
end

figure;
bz_min = min(Bz_DC(:));
bz_max = max(Bz_DC(:));
[X,Y] = meshgrid(x_unique,y_unique);
surf(X, Y, Bz_DC, 'EdgeColor', 'none');
colormap jet;
colorbar;
title(['DC Bz']);
xlabel('X'); ylabel('Y'); zlabel('Bz');
axis tight;
zlim([min(Bz_DC(:)), max(Bz_DC(:))]);
caxis([bz_min, bz_max]);
view(2);

for i = 1:size(Bz_AC,3)
    %Bz_plim = prctile(Bz_AC(:,i), 99.999);
    %Bz_mlim = prctile(Bz_AC(:,i), 0.001);

    Bz_plim = 8e-10;
    Bz_mlim = -8e-10;

    Bz_AC(:,:,i) = min(max(Bz_AC(:,:,i), Bz_mlim), Bz_plim);
end
%%
figure;
bz_min = min(Bz_AC(:));
bz_max = max(Bz_AC(:));

[X,Y] = meshgrid(x_unique,y_unique);
while true
    for i = 1:size(Bz_AC,3)
        surf(X, Y, Bz_AC(:,:,i), 'EdgeColor', 'none');
        colormap jet;
        colorbar;
        title(['AC Bz at frame ', num2str(i)]);
        xlabel('X'); ylabel('Y'); zlabel('Bz');
        axis tight;
        %axis square;
        zlim([min(Bz_AC(:)), max(Bz_AC(:))]);
        caxis([bz_min, bz_max]);
        view(2);
        pause(0.1);
    end
end


%%
t = 0:5e-3:(size(Bz,3)-1)*5e-3;
x = 2;
y = 3;
[~,idx] = min(abs(x_unique-x));
[~,idy] = min(abs(y_unique-y));
plot(t,squeeze(Bz_AC(idx,idy,:)))

%%
Bz_rms = rms(Bz_AC,3);

figure;
bz_min = min(Bz_rms(:));
bz_max = max(Bz_rms(:));
[X,Y] = meshgrid(x_unique,y_unique);
surf(X, Y, Bz_rms, 'EdgeColor', 'none');
colormap jet;
colorbar;
title(['RMS Bz']);
xlabel('X'); ylabel('Y'); zlabel('Bz');
axis tight;
zlim([min(Bz_rms(:)), max(Bz_rms(:))]);
caxis([bz_min, bz_max]);
view(2);

transect = 0;

[~,idt] = min(abs(y_unique-transect))

Bz_transect = Bz_rms(idt,:);

figure
plot(x_unique,Bz_transect)

%%
for i = 1:size(Ez,1)
    for j = 1:size(Ez,2)
        S = squeeze(Ez(i,j,:));
        DC_comp = mean(S);
        S = S - DC_comp;
        Ez_DC(i,j) = DC_comp;
        Ez_AC(i,j,:) = S;

    end
end
%%
transect = 0;

[~,idt] = min(abs(x_unique-transect))

Ez_transect = Ez(:,idt);

Ez_transect(Ez_transect<-4300) = -4300
figure
plot(y_unique,Ez_transect)

%%
temp_min = -1e6;
temp_max = 100;
while true
    for i = 1:size(Ez,3)
        surf(X, Y, Ez(:,:,i), 'EdgeColor', 'none');
        colormap jet;
        colorbar;
        title(['Voltage at frame ', num2str(i)]);
        xlabel('X'); ylabel('Y'); zlabel('V');
        axis tight;
        %axis square;
        zlim([min(Ez(:)), max(Ez(:))]);
        caxis([temp_min, temp_max]);
        view(2);
        pause(0.1);
    end
end
%%
temp_rms = rms(temp,3);

figure;
temp_min = min(temp_rms(:));
temp_max = max(temp_rms(:));
[X,Y] = meshgrid(x_unique,y_unique);
surf(X, Y, temp_rms, 'EdgeColor', 'none');
colormap jet;
colorbar;
title(['RMS Bz']);
xlabel('X'); ylabel('Y'); zlabel('Bz');
axis tight;
zlim([temp_min, temp_max]);
caxis([temp_min, temp_max]);
view(2);

transect = 0;

[~,idt] = min(abs(y_unique-transect))

temp_transect = temp_rms(idt,:);

figure
plot(x_unique,temp_transect)