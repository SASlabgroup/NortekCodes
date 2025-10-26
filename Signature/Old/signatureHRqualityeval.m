% script to explore quality of signature HR data,
% as offloaded from SWIFT v4 and read from binary files using readSWIFT_SIG.m
% J. Thomson, Jan 2018... adapted from Aug 2015 testing
%

% HR cell depths 
hulldepth = 0.2; % meters
z = hulldepth + burst.Blanking + burst.CellSize * [1:size(burst.VelocityData,2)];

% time step (note that burst.time is discretized to the sec, but data is 8 Hz)
dt = range(burst.time)./length(burst.time)*24*3600;

%% despiking (optional, doesn't have large effects)
for i = 1:length(z)
    [newu replacedindex] = func_despike_phasespace3d( burst.VelocityData(:,i), 0, 2 );
    burst.VelocityData(:,i) = newu;
end

%% simple motion correction (bobbing only)
w = gradient(burst.Pressure,dt);
burst.VelocityData = burst.VelocityData + w'*ones(1,length(z));


%% profiles of amp, corr, and velocity
figure(1), clf

subplot(1,3,1)
errorbar(z ,nanmean(burst.AmplitudeData), nanstd(burst.AmplitudeData))%,'linewidth','1.5')
set(gca,'CameraUpVector',[-6 0 0])
set(gca,'Fontsize',16,'fontweight','demi')
ylabel('Amp [-]')
xlabel('z [m]')
grid

subplot(1,3,2)
errorbar(z ,nanmean(burst.CorrelationData), nanstd(burst.CorrelationData))%,'linewidth','1.5')
set(gca,'CameraUpVector',[-6 0 0])
set(gca,'Fontsize',16,'fontweight','demi')
ylabel('Cor [%]')
grid 

subplot(1,3,3)
errorbar(z ,nanmean(burst.VelocityData), nanstd(burst.VelocityData))%,'linewidth','1.5')
set(gca,'CameraUpVector',[-6 0 0])
set(gca,'Fontsize',16,'fontweight','demi')
ylabel('w [m/s]')
grid 


%% histograms, by depth... look for wrapping
figure(2), clf
hist(burst.VelocityData,[-0.6:0.1:0.6])
set(gca,'Fontsize',16,'fontweight','demi')
colorbar
xlabel('w [m/s]')
title('HR velocity histogram, by range bin')
ylabel('Number of observations')
colorbar('peer',gca,'EastOutside','YTickLabel',linspace(min(z),max(z),5),'Ytick',[0:25:100],'ydir','reverse')


%% timeseries of a selected bin

binnum = 128;
figure(3), clf
plot(burst.VelocityData(:,binnum),'-x')
title(['Bin ' num2str(binnum)])
xlabel('Points at 8 Hz')
ylabel('w [m/s]')


%% spectra
clear allpsd

samplingfreq = 8;
window = 64;
noverlap = [];
nfft = [];

figure(4), clf
cmap = colormap;

for i=1:length(z),
    cind = round(64 * i/length(z));
    [psd f ] = pwelch(detrend(burst.VelocityData(:,i)),window,noverlap,nfft,samplingfreq);
    loglog(f,psd,'color',cmap(cind,:)), hold on
    %loglog(f,psd.*f.^(5/3),'color',cmap(cind,:)), hold on
    allpsd(:,i) = psd;    
end

xlabel('Frequency [Hz]')
ylabel('TKE [m^2/s^2/Hz]')
title('Signature HR Spectra at each range bin')

colorbar('peer',gca,'EastOutside','YTickLabel',linspace(min(z),max(z),5),'Ytick',[0:.25:1],'ydir','reverse')


%% structure functions

deltar = zeros(size(z));

figure(4), 
[tke , epsilon , residual, A, Aerror, N, Nerror] = dissipation_debug(burst.VelocityData', z, size(burst.VelocityData,1), 1, deltar);
%epsilon = epsilon./1024;

figure(6), clf
subplot(1,2,1)
semilogx(epsilon,z,'x','linewidth',2)
set(gca,'Fontsize',16,'fontweight','demi')
xlabel('\epsilon [m^2/s^3]')
set(gca,'ydir','reverse','ylim',[0 5.5])
set(gca,'xlim',[1e-6 1e-2])
grid

subplot(1,2,2)
plot(N,z,'x','linewidth',2)
set(gca,'Fontsize',16,'fontweight','demi')
set(gca,'ydir','reverse','ylim',[0 5.5])
xlabel('N [m^2/s^2]')
grid









