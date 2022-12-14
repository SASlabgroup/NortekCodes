% process xIMU data from Vector on TTM
% J. Thomson, 
 
addpath('ximu_matlab_library');     % include library
close all;                          % close all figures
clear;                              % clear all variables
clc;                                % clear the command terminal
 
%% Import and plot data

xIMUdata = xIMUdataClass('/Volumes/Data/ADMIRALTY/ADMIRALTY_Jun2012/TTM_xIMU_Jun2012/TTM_APLUWxIMU_Jun2012');
xIMUdata.Plot();

%% pull out just accelerometer data (already in TTM coordinates)

xIMU.fs = 1 ./ xIMUdata.CalInertialAndMagneticData.SamplePeriod;
accx = xIMUdata.CalInertialAndMagneticData.Accelerometer.X;
accy = xIMUdata.CalInertialAndMagneticData.Accelerometer.Y;
accz = xIMUdata.CalInertialAndMagneticData.Accelerometer.Z;

time = linspace( min(xIMUdata.DateTimeData.Serial), max(xIMUdata.DateTimeData.Serial), length(accx) );


%% trim to ADV only

ADVstart = datenum(2012,6,12,17,30,0);
first = find( time >= ADVstart, 1 );
time(1:first) = [];
accx(1:first) = [];
accy(1:first) = [];
accz(1:first) = [];

ADVend = datenum(2012,6,14,14,30,0); 
last = find( time > ADVend, 1 );
time(last:length(time)) = [];
accx(last:length(time)) = [];
accy(last:length(time)) = [];
accz(last:length(time)) = [];


%% parse into bursts, starting with first ADV time

burstlength = 300 % burst length, secs


bursttime = reshape(time(1:(length(time)-rem(length(time),xIMU.fs*burstlength))),xIMU.fs*burstlength,[]);
burstaccx = reshape(accx(1:(length(time)-rem(length(time),xIMU.fs*burstlength))),xIMU.fs*burstlength,[]);
burstaccy = reshape(accy(1:(length(time)-rem(length(time),xIMU.fs*burstlength))),xIMU.fs*burstlength,[]);
burstaccz = reshape(accz(1:(length(time)-rem(length(time),xIMU.fs*burstlength))),xIMU.fs*burstlength,[]);

xIMU.time = mean(bursttime);

%% burst means

xIMU.axbar = nanmean(burstaccx);
xIMU.aybar = nanmean(burstaccy);
xIMU.azbar = nanmean(burstaccz);


%% spectra

for i=1:length(xIMU.time),
    
    [psd f] = pwelch( detrend( burstaccx(:,i) ), [], [], [], xIMU.fs);
    psd(1) = []; f(1) = []; % remove mean
    xIMU.Saxf(:,i) = psd;  % power spectral density of acceleration mag, m^2 s^-2 Hz^-1
    
    [psd f] = pwelch( detrend( burstaccy(:,i) ), [], [], [], xIMU.fs);
    psd(1) = []; f(1) = []; % remove mean
    xIMU.Sayf(:,i) = psd;  % power spectral density of acc mag, m^2 s^-2 Hz^-1
    
    [psd f] = pwelch( detrend( burstaccz(:,i) ), [], [], [], xIMU.fs);
    psd(1) = []; f(1) = []; % remove mean
    xIMU.Sazf(:,i) = psd;  % power spectral density of acc mag, m^2 s^-2 Hz^-1
    
end

xIMU.f = f;

%% plots
loglog(f,nanmean(xIMU.Saxf'),f,nanmean(xIMU.Sayf'),f,nanmean(xIMU.Sazf'),'linewidth',3),
set(gca,'Fontsize',16,'Fontweight','demi')
legend('x','y','z')
ylabel('Acceleration spectra (means)'),xlabel('Frequency [Hz]')

print -dpng TTM_accspectra.png


%% 

save TTM_APLUWxIMU_Jun2012_processed xIMU

save('/Volumes/Data/ADMIRALTY/ADMIRALTY_Jun2012/TTM_xIMU_Jun2012/TTM_APLUWxIMU_Jun2012_burst.mat','burst*')

