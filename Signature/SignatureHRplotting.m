clear all, %close all, clc

fpath = './'
filename = 'S100882A002_LakeK_10'

load([fpath filename '.mat'])

%% plot raw data (BurstHR data from beam 5)

z = double(Config.BurstHR_BlankingDistance + Config.BurstHR_CellSize .* [1:Config.BurstHR_NCells] );

figure(1), clf

ax(1) = subplot(3,1,1);
pcolor(Data.BurstHR_Time-datenum(2018,0,0), z, double(Data.BurstHR_VelBeam5)')
%pcolor(double(Data.BurstHR_VelBeam5)')
shading flat
colorbar
datetick    
set(gca,'YDir','reverse')
ylabel('range [m]')
legend('Vel [m/s]')


ax(2) = subplot(3,1,2);
pcolor(Data.BurstHR_Time-datenum(2018,0,0),double(Config.BurstHR_BlankingDistance + Config.BurstHR_CellSize .* [1:Config.BurstHR_NCells] ),double(Data.BurstHR_CorBeam5)')
shading flat
colorbar
datetick
set(gca,'YDir','reverse')
ylabel('range [m]')
legend('Cor []')

ax(3) = subplot(3,1,3);
pcolor(Data.BurstHR_Time-datenum(2018,0,0),double(Config.BurstHR_BlankingDistance + Config.BurstHR_CellSize .* [1:Config.BurstHR_NCells] ),double(Data.BurstHR_AmpBeam5)')
shading flat
colorbar
datetick    
set(gca,'YDir','reverse')
ylabel('range [m]')
legend('Amp []')

linkaxes(ax,'x')

print('-dpng',[fpath filename '_raw.png'])

%% find good data
AvgAmp = mean(double(Data.BurstHR_AmpBeam5)');
AvgCor = mean(double(Data.BurstHR_CorBeam5)');
goodinds = find(AvgCor>85);


%% spectra

clear allpsd

samplingfreq = 8;
window = 128;
noverlap = [];
nfft = [];

figure(2), clf
colormap winter
cmap = colormap;

for i=1:length(z),
    cind = ceil(63 * i/length(z));
    [psd f ] = pwelch(detrend(double(Data.BurstHR_VelBeam5( goodinds , i ))),window,noverlap,nfft,samplingfreq);
    loglog(f,psd,'color',cmap(cind,:)), hold on
    allpsd(:,i) = psd;    
end

xlabel('Frequency [Hz]')
ylabel('TKE [m^2/s^2/Hz]')
title('Signature HR Spectra at each range bin [m]')

colorbar('peer',gca,'EastOutside','YTickLabel',linspace(min(z),max(z),7))

print('-dpng',[fpath filename '_spectra.png'])

%% plot other 4 beams also ("averaged data")

AvgCorBeam1 = mean(double(Data.Average_CorBeam1)');
AvgAmpBeam1 = mean(double(Data.Average_AmpBeam1)');
%goodinds = find(AvgCorBeam1>85 & AvgAmpBeam1 > 45);

zavg = double(Config.Average_BlankingDistance + Config.Average_CellSize .* double( [1:Config.Average_NCells] ) );

figure(3), clf

for i = 1:4,
    
    ax(i) = subplot(4,1,i); 
    
    pcolor(Data.Average_Time-datenum(2018,0,0), zavg, double(eval(['Data.Average_VelBeam' num2str(i) '']))' ), 
    shading flat, 
    datetick
    set(gca,'YDir','reverse')
    ylabel(['Vel Beam ' num2str(i)])
    %caxis([-2 2])
    colorbar,
    if i==1, title('Average mode data'), else end
   
end

print('-dpng',[fpath filename '_4beams_Avgvel.png'])



figure(4), clf

for i = 1:4,
    
    ax(i) = subplot(4,1,i); 
    
    pcolor(Data.Average_Time-datenum(2018,0,0), zavg, double(eval(['Data.Average_AmpBeam' num2str(i) '']))' ), 
    shading flat, 
    datetick
    set(gca,'YDir','reverse')
    ylabel(['Amp Beam ' num2str(i)])
    %caxis([-2 2])
    colorbar,
    if i==1, title('Average mode data'), else end

    
end

print('-dpng',[fpath filename '_4beams_Avgamp.png'])

