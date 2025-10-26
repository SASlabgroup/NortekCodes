%function [h1 h2 h3] = plotSigAvg( filename );
% function to plot Nortek Signature "averaged" data, 
% assuming bottom-mounted uplooking
%
%   [h1 h2 h3] = plotSigAvg( filename );
%
% where the outputs are plot handles
%
% J. Thomson, 2015 (original),
%               2019 (modified for new exporting from deployment software)
%               Sep 2020 trim velocity profiles for surface reflections

clear all, close all

filename = 'CODA_S1A1_Sig_2020_avgd_avgd_1.mat';

load([filename])

doff = 0.75;  % height of a seaspider above the seabed
mincorr = 50;

z = double( doff + Config.Average_BlankingDistance + Config.Average_CellSize .* [1: Data.Average_NCells(1)] );


%% trim data to close to surface (or out of water) and low correlations

beamangle = 25;
maxrange = double(Data.Average_Pressure) * cosd(beamangle); % times series of max range
%trimcells = ones(size(maxrange))*(z-zoffset) > .95*maxrange*ones(size(z)) ;
trimcells = ones(size(maxrange))*(z-doff) > .95*maxrange*ones(size(z))  | Data.Average_CorBeam1 < mincorr | ...
        Data.Average_CorBeam2 < mincorr | Data.Average_CorBeam3 < mincorr | Data.Average_CorBeam4 < mincorr;
    
Data.Average_VelEast(trimcells) = NaN;
Data.Average_VelNorth(trimcells) = NaN;
Data.Average_VelUp1(trimcells) = NaN;
Data.Average_VelUp2(trimcells) = NaN;
for i=1:4
    eval(['Data.Average_AmpBeam' num2str(i) '(trimcells) = NaN;'])
    eval(['Data.Average_CorBeam' num2str(i) '(trimcells) = NaN;'])
end

save([filename '_trimmed'])

%% velocities
figure(1), clf

ax(1) = subplot(3,1,1);
pcolor(Data.Average_Time, z, Data.Average_VelEast'), 
shading flat, 
hold on
plot(Data.Average_Time,Data.Average_AltimeterPressure,'k-','linewidth',1)
datetick
ylabel(['Dist. above seabed [m]'])
caxis([-.5 .5])
colormap(gca,nawhimar);
cb = colorbar; cb.Label.String = 'East Vel [m/s]';
title(filename,'interp','none')

ax(2) = subplot(3,1,2);
pcolor(Data.Average_Time, z, Data.Average_VelNorth'),
shading flat,
hold on
plot(Data.Average_Time,Data.Average_AltimeterPressure,'k-','linewidth',1)
datetick
ylabel(['Dist. above seabed [m]'])
caxis([-.5 .5])
colormap(gca,nawhimar);
cb = colorbar; cb.Label.String = 'North Vel [m/s]';

ax(3) = subplot(3,1,3);
pcolor(Data.Average_Time, z, Data.Average_VelUp1'),
shading flat,
hold on
plot(Data.Average_Time,Data.Average_AltimeterPressure,'k-','linewidth',1)
datetick
ylabel(['Dist. above seabed [m]'])
caxis([-.5 .5])
colormap(gca,nawhimar);
cb = colorbar; cb.Label.String = 'Vertical Vel [m/s]';

  
linkaxes(ax,'x')
    
print('-dpng',[filename '_Velocities.png'])


%% Amplitudes
figure(2), clf

for i=1:4
ax(i) = subplot(4,1,i);
pcolor(Data.Average_Time, z, eval(['Data.Average_AmpBeam' num2str(i)])'), 
shading flat, 
hold on
plot(Data.Average_Time,Data.Average_AltimeterPressure,'k-','linewidth',1)
datetick
ylabel(['Dist. above seabed [m]'])
caxis([0 100])
cb = colorbar; cb.Label.String = 'Backscatter [dB]';
if i==1, title(filename,'interp','none'), else end
end

  
linkaxes(ax,'x')
    
print('-dpng',[filename '_Backscatter.png'])


%% Correlation
figure(3), clf

for i=1:4
ax(i) = subplot(4,1,i);
pcolor(Data.Average_Time, z, eval(['Data.Average_CorBeam' num2str(i)])'), 
shading flat, 
hold on
plot(Data.Average_Time,Data.Average_AltimeterPressure,'k-','linewidth',1)
datetick
ylabel(['Dist. above seabed [m]'])
caxis([0 100])
cb = colorbar; cb.Label.String = 'Correlation [%]';
if i==1, title(filename,'interp','none'), else end
end

  
linkaxes(ax,'x')
    
print('-dpng',[filename '_Correlation.png'])