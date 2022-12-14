% matlab script to read and process Nortek Signature data
% collected on drifting Doppcats at ORPC site in Igiugig, AK
% 
%   J. Thomson, March 2015
%   Maricarmen, July 2015

clear all, close all, clc

readraw = 1;

date = '9JUL2015' %'2_22AUG2014'

addpath(genpath('/Users/Maru/Documents/PhD_UW/TurbinaAlaska/Data2015/'))
%cd('/Users/jthomson/Desktop/Iguigig/Signature')
cd('/Users/Maru/Documents/PhD_UW/TurbinaAlaska/Data2015')
fpath = ['ORPC_' date '/']; 
%fname = 'ORPCSIG_2_22AUG2014.234.14.AD2CPRaw.00000.mat';
sigfile = dir([fpath 'SIG_*.mat']);
NSig=length(sigfile);

% Signature file to process
kfile=2;

% GPS files
GPSfile = dir([fpath '*GPS*.csv']);
Nfiles = length(GPSfile);


%% read signature data and QC

load([ fpath sigfile(kfile).name ])

mincor = 30;

Data.Burst_VelBeam1( Data.Burst_CorBeam1 < 30 ) = NaN;
Data.Burst_VelBeam2( Data.Burst_CorBeam2 < 30 ) = NaN;
Data.Burst_VelBeam3( Data.Burst_CorBeam3 < 30 ) = NaN;
Data.Burst_VelBeam4( Data.Burst_CorBeam4 < 30 ) = NaN;
Data.IBurst_VelBeam5( Data.IBurst_CorBeam5 < 30 ) = NaN;

%% read GPS 

GPSmat = dir([ fpath '*GPSdata_all*.mat']);



    
    if isempty(GPSmat) | readraw==1,
        
        for j=1:Nfiles
            
        [gi track utcdate utctime localdate localtime miliseconds valid lat northsouth lon eastwest altitude speed GPSheading distance] = textread([fpath GPSfile(2).name],'%n%s%s%s%s%s%n%s%n%s%n%s%n%n%n%n','delimiter',',','headerlines',1,'bufsize',1600000);
        % note that this will fail if the Qstarz restarted and has a line
        % of header info somewhere in the middle of the file
        
        for gi = 1:length(lat),
            gpstime(gi) = datenum([char(localdate(gi)) ' ' char(localtime(gi)) ] );%datenum([char(utcdate(gi)) ' ' char(utctime(gi)) ] );
        end
        
        gpstime_all(:,j) = gpstime;% + ( UTCtolocal /24);  % adjust to local time stamp
        valid_all{:,j}=valid;
        lat_all(:,j)=lat;
        lon_all(:,j)=lon;
        %lon = -lon;
        speed = speed * .2777; % convert from Km/h to m/s
        
        speed_all(:,j)=speed;
        
        end
        %save([ fpath GPSfile(j).name '.mat'])%,'gpstime','lat','lon','speed','heading','date')
        save([fpath 'GPSdata_all_' date '.mat'])
    
    else
        save temp3
        load([ fpath GPSmat.name ])
        load temp3
    end
    
    
    
    
    %% change to local coordinates
    % positive x is upstream of turbine
    for j = 1:Nfiles
    [ xaux, yaux ] = IgiugigCoorTransform_2014(lat_all(:,j),lon_all(:,j));
    x(:,j)=xaux;
    y(:,j)=yaux;
    end
    
    %% Check valid GPS measurements
    s1='Estimated (dead reckoning)';
    stcomp1=strcmp(s1,valid_all{1});
    stcomp2=strcmp(s1,valid_all{2});
    nanGPS1=find(stcomp1==1);
    nanGPS2=find(stcomp2==1);
    
    x(nanGPS1,1)=NaN;
    y(nanGPS1,1)=NaN;
    x(nanGPS2,2)=NaN;
    y(nanGPS2,2)=NaN;
    speed_all(nanGPS1,1)=NaN;
    speed_all(nanGPS2,2)=NaN;
    %gpstime_all(nanGPS1,1)=NaN;
    %gpstime_all(nanGPS2,2)=NaN;
    
    %% match positions and speed over ground (sog)
    
    pts = length(Data.Burst_MatlabTimeStamp);
    timeshift = 0 % Both are in local time JUL2015;  8/24 - 5/3600/24;  % signature was recorded in local time %in days
    
    for i = 1:pts,
        
        for j=1:Nfiles
        [tdiff matchindex] = min( abs ( Data.Burst_MatlabTimeStamp(i) + timeshift - gpstime_all(:,j) ) );
        
        if tdiff < 1/60/60/24, %1 segundo del dia
            
            Data.x(i,j) = x(matchindex,j);
            Data.y(i,j) = y(matchindex,j);
            Data.sog(i,j) = speed_all(matchindex,j);
            Data.validGPS(i,j)=valid_all{j}(matchindex); %Check 
            
        else
            Data.x(i,j) = NaN;
            Data.y(i,j) = NaN;
            Data.sog(i,j) = NaN;
            
        end
        end
    end




%% Average two GPS data
Data.x_av=nanmean(Data.x,2);
Data.y_av=nanmean(Data.y,2);
Data.sog_av=nanmean(Data.sog,2);


%% Get distance between two measurements
Data.dist=sqrt((Data.x(:,1)-Data.x(:,2)).^2+(Data.y(:,1)-Data.y(:,2)).^2);


%% plot drift tracks

figure(1), clf

plot([-6 10],[0 0],'k-','linewidth',6),  hold on
scatter(Data.y_av, Data.x_av,  5*ones(pts,1), Data.sog_av),
axis([ -50 50 -200 200])
caxis([0 3])
colorbar


%% plot raw data as tracks

figure(2), clf

hubindex = 5;

relativevel = sqrt( Data.Burst_VelBeam1(:,hubindex).^2 + Data.Burst_VelBeam2(:,hubindex).^2 + Data.Burst_VelBeam3(:,hubindex).^2 + Data.Burst_VelBeam4(:,hubindex).^2 );
%relativevel = Data.IBurst_VelBeam5(:,hubindex);


     plot([-6 10],[0 0],'k-','linewidth',6),      hold on
     scatter(Data.y_av, Data.x_av,  5*ones(pts,1), relativevel ./ Data.sog_av ), 
     %scatter(Data.x' * ones(1,10), ones(pts,1) * double(Data.Burst_Range), ones(pts,10), double(eval(['Data.Burst_VelBeam' num2str(i)])) ), 
     axis([ -50 50 -200 200])
     caxis([0 0.5])
     colorbar

print('-dpng',[fpath 'RelativeVelocityMap_' date '_' int2str(kfile) '.png']) 
     
%%

save([fpath 'Signature_merged_' date '_' int2str(kfile) '.mat'])

%cd('/Users/jthomson/Dropbox/ORPC')

%% Battery

figure(3)
plot(Data.Burst_MatlabTimeStamp-datenum(2015,0,0),Data.Burst_Battery)


%% Get data only from transects


%% old stuff

% 
% figure(1), clf
% 
% for i = 1:4,
%     
%     ax(i) = subplot(5,1,i); 
%     
%     pcolor(Data.Burst_MatlabTimeStamp-datenum(2014,0,0), double(Data.Burst_Range), double(eval(['Data.Burst_VelBeam' num2str(i)]))' ), 
%     shading flat, 
%     datetick
%     set(gca,'YDir','reverse')
%     ylabel(['Vel Beam ' num2str(i)])
%     caxis([-2 2])
%     colorbar,
%    
%     
% end
% 
% ax(5) = subplot(5,1,5);
% pcolor(Data.IBurst_MatlabTimeStamp-datenum(2014,0,0), double(Data.IBurst_Range), double(Data.IBurst_VelBeam5)' ), 
%     shading flat, 
%     datetick
%     set(gca,'YDir','reverse')
%     ylabel(['Vel Beam 5'])
%     caxis([-2 2])
%     colorbar,
% 
% linkaxes(ax,'x')

% %% integral length or time scales
% 
% w = double(Data.IBurst_VelBeam5);
% 
% n = 3824; % = Config.AD2CP_burst_nSamples; % 4096 was intent but is it 3824?
% m = floor(length(w)./n);
% dz = median(diff(double(Data.IBurst_Range)));
% 
% for wi=1:(m*n),
%     lz(wi) = sonic_integral_timescale(w(wi,:),1./dz);
% end
