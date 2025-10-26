% script for processing Signature1000 HR on SWIFT v4beta
% assuming beam 5 is in HR mode and beams 1-4 are in broadband
% J. Thomson Aug 2015, beam 5 testing
%   added beams 1-4 (Jun 2016)
%   changed to 'BurstHRHR' notation to match new firmware (Aug 2016)
%   and skip transoform of average profiles (settings now ENU on instrument)


clear all, close all

si = 0;  % initialize indexing for SWIFT structure

plotflag = true;  % turn on and off raw data plots
strplotflag = false; % turn on and off structure function plots
timestep = 600; % BurstHR length in seconds (all results based on this)
maxQCratio = 0.3; % max ratio of removed points to total points allowed
mincor = 75;  % minimum correlation allowed, default is 75
minamp = 30;  % minimum amplitude allowed, default is 30

parentdir = './'%/Users/jthomson/Dropbox/SWIFT_v4.x/NortekSignature/SWIFTv4b_Sig_04Aug2016';
cd([parentdir])
wd = pwd;
wdi = find(wd == '/',1,'last');
wd = wd((wdi+1):length(wd));

flist = dir(['*ad2cp*.mat']);

for fi = 1:length(flist),
    
load([ flist(fi).name ])


%% Quality control the HR (beam 5) data... and the slant beams too

Data.BurstHR_VelBeam5(Data.BurstHR_CorBeam5<mincor) = NaN;

Data.BurstHR_VelBeam5(Data.BurstHR_AmpBeam5<minamp) = NaN;

bad = (Data.Average_CorBeam1<mincor | Data.Average_CorBeam2<mincor | Data.Average_CorBeam3<mincor | Data.Average_CorBeam4<mincor);
Data.Average_VelEast(bad) = NaN;
Data.Average_VelNorth(bad) = NaN;
Data.Average_VelUp1(bad) = NaN;
Data.Average_VelUp2(bad) = NaN;

bad = (Data.Average_AmpBeam1<minamp | Data.Average_AmpBeam2<minamp | Data.Average_AmpBeam3<minamp | Data.Average_AmpBeam4<minamp);
Data.Average_VelEast(bad) = NaN;
Data.Average_VelNorth(bad) = NaN;
Data.Average_VelUp1(bad) = NaN;
Data.Average_VelUp2(bad) = NaN;


%% geometry and motion correction options for HR processing,

%z5 = Data.BurstHR_Blanking(1)+Data.BurstHR_CellSize(1).*[1:Config.bursthr_nCells]; %-Data.BurstHR_CellSize(1);
z5 = Data.BurstHR_Range; % prefer to use range field (Aug 2016)

% vertical motion from pressure
platformw = gradient(double(Data.BurstHR_Pressure), 1./Config.burst_sampleRate);
%Data.BurstHR_VelBeam5 = Data.BurstHR_VelBeam5 - ( platformw * ones(1,length(z)) );  % doesn't make much difference

% range shifting (smearing) to be applied in call to dissipation.m
%deltar = zeros(length(z5),1); % null option
deltar = std(double(Data.BurstHR_Pressure)) * ones(1,length(z5)); % doesn't make much difference

%% geometry and coordinate transform for slant beams (skip as of Aug 2016, set ENU on instrument)

%[ Data, Config, T_beam2xyz ] = signatureAD2CP_beam2xyz_enu( Data, Config, 'avg', 0 );

%z1234 = Data.Average_Blanking(1)+Data.Average_CellSize(1).*[1:Config.avg_nCells]; %-Data.Average_CellSize(1);
z1234 = Data.Average_Range; % prefer as of Aug 2016;

%% loop thru BurstHRs, processing only if enough good data

%BurstHRtimes = [ Data.BurstHR_MatlabTimeStamp(1) : timestep/(24*3600) : Data.BurstHR_MatlabTimeStamp(end) ]; % could round the first one to make this cleaner
BurstHRtimes = Data.BurstHR_MatlabTimeStamp( Data.BurstHR_EnsembleCount == 1); 

for bi = 1:length(BurstHRtimes);
    
    si = si + 1;
    
    BurstHRindices = find( Data.BurstHR_MatlabTimeStamp > BurstHRtimes(bi) & Data.BurstHR_MatlabTimeStamp < (BurstHRtimes(bi) + timestep/(24*3600)) );
    avgindices = find( Data.Average_MatlabTimeStamp > BurstHRtimes(bi) & Data.Average_MatlabTimeStamp < (BurstHRtimes(bi) + timestep/(24*3600)) );
    
    SWIFT(si).time = BurstHRtimes(bi);
    
    QCratio = sum(sum(isnan(Data.BurstHR_VelBeam5(BurstHRindices,:)))) ./ numel(Data.BurstHR_VelBeam5(BurstHRindices,:));
    
    if QCratio < maxQCratio,  % proceed
        
        good(si) = true;
        
        % dissipation along beam5 from structure function
        [ tke epsilon ] = dissipation(double(Data.BurstHR_VelBeam5(BurstHRindices,:))', z5, length(BurstHRindices), strplotflag, deltar);
        epsilon( epsilon ==0 ) = NaN;
        SWIFT(si).signature.HRprofile.tkedissipationrate = epsilon';
        SWIFT(si).signature.HRprofile.z = z5;  
        
        % could impliment dissipation along the slant beams too, but must be careful of applications with limiting lengths scales
        
        % current profiles from the converted slant beams
        % these relative components will still need to be adjusted to absolute using the GPS drift components
        SWIFT(si).signature.profile.east = nanmean( double(Data.Average_VelEast(avgindices,:) ) );     
        SWIFT(si).signature.profile.north = nanmean( double(Data.Average_VelNorth(avgindices,:) ) );
        SWIFT(si).signature.profile.up = nanmean( double(Data.Average_VelUp1(avgindices,:) ) );
        SWIFT(si).signature.profile.z = z1234;

    else
        
        good(si) = false;
        % assign NaNs if bad BurstHR (usually out of the water)
        SWIFT(si).signature.HRprofile.tkedissipationrate = NaN(1,length(z5));
        SWIFT(si).signature.HRprofile.z = NaN(1,length(z5));
        
        SWIFT(si).signature.profile.east = NaN(1,length(z1234));     
        SWIFT(si).signature.profile.north = NaN(1,length(z1234));
        SWIFT(si).signature.profile.up = NaN(1,length(z1234));
        SWIFT(si).signature.profile.z = NaN(1,length(z1234));
    
    end 

end % close BurstHR loop


%%%%%% PLOTTING %%%%%%%

if plotflag==true, 

    
%% beam 5 HR data

figure(1), 

ax(1) = subplot(3,1,1);
pcolor(Data.BurstHR_MatlabTimeStamp-datenum(2015,0,0),double(Data.BurstHR_Range),double(Data.BurstHR_VelBeam5)')
hold on
if sum(good)~= 0, plot([SWIFT(good).time]-datenum(2015,0,0),min(z5),'ks','linewidth',5), else end
%pcolor(double(Data.BurstHR_VelBeam5)')
set(gca,'YDir','reverse')
shading flat
colorbar
datetick
ylabel('range [m]')
legend('Vel [m/s]')

ax(2) = subplot(3,1,2);
pcolor(Data.BurstHR_MatlabTimeStamp-datenum(2015,0,0),double(Data.BurstHR_Range),double(Data.BurstHR_CorBeam5)')
hold on
shading flat
set(gca,'YDir','reverse')
colorbar
datetick
ylabel('range [m]')
legend('Cor []')

ax(3) = subplot(3,1,3);
pcolor(Data.BurstHR_MatlabTimeStamp-datenum(2015,0,0),double(Data.BurstHR_Range),double(Data.BurstHR_AmpBeam5)')
hold on
shading flat
set(gca,'YDir','reverse')
colorbar
datetick
ylabel('range [m]')
legend('Amp []')

linkaxes(ax,'x')

print('-dpng',[ wd '_rawHRdata.png'])

%% beams 1-4 raw data

figure(2),  % broadband velocities (raw from slant beams)
ax(1) = subplot(4,1,1);
pcolor(Data.Average_MatlabTimeStamp-datenum(2015,0,0), double(Data.Average_Range), double(Data.Average_VelEast)' ),
hold on
if sum(good)~= 0, plot([SWIFT(good).time]-datenum(2015,0,0),min(z1234),'ks','linewidth',5), else end
shading flat,
datetick
set(gca,'YDir','reverse')
ylabel('East [m/s]')
caxis([-2 2])
colorbar,

ax(2) = subplot(4,1,2);
pcolor(Data.Average_MatlabTimeStamp-datenum(2015,0,0), double(Data.Average_Range), double(Data.Average_VelNorth)' ),
hold on
if sum(good)~= 0, plot([SWIFT(good).time]-datenum(2015,0,0),min(z1234),'ks','linewidth',5), else end
shading flat,
datetick
set(gca,'YDir','reverse')
ylabel('North [m/s]')
caxis([-2 2])
colorbar,

ax(3) = subplot(4,1,3);
pcolor(Data.Average_MatlabTimeStamp-datenum(2015,0,0), double(Data.Average_Range), double(Data.Average_VelUp1)' ),
hold on
if sum(good)~= 0, plot([SWIFT(good).time]-datenum(2015,0,0),min(z1234),'ks','linewidth',5), else end
shading flat,
datetick
set(gca,'YDir','reverse')
ylabel('Up1 [m/s]')
caxis([-2 2])
colorbar,

ax(4) = subplot(4,1,4);
pcolor(Data.Average_MatlabTimeStamp-datenum(2015,0,0), double(Data.Average_Range), double(Data.Average_VelUp2)' ),
hold on
if sum(good)~= 0, plot([SWIFT(good).time]-datenum(2015,0,0),min(z1234),'ks','linewidth',5), else end
shading flat,
datetick
set(gca,'YDir','reverse')
ylabel('Up2 [m/s]')
caxis([-2 2])
colorbar,

print('-dpng',[ wd '_BBvels.png'])


figure(3), % broadband amplitudes
for i = 1:4,
    ax(i) = subplot(4,1,i);   
    pcolor(Data.Average_MatlabTimeStamp-datenum(2015,0,0), double(Data.Average_Range), double(eval(['Data.Average_AmpBeam' num2str(i)]))' ), 
    hold on
    shading flat, 
    datetick
    set(gca,'YDir','reverse')
    ylabel(['Amp Beam ' num2str(i)])
    caxis([0 100])
    colorbar,
end
print('-dpng',[ wd '_rawBBamps.png'])

figure(4), % broadband correlations
for i = 1:4,
    ax(i) = subplot(4,1,i);   
    pcolor(Data.Average_MatlabTimeStamp-datenum(2015,0,0), double(Data.Average_Range), double(eval(['Data.Average_CorBeam' num2str(i)]))' ), 
    hold on
    shading flat, 
    datetick
    set(gca,'YDir','reverse')
    ylabel(['Cor Beam ' num2str(i)])
    caxis([0 100])
    colorbar,
end
print('-dpng',[ wd '_rawBBcors.png'])


end % close plot loop


end % close file loop


%% profiles 

figure(5), clf
figure(6), clf

for si = 1:length(SWIFT),
    
% dissipation profiles 
figure(5), 
plot(SWIFT(si).signature.HRprofile.tkedissipationrate,SWIFT(si).signature.HRprofile.z,'linewidth',2)
%scatter(ones(length(z5),1)*SWIFT(bi).time,SWIFT(bi).signature.HRprofile.z,SWIFT(bi).signature.HRprofile.tkedissipationrate,SWIFT(bi).signature.HRprofile.tkedissipationrate)
hold on
set(gca,'YDir','reverse')
set(gca,'FontWeight','demi','Fontsize',16)
xlabel('epsilon [m^2/s^3]')
ylabel('z [m]')

% current profiles
figure(6)
magnitude = ( SWIFT(si).signature.profile.east.^2 + SWIFT(si).signature.profile.north.^2 ).^0.5;
plot(magnitude,SWIFT(si).signature.profile.z,'linewidth',2)
hold on
set(gca,'YDir','reverse')
set(gca,'FontWeight','demi','Fontsize',16)
xlabel('relative current [m/s]')
ylabel('z [m]')

end % close BurstHR plotting loop

figure(5), print('-dpng',[ wd '_dissipationprofiles.png'])
figure(6), print('-dpng',[ wd '_relativecurrentprofiles.png'])



%%

save([ wd '.mat'],'SWIFT')
