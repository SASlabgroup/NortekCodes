% Matlab script to read and plot Nortek AWAC ASCII files 
% (extracted from .prf binary files using Nortek software)
% with some basic quality control
%
% J. Thomson, 8/2010 
%
% 

clear


%fpath = '/Volumes/THOMSONDATA/ADMIRALTY_Nov2010/ADCPdata_Nov2010/';
fpath = './';
prefix = 'SS04_AWAC_May-Jun2011';
%prefix = 'Admiralty_AWAC_Aug-Nov2010_';

filebase = [ fpath prefix ];

doff = 0.69;  % distance above seabed [m]

% hardwire from hdr file:
res = 1; % m
blanking = 0.4; % m

sen = load( [ filebase '.sen' ]);

% timestamps
month = sen(:, 1);
day = sen(:, 2);
year = sen(:, 3);
hour = sen(:, 4);
minute = sen(:, 5);
second = sen(:, 6);
time = datenum( year, month, day, hour, minute, second);

% ancillary data
voltage = sen( :, 9);
soundspd = sen( :, 10); % m/s
heading = sen( :, 11); % deg M
pitch = sen( :, 12); % deg
roll = sen( :, 13); % deg
pres = sen( :, 14); % dbar *** GAGE VALUE ***
temp = sen( :, 15); % deg C

% doppler data as arrays of time x cell: velocity (m/s), amplitude (cts), correlation (%)
v1 = load( [ filebase '.v1' ]);  
v2 = load( [ filebase '.v2' ]);  
v3 = load( [ filebase '.v3' ]);  
%a1 = load( [ filebase '.a1' ]);  
%a2 = load( [ filebase '.a2' ]);  
%a3 = load( [ filebase '.a3' ]);  
%%

[pts cells] = size(v1);

% quality control (based on signal amplitude)
% exclude = find( a1 < 50 | a2 < 50 | a3 < 50);
% v1 (exclude)  = NaN;
% v2 (exclude) = NaN;
% v3 (exclude) = NaN;
% ratio = size(exclude)./prod(size(v1));

% profile
z = doff+blanking+[0:(cells-1)]*res;  % actual data

% remove surface pts (pressure cutoff)
for i=1:pts,
    air = find( z > pres(i), 1);
    v1(i,air:cells) = NaN;
    v2(i,air:cells) = NaN;
    v3(i,air:cells) = NaN;
end
    

% option to prescribe no-slip at bottom (for pretty plots)
%z = [ 0 z ];  
%v1 = [ zeros(pts,1) v1 ];
%v2 = [ zeros(pts,1) v2 ];
%v3 = [ zeros(pts,1) v3 ];


%% wave bursts
% call subroutine AWACanalysis.m for wave directional spectra, moments, etc

% wavefiles = dir([fpath '*.wad']);
% nbursts = length(wavefiles);
% samplingrate = 1; % Hz
% window = 512;
% noverlap = 128;
% nfft = 512;
% 
% for i=1:nbursts,
%    wad = load([fpath wavefiles(i).name]);
%    distance1 = wad(:,8); % m
%    distance2 = wad(:,9); % m
%    Hsig(i) = mean( [ 4*std(distance1) 4*std(distance2) ] ); % m
%    [psd f ] = pwelch(distance1, window, noverlap, nfft, samplingrate); 
%    clf, loglog(f, psd), drawnow
%    Hpsd = psd(i,:);
% end

%%

save([prefix])

%% deal with large data files

clear all
bsecs = 600;

%load SS04_AWAC_May-Jun2011 v1 time z prefix
% for i = 1:length(time)./bsecs,
%     east(i,:) = mean(v1(((i-1)*bsecs+1):i*bsecs,:) ,1);
% end
% clear v1
%load SS04_AWAC_May-Jun2011 v2 
% for i = 1:length(time)./bsecs,
%     north(i,:) = mean(v2(((i-1)*bsecs+1):i*bsecs,:) ,1);
% end
% clear v2
%vmag = sqrt(east.^2 + north.^2);

%load SS04_AWAC_May-Jun2011 v1 time z prefix
%v1 = v1([1:bsecs:length(time)],:);
%load SS04_AWAC_May-Jun2011 v2 
%v2 = v2([1:bsecs:length(time)],:);
%vmag = sqrt(v1.^2 + v2.^2);

load SS04_AWAC_May-Jun2011 v1 v2 time z prefix
for i = 1:length(time)./bsecs,
     vmean(i,:) = mean( sqrt( v1(((i-1)*bsecs+1):i*bsecs,:).^2 + v2(((i-1)*bsecs+1):i*bsecs,:).^2 ) ,1);
     vstd(i,:) = std( sqrt( v1(((i-1)*bsecs+1):i*bsecs,:).^2 + v2(((i-1)*bsecs+1):i*bsecs,:).^2 ) ,1);
end
clear v1 v2

figure(1), clf
pcolor(time([1:bsecs:length(time)-bsecs])-datenum(2011,0,0), z, vmag'), shading flat
datetick,
ylabel('z (m)')
cb = colorbar;
axes(cb), ylabel('m/s')
print('-dpng',[prefix '_velocity.png'])

%% plotting

figure(1), clf

for i=1:3,
    
    ax(i) = subplot(3,1,i);
    pcolor(time([1:300:length(time)])-datenum(2010,0,0),z,eval(['v' num2str(i)])'), 
    shading interp,
    %caxis([-1 1]),
    hc(i) = colorbar;
    hold on,
    plot(time-datenum(2010,0,0), pres, 'k.', 'MarkerSize', 2),
    %axis([ -inf inf 0 max(z) ]),
    datetick,
    axis tight
    ylabel('z (m)')
    if i==1, 
        title(prefix,'interpreter','none'), 
        axes(hc(i)), ylabel('east (m/s)')
    elseif i==2,
        axes(hc(i)), ylabel('north (m/s)')
    elseif i==3,
        axes(hc(i)), ylabel('up (m/s)')
    else end
        
end
linkaxes(ax)
print('-dpng',[prefix '_velocity.png'])


figure(2), clf

for i=1:3,
    
    ax(i) = subplot(3,1,i);
    pcolor(time-datenum(2010,0,0),z,eval(['a' num2str(i)])'), 
    shading interp,
    %caxis([-1 1]),
    hc(i) = colorbar;
    hold on,
    plot(time-datenum(2010,0,0), pres, 'k.', 'MarkerSize', 2),
    %axis([ -inf inf 0 max(z) ]),
    datetick,
    axis tight
    ylabel('z (m)')
    if i==1, 
        title(prefix,'interpreter','none'), 
    else end
    axes(hc(i)), ylabel(['Amp ' num2str(i)])
        
end
linkaxes(ax)
print('-dpng',[prefix '_amplitude.png'])


figure(3), clf

ax(1) = subplot(3,1,1);
plot(time-datenum(2010,0,0), heading./100, time-datenum(2010,0,0), pitch, time-datenum(2010,0,0), roll),
ylabel('deg')
datetick
legend('heading./100','pitch','roll')
title(prefix,'interpreter','none'), 

ax(2) = subplot(3,1,2);
plot(time-datenum(2010,0,0), temp), 
ylabel('Temp (deg C)')
datetick

ax(3) = subplot(3,1,3);
plot(time-datenum(2010,0,0), voltage), 
ylabel('Battery Voltage')
datetick

linkaxes(ax,'x')
print('-dpng',[prefix '_ancillary.png'])