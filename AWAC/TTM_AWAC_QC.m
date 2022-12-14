% Matlab script to read and plot Nortek AWAC ASCII files 
% (extracted from .prf binary files using Nortek software)
% with some basic quality control
%
% J. Thomson, 8/2010 
%
% 

clear

minamp = 60;  % removes low signal returns 
maxamp = 160;  % removes mooring interference (high signal returns)

fpath = '/Volumes/Data/ADMIRALTY/ADMIRALTY_Jun2012/TTM_AWAC/';
prefix = 'TTM_AWAC_Jun2012';

filebase = [ fpath prefix ];

doff = 0.74;  % distance above seabed [m]

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

mean(pres)

% doppler data as arrays of time x cell: velocity (m/s), amplitude (cts), correlation (%)
v1 = load( [ filebase '.v1' ]);  
v2 = load( [ filebase '.v2' ]);  
v3 = load( [ filebase '.v3' ]);  
a1 = load( [ filebase '.a1' ]);  
a2 = load( [ filebase '.a2' ]);  
a3 = load( [ filebase '.a3' ]);  

%% basic QC

[pts cells] = size(v1);

%quality control (based on signal amplitude)
exclude = find( a1 < minamp | a2 < minamp | a3 < minamp  |  a1 > maxamp | a2 > maxamp | a3 > maxamp);
v1 (exclude)  = NaN;
v2 (exclude) = NaN;
v3 (exclude) = NaN;
QCratio = size(exclude)./prod(size(v1))

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



%%
save([prefix])%,'u','v','w','z','time','pitch','roll','heading','pres')


%% plotting

figure(1), clf

for i=1:3,
    
    ax(i) = subplot(3,1,i);
    pcolor(time-datenum(2010,0,0),z,eval(['v' num2str(i)])'), 
    shading interp,
    caxis([-2 2]),
    hc(i) = colorbar;
    hold on,
    %plot(time-datenum(2010,0,0), pres, 'k.', 'MarkerSize', 5),
    %axis([ -inf inf 0 max(z) ]),
    datetick,
    axis tight
    ylabel('z (m)')
    if i==1, 
        title([prefix ', Depth = ' num2str(mean(pres),3) ' m'],'interpreter','none'), 
        axes(hc(i)), ylabel('east (m/s)')
    elseif i==2,
        axes(hc(i)), ylabel('north (m/s)')
    elseif i==3,
        caxis([-0.5 0.5]),
        axes(hc(i)), ylabel('up (m/s)'),
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
    %plot(time-datenum(2010,0,0), pres, 'k.', 'MarkerSize', 5),
    %axis([ -inf inf 0 max(z) ]),
    datetick,
    axis tight
    ylabel('z (m)')
    if i==1, 
        title([prefix ', Depth = ' num2str(mean(pres),3) ' m'],'interpreter','none'), 
    else end
    axes(hc(i)), ylabel(['Amp ' num2str(i)])
        
end
linkaxes(ax)
print('-dpng',[prefix '_amplitude.png'])


figure(3), clf
ax(1) = subplot(4,1,1);
plot(time-datenum(2010,0,0), heading),
ylabel('deg')
datetick
legend('heading')
title([prefix ', Depth = ' num2str(mean(pres),3) ' m'],'interpreter','none'), 
        
ax(2) = subplot(4,1,2);
plot(time-datenum(2010,0,0), pitch, time-datenum(2010,0,0), roll),
ylabel('deg')
datetick
legend('pitch','roll')

ax(3) = subplot(4,1,3);
plot(time-datenum(2010,0,0), temp), 
ylabel('Temp (deg C)')
datetick

ax(4) = subplot(4,1,4);
plot(time-datenum(2010,0,0), voltage), 
ylabel('Battery Voltage')
datetick

linkaxes(ax,'x')
print('-dpng',[prefix '_ancillary.png'])