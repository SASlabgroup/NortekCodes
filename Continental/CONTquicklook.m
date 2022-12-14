% Matlab script to read and plot Nortek Continental ASCII files 
% (extracted from .prf binary files using Nortek software)
% with some basic quality control
%
% J. Thomson, 8/2010 
%
% 

clear

%fpath = '/Volumes/THOMSONDATA/ADMIRALTY_Feb2011/ADCPdata_Feb2011/';
fpath = './';
prefix = 'SS02_Cont_Apr-Jul2013';

filebase = [ fpath prefix ];


% hardwire from hdr file:
res = 1.0; % m
blanking = 1.0; % m

sen = load( [ fpath prefix '.sen' ]);

% timestamps
month = sen(:, 1);
day = sen(:, 2);
year = sen(:, 3);
hour = sen(:, 4);
minute = sen(:, 5);
second = sen(:, 6);
time = datenum( year, month, day, hour, minute, second);

% ancillary data
error = sen(:,7);
status = sen(:,8);
voltage = sen( :, 9);
soundspd = sen( :, 10); % m/s
heading = sen( :, 11); % deg M
pitch = sen( :, 12); % deg
roll = sen( :, 13); % deg
pres = sen( :, 14); % dbar *** GAGE VALUE ***
temp = sen( :, 15); % deg C

%% doppler data as arrays of time x cell: velocity (m/s), amplitude (cts), correlation (%)
v1 = load( [ filebase '.v1' ]);  
v2 = load( [ filebase '.v2' ]);  
v3 = load( [ filebase '.v3' ]);  
a1 = load( [ filebase '.a1' ]);  
a2 = load( [ filebase '.a2' ]);  
a3 = load( [ filebase '.a3' ]);  

[pts cells] = size(v1);

% quality control (based on signal amplitude)
exclude = find( a1 < 50 | a2 < 50 | a3 < 50);
v1 (exclude)  = NaN;
v2 (exclude) = NaN;
v3 (exclude) = NaN;
ratio = size(exclude)./prod(size(v1));

% profile
z = blanking+[0:(cells-1)]*res;  % actual data

%% remove surface pts (pressure cutoff)
for i=1:pts,
    air = find( z > 0.9*pres(i), 1);
    v1(i,air:cells) = NaN;
    v2(i,air:cells) = NaN;
    v3(i,air:cells) = NaN;
end
    

% option to prescribe no-slip at bottom (for pretty plots)
%z = [ 0 z ];  
%v1 = [ zeros(pts,1) v1 ];
%v2 = [ zeros(pts,1) v2 ];
%v3 = [ zeros(pts,1) v3 ];


%% plotting

figure(1), clf

for i=1:3,
    
    ax(i) = subplot(3,1,i);
    pcolor(time-datenum(2010,0,0),z,eval(['v' num2str(i)])'), 
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
        caxis([-.5 .5])
    else end
        
end
linkaxes(ax)
print('-dpng',[prefix '_velocity.png'])

%%


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

%%

figure(3), clf

ax(1) = subplot(3,1,1);
plot(time-datenum(2010,0,0), heading./10, time-datenum(2010,0,0), pitch, time-datenum(2010,0,0), roll),
ylabel('deg')
datetick
legend('heading./10','pitch','roll')
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

%%

save([ prefix '.mat'])
