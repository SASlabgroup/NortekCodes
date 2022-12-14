% Matlab script to read and plot Nortek Aquadopp Profiler ASCII files 
% (extracted from .prf files using AquaPro software)
% with some basic quality control
%
% J. Thomson, 2/2008
%           rev. 11/2011 for doppcats 
%           rev. 14/2014 for GPS integration and statistics
%

clear

% files to process *** CHECK HDR SETTINGS FOR BIN DETAILS **
date = 'Jun2022', % date (and directory) to process
ID ='', % Doppcat to process
readraw = 0;  % binary flag (0 or 1) to force read in of raw data
fpath = './';% ['./' date '/']; % path to directory with that days data
fname = ['Doppcat' ID '_AQD_' date]; % filenaming convention (must be exact)

% minimum echo amplitude (for quality control)
minamp = 30;

% processing parameters
timestep = 1; % time step (in minutes) for processing statistics and positions
zhub = 1; % hub height (distance in m below surface)

% profile (settings from .hdr file)
res = 0.50; % m
blanking = 0.2; % m
cells = 20;
z = blanking+[0:(cells-1)]*res;  % actual data
[ zoffset hubindex ] = min(abs(z-zhub)); 

% Doppler noise 
if str2num(ID) <= 4; 
    noise = 0.06; % std dev of horizontal velocity
else
    noise = 0.22; % std dev of horizontal velocity
end

%% load (or read in) Aquadopps

if isempty( dir([fpath fname '.mat']) ) | readraw==1,

% read sensor data
sen = load( [ fpath fname '.sen' ]);

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
pres = sen( :, 14); % dbar *** GAGE VALUE, NEEDS ATMOSPHERIC CORRECTION ***
temp = sen( :, 15); % deg C

% doppler data as arrays of time x cell: velocity (m/s), amplitude (cts), correlation (%)
u = load( [ fpath fname '.v1' ]);  
v = load( [ fpath fname '.v2' ]);  
w = load( [ fpath fname '.v3' ]);  
a1 = load( [ fpath fname '.a1' ]);  
a2 = load( [ fpath fname '.a2' ]);  
a3 = load( [ fpath fname '.a3' ]);  

[pts cells] = size(u);

else 
    save temp
    load([fpath fname '.mat'])
    load temp
end


%% quality control Aquadopps (based on signal amplitude)
exclude = find( a1 < minamp | a2 < minamp | a3 < minamp );% 
u (exclude)  = NaN;
v (exclude) = NaN;
w (exclude) = NaN;
QCratio = size(exclude)./prod(size(u));

% mean amplitude (all beams)
a = (a1+a2+a3)./3;

%% remove seabed bins
% if ID=='01',
%     seabedpts = find(z>5);
% elseif ID =='02',
%     seabedpts = find(z>2);
% else
%     seabedpts = find(z==max(z));
% end
% 
% u(:,seabedpts)=NaN;
% v(:,seabedpts)=NaN;
% w(:,seabedpts)=NaN;

%% make means and stats

timemeans =  [min(time):timestep/(24*60):max(time)]';

for ti=1:(length(timemeans)-1),
    thisminute = find(time>timemeans(ti) & time<timemeans(ti+1));
    meanu(ti,:) = nanmean(u(thisminute,:));
    stdu(ti,:) = nanstd(u(thisminute,:));
    meanv(ti,:) = nanmean(v(thisminute,:));
    stdv(ti,:) = nanstd(v(thisminute,:));
    meanw(ti,:) = nanmean(w(thisminute,:));
    stdw(ti,:) = nanstd(w(thisminute,:));
end
timemeans(length(timemeans))=[];

TI = real( sqrt( stdu.^2 - noise.^2) ./ abs(meanu) );


%% plot data

figure(1), clf
    
s(1)=subplot(4,1,1);
pcolor(timemeans-datenum(2011,0,0),z',meanu'), hold on
shading interp,
datetick
caxis([-3 3]),
ylabel('z (m)')
title(fname,'interpreter','none'), 
set(gca,'YDir','reverse')
hc = colorbar; %axes(hc), 
title('u (m/s)')


s(2)=subplot(4,1,2);
pcolor(timemeans-datenum(2011,0,0),z',meanv'), hold on
shading interp,
datetick
caxis([-3 3]),
ylabel('z (m)')
%title(fname,'interpreter','none'), 
set(gca,'YDir','reverse')
hc = colorbar; %axes(hc), 
title('v (m/s)')

s(3)=subplot(4,1,3);
pcolor(timemeans-datenum(2011,0,0),z',meanw'), hold on
shading interp,
datetick
caxis([-3 3]),
ylabel('z (m)')
%title(fname,'interpreter','none'), 
set(gca,'YDir','reverse')
hc = colorbar; %axes(hc), 
title('w (m/s)')
    
        
s(4)=subplot(4,1,4);
pcolor(time-datenum(2011,0,0),z',a'), hold on
shading interp,
datetick
caxis([0 200]),
ylabel('z (m)')
%title(fname,'interpreter','none'), 
set(gca,'YDir','reverse')
hc = colorbar; %axes(hc), 
title('Echo amp (dB)')

linkaxes(s,'x')  % allows simultaneuos scrolling 

print('-dpng',[fpath fname '_velocityprofiles.png'])

%%
figure(2), clf, 

t(1) = subplot(2,1,1);
plot(time,u(:,hubindex),timemeans,meanu(:,hubindex),'r')
datetick
ylabel('Surface Velocity [m/s]')
title(fname)
title(fname,'interpreter','none')

t(2) = subplot(2,1,2);
plot(timemeans,100*TI(:,hubindex),'r.')
datetick
ylabel('TI [%]')
set(gca,'YLim',[0 30])

linkaxes(t,'x')  % allows simultaneuos scrolling 

print('-dpng',[fpath fname '_hubvelocity.png'])

%%

save([fpath fname])
save([fpath fname '_means'],'*mean*','z')
