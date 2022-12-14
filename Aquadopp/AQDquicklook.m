% Matlab script to read and plot Nortek Aquadopp Profiler ASCII files 
% (extracted from .prf files using AquaPro software)
% with some basic quality control
%
% J. Thomson, 2/2008 
%
% this version specific to fixed deployment at Duck FRF 
% *** WHICH IS NON-HR MODE ***
%

clear

fpath = '/Volumes/THOMSONDATA/Duck-SurfZoneOptics/FixedAquadopp/';
%prefix = 'SurfzoneOptics_Aquadopp_011-15Sep2010';
prefix = 'FRFbar_11-15Sep2010';


filebase = [ fpath prefix ];


% hardwire from hdr file:
res = 0.10; % m
blanking = 0.2; % m

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
pres = sen( :, 14); % dbar *** GAGE VALUE, NEEDS ATMOSPHERIC CORRECTION ***
temp = sen( :, 15); % deg C

% doppler data as arrays of time x cell: velocity (m/s), amplitude (cts), correlation (%)
u = load( [ filebase '.v1' ]);  
v = load( [ filebase '.v2' ]);  
w = load( [ filebase '.v3' ]);  
a1 = load( [ filebase '.a1' ]);  
a2 = load( [ filebase '.a2' ]);  
a3 = load( [ filebase '.a3' ]);  

[pts cells] = size(u);

% quality control (based on signal amplitude)
exclude = find( a1 < 18 | a2 < 18 | a3 < 18 );% & (time*ones(1,40)) > time(700));
u (exclude)  = NaN;
v (exclude) = NaN;
w (exclude) = NaN;
%ratio = size(exclude)./prod(size(u));

% profile
z = blanking+[0:(cells-1)]*res;  % actual data

% mean amplitude (all beams)
a = (a1+a2+a3)./3;

% remove low tide pts (pressure cutoff)
for i=1:pts,
    air = find( z > pres(i), 1);
    u(i,air:cells) = NaN;
    v(i,air:cells) = NaN;
    w(i,air:cells) = NaN;
    a(i,air:cells) = NaN;
end

 
%% plot burst

figure(1), clf

houri = 75;
burstlength = 300;
bursti = [(houri*3600+1):(houri*3600+burstlength)];
    
s(1)=subplot(4,1,1);
pcolor(1:burstlength,z',u(bursti,:)'), 
shading interp,
caxis([-1 1]),
hc = colorbar;
hold on,
plot(1:burstlength, pres(bursti), 'k', 'MarkerSize', 2),
%axis([ -inf inf 0 max(z) ]),
axis tight
ylabel('z (m)')
title(datestr(time(bursti(1))),'interpreter','none'), 
set(gca,'YLim',[0 3])
axes(hc), ylabel('u (m/s)')


s(2)=subplot(4,1,2);
pcolor(1:burstlength,z',v(bursti,:)'), 
shading interp,
caxis([-1 1]),
hc = colorbar;
hold on,
plot(1:burstlength, pres(bursti), 'k', 'MarkerSize', 2),
%axis([ -inf inf 0 max(z) ]),
axis tight
ylabel('z (m)')
set(gca,'YLim',[0 3])
axes(hc), ylabel('v (m/s)')


s(3)=subplot(4,1,3);
pcolor(1:burstlength,z',w(bursti,:)'), 
shading interp,
caxis([-1 1]),
hc = colorbar;
hold on,
plot(1:burstlength, pres(bursti), 'k', 'MarkerSize', 2),
%axis([ -inf inf 0 max(z) ]),
axis tight
ylabel('z (m)')
set(gca,'YLim',[0 3])
axes(hc), ylabel('w (m/s)')
    
        
s(4)=subplot(4,1,4);
pcolor(1:burstlength,z',a(bursti,:)'), hold on
shading flat
plot(1:burstlength, pres(bursti), 'k', 'MarkerSize', 2),
caxis([0 200])
axis tight
ylabel('z (m)'),
xlabel('t (s)'),
set(gca,'YLim',[0 3])
hc = colorbar; 
axes(hc), ylabel('Amp (dB)')

linkaxes(s,'x')

 print('-dpng',['SurfzoneOptics_Aquadopp_' datestr(time(bursti(1))) '_' num2str(burstlength) 'secburst.png'])



%% means (for plotting)

maxpts = pts; % use "pts" for whole record

for i=[1:60:(maxpts-60)],
    inds = [i:(i+59)];
    j=((i-1)/60+1);
    umeans(j,:) = nanmean(u(inds,:),1);
    vmeans(j,:) = nanmean(v(inds,:),1);
    wmeans(j,:) = nanmean(w(inds,:),1);
    a1means(j,:) = nanmean(a1(inds,:),1);
    a2means(j,:) = nanmean(a2(inds,:),1);
    a3means(j,:) = nanmean(a3(inds,:),1);
    timemeans(j) = nanmean(time(inds));    
    presmeans(j) = nanmean(pres(inds));
    presvar(j) = var(pres(inds));
    ameans(j,:) = nanmean(a(inds,:),1);
end


%ampmeans = (a1means + a2means + a3means )./3;
%for i=1:length(ampmeans),
%    air = find( z > presmeans(i), 1);
%    ampmeans(i,air:cells) = NaN;
%end


%% plot means

figure(2), clf
    
    s(1)=subplot(4,1,1);
    pcolor(timemeans,z',umeans'), 
    shading interp,
    caxis([-1 1]),
    hc = colorbar;
    hold on,
    plot(timemeans, presmeans, 'k.', 'MarkerSize', 2),
    %axis([ -inf inf 0 max(z) ]),
    datetick,
    axis tight
    ylabel('z (m)')
    title(prefix,'interpreter','none'), 
    axes(hc), ylabel('u (m/s)')
         
s(3)=subplot(4,1,2);
    pcolor(timemeans,z',vmeans'), 
    shading interp,
    caxis([-1 1]),
    hc = colorbar;
    hold on,
    plot(timemeans, presmeans, 'k.', 'MarkerSize', 2),
    %axis([ -inf inf 0 max(z) ]),
    datetick,
    axis tight
    ylabel('z (m)')
    axes(hc), ylabel('v (m/s)')

    s(3)=subplot(4,1,3);
    pcolor(timemeans,z',wmeans'), 
    shading interp,
    caxis([-1 1]),
    hc = colorbar;
    hold on,
    plot(timemeans, presmeans, 'k.', 'MarkerSize', 2),
    %axis([ -inf inf 0 max(z) ]),
    datetick,
    axis tight
    ylabel('z (m)')
    axes(hc), ylabel('w (m/s)')
    
s(4) = subplot(4,1,4);
pcolor(timemeans,z',ameans'), hold on
shading flat
plot(timemeans, presmeans, 'k.', 'MarkerSize', 2),
datetick,
axis tight
ylabel('z (m)')
hc= colorbar; 
axes(hc), ylabel('Amp (dB)')

linkaxes(s,'x')

 print('-dpng',[prefix '_60secmeans.png'])

%%

README = 'Bottom-mount Aquadopp deployed by Jim Thomson (jthomson@apl.uw.edu) at [N 36.18425, W 075.75014] during the Surfzone Optics Experiment at the FRF (Duck).  Timestamps are UTC, and vertical (z) coordinates are relative to local sand level.  Velocities are u (positive onshore), v (postive alongshore/south), w (positive upwards).  Near surface velocities are to be used with caution.  Raw 1 Hz data have velocity uncertainty of 0.25 m/s.';
 
 save(prefix)
