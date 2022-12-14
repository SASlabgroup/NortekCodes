% Matlab script to read and plot Nortek Aquadopp Profiler ASCII files 
% (extracted from .prf files using AquaPro software)
% with some basic quality control
%
% J. Thomson, 2/2008, rev 10/2010 
%
% this version specific to fixed deployments
% *** WHICH IS NON-HR MODE AND ENU COORDINATES***
%

clear

fpath = './';

wd = pwd;
wdi = find(wd == '/',1,'last');
prefix  =  wd((wdi+1):length(wd));

filebase = [ fpath prefix ];


% hardwire from hdr file:
res = 0.50; % m
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
v1 = load( [ filebase '.v1' ]);  
v2 = load( [ filebase '.v2' ]);  
v3 = load( [ filebase '.v3' ]);  
a1 = load( [ filebase '.a1' ]);  
a2 = load( [ filebase '.a2' ]);  
a3 = load( [ filebase '.a3' ]);  

[pts cells] = size(v1);

% profile
z = blanking+[1:cells]*res;  % actual data

% mean amplitude (all beams)
a = (a1+a2+a3)./3;



%% Quality control

minamp = 20;

% quality control (based on signal amplitude)
exclude = find( a1 < minamp | a2 < minamp | a3 < minamp );% & (time*ones(1,40)) > time(700));
v1 (exclude)  = NaN;
v2 (exclude) = NaN;
v3 (exclude) = NaN;
%ratio = size(exclude)./prod(size(u));


% remove low tide pts (pressure cutoff)
for i=1:pts,
    air = find( z > (pres(i)-.2), 1);
    v1(i,air:cells) = NaN;
    v2(i,air:cells) = NaN;
    v3(i,air:cells) = NaN;
    u(i,air:cells) = NaN;
    v(i,air:cells) = NaN;
    w(i,air:cells) = NaN;
    a(i,air:cells) = NaN;
end

%% ENU
east = v1;
north = v2;
up = v3;


%%
 
 save(prefix)
 
 %% plot
 
figure(2), clf, clear s
    
    
    s(1)=subplot(4,1,1);
    pcolor(time,z',east'), 
    shading interp,
    caxis([-1 1]),
    colorbar;
    hold on,
    plot(time, pres, 'k.', 'MarkerSize', 2),
    %axis([ -inf inf 0 max(z) ]),
    datetick,
    axis tight
    ylabel('z (m)')
    title(prefix,'interpreter','none'), 
         
    s(2)=subplot(4,1,2);
    pcolor(time,z',north'), 
    shading interp,
    caxis([-1 1]),
    colorbar;
    hold on,
    plot(time, pres, 'k.', 'MarkerSize', 2),
    %axis([ -inf inf 0 max(z) ]),
    datetick,
    axis tight
    ylabel('z (m)')

    s(3)=subplot(4,1,3);
    pcolor(time,z',up'), 
    shading interp,
    caxis([-1 1]),
    colorbar;
    hold on,
    plot(time, pres, 'k.', 'MarkerSize', 2),
    %axis([ -inf inf 0 max(z) ]),
    datetick,
    axis tight
    ylabel('z (m)')
    
    s(4) = subplot(4,1,4);
    pcolor(time,z',a'), hold on
    shading flat
    colorbar;
    plot(time, pres, 'k.', 'MarkerSize', 2),
    datetick,
    axis tight
    ylabel('z (m)')

linkaxes(s,'x')

 print('-dpng',[prefix '.png'])
 
%% plot burst

% figure(1), clf
% 
% houri = 3;
% burstlength = 300;
% bursti = [(houri*3600+1):(houri*3600+burstlength)];
%     
% b(1)=subplot(4,1,1);
% pcolor(1:burstlength,z',east(bursti,:)'), set(gca,'XTickLabel',[]),
% shading interp,set(gca,'Fontsize',14,'fontweight','demi'),
% caxis([-1 1]),
% hc = colorbar;
% hold on,
% plot(1:burstlength, pres(bursti), 'k', 'MarkerSize', 2),
% %axis([ -inf inf 0 max(z) ]),
% axis tight
% ylabel('z (m)')
% title(datestr(time(bursti(1))),'interpreter','none'), 
% set(gca,'YLim',[0 3])
% axes(hc), set(gca,'Fontsize',12,'fontweight','demi'),ylabel('u (m/s)')
% 
% 
% b(2)=subplot(4,1,2);
% pcolor(1:burstlength,z',north(bursti,:)'), set(gca,'XTickLabel',[]),
% shading interp,set(gca,'Fontsize',14,'fontweight','demi'),
% caxis([-1 1]),
% hc = colorbar;
% hold on,
% plot(1:burstlength, pres(bursti), 'k', 'MarkerSize', 2),
% %axis([ -inf inf 0 max(z) ]),
% axis tight
% ylabel('z (m)')
% set(gca,'YLim',[0 3])
% axes(hc), set(gca,'Fontsize',12,'fontweight','demi'),ylabel('v (m/s)')
% 
% 
% b(3)=subplot(4,1,3);
% pcolor(1:burstlength,z',up(bursti,:)'), set(gca,'XTickLabel',[]),
% shading interp,set(gca,'Fontsize',14,'fontweight','demi'),
% caxis([-1 1]),
% hc = colorbar;
% hold on,
% plot(1:burstlength, pres(bursti), 'k', 'MarkerSize', 2),
% %axis([ -inf inf 0 max(z) ]),
% axis tight
% ylabel('z (m)')
% set(gca,'YLim',[0 3])
% axes(hc), set(gca,'Fontsize',12,'fontweight','demi'),ylabel('w (m/s)')
%     
%         
% b(4)=subplot(4,1,4);
% pcolor(1:burstlength,z',a(bursti,:)'), hold on
% shading flat
% plot(1:burstlength, pres(bursti), 'k', 'MarkerSize', 2),set(gca,'Fontsize',14,'fontweight','demi'),
% caxis([0 200])
% axis tight
% ylabel('z (m)'),
% xlabel('t (s)'),
% set(gca,'YLim',[0 3])
% hc = colorbar; 
% axes(hc), set(gca,'Fontsize',12,'fontweight','demi'),ylabel('Amp (dB)')
% 
% linkaxes(b,'x')

 %print('-dpng',['RIVET_' datestr(time(bursti(1))) '_' num2str(burstlength) 'secburst.png'])



%% means (for plotting)
% % note that using nanmean would bias by waves, because they are only onshore velocities above the trough
% % whereas a straight mean will ignore these z levels since there are some NaNs there in a given burst
% 
% maxpts = pts; % use "pts" for whole record
% 
% for i=[1:60:(maxpts-60)],
%     inds = [i:(i+59)];
%     j=((i-1)/60+1);
%     eastmeans(j,:) = mean(east(inds,:),1);
%     northmeans(j,:) = mean(north(inds,:),1);
%     upmeans(j,:) = mean(up(inds,:),1);
%     a1means(j,:) = mean(a1(inds,:),1);
%     a2means(j,:) = mean(a2(inds,:),1);
%     a3means(j,:) = mean(a3(inds,:),1);
%     timemeans(j) = mean(time(inds));    
%     presmeans(j) = mean(pres(inds));
%     presvar(j) = var(pres(inds));
%     ameans(j,:) = mean(a(inds,:),1);
% end
% 

%ampmeans = (a1means + a2means + a3means )./3;
%for i=1:length(ampmeans),
%    air = find( z > presmeans(i), 1);
%    ampmeans(i,air:cells) = NaN;
%end


%% plot means
% 
% figure(2), clf, clear s
%     
%     
%     s(1)=subplot(4,1,1);
%     pcolor(timemeans,z',eastmeans'), 
%     shading interp,
%     caxis([-1 1]),
%     hc = colorbar;
%     hold on,
%     plot(timemeans, presmeans, 'k.', 'MarkerSize', 2),
%     %axis([ -inf inf 0 max(z) ]),
%     datetick,
%     axis tight
%     ylabel('z (m)')
%     title(prefix,'interpreter','none'), 
%     axes(hc), ylabel('u (m/s)')
%          
%     s(2)=subplot(4,1,2);
%     pcolor(timemeans,z',northmeans'), 
%     shading interp,
%     caxis([-1 1]),
%     hc = colorbar;
%     hold on,
%     plot(timemeans, presmeans, 'k.', 'MarkerSize', 2),
%     %axis([ -inf inf 0 max(z) ]),
%     datetick,
%     axis tight
%     ylabel('z (m)')
%     axes(hc), ylabel('v (m/s)')
% 
%     s(3)=subplot(4,1,3);
%     pcolor(timemeans,z',upmeans'), 
%     shading interp,
%     caxis([-1 1]),
%     hc = colorbar;
%     hold on,
%     plot(timemeans, presmeans, 'k.', 'MarkerSize', 2),
%     %axis([ -inf inf 0 max(z) ]),
%     datetick,
%     axis tight
%     ylabel('z (m)')
%     axes(hc), ylabel('w (m/s)')
%     
%     s(4) = subplot(4,1,4);
%     pcolor(timemeans,z',ameans'), hold on
%     shading flat
%     hca= colorbar;
%     plot(timemeans, presmeans, 'k.', 'MarkerSize', 2),
%     datetick,
%     axis tight
%     ylabel('z (m)')
%     axes(hca), ylabel('Amp (dB)')
% 
% linkaxes(s,'x')
% 
%  print('-dpng',[prefix '_60secmeans.png'])
% 
