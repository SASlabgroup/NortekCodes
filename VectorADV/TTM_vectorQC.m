% process Vector data from TTM in daily blocks (continuous sampling)
% 16 Hz sampling
%
%   J. Thomson, 6/2012, modified from 2011 TTT processing

clear all, close all

filebase = 'Killarney Deployment 1 Vector_all';
fpath = './';


%% load raw velocity data (16 Hz)
data = load([fpath filebase '.dat']);
burstnum = data(:,1); 
u = data(:,3);
v = data(:,4);
w = data(:,5);
cor1 = data(:,12);
cor2 = data(:,13);
cor3 = data(:,14);
pres = data(:,15);

%% apply correlation cutoff
mincor = 60;
bad = cor1 < mincor | cor2 < mincor | cor3 < mincor;
u(bad) = NaN;
v(bad) = NaN;
w(bad) = NaN;

percentbad = sum(bad)./length(u)

%% load ancillary data (1 Hz)

sen = load([fpath filebase '.sen']);
mo = sen(:,1);
day = sen(:,2);
year = sen(:,3);
hour = sen(:,4);
min = sen(:,5);
sec = sen(:,6);
time = datenum(year,mo,day,hour,min,sec);

heading = sen(:,11);
pitch = sen(:,12);
roll = sen(:,13);

% if filebase(5:7)=='APL', % cannot use prh data when vertically mounted
%     heading = NaN;
%     pitch = NaN;
%     roll = NaN;
% elseif filebase(5:8)=='NREL', % not prf data in INS mode
%     heading = NaN;
%     pitch = NaN;
%     roll = NaN;
%     trim = 1:(32*3600*5.5);  % ASCII conversion did not apply date limits
%     u(trim) = []; v(trim)=[]; w(trim) = [];
%     pres = NaN (size(pres) ); % NREL pres was out of range
% end

%% phase spacce despiking 
%[ u v w spikeinds] = func_despike_phasespace3d_3var( u, v, w, 2 );

%%

save([fpath filebase '.mat'])

%% plots
% 
pts = length(u);

figure(1), clf
subplot(2,1,1)
plot(1:pts,u,1:pts,v,1:pts,w)
set(gca,'Ylim',[-3 1])
title(filebase,'interpreter','none'),
ylabel('m/s')
legend('u','v','w')
subplot(2,1,2),
plot(1:pts,pres,'k'), hold on
ylabel('m')
print('-dpng',[filebase '_uvwp.png'])

figure(2), clf
plot(time,pitch,time,roll,time,heading), hold on
datetick, axis tight
legend('pitch','roll','heading')
title(filebase,'interpreter','none')
print('-dpng',[filebase '_prh.png'])


%% mooring behavior

% speed = resample(sqrt(u.^2 + w.^2),1,16);
% figure(3), clf, 
% plot(speed,pitch,'k.')
% set(gca,'fontweight','demi','fontsize',14)
% xlabel('speed [m/s]')
% ylabel('pitch [deg]')
% title(filebase,'interpreter','none')
% axis([0 2.5 0 25])
% print('-dpng',[filebase '_pitchvsspeed.png'])
