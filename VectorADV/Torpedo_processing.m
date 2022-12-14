% process Vector ADV data from "Turbulence torpedo" 
% assumes continuous sampling
%
% J. Thomson, 8/2014... modified from TTT codes, 5/2011

clear all, close all

date = '13May2015';
fname = ['TTT_AFrame_ADV_' date ];
fpath = ['./']; 
%fpath = ['./' date '/'];

samplingrate = 16;

hourstotrim = 0; % hours of junk data before deployed (export times don't work in Nortek software)


%% load raw velocity data (samplingrate Hz)

trim = 1:(samplingrate*3600*hourstotrim);  % ASCII conversion did not apply date limits
data = load([fpath fname '.dat']);
data(trim,:)=[];

burstnum = data(:,1); 
u = data(:,3);
v = data(:,4);
w = data(:,5);
amp1 = data(:,6);
amp2 = data(:,7);
amp3 = data(:,8);
snr1 = data(:,9);
snr2 = data(:,10);
snr3 = data(:,11);
cor1 = data(:,12);
cor2 = data(:,13);
cor3 = data(:,14);
pres = data(:,15);


%% load ancillary data (1 Hz)

sen = load([fpath fname '.sen']); % this does get trimmed by Nortek software
mo = sen(:,1);
day = sen(:,2);
year = sen(:,3);
hour = sen(:,4);
min = sen(:,5);
sec = sen(:,6);
time = datenum(year,mo,day,hour,min,sec);

voltage = sen(:,9);

heading = NaN;
pitch = NaN;
roll = NaN;


%% QC for out of water

u( find(pres<1) ) = NaN;
v( find(pres<1) ) = NaN;
w( find(pres<1) ) = NaN;


%% QC (remove low cor points)
corcutoff = 50;
u( find (cor1 < corcutoff | cor2 < corcutoff | cor3 < corcutoff ) ) = NaN;
v( find (cor1 < corcutoff | cor2 < corcutoff | cor3 < corcutoff ) ) = NaN;
w( find (cor1 < corcutoff | cor2 < corcutoff | cor3 < corcutoff ) ) = NaN;


%% spike cuttoff

wrap = 3.5; % wrap velocity from settings (horizontal range)
u( abs(u) > 3.5 ) = NaN;
v( abs(v) > 3.5 ) = NaN;
w( abs(w) > 3.5 ) = NaN;

%% forward cutoff

u(u>-1) = NaN;


%% phase spacce despiking (optional)
%[ u v w spikeinds] = func_despike_phasespace3d_3var( u, v, w, 2 );
%[ u spikeinds] = func_despike_phasespace3d( u, 0 );


%%

save([fpath fname '.mat'])

%% plots
% 
pts = length(u);

figure(1), clf
s(1) = subplot(2,1,1);
plot([1:pts]./samplingrate,-u,'.')%1:pts,v,1:pts,w)
title(fname,'interpreter','none'),
ylabel('u [m/s]')
xlabel('time [s]')
%legend('u','v','w')

s(2) = subplot(2,1,2);
plot([1:pts]./samplingrate,pres,'k'), hold on
ylabel('depth [m]')
xlabel('time [s]')
print('-dpng',[fpath fname '_raw.png'])

linkaxes(s,'x')

%% spectra

% figure(2),clf
% 
% cmap = colormap;
%     
% for i = 1:(pts/(samplingrate*60)),
%     burstu = u( (i-1)*samplingrate*60 + [1:(samplingrate*60)] );
%     meanu = nanmean(abs(burstu));
%     burstp = pres( (i-1)*samplingrate*60 + [1:(samplingrate*60)] );
%     good = sum(~isnan(burstu));
%     if good/(samplingrate*60) > .8,
%     [upsd f] = pwelch(detrend(burstu(~isnan(burstu))),[],[],[],samplingrate);
%     %[ppsd f] = pwelch(detrend(burstp),[],[],[],samplingrate);
%     loglog(f,upsd,'color',cmap(ceil(meanu/3*64),:) ), hold on
%     %loglog(f,ppsd,'r'), hold on
%     else
%     end
% end
% %colorbar
% xlabel('freq [Hz]'),
% ylabel('TKE')
% title(fname)
% 
% print('-dpng',[fpath fname '_spectra.png'])
% 
%     
