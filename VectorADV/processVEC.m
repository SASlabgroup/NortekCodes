% process Vector data (continuous sampling)
% 32 Hz sampling
%
% J. Thomson, 2011 (TTT version)
%               2022 (generalize)

clear all, close all

filebase = 'VEC_402';
fpath = './';
fs = 32; % Hz 
starttime = datenum(2022,8,9,16,0,0); % 8/9/2022 4:00:00 PM

data = load([fpath filebase '.dat']);

u = data(:,3);
v = data(:,4);
w = data(:,5);

%% scalar magnitude

fastmag = sqrt(u.^2 + v.^2);

%% full time
fastt = [1/fs:1/fs:length(fastmag)./fs];
%fastt = starttime + fastt ./ (24*3600);

%% make a smoothed magnitude product 

dt = 600; % s
slowt = decimate( fastt, dt*fs);
slowmag = decimate( fastmag, dt*fs);
slowu = decimate( u, dt*fs);
slowv = decimate( v, dt*fs);
slowmag = sqrt(slowu.^2 + slowv.^2);

%% save as mat 

save([fpath filebase '.mat'],'u','v','w','fastt','fastmag','slowt','slowmag','starttime')
save([fpath filebase '_smoothed.mat'],'slowt','slowmag','starttime')

%%

figure(1), clf

subplot(2,1,1),
plot(fastt,fastmag,'b.',slowt,slowmag,'r-'), hold on,
set(gca,'FontSize',16,'FontWeight','demi')
%axis([0 24 0 2]),
ylabel('[m/s]'), 
xlabel('Time [s]'),
legend('raw',['dt=' num2str(dt,3)])

print('-dpng',[fpath filebase '.png'])

