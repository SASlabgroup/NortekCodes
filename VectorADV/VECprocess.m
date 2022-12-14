% Matlab script to run processing on multiple burst files of Nortek Vector ADV data
%
% J. Thomson, 10/2002, rev 9/2008, 6/2010, 11/2011 
%
%

clear

filepath = './';
fileprefix = 'SouthPiling_ADV_Apr2018_';

files= dir([ filepath fileprefix '*.dat' ]) % note that this file list may not be sequential because of digit changes in naming

nf = length(files);

samples = 16384;
sample_rate = 16; % Hz


%% timestamps
vhd = load( [ filepath fileprefix '.vhd' ]);  % one pt per burst, use sen file for 1 Hz  instead
month = vhd(:, 1);
day = vhd(:, 2);
year = vhd(:, 3);
hour = vhd(:, 4);
minute = vhd(:, 5);
second = vhd(:, 6);
time = datenum( year, month, day, hour, minute, second);


%% analysis loop

for i = 2:length(time),
    
    % read data

data = load([filepath files(i).name]);
uraw = data(:,3);
vraw = data(:,4); 
wraw = data(:,5); 
snr1 = data(:,6);
snr2 = data(:,7); 
snr3 = data(:,8); 
cor1 = data(:,12);
cor2 = data(:,13); 
cor3 = data(:,14); 
p = data(:, 15);

% quality control 
u = uraw; v=vraw; w=wraw;

% low beam correlations
corcutoff = 30 + 40*sqrt(sample_rate/25);  % Elgar et al 2001
lowcor = find( cor1 < corcutoff | cor2 < corcutoff | cor3 < corcutoff); 
u(lowcor) = NaN;
v(lowcor) = NaN;
w(lowcor) = NaN;
    

    filenum = str2num(files(i).name((length(fileprefix)+1):(length(files(i).name)-4)));  % determine place in squence
if length(lowcor) < 0.1*length(u),
    %[ubari, vbari, wbari, uprimei, vprimei, wprimei, upsdi, vpsdi, wpsdi, f, uwi, vwi, depthi, percentbadi] = VECanalysis(filepath, files(i).name, time(filenum+1));
    meanu(filenum-1) = nanmean(u);
    meanv(filenum-1) = nanmean(v);
    meanw(filenum-1) = nanmean(w);
    depth(filenum-1) = nanmean(p)+0.6;
else
    meanu(filenum-1) = NaN;
    meanv(filenum-1) = NaN;
    meanw(filenum-1) = NaN;
    depth(filenum-1) = nanmean(p)+0.6;
end
   
end

time(1)=[];

% %% TKE
% 
% slackcutoff = 0.0;
% slack = find( sqrt(ubar.^2 + vbar.^2) < slackcutoff );
% notslack = find( sqrt(ubar.^2 + vbar.^2) > slackcutoff );
% 
% % turbulent intensity
% n = 0.04;  % Doppler noise from Nortek specs... 1% of 4 m/s range
% I = 100* (sqrt(uprime.^2 + vprime.^2) - n) ./ sqrt(ubar.^2 + vbar.^2);
% I( I <=0 ) = NaN; % remove any negative values resulting from noise correction
% 
% % dissipation rate 
% finds = find(f>.5 & f < 5);
% [wfit stats] = robustfit(f(finds).^(-5/3),wpsd(finds),'bisquare',1,'off');
% epsilon(fi) =  ( wfit ./ ( ( sqrt(ubar(fi)) ./ (2*pi) ).^(2/3)  .* K ) )^(3/2); % W/m^3
%     slopeerror(fi) = stats.se;
%     epsilonerror(fi) = rho* ( (wfit+stats.se) ./ ( ( meanspeed(fi) ./ (2*pi) ).^(2/3)  .* K ) )^(3/2); % W/m^3
% 
% 
% 
% fids = find( f < 2 & f > .5);
% for i=1:nf,
%     fit = polyfit( f(fids).^(-5/3), squeeze(wpsd(fids,i)), 1 );
%     if fit(1)>0,
%         epsilon(i) = fit(1).^(3/2) * 2;  % W kg^{-1}
%     else
%         epsilon(i)=NaN;
%     end
% end
% density = 1030; % for converson to W/m^3
% epsilon = density*epsilon;

%% save output

% save([fileprefix '_raw.mat'])
save([fileprefix '_means.mat'],'mean*','time','depth')
% 
% %% plots
% 
figure(1), clf
subplot(2,1,1)
plot(time,depth),
datetick
ylabel('Depth [m]')
title(fileprefix,'interpreter','none')

subplot(2,1,2)
plot(time,meanu,time,meanv,time,meanw),
datetick,
ylabel('Currents [m/s]')
legend('u','v','w')

print('-dpng',[fileprefix '_means.png'])

% figure(3), clf % time series
% subplot(4,1,1)
% plot(time,ubar,'r',time,vbar,'g',time,wbar,'b'),
% %set(gca,'FontSize',14,'Fontweight','demi'),
% datetick, 
% %axis tight
% l1 = legend('$\bar{u}$','$\bar{v}$','$\bar{w}$');
% set(l1,'Interpreter','latex')
% ylabel('m s^{-1}'),title([fileprefix],'interpreter','none')
% 
% subplot(4,1,2)
% plot(time,uprime,'r',time,vprime,'g',time,wprime,'b'),
% %set(gca,'FontSize',14,'Fontweight','demi'),
% datetick, 
% %axis tight
% l2=legend('$u\prime$','$v\prime$','$w\prime$');
% set(l2,'Interpreter','latex')
% ylabel('m s^{-1}')
% 
% subplot(4,1,3)
% plot(time,I,'k'), hold on
% plot(time(slack), I(slack),'color',[.8 .8 .8])
% %set(gca,'FontSize',14,'Fontweight','demi'),
% datetick, 
% %axis tight
% l3 = legend('$I_{u,v}$');
% set(l3,'Interpreter','latex')
% ylabel('%'), set(gca,'YLim',[0 50])
% 
% subplot(4,1,4)
% semilogy(time,epsilon,'k'), hold on
% %plot(time(slack), epsilon(slack),'color',[.8 .8 .8])
% datetick, 
% %axis tight
% l4 = legend('$\epsilon$');
% set(l4,'Interpreter','latex'),
% ylabel('W kg^{-1}'),
% set(gca,'YTick',[1e-10, 1e-5],'YLim',[1e-11,1e-3])
% 
% print('-dpng',[fileprefix 'results.png'])
% print('-depsc',[fileprefix 'results.eps'])
% 
%  
% figure(4), clf % spectra
% lowf = 5;
% 
% colormap(flipud([0:.1:1]'*[1 1 1])),
% 
% subplot(3,1,1)
% pcolor(time,f(lowf:length(f)),log(upsd(lowf:length(f),:))), 
% set(gca,'FontSize',14,'Fontweight','demi'),
% shading flat, title('u'),
% datetick, axis tight
% caxis([-10 -4]),
% ylabel('Hz')
% 
% subplot(3,1,2)
% pcolor(time,f(lowf:length(f)),log(vpsd(lowf:length(f),:))), 
% set(gca,'FontSize',14,'Fontweight','demi'),
% shading flat, title('v')
% datetick, axis tight
% caxis([-10 -4]),
% ylabel('Hz')
% 
% subplot(3,1,3)
% pcolor(time,f(lowf:length(f)),log(wpsd(lowf:length(f),:))), 
% set(gca,'FontSize',14,'Fontweight','demi'),
% shading flat, title('w')
% datetick, axis tight
% caxis([-10 -4]),
% ylabel('Hz')
% 
% print('-dpng',[fileprefix 'spectra.png'])
% print('-depsc',[fileprefix 'spectra.eps'])
% 
% 
% figure(5), clf
% 
% subplot(2,1,1),
% semilogy(sqrt(ubar.^2 + vbar.^2),I,'kx'),
% xlabel('$\sqrt{ u^2 + v^2}, m s^{-1}$','Interpreter','latex'),
% ylabel('I, %'),
% set(gca,'YLim',[5e0 1e3])
% 
% subplot(2,1,2),
% semilogy(sqrt(ubar.^2 + vbar.^2),epsilon,'k+'),
% xlabel('$\sqrt{ u^2 + v^2}, m s^{-1}$','Interpreter','latex'),
% ylabel('\epsilon, W kg^{-1}'),
% set(gca,'YLim',[1e-10 1e-3])
% 
% print('-dpng',[fileprefix 'I_epsilon.png'])
% print('-depsc',[fileprefix 'I_epsilon.eps'])
