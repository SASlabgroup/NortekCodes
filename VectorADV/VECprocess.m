% Matlab script to run processing on multiple burst files of Nortek Vector ADV data
%
% J. Thomson, 10/2002, rev 9/2008, 6/2010, 11/2011
%
%

clear, close all

filepath = './';
fileprefix = 'AgatePass_Vector_16Jul2024_VR0_3_';

files = dir([ filepath fileprefix '*.dat' ]); % note that this file list may not be sequential because of digit changes in naming

nf = length(files);

samples = 19200;
sample_rate = 64; % Hz
fftwindow = 4096;


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

for i = 2:length(time)

    % read data

    data = load([filepath files(i).name]);
    uraw = data(:,3);
    vraw = data(:,4);
    wraw = data(:,5);
    amp1 = data(:,6);
    amp2 = data(:,7);
    amp3 = data(:,8);
    cor1 = data(:,12);
    cor2 = data(:,13);
    cor3 = data(:,14);
    p = data(:, 15);
    t = [1:length(uraw)]./sample_rate; % seconds

    filenum = str2num(files(i).name((length(fileprefix)+1):(length(files(i).name)-4)));  % determine place in squence

    figure(1),clf % raw data
    subplot(4,1,1)
    plot(t,p, 'k.')
    legend, xlabel('t [secs]'), ylabel('m')
    title(files(i).name,'Interpreter','none')

    subplot(4,1,2)
    plot(t,uraw,'r.',t,vraw,'g.',t,wraw,'b.')
    legend, xlabel('t [secs]'), ylabel('m/s')

    subplot(4,1,3)
    plot(t,amp1,'r.',t,amp2,'g.',t,amp3,'b.')
    legend, xlabel('t [secs]'), ylabel('Amp')

    subplot(4,1,4)
    plot(t,cor1,'r.',t,cor2,'g.',t,cor3,'b.')
    legend, xlabel('t [secs]'), ylabel('%')

    print('-dpng',[filepath files(i).name '_rawtimeseries.png'])

    % quality control
    u = uraw; v=vraw; w=wraw;

    % low beam correlations
    %corcutoff = 30 + 40*sqrt(sample_rate/25);  % Elgar et al 2001
    corcutoff = 50;
    lowcor = find( cor1 < corcutoff | cor2 < corcutoff | cor3 < corcutoff);
    badpoints(filenum-1) = length(lowcor);
    %u(lowcor) = NaN;
    %v(lowcor) = NaN;
    %w(lowcor) = NaN;


    if length(lowcor) < 0.1*length(u)
        %[ubari, vbari, wbari, uprimei, vprimei, wprimei, upsdi, vpsdi, wpsdi, f, uwi, vwi, depthi, percentbadi] = VECanalysis(filepath, files(i).name, time(filenum+1));

        % spectra
        [Suu f] = pwelch(u, fftwindow, [], [], sample_rate);
        [Svv f] = pwelch(v, fftwindow, [], [], sample_rate);
        [Sww f] = pwelch(w, fftwindow, [], [], sample_rate);

        figure(2)
        loglog(f,Suu,'r', f,Svv,'g', f,Sww,'b')
        hold on

        psdu(filenum-1,:) = Suu;
        psdv(filenum-1,:) = Sww;
        psdw(filenum-1,:) = Sww;

        % means
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
save([fileprefix '_means_and_spectra.mat'],'mean*','time','depth','f','psd*')
%

%% summary plots

% all spectra 
figure(2),
set(gca,'fontsize',16,'fontweight','demi')
ylabel('TKE [m^2/s^2/Hz]')
xlabel('Frequency [Hz]')
legend('b1','b2','b3')
title(fileprefix,'interpreter','none')
print('-dpng',[fileprefix '_spectra.png'])

% means
figure(3), clf
subplot(2,1,1)
plot(time,depth),
datetick
ylabel('Depth [m]')
title(fileprefix,'interpreter','none')
set(gca,'fontsize',16,'fontweight','demi')

subplot(2,1,2)
plot(time,meanu,'r',time,meanv,'g',time,meanw,'b'),
datetick,
ylabel('Currents [m/s]')
legend('b1','b2','b3')
set(gca,'fontsize',16,'fontweight','demi')
print('-dpng',[fileprefix '_means.png'])


%% spectra sorted by mean flow
ds = 0.1;
spdgrid = [0:ds:0.5];
figure(4),clf
cmap = colormap;
for si=1:length(spdgrid)
    si
    inds = find( abs(meanu-spdgrid(si)) < ds/2 )
    meanpsd = nanmean( psdu(inds,:) );
    ci = round(si / length(spdgrid) * 256); 
    loglog(f, meanpsd,'color',cmap(ci,:) ,'linewidth',2), hold on

        inds = find( abs(meanv-spdgrid(si)) < ds/2 )
    meanpsd = nanmean( psdu(inds,:) );
    ci = round(si / length(spdgrid) * 256); 
    loglog(f, meanpsd,'color',cmap(ci,:) ,'linewidth',2), hold on

        inds = find( abs(meanw-spdgrid(si)) < ds/2 )
    meanpsd = nanmean( psdu(inds,:) );
    ci = round(si / length(spdgrid) * 256); 
    loglog(f, meanpsd,'color',cmap(ci,:) ,'linewidth',2), hold on

end
set(gca,'fontsize',16,'fontweight','demi')
title(fileprefix,'interpreter','none')
ylabel('TKE [m^2/s^2/Hz]'), xlabel('Frequency [Hz]')
cb = colorbar('location','north');cb.Limits = [spdgrid(1) spdgrid(end)];
caxis([spdgrid(1) spdgrid(end)])
cb.Label.String = 'm/s';
print('-dpng',[fileprefix '_sortedspectra.png'])




%% OLD
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
