% script to process TTT ADV data from Feb 2011
%
% J. Thomson, 5/2011
tic

clear all
% load QC data from Vibhav
load /Volumes/Data/Marrowstone/TTT_Vector_Feb2011/QC_ADV_1.mat
ADV.n = 0.02 % Doppler Noise (1% of velocity scale)
ADV.fs = 32; % sampling frequency (Hz)

%% phase spacce despiking (optional, data already despiked from Vibhav)
[ u v w spikeinds] = func_despike_phasespace3d_3var( u, v, w, 2 );

%% parse into bursts 
burstlength = 300 % burst length, secs

bursttime = reshape(time(1:(length(time)-rem(length(time),ADV.fs*burstlength))),ADV.fs*burstlength,[]);
burstu = reshape(u(1:(length(time)-rem(length(time),ADV.fs*burstlength))),ADV.fs*burstlength,[]);
burstv = reshape(v(1:(length(time)-rem(length(time),ADV.fs*burstlength))),ADV.fs*burstlength,[]);
burstw = reshape(w(1:(length(time)-rem(length(time),ADV.fs*burstlength))),ADV.fs*burstlength,[]);

ADV.time = mean(bursttime);

%% remove trends (tidal signal)
for bi=1:length(ADV.time),
    ufit = polyfit([1:ADV.fs*burstlength]',sqrt( burstu(:,bi).^2 + burstv(:,bi).^2 ),2);
    utrend(:,bi) = ufit(1)*[1:ADV.fs*burstlength]'.^2 + ufit(2)*[1:ADV.fs*burstlength]' + ufit(3);
    thetafitrads = polyfit([1:ADV.fs*burstlength]',abs(atan2(burstv(:,bi),burstu(:,bi))),2);
    thetatrendrads(:,bi) = thetafitrads(1)*[1:ADV.fs*burstlength]'.^2 + thetafitrads(2)*[1:ADV.fs*burstlength]' + thetafitrads(3);
end

%% determine SPEED: mean, standard deviation, and intensity for each burst
% ** scalar speed, assumed along principle axis **
ADV.ubar = nanmean( sqrt( burstu.^2 + burstv.^2 ) ); 
ADV.uprime = nanstd( sqrt( burstu.^2 + burstv.^2 ) - utrend );
ADV.uprimemax = max( sqrt( burstu.^2 + burstv.^2 ) - utrend );
ADV.ustddev = nanstd(burstu);
ADV.vstddev = nanstd(burstv);
%ADV.Iu = 100 * sqrt(ADV.uprime.^2 - n^2 )./ADV.ubar;


%% determine DIRECTION: mean, standard deviation, and intensity for each burst
% ** still need to correct for noise **
ADV.thetabar = rad2deg( nanmean( atan2(burstv,burstu) ) );
ADV.thetaprime = rad2deg( nanstd( abs(atan2(burstv,burstu)) - thetatrendrads ) );  % abs gets rid of wrapping at 180 deg
%ADV.Itheta = 100 * abs(ADV.thetaprime)./90 ; % normalize by 90, which would be complete off-axis flow


%% remove slack periods
cutin = 0.8
slack = find(ADV.ubar < cutin);
%ADV.Iu(slack) = NaN;
%ADV.Itheta(slack) = NaN;

%% frequency spectra and dissipation 
tic

K = 0.69 ; % Kolmogorov const (vertical)
rho = 1024; % water density [kg/m^3]


for i=1:length(ADV.time),
    
    %horizontal
    [psd f] = pwelch( sqrt( burstu(:,i).^2 + burstv(:,i).^2 ) - ADV.ubar(i), [], [], [], 32);
    psd(1) = []; f(1) = []; % remove mean
    ADV.Suuf(:,i) = psd;  % power spectral density of velocity mag, m^2 s^-2 Hz^-1
    
    [upsd f] = pwelch(  detrend(burstu(:,i)) , [], [], [], 32);
    upsd(1) = []; f(1) = []; % remove mean
    ADV.Surawf(:,i) = upsd;  % power spectral density of velocity mag, m^2 s^-2 Hz^-1
    
    [vpsd f] = pwelch(  detrend(burstv(:,i)) , [], [], [], 32);
    vpsd(1) = []; f(1) = []; % remove mean
    ADV.Svrawf(:,i) = vpsd;  % power spectral density of velocity mag, m^2 s^-2 Hz^-1
    
    
    % vertical
    [wpsd f] = pwelch( burstw(:,i) - mean(burstw(:,i)), [], [], [], 32);
    wpsd(1) = []; f(1) = []; % remove mean
    ADV.Swwf(:,i) = wpsd;  % power spectral density of velocity mag, m^2 s^-2 Hz^-1
    
    % dissipation
    finds = find( f > 0.1 & f < 1 );
    warning('off','stats:statrobustfit:IterationLimit')
    [speedfit stats] = robustfit(f(finds).^(-5/3),wpsd(finds),'bisquare',1,'off');
    ADV.epsilon_old(i) = speedfit.^(3/2) * 2;  % m^2 /s^3
    ADV.epsilon(i) = ( speedfit ./ ( ( ADV.ubar(i) ./ (2*pi) ).^(2/3)  .* K ) )^(3/2); % m^2/s^3
end

%ADV.Suuf(:,slack)= NaN;
%ADV.Swwf(:,slack)= NaN;
ADV.df = median(diff(f));
ADV.f = f;

toc

%% wavenumber spectra and fractional turbulent intensity
l = (ones(length(f),1)*ADV.ubar) ./ (f*ones(1,length(ADV.ubar))) ; % length scale, m, from Taylor hypothesis
psdl = ADV.Suuf ./ (ones(length(f),1)*ADV.ubar) .* (f*ones(1,length(ADV.ubar))).^2; % power spectral density, m^2 s^-2 m^-1
[da dl] = gradient(l); dl = abs(dl); % length bandwidths (non-uniform)

% integrate variance and map to uniform length bands
ADV.dL = 1;
ADV.L = [ADV.dL:ADV.dL:(max(ADV.ubar)./min(f))]; % uniform length scale bins, m

for i=find(ADV.ubar>cutin),
    ADV.SuuL(:,i) = interp1(l(:,i), psdl(:,i).*dl(:,i), ADV.L); % integrated velocity variance in each uniform length band, m^2 s^-2
end
ADV.SuuL(:,slack)= NaN;

ADV.iul = 100 * sqrt(ADV.SuuL) ./ (ones(length(ADV.L),1)*ADV.ubar);


%% Reynolds stresses

% componet based
for bi=1:length(ADV.time),
    rotation = rad2deg(atan2( mean( burstu(:,bi)), mean(burstv(:,bi)) ) ) - 90;
    along = burstu(:,bi) .* cos(deg2rad(rotation))   -   burstv(:,bi) .* sin (deg2rad(rotation));;
    across = burstu(:,bi) .* sin(deg2rad(rotation))   +   burstv(:,bi) .* cos (deg2rad(rotation));
    ADV.along(bi) = mean(along);
    ADV.across(bi) = mean(across);
    ADV.Reynolds(bi) = mean( detrend(along) .* detrend(burstw(:,bi)) );
end

% magnitude based
ADV.Reynolds_alt = mean( (burstu - ones(ADV.fs*burstlength,1)*mean(burstu)).* (burstw - ones(ADV.fs*burstlength,1)*mean(burstw)) );


%% save results 

save(['TTT_ADV_Feb2011processed_' num2str(burstlength) 's.mat'],'ADV') 

%% plots

%clear all
%load(['TTT_' num2str(burstlength) 'sprocessedADV_Feb2011.mat'],'ADV')

% figure(1),clf
% hist([ADV.Iu; ADV.Itheta;]',[0:1:100]),
% axis([0 20 0 inf]),
% legend('I_u [%]','I_\theta [%]')
% print -dpng ADV_intensity.png

figure(1), clf
plot(ADV.ubar,ADV.uprime,'.')
set(gca,'Fontsize',14,'Fontweight','demi')
xlabel('<u> [m/s]'),ylabel('u\prime [m/s]')
title('Nodule Point ADV')
%print -dpng ADV_uprime_vs_ubar.png

%%

ADV.Swwf(:, find(ADV.ubar<0.8) )= NaN;
ADV.Suuf(:, find(ADV.ubar<0.8) )= NaN;

figure(2), clf

loglog(ADV.f,nanmean(ADV.Suuf'),'r','linewidth',5), hold on
loglog(ADV.f,nanmean(ADV.Swwf'),'--','linewidth',3,'color',[.4 .1 .1]), hold on

loglog(ADV.f,ADV.Swwf(:,10:1:1250),'color',[.7 .5 .5],'linewidth',.5), hold on
loglog(ADV.f,ADV.Suuf(:,10:1:1250),'color',[1 0.5 .5],'linewidth',.5), hold on

loglog(ADV.f,nanmean(ADV.Suuf'),'r','linewidth',5), hold on
loglog(ADV.f,nanmean(ADV.Swwf'),'--','linewidth',3,'color',[.4 .1 .1]), hold on

loglog([.2 5],.01*[.2 5].^(-5/3),'k:','linewidth',2), hold on
text(.6,.06,'f^{-5/3}','Fontsize',14,'Fontweight','bold')
%text(0.02,3e-5,['df = ' num2str(ADV.df,2) ' Hz'],'Fontsize',14,'Fontweight','demi')

set(gca,'Fontsize',14,'Fontweight','demi')
xlabel('Frequency [Hz]'),ylabel('Turbulent Kinetic Energy [m^2/s^2 Hz^{-1}]'),
title('Nodule Point ADV, z_{hub} = 4.7 m'),
axis([min(ADV.f) 10 1e-5 1]) 
legend('horizonal','vertical')

dn = text(6,8e-5,'Doppler noise','Fontsize',12,'Fontweight','bold');
iso = text(.2,1e-3,'isotropic turbulence','Fontsize',12,'Fontweight','bold');
aniso = text(.01,5e-2,'anisotropic eddies','Fontsize',12,'Fontweight','bold');

%print -dpdf TTT-ADV-tkespectra-freq.pdf
%print -depsc TTT-ADV-tkespectra-freq.eps


% figure(2), clf
% loglog(l,tkepsd,'color',[.7 .7 .7]), hold on
% set(gca,'Fontsize',20,'Fontweight','demi')
% xlabel('Length scale [m]'),ylabel('Turbulent Kinetic Energy [m^2/s^2 Hz^{-1}]'),title('Marrowstone TTT-ADV'),
% loglog(nanmean(tkepsd.*l,2)./nanmean(tkepsd,2),nanmean(tkepsd,2),'k','linewidth',5), hold on
% loglog([2 20],.0002*[2 20].^(5/3),'r--','linewidth',2), hold on
% text(20,.01,'L^{5/3}','color',[1 0 0],'Fontsize',20,'Fontweight','bold'),
% text(2,.5,['dL = <v> f^{-2} df'],'Fontsize',16,'Fontweight','demi')
% axis([1e-2 1e2 1e-5 1])
% print -dpng TTT-ADV-tkespectra-length.png


%%


figure(3), clf
%Linds = find(ADV.L<150);
%Lfit = polyfit(ADV.L(Linds).^(2/3),nanmean(ADV.iul(Linds,:),2)',1);
plot(ADV.L,ADV.iul,'color',[1 0.5 .5]), hold on, grid
set(gca,'Fontsize',14,'Fontweight','demi')
xlabel('Length scale, L [m]'),ylabel('Fractional turbulent intensity, i_u(L) [%]'),
title('Nodule Point ADV, z_{hub} = 4.7 m'),
plot(ADV.L,nanmean(ADV.iul,2),'r','linewidth',5), hold on
%plot(ADV.L(Linds),Lfit(1)*ADV.L(Linds).^(2/3),'r--','linewidth',2), hold on
%text(9,1.8,'L^{2/3}','color',[1 0 0],'Fontsize',14,'Fontweight','bold')%,'fontname','times')
%text(8,7,['dL = ' num2str(ADV.dL,2) ' m'],'Fontsize',14,'Fontweight','demi')
axis([0 200 0 10])
text(50,3.5,'anisotropic eddies','fontsize',12,'fontweight','bold'),
text(5,1,'isotropic turbulence','fontsize',12,'fontweight','bold')
%print -dpdf TTT-ADV-intensity-length.pdf
%print -depsc TTT-ADV-intensity-length.eps


%%
figure(4),clf
plotclr(ADV.thetabar,ADV.thetaprime,ADV.ubar)
xlabel('<\theta> [deg]')
ylabel('\theta\prime [deg]')
%print -dpng ADV_directions.png

toc