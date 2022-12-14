function [ubar, vbar, wbar, uprime, vprime, wprime, upsd, vpsd, wpsd, f, uw, vw, depth, percentbad] = VECanalysis(filepath, filename, time);
% 
% Matlab function to analysize Nortek Vector ADV data
% (using ASCII files extracted from .vec files with AquaPro software)
% includes basic quality control, ASSUMING WIDE STATIONARITY (steady mean)
%
% J. Thomson, 10/2002, rev. 10/2008, 6/2010
%
% this version specific to Tidal Turbulence Tripod
% 


%filepath = '/Volumes/Data/MarrowstoneT3/T3_May2010_NodulePoint-VectorADV/';
%filename = 'T3_NodulePoint_May20102200';

samples = 2048;
pts = [1:samples];
sample_rate = 32; % Hz


%% read data

data = load([filepath filename]);
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



%% quality control 
u = uraw; v=vraw; w=wraw;

% low beam correlations
corcutoff = 30 + 40*sqrt(sample_rate/25);  % Elgar et al 2001
lowcor = find( cor1 < corcutoff | cor2 < corcutoff | cor3 < corcutoff); 
u(lowcor) = NaN;
v(lowcor) = NaN;
w(lowcor) = NaN;

% jumps
jumpstd = 3;  % set higher for weaker despiking
ujump = find( abs(diff(u)) > (jumpstd*nanstd(u)) );
vjump = find( abs(diff(v)) > (jumpstd*nanstd(v)) );
wjump = find( abs(diff(w)) > (jumpstd*nanstd(w)) );
u(ujump) = NaN;
v(vjump) = NaN;
w(wjump) = NaN;

% spikes
spikestd = 6;  % set higher for weaker despiking
uspike = find( abs(u) > (spikestd*nanstd(u) + abs(nanmean(u))) );
vspike = find( abs(v) > (spikestd*nanstd(v) + abs(nanmean(v))) );
wspike = find( abs(w) > (spikestd*nanstd(w) + abs(nanmean(w))) );
u(uspike) = NaN;
v(vspike) = NaN;
w(wspike) = NaN;

% reaplce bad points with burst means
u(isnan(u)) = nanmean(u);
v(isnan(v)) = nanmean(v);
w(isnan(w)) = nanmean(w);


% record how many bad points
percentbad = 100* ( length(lowcor) + mean( [ length(uspike), length(vspike), length(wspike) ]) ) ./ samples;

% save cleaned data, using timestamp as filename
if percentbad < 5, 
    save(['./CleanedData/NodulePoint ADV ' datestr(time) ' 32 Hz.mat'], 'u', 'v', 'w', 'p')
else
end

%% stats 

% means
ubar = mean(u);
vbar = mean(v);
wbar = mean(w);
depth = mean(p);

% fluctuations
uprime = std(u);
vprime = std(v);
wprime = std(w);
pprime = std(p);

% eddy covariance
uw = mean( (u-mean(u)) .* (w-mean(w)) );
vw = mean( (v-mean(v)) .* (w-mean(w)) );


%% quick and dirty spectra
window = 512; % window length in points
[upsd f ]= pwelch(u-mean(u),window,[],[],sample_rate);
[vpsd f ]= pwelch(v-mean(v),window,[],[],sample_rate);
[wpsd f ]= pwelch(w-mean(w),window,[],[],sample_rate);
df = f(2)-f(1);
dof = samples./window * 2;


%% plotting

secs = pts./sample_rate;

figure(1), clf

    subplot(3,1,1), hold on
    plot(secs,uraw,'r.',secs, vraw,'g.',secs,wraw,'b.'), legend('V_x','V_y','V_z')
    plot(secs(lowcor),uraw(lowcor),'kx',secs(lowcor),vraw(lowcor),'kx',secs(lowcor),wraw(lowcor),'kx')
    plot(secs(uspike),uraw(uspike),'kx',secs(vspike),vraw(vspike),'kx',secs(wspike),wraw(wspike),'kx')
    plot(secs(ujump),uraw(ujump),'kx',secs(vjump),vraw(vjump),'kx',secs(wjump),wraw(wjump),'kx')
    set(gca,'FontSize',14,'FontWeight','demi'), 
    ylabel('m/s'),title([filename '   ' datestr(time)],'interpreter','none'),

    subplot(3,1,2), hold on
    plot(secs,cor1,'r.'),plot(secs,cor2,'g.'),plot(secs,cor3,'b.'),
    plot([0 max(secs)],[corcutoff corcutoff],'k--')
    set(gca,'FontSize',14,'FontWeight','demi'), 
    ylabel('cor')

    subplot(3,1,3), hold on
    plot(secs,snr1,'r.'),plot(secs,snr2,'g.'),plot(secs,snr3,'b.'),
    set(gca,'FontSize',14,'FontWeight','demi'), 
    ylabel('snr')

    print('-dpng',['./BurstPlots/' filename '_' datestr(time) '_raw.png'])




figure(2), clf

    loglog(f,upsd,'r',f,vpsd,'g',f,wpsd,'b','linewidth',2),legend('V_x','V_y','V_z')
    hold on
    loglog([.4 1.5],1e-5*[.3 1.5].^(-5/3),'k--'),
    set(gca,'FontSize',14,'FontWeight','demi'), 
    text(.4,1e-5,'f^{-5/3}','FontSize',14,'FontWeight','demi'),
    xlabel('f [Hz]'), ylabel('S [m^2/s^2/Hz]'), 
    title([filename '  ' datestr(time)],'interpreter','none'),

    print('-dpng',['./BurstPlots/' filename '_' datestr(time) '_spectra.png'])
    
    
