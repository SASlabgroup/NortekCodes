% script to process read and process AQD wave data
% stored as individual .wad burst files
% assuming velocites are already east-north-up
%
% calls a sub routine from the SWIFT codes for the actual processing
%   assuming deep water
% and stores results in a SWIFT compatible structure
%   actually three structures, for ast1, ast2, and pressure
%
% ancillary processing includes range - pressure difference to determine ice thickness
%
% J. Thomson,  4/2017  ** PRESSURE ATTENTUATION CORRECTION STILL NOT APPLIED ***
%

clear all, close all

astqualitycutoff = 120;
burstplots = false;  % binary flag for plotting individual burst results (add a lot of time)
readraw = false; % binary flag to force reading raw .wad files, even if .mat exist
lat = []; lon = []; 
%lat = 34.8980; lon = -120.6502; % APL04 at Inner Shelf
%lat = 34.8862; lon = -120.6444; % APL05 at Inner Shelf


doff = 0.7;  % distance instrument is off the bottom

fpath = './'

flist = dir([fpath '*.wad'])


% loop thru burst files, reading raw or loading existing mat, then processing

for fi=1:length(flist),
    
    disp([ num2str(fi) ' of ' num2str(length(flist)) ])
    
    matfile = dir([ flist(fi).name '.mat']);
    
    %% read in data
    
    if readraw | isempty(matfile),
        
        % fixed mode
        [ mo day yr hr minute sec pres spare analog1 u v w vel s1 s2 s3 s4x] = textread([fpath flist(fi).name]); heading = []; pitch = []; roll = [];
        
        % SUV mode (with raw pitch roll heading)
        %[ mo day yr hr minute sec pres ast1 ast2 astq heading pitch roll u v w amp1 amp2 amp3] = textread([fpath flist(fi).name]);
        
        save([flist(fi).name '.mat'], 'mo', 'day', 'yr', 'hr', 'minute', 'sec', 'pres',  'heading', 'pitch', 'roll', 'u', 'v', 'w')
        
    else
        
        load(matfile.name);
        
    end
    

    
    
    %% calc basic stats
    
    time = datenum(yr(1),mo(1),day(1),hr(1),minute(1),sec(1));
    depth = mean(pres);
    waterdepth = depth + doff;

    avgheading = mean(heading);
    avgpitch = mean(pitch);
    avgroll = mean(roll);
    
 
    rate = 1./median(diff(sec));
    
    
    %% pressure wave processing
    [ Hs, Tp, Dp, E, f, a1, b1, a2, b2, check] = UVZwaves(-u, -v, pres, rate); % assumes v1 and v2 are already east and north
    
    AQDpres(fi).time = time; AQDpres(fi).lat = lat; AQDpres(fi).lon = lon;
    AQDpres(fi).sigwaveheight = Hs;
    AQDpres(fi).peakwaveperiod = Tp;
    AQDpres(fi).peakwavedirT = Dp;
    AQDpres(fi).wavespectra.energy = E';
    AQDpres(fi).wavespectra.freq = f';
    AQDpres(fi).wavespectra.a1 = a1';
    AQDpres(fi).wavespectra.b1 = b1';
    AQDpres(fi).wavespectra.a2 = a2';
    AQDpres(fi).wavespectra.b2 = b2';
    AQDpres(fi).wavespectra.check = check';
    
    
    %% should need correct pressure spectra for attenutation, but this over-corrects
    for j=1:length(f), 
        k(j) = wavenumber( f(j), waterdepth );
    end
    attenuation = cosh( k .* waterdepth ) ./ cosh( k.*(-depth+waterdepth) ) ;
        %attenuation = exp(k*depth); % should be same as cosh for deep water
    attenuation = attenuation.^2; % square for energy
    attenuation( attenuation > 100 ) = NaN; % cut if off, don't amplify noise
    E = E.*attenuation;
    Hs = 4 * sqrt( nansum(E) * median(diff(f)) );
    AQDpres(fi).sigwaveheight = Hs;
    AQDpres(fi).wavespectra.energy = E';
    
    
    
end

%% sort

[time sortinds] = sort([AQDpres.time]);
AQDpres = AQDpres(sortinds);


%% QC and make a final product

AQD = AQDpres;  % buld final product from 


%% save final

save AQDwaves.mat AQD

%% plot final
figure(3), clf, n = 3;


ax(1) = subplot(n,1,1);
plot( [AQD.time],[AQD.sigwaveheight],'k+','linewidth',2)
set(gca,'Fontsize',16,'fontweight','demi')
datetick
ylabel([ 'H_s [m]'])
set(gca,'Ylim',[0 ceil(max([AQD.sigwaveheight])) ])
title(pwd,'interpreter','none','fontsize',10)

ax(2) = subplot(n,1,2);
plot( [AQD.time],[AQD.peakwaveperiod],'k+','linewidth',2)
set(gca,'Fontsize',16,'fontweight','demi')
datetick
ylabel([ 'T_p [s]'])
set(gca,'Ylim',[0 20])

ax(3) = subplot(n,1,3);
plot([AQD.time],[AQD.peakwavedirT],'k+','linewidth',2), hold on
set(gca,'Fontsize',16,'fontweight','demi')
datetick
ylabel('D_p [^\circ T]')
set(gca,'Ylim',[0 360])
set(gca,'YTick',[0 180 360])

linkaxes(ax,'x')
set(gca,'XLim',[min([AQD.time]) max([AQD.time])] )
print('-dpng',['wavestats.png'])

figure(4), clf
for ai=1:length(AQD),
    loglog(AQD(ai).wavespectra.freq,AQD(ai).wavespectra.energy,'k-'), hold on
end
set(gca,'Fontsize',16,'fontweight','demi')
xlabel('freq [Hz]')
ylabel('Energy [m^2/Hz')
title('Scalar wave spectra')
title(pwd,'interpreter','none','fontsize',10)
print -dpng spectra.png