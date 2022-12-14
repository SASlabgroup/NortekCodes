% script to process read and process AWAC wave (and ice) data
% stored as individual .wad burst files
% assuming velocites are already east-north-up (**this may not be true**)
%
% calls a sub routine from the SWIFT codes for the actual processing
% and stores results in a SWIFT compatible structure
%   actually three structures, for ast1, ast2, and pressure
%
% ancillary processing includes range - pressure difference to determine ice thickness
%
% J. Thomson,  4/2017   adapted from PUVspectra (2001)
%              8/2019   updated to make a hybrid product,
%                       using pressure for low freq and acoustic surface
%                       tracking (ast) for high freq
%               9/2020  back to using ast exclusively, with simpler despiking
%                       and spectral shapes as quality control
%
% ** DIRECTIONS STILL NOT RELIABLE FROM ARCTIC DEPLOYMENTS**
% FOR FIXED DEPLOYMENTS, NEED TO TRY LOADING HEADING FROM WHD FILE AND MAKING ROTAIONS
% BASED ON THAT.  QUINUALT HEADING WAS 260 DEG, AND RESULTS LOOK REASONABLE
% WITH -U, V.  THAT'S CONSISENT WITH U POSITIVE WEST(ISH) AND V POSITIVE
% SOUTH(ISH).  INNER SHELF APL01 HEADING WAS 290 DEG, SO CAN CHECK WITH
% THAT... SHOULD ALSO CONSIDER AS BEAM VELOCITIES (THOUGH CHECK FACTOR WOULD BE OFF)

clear all

tic

%% processing options
headingoffset = 20;  % AWAC heading relative to 270 (for experimentation with fixed heading rotations)
despike = true; % despike raw data
%astqualitycutoff = 120;
burstplots = false;  % binary flag for plotting individual burst results (adds a lot of time)
readraw = false; % binary flag to force reading raw .wad files, even if .mat exist
extrapequilibriumrange = false; % extrapolate the pressure spectra when beyond noise floor
declination = 15.6 + headingoffset; % deg (positive east)
finalscreening = true;  % optional final screening of results
maxdepth = 50; % max depth (including mooring blow-down) from which to use data
maxwaveperiod = 20; % max wave period allowed during final screening
minwaveperiod = 3; % min wave period allowed during final screening
minicethickness = 0.1; % min ice thickness considered valid during final screening (10 cm is barometric noise)
maxicethickness = 10; % min ice thickness considered valid during final screening (10 cm is barometric noise)
minwaveheight = 0.1; % smallest wave height observable
maxtailshapeexponent = -2;  % max value for f^q in the tail (f> 0.3 Hz)
icebincenters = [0:.5:10]; % meters 

%% Deployment details
%lat = []; lon = [];
lat = 74; lon = -140; % BGEP-D
%lat = 75; lon = -150; % BGEP-A
%lat = 35.0593; lon = -120.6484; % APL01 at Inner Shelf
%lat = 47 + 20.450/60; lon = -( 124 + 18.560);  % AWAC at Quinault

doff = 3000;  % distance instrument is off the bottom (0.7 for a sea spider tripod)

SUV = true;  % SUV mode, set to false for fixed deployments (changes the raw read format)

fpath = './burstfiles/';

wd = pwd;
lastslash = find(wd=='/',1,'last') + 1;
wd = wd( lastslash : end );

flist = dir([fpath '*.wad']);


%% loop thru burst files, reading raw or loading existing mat, then processing

for fi=1:length(flist),
    
    if burstplots,
        figure(1), clf
    end
    
    disp([ num2str(fi) ' of ' num2str(length(flist)) ])
    
    matfile = dir([ fpath flist(fi).name '.mat']);
    
    %% read in data
    
    if readraw | isempty(matfile),
        
        if SUV==false, % fixed mode
            [ mo day yr hr minute sec pres ast1 ast2 astq analog1 u v w amp1 amp2 amp3] = ...
                textread([fpath flist(fi).name]);
            heading = []; pitch = []; roll = [];
            
        elseif SUV==true, % SUV mode (with raw pitch roll heading)
            [ mo day yr hr minute sec pres ast1 ast2 astq heading pitch roll u v w amp1 amp2 amp3] = ...
                textread([fpath flist(fi).name]);
        else
            
        end
        
        save([fpath flist(fi).name '.mat'], 'mo', 'day', 'yr', 'hr', 'minute', 'sec', 'pres', 'ast1', 'ast2', 'astq', 'heading',...
            'pitch', 'roll', 'u', 'v', 'w', 'amp1', 'amp2', 'amp3')
        
    else
        
        load([fpath matfile.name]);
        
    end
    
    %% despike
    % plot before despike
    if burstplots
        figure(1),
        subplot(2,1,1)
        hold on
        plot(ast1,'b.')
        plot(ast2,'r.')
        plot(pres,'y.')
    end
    
    if despike
        ast1 = filloutliers(ast1,'linear');
        ast2 = filloutliers(ast2,'linear');
        u = filloutliers(u,'linear');
        v = filloutliers(v,'linear');
    end
    % plot after despike
    if burstplots
        figure(1),
        subplot(2,1,1)
        hold on
        plot(ast1,'b-')
        plot(ast2,'r-')
        plot(pres,'y-')
    end
    
    %% rotate velocities for declination
    east = u;
    north = v;
    
    u = east .* cos(deg2rad(declination))   -   north .* sin (deg2rad(declination));;
    
    v = east .* sin(deg2rad(declination))   +   north .* cos (deg2rad(declination));
    
    %% calc basic stats
    
    bursttime = datenum(yr(1),mo(1),day(1),hr(1),minute(1),sec(1));
    if burstplots, title(datestr(bursttime)), end
    depth = mean(pres);
    waterdepth = depth + doff;
    range1 = mean(ast1);
    range2 = mean(ast2);
    
    ice1 = depth - range1;
    ice2 = depth - range2;
    
    avgheading(fi) = mean(heading);
    avgpitch(fi) = mean(pitch);
    avgroll(fi) = mean(roll);
    
    alldepth(fi) = depth;
    
    
    %% ast wave processing
    rate = 1./median(diff(sec));
    
    % ast 1
    [ Hs, Tp, Dp, E, f, a1, b1, a2, b2, check] = UVZwaves(-u, v, ast1, rate); % assumes v1 and v2 are already east and north
    
    AWACast1(fi).time = bursttime; AWACast1(fi).lat = lat; AWACast1(fi).lon = lon;
    AWACast1(fi).icethickness = ice1;
    AWACast1(fi).icehistogram.Nobs = hist( depth - ast1 , icebincenters);
    AWACast1(fi).icehistogram.bincenters = icebincenters;
    AWACast1(fi).sigwaveheight = Hs;
    AWACast1(fi).peakwaveperiod = Tp;
    AWACast1(fi).peakwavedirT = Dp;
    AWACast1(fi).wavespectra.energy = E';
    AWACast1(fi).wavespectra.freq = f';
    AWACast1(fi).wavespectra.a1 = a1';
    AWACast1(fi).wavespectra.b1 = b1';
    AWACast1(fi).wavespectra.a2 = a2';
    AWACast1(fi).wavespectra.b2 = b2';
    AWACast1(fi).wavespectra.check = check';
    bandwidthast1(fi) = ( nansum(E) * nansum(E.*f.^2) ./ nansum(E.*f).^2 - 1 ).^.5';
    fit = polyfit(log(f(f>0.3)),log(E(f>0.3)),1);
    tailshapeast1(fi) = fit(1);
    
    if burstplots && ~all(isnan(E))
        figure(1),
        subplot(2,1,2)
        set(gca,'xscale','log','yscale','log')
        hold on
        loglog(f,E)
    end
    
    % ast 2
    [ Hs, Tp, Dp, E, f, a1, b1, a2, b2, check] = UVZwaves(-u, v, ast2, rate); % assumes v1 and v2 are already east and north
    
    AWACast2(fi).time = bursttime; AWACast2(fi).lat = lat; AWACast2(fi).lon = lon;
    AWACast2(fi).icethickness = ice2;
    AWACast2(fi).icehistogram.Nobs = hist( [depth - ast2] , icebincenters);
    AWACast2(fi).icehistogram.bincenters = icebincenters;
    AWACast2(fi).sigwaveheight = Hs;
    AWACast2(fi).peakwaveperiod = Tp;
    AWACast2(fi).peakwavedirT = Dp;
    AWACast2(fi).wavespectra.energy = E';
    AWACast2(fi).wavespectra.freq = f';
    AWACast2(fi).wavespectra.a1 = a1';
    AWACast2(fi).wavespectra.b1 = b1';
    AWACast2(fi).wavespectra.a2 = a2';
    AWACast2(fi).wavespectra.b2 = b2';
    AWACast2(fi).wavespectra.check = check';
    bandwidthast2(fi) = ( nansum(E) * nansum(E.*f.^2) ./ nansum(E.*f).^2 - 1 ).^.5';
    fit = polyfit(log(f(f>0.3)),log(E(f>0.3)),1);
    tailshapeast2(fi) = fit(1);
    
    if burstplots && ~all(isnan(E))
        figure(1),
        subplot(2,1,2)
        hold on
        loglog(f,E)
    end
    
    %% pressure wave processing
    [ Hs, Tp, Dp, E, f, a1, b1, a2, b2, check] = UVZwaves(u, v, pres, rate); % assumes v1 and v2 are already east and north
    
    AWACpres(fi).time = bursttime; AWACpres(fi).lat = lat; AWACpres(fi).lon = lon;
    AWACpres(fi).sigwaveheight = Hs;
    AWACpres(fi).peakwaveperiod = Tp;
    AWACpres(fi).peakwavedirT = Dp;
    AWACpres(fi).wavespectra.energy = E';
    AWACpres(fi).wavespectra.freq = f';
    AWACpres(fi).wavespectra.a1 = a1';
    AWACpres(fi).wavespectra.b1 = b1';
    AWACpres(fi).wavespectra.a2 = a2';
    AWACpres(fi).wavespectra.b2 = b2';
    AWACpres(fi).wavespectra.check = check';
    bandwidthpres(fi) = ( nansum(E) * nansum(E.*f.^2) ./ nansum(E.*f).^2 - 1 ).^.5';
    fit = polyfit(log(f(f>0.3)),log(E(f>0.3)),1);
    tailshapepres(fi) = fit(1);
    
    if burstplots && ~all(isnan(E))
        figure(1),
        subplot(2,1,2)
        hold on
        loglog(f,E)
    end
    
    %% correct pressure spectra for Bernouli?
    % equiparition of energy suggests velocity double pressure amplitude, so 4 times energy
    E = E ./4;
    
    %% correct pressure spectra for depth attenutation,
    %     for j=1:length(f),
    %         k(j) = wavenumber( f(j), waterdepth );
    %     end
    k = (2*pi*f).^2 ./ 9.8;
    attenuation = cosh( k .* waterdepth ) ./ cosh( k.*(waterdepth-depth) ) ;
    %attenuation = exp(k*depth); % should be same as cosh for deep water
    attenuation = attenuation.^2; % square for energy
    noise = attenuation > 100 | isnan(attenuation); % limit the size of the attenuation correction
    attenuation( noise ) = NaN; % cut it off when correction too big, don't amplify noise
    E = E.*attenuation;
    if extrapequilibriumrange % extrapolate the equilibruim range
        E( noise ) = ( E( min(noise) - 1 ) .* f( min(noise) - 1 ).^4 ) .* f(noise).^-4; % extrapolate equilibrium range
    end
    Hs = 4 * sqrt( nansum(E) * median(diff(f)) );
    AWACpres(fi).sigwaveheight = Hs;
    AWACpres(fi).wavespectra.energy = E';
    
    
    if burstplots
        figure(1),
        hold on
        subplot(2,1,2)
        loglog(f,E,'--')
        legend('ast1','ast2','pres','pres, corrected','Location','NorthWest')
        ylabel('E'),xlabel('f')
        title([flist(fi).name ', Hs = ' num2str(AWACast1(fi).sigwaveheight,2) ' m, Ice = ' num2str(nanmean([ice1 ice2]),2) ' m' ]...
            ,'interpreter','none')
        print('-dpng',[fpath flist(fi).name '_spectra.png' ])
    end
    
    
end



%% sort, save, and plot raw
[time sortinds] = sort([AWACast1.time]);
AWACast1 = AWACast1(sortinds);
tailshapeast1 = tailshapeast1(sortinds);
bandwidthast1 = bandwidthast1(sortinds);
save([wd '_ast1_raw.mat'], 'AWACast1')

[time sortinds] = sort([AWACast2.time]);
AWACast2 = AWACast2(sortinds);
tailshapeast2 = tailshapeast2(sortinds);
bandwidthast2 = bandwidthast2(sortinds);
save([wd '_ast2_raw.mat'], 'AWACast2')

[time sortinds] = sort([AWACpres.time]);
AWACpres = AWACpres(sortinds);
tailshapepres = tailshapepres(sortinds);
bandwidthpres = bandwidthpres(sortinds);
save([wd '_pres_raw.mat'], 'AWACpres')

avgpitch =avgpitch(sortinds);
avgroll = avgroll(sortinds);
avgheading = avgheading(sortinds);
alldepth = alldepth(sortinds);

figure(1), clf
plot([AWACast1.time],[AWACast1.sigwaveheight],'bx'), hold on
plot([AWACast2.time],[AWACast2.sigwaveheight],'g+')
plot([AWACpres.time],[AWACpres.sigwaveheight],'ro')
legend('ast1','ast2','pres')
datetick, ylabel('Hs [m]')
print('-dpng',[wd '_rawwaves.png'])

figure(2), clf
subplot(2,1,1)
plot([AWACast1.time],[AWACast1.icethickness],'bx'), hold on
plot([AWACast2.time],[AWACast2.icethickness],'g+')
datetick, ylabel('Ice [m]')
subplot(2,1,2)
plot(time,alldepth,'ks')
set(gca,'YDir','reverse')
datetick, ylabel('Depth [m]')
print('-dpng',[wd '_rawice.png'])


%% use ast1 as final product

AWAC = AWACast1;


%% always remove directions for SUV mode, not reliable yet

if SUV
for ai = 1:length(AWAC),
    AWAC(ai).peakwavedirT = NaN;
end
end

%% final screening of results

if finalscreening,
    
    for ai = 1:length(AWAC),
        
        if  [AWAC(ai).peakwaveperiod] > maxwaveperiod | [AWAC(ai).peakwaveperiod] < minwaveperiod | ...
                [AWAC(ai).sigwaveheight] < minwaveheight | tailshapeast1(ai) > maxtailshapeexponent |...
                isnan([AWAC(ai).sigwaveheight]) | alldepth(ai) > maxdepth,
            
            AWAC(ai).sigwaveheight = NaN;
            AWAC(ai).peakwaveperiod = NaN;
            AWAC(ai).peakwavedirT = NaN;
            AWAC(ai).wavespectra.energy = NaN(size(AWAC(ai).wavespectra.energy));
            AWAC(ai).wavespectra.freq = NaN(size(AWAC(ai).wavespectra.freq));
            AWAC(ai).wavespectra.a1 = NaN(size(AWAC(ai).wavespectra.a1));
            AWAC(ai).wavespectra.b1 = NaN(size(AWAC(ai).wavespectra.b1));
            AWAC(ai).wavespectra.a2 = NaN(size(AWAC(ai).wavespectra.a2));
            AWAC(ai).wavespectra.b2 = NaN(size(AWAC(ai).wavespectra.b2));
            AWAC(ai).wavespectra.check = NaN(size(AWAC(ai).wavespectra.check));
        end
        
        if AWAC(ai).icethickness < minicethickness | AWAC(ai).icethickness > maxicethickness | alldepth(ai) > maxdepth,
            AWAC(ai).icethickness = NaN; %  non-physical ice results
            AWAC(ai).icehistogram.Nobs = NaN(size(icebincenters)); %  non-physical ice results
            AWAC(ai).icehistogram.bincenters = NaN(size(icebincenters)); %  non-physical ice results

        end
    end
end

%% save final

save([wd '.mat'], 'AWAC')


%% plot final
figure(3), clf, n = 3;


ax(1) = subplot(n,1,1);
plot( [AWAC.time],[AWAC.sigwaveheight],'k+','linewidth',2)
set(gca,'Fontsize',16,'fontweight','demi')
datetick
ylabel([ 'H_s [m]'])
%set(gca,'Ylim',[0 ceil(max([AWAC.sigwaveheight])) ])
title(wd,'interpreter','none','fontsize',10)

ax(2) = subplot(n,1,2);
plot( [AWAC.time],[AWAC.peakwaveperiod],'k+','linewidth',2)
set(gca,'Fontsize',16,'fontweight','demi')
datetick
ylabel([ 'T_p [s]'])
set(gca,'Ylim',[0 12])

ax(3) = subplot(n,1,3);
plot([AWAC.time],[AWAC.icethickness],'k+','linewidth',2), hold on
set(gca,'Fontsize',16,'fontweight','demi')
datetick
ylabel('Ice [m]')
set(gca,'Ylim',[0 5])

linkaxes(ax,'x')
set(gca,'XLim',[min([AWAC.time]) max([AWAC.time])] )
print('-dpng',[wd '_waveicestats.png'])

figure(4), clf
for ai=1:length(AWAC),
    loglog(AWAC(ai).wavespectra.freq,AWAC(ai).wavespectra.energy,'k-'), hold on
    allE(:,ai) = AWAC(ai).wavespectra.energy;
    allice(:,ai) = AWAC(ai).icehistogram.Nobs';
end
set(gca,'Fontsize',16,'fontweight','demi')
xlabel('freq [Hz]')
ylabel('Energy [m^2/Hz')
title('Scalar wave spectra')
title(wd,'interpreter','none','fontsize',10)
print('-dpng',[wd '_wavespectra.png'])

figure(5), clf
pcolor([AWAC.time],f,log10(allE)), shading flat
datetick
ylabel('freq [Hz]')
print('-dpng',[wd '_spectrogram.png'])


figure(6), clf
subplot(3,1,1)
plot(time,avgpitch,'mx')
datetick
ylabel('roll [deg]')
subplot(3,1,2)
plot(time,avgroll,'cx')
datetick
ylabel('roll [deg]')
subplot(3,1,3)
plot(time,avgheading,'yx')
datetick
ylabel('heading [deg]')
print('-dpng',[wd '_pitchrollheading.png'])

figure(7), clf
subplot(2,1,1)
plot(time,bandwidthast1,'bx',time, bandwidthast2,'gx',time,bandwidthpres,'rx')
datetick
ylabel('\nu')
subplot(2,1,2)
plot(time,tailshapeast1,'bx',time, tailshapeast2,'gx',time,tailshapepres,'rx')
datetick
ylabel('f^q')
print('-dpng',[wd '_spectralshape.png'])

figure(8), clf
pcolor([AWAC.time],icebincenters,allice), shading flat
datetick
ylabel('ice bin [m]')
colormap gray
cmap = colormap;
colormap(flipud(cmap))
cb = colorbar; cb.Label.String = 'Number of Observations';
print('-dpng',[wd '_icehistogram.png'])

toc

