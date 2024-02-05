% script to process raw Nortek Signature data, both burst and average cycles
% assuming an export from Signature software (not MIDAS)
% with coordinate rotations to ENU already included in the export
%
% Note that Burst and Average data usually are on different time steps
%
% J. Thomson, Sep 2020
%       modibed from original AWAC processing (c. 2012)
%       using SWIFT codes for wave processing
%       and parts of Sam Brenner signature processing suite
%

clear all
tic


%% Set constants and QC criteria
%site = 'S1A1'; lat = 70.48695; lon = -162.28278; doff=0.75; % CODA S1-A1, a sea spider 0.75 m off the seabed
%site = 'S2A1'; lat = 70.77422; lon = -149.47707; doff=0.75; % CODA S1-A1, a sea spider 0.75 m off the seabed
%site = 'S3A1'; lat = 70.39953; lon = -145.85623; doff=0.75; % CODA S1-A1, a sea spider 0.75 m off the seabed
site = 'NORSE'; lat = 70.831946; lon = -6.399105; doff=385; % CODA S1-A1, a sea spider 0.75 m off the seabed

despike = true;
mindepth = 4;
mincorr = 40; % broadband Doppler correlation
maxwaveperiod = 20; % max wave period allowed during final screening
minwaveperiod = 2; % min wave period allowed during final screening
minicethickness = 0.6; % min ice thickness considered valid during final screening (10 cm is barometric noise)
maxicethickness = 20; % min ice thickness considered valid during final screening (10 cm is barometric noise)
minwaveheight = 0.1; % smallest wave height observable
minwaveheightfordir = 1; % smallest wave height for directional estimates
maxtailshapeexponent = -2;  % max value for f^q in the tail (f> 0.3 Hz)
icebincenters = [0:.5:10]; % meters

%% Load data

dataDir = ['./'];

fNameSigBase = ['S100793A014*.mat'];

flist = dir([dataDir fNameSigBase ]);

%% Loop through files

bcounter = 1; % burst counter (across multiple files)
acounter = 1; % Avg counter (across multiple files)

for fi=1:length(flist)
    
    disp(['file ' num2str(fi) ' of ' num2str(length(flist))])
    load([dataDir flist(fi).name]) % this ble will have N bursts and M averages
    
    %% deal with MIDAS export by adding fieldnames compliant with Sig export
    if ~isfield(Data,'Burst_Time'),
       disp('File from MIDAS')
       Data.Burst_Time = Data.Burst_MatlabTimeStamp;
       Data.Burst_AltimeterDistanceAST = Data.Burst_AltimeterAST;
       Config.Burst_BlankingDistance = Config.burst_blankingDistance;
       Config.Burst_CellSize = Config.burst_cellSize;
       Config.Burst_NCells = Config.burst_nCells;
       Config.Burst_SamplingRate = Config.burst_sampleRate;
       Data.Burst_Temperature = Data.Burst_WaterTemperature;
       Config.Average_BlankingDistance = Config.avg_blankingDistance;
       Config.Average_CellSize = Config.avg_cellSize;
       Data.Average_NCells = Config.avg_nCells;
       Data.Average_Time = Data.Average_MatlabTimeStamp;
       Data.Average_Temperature = Data.Average_WaterTemperature;
       Config.Average_IceTrack = Config.avg_ice;
       Data = sigBeamMapping(Data,Config,'burst','enu');  % use Sam Brenner's code for ENU of the burst data (or assign NaN)
       %Data.Burst_VelEast = NaN(size(Data.Burst_VelBeam1));  % ENU are not part of MIDAS export!
       %Data.Burst_VelNorth = NaN(size(Data.Burst_VelBeam1));  % ENU are not part of MIDAS export!
    end
    
    %% Loop through bursts (of that file) making wave and ice products
    
    % Find indices for each burst
    firstindex = find(Data.Burst_EnsembleCount==1);
    spacing = mean(diff(firstindex));
    if isnan(spacing), spacing=length(Data.Burst_EnsembleCount)-firstindex; end % exception for only one burst in file
    nb = length(firstindex);
    
    % Loop through bursts and process
    for b = 1:nb
        
        %% select burst data
        burstInd = firstindex(b) + [1:spacing] - 1;
        burstInd( burstInd > length(Data.Burst_Time) ) = [];
        ast = double( Data.Burst_AltimeterDistanceAST(burstInd) );
        range = Config.Burst_BlankingDistance + Config.Burst_CellSize*[1:Config.Burst_NCells];
        velindex = find(range < min(ast),1,'last'); % bin of burst velocity profile to use in wave processing
        if isempty(velindex), velindex = 1; end
        if isfield(Data,'Burst_VelEast')
            u = double( Data.Burst_VelEast(burstInd, velindex ) );
            v = double( Data.Burst_VelNorth(burstInd, velindex ) );
        else
            u = NaN* double( Data.Burst_VelBeam1(burstInd, velindex ) );
            v = NaN* double( Data.Burst_VelBeam2(burstInd, velindex ) );
            disp('** No burst ENU transform, wave directios will be spurious')
        end
        bursttime = Data.Burst_Time(burstInd);
        rate = double( Config.Burst_SamplingRate );
        
        %% meta data
        sigBurst(bcounter).time = bursttime(1); 
        sigBurst(bcounter).lat = lat; 
        sigBurst(bcounter).lon = lon;
        sigBurst(bcounter).watertemp = mean( Data.Burst_Temperature(burstInd) );
        sigBurst(bcounter).depth = mean( Data.Burst_Pressure(burstInd) );
        waterdepth = sigBurst(bcounter).depth + doff;
        
        
        %% QC for out of water
        if sigBurst(bcounter).depth < mindepth,
            badburst(bcounter)=true;
        else
            badburst(bcounter)=false;
        end
        
        %% despike raw before wave and ice processing
        if despike
            ast = filloutliers(ast,'linear');
            u = filloutliers(u,'linear');
            v = filloutliers(v,'linear');
        end
        
        %% wave processing
        [ Hs, Tp, Dp, E, f, a1, b1, a2, b2, check] = UVZwaves(u, v, ast, rate); % assumes already east and north
        if ~isnan(E), 
            bt = polyfit(log(f(f>0.3)),log(E(f>0.3)),1);
            tailshape = bt(1);
        else
            tailshape = inf;
        end
        if tailshape <= maxtailshapeexponent & Hs > minwaveheight & Tp > minwaveperiod & Tp < maxwaveperiod,
            sigBurst(bcounter).sigwaveheight = Hs;
            sigBurst(bcounter).peakwaveperiod = Tp;
            sigBurst(bcounter).peakwavedirT = Dp;
            sigBurst(bcounter).wavespectra.energy = E';
            sigBurst(bcounter).wavespectra.freq = f';
            sigBurst(bcounter).wavespectra.a1 = a1';
            sigBurst(bcounter).wavespectra.b1 = b1';
            sigBurst(bcounter).wavespectra.a2 = a2';
            sigBurst(bcounter).wavespectra.b2 = b2';
            sigBurst(bcounter).wavespectra.check = check';
        else
            sigBurst(bcounter).sigwaveheight = NaN;
            sigBurst(bcounter).peakwaveperiod = NaN;
            sigBurst(bcounter).peakwavedirT = NaN;
            sigBurst(bcounter).wavespectra.energy = NaN(size(E'));
            sigBurst(bcounter).wavespectra.freq = NaN(size(f'));
            sigBurst(bcounter).wavespectra.a1 = NaN(size(a1'));
            sigBurst(bcounter).wavespectra.b1 = NaN(size(b1'));
            sigBurst(bcounter).wavespectra.a2 = NaN(size(a2'));
            sigBurst(bcounter).wavespectra.b2 = NaN(size(b2'));
            sigBurst(bcounter).wavespectra.check = NaN(size(check'));
        end
        
        %% ice stats
        draft = sigBurst(bcounter).depth - mean(ast);
        if draft > minicethickness & draft < maxicethickness,
            sigBurst(bcounter).icethickness = draft;
            sigBurst(bcounter).icehistogram.Nobs = hist( sigBurst(bcounter).depth - ast , icebincenters);
            sigBurst(bcounter).icehistogram.bincenters = icebincenters;
        else
            sigBurst(bcounter).icethickness = NaN;
            sigBurst(bcounter).icehistogram.Nobs = NaN(size(icebincenters));
            sigBurst(bcounter).icehistogram.bincenters = NaN(size(icebincenters));
        end
        
        %% increment overall burst counter (keeps counting across files)
        bcounter = bcounter+1;
        
    end % close burst loop for that file
    
    
    %% loop through Avgd ensembles (of that file) making profiles and ice

    %% QC data to close to surface (or out of water) and with low correlations
    
    z = double( doff + Config.Average_BlankingDistance + Config.Average_CellSize .* [1: Data.Average_NCells(1)] );
    beamangle = 25;
    maxrange = double(Data.Average_AltimeterPressure) * cosd(beamangle); % times series of max range
    trimcells = ones(size(maxrange))*(z-doff) > 0.95*maxrange*ones(size(z))  | Data.Average_CorBeam1 < mincorr | ...
        Data.Average_CorBeam2 < mincorr | Data.Average_CorBeam3 < mincorr | Data.Average_CorBeam4 < mincorr;
    Data.Average_VelEast(trimcells) = NaN;
    Data.Average_VelNorth(trimcells) = NaN;
    Data.Average_VelUp1(trimcells) = NaN;
    Data.Average_VelUp2(trimcells) = NaN;
    Data.Average_AmpBeam1(trimcells) = NaN;
    Data.Average_AmpBeam2(trimcells) = NaN;
    Data.Average_AmpBeam3(trimcells) = NaN;
    Data.Average_AmpBeam4(trimcells) = NaN;
    
    % Find indices for each ensemble
    firstindex = find(Data.Average_EnsembleCount==1);
    spacing = mean(diff(firstindex));
    na = length(firstindex);
    
    for a=1:na
        
        AvgInd = firstindex(a) + [1:spacing] - 1;
        AvgInd( AvgInd > length(Data.Average_Time) ) = [];
        
        if ~isnan(AvgInd)

        %% meta data
        sigAverage(acounter).time = Data.Average_Time(AvgInd(1));   
        sigAverage(acounter).lat = lat;  
        sigAverage(acounter).lon = lon;
        sigAverage(acounter).watertemp = median(Data.Average_Temperature(AvgInd));
        sigAverage(acounter).depth = median(Data.Average_Pressure(AvgInd));
        
                %% QC for out of water
        if sigAverage(acounter).depth < mindepth,
            badavg(acounter)=true;
        else
            badavg(acounter)=false;
        end
        
        %% profiles
        sigAverage(acounter).z = z;
        if size(Data.Average_VelEast,1) == 1 % check if ENU needs to be reshaped (as in NORSE, but not CODA)
            checklength = length(Data.Average_Time) * length(z) - size(Data.Average_VelEast,2);
            if checklength ~=0
                Data.Average_VelEast(end + [1:checklength]) = NaN;
                Data.Average_VelNorth(end + [1:checklength]) = NaN;
                Data.Average_VelUp1(end + [1:checklength]) = NaN;
                Data.Average_VelUp2(end + [1:checklength]) = NaN;
            end
            Data.Average_VelEast = reshape(Data.Average_VelEast, length(Data.Average_Time), length(z));
            Data.Average_VelNorth = reshape(Data.Average_VelNorth, length(Data.Average_Time), length(z));
            Data.Average_VelUp1 = reshape(Data.Average_VelUp1, length(Data.Average_Time), length(z));
            Data.Average_VelUp2 = reshape(Data.Average_VelUp2, length(Data.Average_Time), length(z));
        end
        sigAverage(acounter).east = nanmean( Data.Average_VelEast(AvgInd,:) );
        sigAverage(acounter).north = nanmean( Data.Average_VelNorth(AvgInd,:) );
        sigAverage(acounter).up = ( nanmean(Data.Average_VelUp1(AvgInd,:)) + nanmean(Data.Average_VelUp2(a,:)) ) ./ 2;
        sigAverage(acounter).backscatter1 = nanmean( Data.Average_AmpBeam1(AvgInd,:) );
        sigAverage(acounter).backscatter2 = nanmean( Data.Average_AmpBeam2(AvgInd,:) );
        sigAverage(acounter).backscatter3 = nanmean( Data.Average_AmpBeam3(AvgInd,:) );
        sigAverage(acounter).backscatter4 = nanmean( Data.Average_AmpBeam4(AvgInd,:) );
        
        %% ice products (if tracking enabled)
        sigAverage(acounter).icethickness = NaN;
        
        if Config.Average_IceTrack == 'True' & Config.Average_Altimeter == 'True' & ...
                Config.Average_AltimeterEnd > sigAverage(acounter).depth
            distance = mean( [Data.Average_AltimeterDistanceAST(AvgInd)' ]);
            %distance = mean( [Data.AverageIce_DistanceBeam1(AvgInd)' Data.AverageIce_DistanceBeam2(AvgInd)' ...
                %Data.AverageIce_DistanceBeam3(AvgInd)' Data.AverageIce_DistanceBeam4(AvgInd)' ]); % if no altimeter
            draft = sigAverage(acounter).depth - distance;
            if draft>minicethickness & draft<maxicethickness & isfield(Data,'AverageIce_VelSpeed'),
                sigAverage(acounter).icethickness = draft;
                sigAverage(acounter).icespeed = mean( Data.AverageIce_VelSpeed(AvgInd) );
                sigAverage(acounter).icedir = mean( Data.AverageIce_VelDirection(AvgInd) );
            else
                sigAverage(acounter).icethickness = NaN;
                sigAverage(acounter).icespeed = NaN;
                sigAverage(acounter).icedir = NaN;
            end
        else
        end
        
        %% increment the overall average counter
        acounter = acounter+1;
        
        else
        end
        
    end % close Avgd loop
    
    
    
    
end % close file loop

%% apply burst QC
sigBurst(badburst) = [];
sigAverage(badavg) = [];

%% trim velocity profiles for ice

for sa = 1:length(sigAverage),
    icecells = find( z > 0.6*(sigAverage(sa).depth - sigAverage(sa).icethickness) );
    sigAverage(sa).east( icecells  ) = NaN;
    sigAverage(sa).north( icecells  ) = NaN;
    sigAverage(sa).up( icecells  ) = NaN;
    sigAverage(sa).backscatter1( icecells  ) = NaN;
    sigAverage(sa).backscatter2( icecells  ) = NaN;
    sigAverage(sa).backscatter3( icecells  ) = NaN;
    sigAverage(sa).backscatter4( icecells  ) = NaN;
end

%% remove peak direction if spread too large or waves too small (SNR problem)


for si=1:length(sigBurst), 
    [m fpindex] = min( abs( sigBurst(si).peakwaveperiod - sigBurst(si).wavespectra.freq ) );
    a1 = sigBurst(si).wavespectra.a1(fpindex);
    b1 = sigBurst(si).wavespectra.b1(fpindex);
    a2 = sigBurst(si).wavespectra.a2(fpindex);
    b2 = sigBurst(si).wavespectra.b2(fpindex);
    dir1(si) = atan2(b1,a1) ;  % [rad], 4 quadrant
    dir2(si) = atan2(b2,a2)/2 ; % [rad], only 2 quadrant
    spread1(si) = sqrt( 2 * ( 1 - sqrt(a1.^2 + b1.^2) ) );
    spread2(si) = sqrt( abs( 0.5 - 0.5 .* ( a2.*cos(2.*dir2(si)) + b2.*cos(2.*dir2(si)) )  ));
    if 180/3.14*spread1(si) > 85 | sigBurst(si).sigwaveheight < minwaveheightfordir | isnan( sigBurst(si).sigwaveheight )
        sigBurst(si).peakwavedirT=NaN; 
    end, 
end

%% sort results
[time sortinds] = sort([sigBurst.time]);
sigBurst = sigBurst(sortinds);

[time sortinds] = sort([sigAverage.time]);
sigAverage = sigAverage(sortinds);

%% save output

save( [site '_Signature_all_processed.mat'], 'sigBurst','sigAverage','site');

toc

%% plots

panels = 6;
chunks = find(diff([sigAverage.time])>1/24); % indices of contiguous data chunks (for pcolors)
nc = length(chunks) + 1; % number of chunks 
chunks = [0 chunks length([sigAverage.time])];

figure(1), clf

ax(1) = subplot(panels,1,1); % wave heights
plot([sigBurst.time],[sigBurst.sigwaveheight],'.')
datetick
set(gca,'YLim',[0 1.1*max([sigBurst.sigwaveheight])])
cb = colorbar; set(cb,'Visible','off')
ylabel('H_s [m]')
title(site)

ax(2) = subplot(panels,1,2); % ice draft
plot([sigBurst.time],[sigBurst.icethickness],'bo',[sigAverage.time],[sigAverage.icethickness],'g.')
datetick
cb = colorbar; set(cb,'Visible','off')
ylabel('Ice [m]')


ax(3) = subplot(panels,1,3); % temperature
plot([sigBurst.time],[sigBurst.watertemp],'bo',[sigAverage.time],[sigAverage.watertemp],'g.')
datetick
set(gca,'YLim',[-2 10])
cb = colorbar; set(cb,'Visible','off')
ylabel('T [C]')

for ci = 1:nc
    thisc = [ (chunks(ci)+1) : chunks(ci+1) ];

if ci==1, ax(4) = subplot(panels,1,4); else axes(ax(4)), end % east profiles
pcolor([sigAverage(thisc).time], sigAverage(thisc(1)).z,reshape([sigAverage(thisc).east],... 
    length(sigAverage(thisc(1)).z),length([sigAverage(thisc).time])))
hold on
plot([sigAverage(thisc).time],[sigAverage(thisc).depth],'k-')
set(gca,'YLim',[0 1.1*max([sigAverage.depth])+1])
shading flat
datetick
ylabel(['z [m]'])
caxis([-.5 .5])
colormap(gca,nawhimar);
cb = colorbar; cb.Label.String = 'East [m/s]';

if ci==1, ax(5) = subplot(panels,1,5); else axes(ax(5)), end % north profiles
pcolor([sigAverage(thisc).time], sigAverage(thisc(1)).z,reshape([sigAverage(thisc).north],...
    length(sigAverage(thisc(1)).z),length([sigAverage(thisc).time])))
hold on
plot([sigAverage(thisc).time],[sigAverage(thisc).depth],'k-')
set(gca,'YLim',[0 1.1*max([sigAverage.depth])+1])
shading flat
datetick
ylabel(['z [m]'])
caxis([-.5 .5])
colormap(gca,nawhimar);
cb = colorbar; cb.Label.String = 'North [m/s]';

if ci==1, ax(6) = subplot(panels,1,6); else axes(ax(6)), end % backscatter profiles
pcolor([sigAverage(thisc).time], sigAverage(thisc(1)).z,reshape([sigAverage(thisc).backscatter1],...
    length(sigAverage(thisc(1)).z),length([sigAverage(thisc).time])))
hold on
plot([sigAverage(thisc).time],[sigAverage(thisc).depth],'k-')
set(gca,'YLim',[0 1.1*max([sigAverage.depth])+1])
shading flat
datetick
ylabel(['z [m]'])
%caxis([-.5 .5])
cb = colorbar; cb.Label.String = 'Echo [dB]';

end

linkaxes(ax,'x')

print('-dpng',[site '_Signature_all_processed.png'])
savefig([site '_Signature_all_processed.fig'])


