% script to process Nortek Signature average ouput (avgd_avgd)
% assuming an export from Signature software (not MIDAS)
% with coordinate rotations to ENU already included in the export
%
% Note that Burst and Average data usually are on different time steps
%   this code does only the avgd_avgd data (ensembles already made by
%   Nortek)
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
site = 'S3A1'; lat = 70.39953; lon = -145.85623; doff=0.75; % CODA S1-A1, a sea spider 0.75 m off the seabed


despike = true;
mindepth = 4;
mincorr = 40; % broadband Doppler correlation
minicethickness = 0.1; % min ice thickness considered valid during final screening (10 cm is barometric noise)
maxicethickness = 10; % min ice thickness considered valid during final screening (10 cm is barometric noise)


%% Load data

dataDir = ['./'];

fNameSigBase = ['*avgd_avgd*.mat'];

flist = dir([dataDir fNameSigBase ]);

%% Loop through files

acounter = 1; % Avg counter (across multiple files)

for fi=1:length(flist)
    
    disp(['file ' num2str(fi) ' of ' num2str(length(flist))])
    load([dataDir flist(fi).name]) % this ble will have N bursts and M averages
    
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
    
    na = length(Data.Average_Time);
    
    for a=1:na
        
        AvgInd = a; 

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
        sigAverage(acounter).east = Data.Average_VelEast(AvgInd,:) ;
        sigAverage(acounter).north = Data.Average_VelNorth(AvgInd,:) ;
        sigAverage(acounter).up = ( Data.Average_VelUp1(AvgInd,:)) + mean(Data.Average_VelUp2(a,:) ) ./ 2;
        sigAverage(acounter).backscatter1 = Data.Average_AmpBeam1(AvgInd,:) ;
        sigAverage(acounter).backscatter2 =  Data.Average_AmpBeam2(AvgInd,:) ;
        sigAverage(acounter).backscatter3 =  Data.Average_AmpBeam3(AvgInd,:) ;
        sigAverage(acounter).backscatter4 =  Data.Average_AmpBeam4(AvgInd,:) ;
        
        %% ice products (if tracking enabled)
        sigAverage(acounter).icethickness = NaN;
        
        if Config.Average_Altimeter == 'True' & Config.Average_AltimeterEnd > sigAverage(acounter).depth,
            distance = mean( [Data.Average_AltimeterDistanceAST(AvgInd)' ]); 
            draft = sigAverage(acounter).depth - distance;
            if draft>minicethickness & draft<maxicethickness,
                sigAverage(acounter).icethickness = draft;
            else
                sigAverage(acounter).icethickness = NaN;
            end
        else
            sigAverage(acounter).icethickness = NaN;
        end
        
        %% increment the overall average counter
        acounter = acounter+1;
        
    end % close Avgd loop
    
    
    
    
end % close file loop

%% apply burst QC
sigAverage(badavg) = [];

%% trim velocity profiles for ice

for sa = 1:length(sigAverage),
    icecells = find( z > 0.8*(sigAverage(sa).depth - sigAverage(sa).icethickness) );
    sigAverage(sa).east(icecells  ) = NaN;
    sigAverage(sa).north(icecells  ) = NaN;
    sigAverage(sa).up(icecells  ) = NaN;
end

%% sort results

[time sortinds] = sort([sigAverage.time]);
sigAverage = sigAverage(sortinds);

%% save output

save( [site '_Signature_avgd_processed.mat'], 'sig*');

toc

%% plots

panels = 6;

figure(1), clf

%ax(1) = subplot(panels,1,1); % wave heights
% plot([sigBurst.time],[sigBurst.sigwaveheight],'x')
% datetick
% cb = colorbar; set(cb,'Visible','off')
% ylabel('H_s [m]')
title([site])

ax(2) = subplot(panels,1,2); % ice draft
plot([sigAverage.time],[sigAverage.icethickness],'g.')
datetick
cb = colorbar; set(cb,'Visible','off')
ylabel('Ice [m]')
title([site])

ax(3) = subplot(panels,1,3); % temperature
plot([sigAverage.time],[sigAverage.watertemp],'g.')
datetick
cb = colorbar; set(cb,'Visible','off')
ylabel('T [C]')


ax(4) = subplot(panels,1,4); % east profiles
pcolor([sigAverage.time], sigAverage(1).z,reshape([sigAverage.east],length(sigAverage(1).z),length([sigAverage.time])))
hold on
plot([sigAverage.time],[sigAverage.depth],'k-')
set(gca,'YLim',[0 1.1*max([sigAverage.depth])+1])
shading flat
datetick
ylabel(['z [m]'])
caxis([-.5 .5])
colormap(gca,nawhimar);
cb = colorbar; cb.Label.String = 'East [m/s]';

ax(5) = subplot(panels,1,5); % east profiles
pcolor([sigAverage.time], sigAverage(1).z,reshape([sigAverage.north],length(sigAverage(1).z),length([sigAverage.time])))
hold on
plot([sigAverage.time],[sigAverage.depth],'k-')
set(gca,'YLim',[0 1.1*max([sigAverage.depth])+1])
shading flat
datetick
ylabel(['z [m]'])
caxis([-.5 .5])
colormap(gca,nawhimar);
cb = colorbar; cb.Label.String = 'North [m/s]';

ax(6) = subplot(panels,1,6); % east profiles
pcolor([sigAverage.time], sigAverage(1).z,reshape([sigAverage.backscatter1],length(sigAverage(1).z),length([sigAverage.time])))
hold on
plot([sigAverage.time],[sigAverage.depth],'k-')
set(gca,'YLim',[0 1.1*max([sigAverage.depth])])
shading flat
datetick
ylabel(['z [m]'])
%caxis([-.5 .5])
cb = colorbar; cb.Label.String = 'Echo [dB]';

linkaxes(ax,'x')

print('-dpng',[site '_Signature_avgd_processed.png'])


