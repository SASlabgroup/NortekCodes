function sodaProcessSigAverage( sodaLabel )
% SODAPROCESSSIGAVERAGE Process "Average" data for SODA Signature500 ADCP's
%
%   sodaProcessSigAverage( sodaLabel ), for 'sodaLabel' one of: 'A', 'B',
%   or 'C'.
%
%   Note: this function requires a specific relative folder-path setup in
%   order to successfully find, load, and then save relevant data.  The
%   function has no direct output, but saves processed data to be used for
%   further analysis.
%
%   S.D.Brenner, 2019


%% Load data

dataDir = ['../data/SODA_',sodaLabel,'/Concatenated_Matfiles/'];
fNameSig = ['SODA_',sodaLabel,'_Average_raw.mat'];
fNameConfig = ['SODA_',sodaLabel,'_Config.mat'];
fNameMcat = ['soda_',lower(sodaLabel),'_microcats.mat'];

load('sodaConstants.mat');
load( [dataDir,fNameSig] );
load( [dataDir,fNameConfig] );
mc = load( ['../data/',fNameMcat] );


%% For SODA-C: re-create beam velocity record (inverse beam mapping)

if strcmpi( sodaLabel, 'C')
    averageData = sigBeamMappingInv(averageData,Config,'avg');
end

%% Remove "in-air" measurements

% find incides of "in-air" data (based on predefined date limits)
datelim = sodaConsts.(['SODA_',sodaLabel]).datelim;
inAirInd = findbetween( averageData.Average_MatlabTimeStamp, datelim);

% Remove data
averageData = subset(averageData,inAirInd);


%% Apply magnetic and declination offset compass corrections

magOffset = sodaConsts.(['SODA_',sodaLabel]).magOffset;
decOffset = sodaConsts.(['SODA_',sodaLabel]).declinationOffset;
decOffsetTime = sodaConsts.(['SODA_',sodaLabel]).declinationOffsetTime;

averageData = sigMagCorrection2D( averageData, 'avg', magOffset);
averageData = sigDeclinationCorrection( averageData, 'avg', decOffset,decOffsetTime);


%% Remove data with likely side-lobe interaction
% Sidelob interaction occurs at hard reflectors, which include both the
% water surface and the underside of the ice. Using altimeter data (instead
% of pressure data) to define sidelobe interaction range would account for
% both, however these data are much noiser than the pressure data.  Just
% using pressure data seems to do a good job, so I'll do that.

cutoffFactor = 0.9;
pCutoff = cutoffFactor*averageData.Average_Pressure;
% Round to nearest bin
pCutoff = 2*round( 0.5*(pCutoff-0.5) ) + 0.5;

% Define bad Indices (ranges above cutoff values):
badInds = ( averageData.Average_Range > pCutoff );

% Loop through beams and remove bins above the cutoff threshold
numBeams = 4;
for n = 1:numBeams    
    averageData.(['Average_VelBeam',num2str(n)])( badInds ) = NaN;
end

%% Apply quality metrics to velocity data

% Define quality thresholds
correlationThreshold1 = 40;
correlationThreshold2 = 15;

% Loop through beams
numBeams = 4;
for n = 1:numBeams
    % Get beam data
    beamVel = averageData.( ['Average_VelBeam',num2str(n)] );
    beamCor = averageData.( ['Average_CorBeam',num2str(n)] );
    
    % Identify bad data:
    % Using Matlab's outlier detection method (well described in the
    % documentation), I identify statical outliers compared to the full 
    % (depth-varying) data.  "Bad data" are those data identified as
    % outliers AND with low correlation values.  Additionally, if
    % correlation value is too low, then it is identified as bad
    % reglardless of outlier status
    [numSamples,numBins]  = size(beamVel);
    beamVelVect = reshape( beamVel, 1, []);
    outlierInds = isoutlier( beamVelVect,'median'); % 'ThresholdFactor',4);
    outlierInds= reshape( outlierInds, size(beamVel) );
    lowCorInds1 = beamCor < correlationThreshold1;
    lowCorInds2 = beamCor < correlationThreshold2;
    badInds = (outlierInds & lowCorInds1) | lowCorInds2;
    
    % Remove bad data (replace with NaN for now)
    beamVel(badInds) = NaN;
    
    % Update data structure
    averageData.(['Average_VelBeam',num2str(n)]) = beamVel;
    
    % Plot histogram results of beam quality corrections for one of the
    % four beams:
    if n==4
        figure(1); clf;
        binEdges = linspace(-6,6,161);
        grey = [0.25,0.25,0.5] ;%0.25*[1,1,1];
        subplot(2,1,1);
        histogram( beamVelVect , binEdges, 'edgecolor',grey);
        set(gca,'yscale','log','ylim',10.^[0,8],'ytick',10.^[0:2:8],'fontsize',10);
        xlabel('Beam 4 velocity [m/s]');
        ylabel('Count');
        text(-6.2,10^7.5,'(a)','fontsize',12);
        
        subplot(2,1,2)
        hold on;
        histogram( beamVelVect(~badInds) , binEdges, 'edgecolor',grey);
        histogram( beamVelVect(~outlierInds) , binEdges, 'edgecolor',grey);
        set(gca,'yscale','log','ylim',10.^[0,8],'ytick',10.^[0:2:8],'fontsize',10);
        xlabel('Beam 4 velocity [m/s]');
        ylabel('Count');
        text(-6.2,10^7.5,'(b)','fontsize',12);
        
        set(gcf,'units','inches','position',[1,1,5.5,5])
    end
    
end

%% Apply quality metrics to altimeter data

% Define quality thresholds (determined semi-empirically) 
minAltimeterCutoff = 20;
altimQualityThresholdAST = 8000;
altimQualityThresholdLE = 30000;


% Apply them:
averageData.Average_AltimeterIceAST(...
    averageData.Average_AltimeterIceAST < minAltimeterCutoff | ...
    averageData.Average_AltimeterQualityAST < altimQualityThresholdAST ) = NaN;
averageData.Average_AltimeterIceLE(...
    averageData.Average_AltimeterIceLE < minAltimeterCutoff | ...
    averageData.Average_AltimeterQualityLE < altimQualityThresholdLE ) = NaN;

%% Calculate altimeter distances
% Appy corrections to account for both sound speed and tilt

% Extract microcata data 
% (the microcat closest to the Stablemoor is the first one in the struct)
mcFlds = fields(mc); 
mcName = mcFlds{1};
mcData = mc.(mcName);
clear mc mcFlds mcName;

% Extract constants data:
lon = sodaConsts.(['SODA_',sodaLabel]).lon;
lat = sodaConsts.(['SODA_',sodaLabel]).lat;

% Get corrected sound speed
averageData = sodaSoundSpeedCorrection(averageData,mcData,'avg',lon,lat,datelim);

% Calculate altimeter correction factors
soundSpeedCor = averageData.Average_SpeedOfSound_Corrected./averageData.Average_SpeedOfSound;
tiltCor = cosd(averageData.Average_Pitch) .* cosd(averageData.Average_Roll);
distCor = soundSpeedCor.*tiltCor;

averageData.Average_AltimeterIceAST_Corrected = distCor.*averageData.Average_AltimeterIceAST;
averageData.Average_AltimeterIceLE_Corrected  = distCor.*averageData.Average_AltimeterIceLE;

%% Perform ensemble averaging

nancount = @(x) sum( ~isnan(x) );

averageDataCnt = sigEnsembleAvg(averageData,'avg',[],nancount);
averageDataStd = sigEnsembleAvg(averageData,'avg',[],@nanstd);
averageData = sigEnsembleAvg(averageData,'avg');

%% Transform BEAM to ENU coordinates

averageData = sigBeamMapping(averageData,Config,'avg','enu');

% %% Make a plot
% 
% figure(10); clf;
% hold on;
% pcolor( averageData.Average_MatlabTimeStamp,...
%         averageData.Average_Range,...
%         averageData.Average_VelEast.');
% shading interp;
% 
% set(gca,'clim',0.25*[-1,1]);
% colormap( cbrewer('div','RdBu') );


%% Save output

saveDir = ['../data/SODA_',sodaLabel,'/Processed_Data/'];
fName = ['SODA_',sodaLabel,'_Average.mat'];
fNameStd = ['SODA_',sodaLabel,'_AverageStd.mat'];
fNameCnt = ['SODA_',sodaLabel,'_AverageCnt.mat'];

save([saveDir,fName],'averageData');
save([saveDir,fNameStd],'averageDataStd');
save([saveDir,fNameCnt],'averageDataCnt');

msgbox('All files saved') 

end

%% EMBEDDED FUNCTIONS %% ==================================================

 function [ Data ] = subset( Data, ind )
    % list all fields in data structure
    flds = fields(Data);
       
    % Loop and subset
    for n = 1:length(flds)
        fldName = flds{n};
        [numSamples,~] = size( Data.(fldName) );
        if numSamples > 1
            Data.(fldName) = Data.(fldName)(ind,:);
        end
    end
    
end