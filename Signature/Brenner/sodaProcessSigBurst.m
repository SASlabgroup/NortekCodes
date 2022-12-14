function sodaProcessSigBurst( sodaLabel )
% SODAPROCESSSIGBURST Process "Burst" data for SODA Signature500
% ADCP's
%
%   sodaProcessSigBurst( sodaLabel ), for 'sodaLabel' one of: 'A', 'B', or
%   'C'.
%
%   Note: this function requires a specific relative folder-path setup in
%   order to successfully find, load, and then save relevant data.  The
%   function has no direct output, but saves processed data to be used for
%   further analysis.
%
%   S.D.Brenner, 2019

%% Load data

dataDir = ['../data/SODA_',sodaLabel,'/Concatenated_Matfiles/'];
fNameSig = ['SODA_',sodaLabel,'_Burst_raw.mat'];
fNameConfig = ['SODA_',sodaLabel,'_Config.mat'];
fNameMcat = ['soda_',lower(sodaLabel),'_microcats.mat'];

load('sodaConstants.mat');
load( [dataDir,fNameSig] );
load( [dataDir,fNameConfig] );
mc = load( ['../data/',fNameMcat] );



%% Remove "in-air" measurements

% find incides of "in-air" data (based on predefined date limits)
datelim = sodaConsts.(['SODA_',sodaLabel]).datelim;
inAirInd = findbetween( burstData.Burst_MatlabTimeStamp, datelim);

% Remove data
burstData = subset(burstData,inAirInd);


%% Apply magnetic and declination offset compass corrections

magOffset = sodaConsts.(['SODA_',sodaLabel]).magOffset;
decOffset = sodaConsts.(['SODA_',sodaLabel]).declinationOffset;
decOffsetTime = sodaConsts.(['SODA_',sodaLabel]).declinationOffsetTime;

burstData = sigMagCorrection2D( burstData, 'burst', magOffset);
burstData = sigDeclinationCorrection( burstData, 'burst', decOffset,decOffsetTime);

%% Remove data with likely side-lobe interaction
% Sidelob interaction occurs at hard reflectors, which include both the
% water surface and the underside of the ice. Using altimeter data (instead
% of pressure data) to define sidelobe interaction range would account for
% both, however these data are much noiser than the pressure data.  Just
% using pressure data seems to do a good job, so I'll do that.

cutoffFactor = 0.9;
pCutoff = cutoffFactor*burstData.Burst_Pressure;
% Round to nearest bin
pCutoff = 2*round( 0.5*(pCutoff-0.5) ) + 0.5;

% Define bad Indices (ranges above cutoff values):
badInds = ( burstData.Burst_Range > pCutoff );

% Loop through beams and remove bins above the cutoff threshold
numBeams = 4;
for n = 1:numBeams    
    burstData.(['Burst_VelBeam',num2str(n)])( badInds ) = NaN;
end

%% Apply quality metrics to velocity data

% Define quality thresholds
correlationThreshold1 = 40;
correlationThreshold2 = 15;

% Loop through beams
numBeams = 4;
for n = 1:numBeams
    % Get beam data
    beamVel = burstData.( ['Burst_VelBeam',num2str(n)] );
    beamCor = burstData.( ['Burst_CorBeam',num2str(n)] );
    
    % Identify bad data:
    % Using Matlab's outlier detection method (well described in the
    % documentation), I identify statical outliers compared to the full 
    % (depth-varying) data.  "Bad data" are those data identified as
    % outliers AND with low correlation values.  Additionally, if
    % correlation value is too low, then it is identified as bad
    % reglardless of outlier status
    [numSamples,numBins]  = size(beamVel);
    beamVelVect = reshape( beamVel, 1, []);
    outlierInds = isoutlier( beamVelVect,'median');
    outlierInds= reshape( outlierInds, size(beamVel) );
    lowCorInds1 = beamCor < correlationThreshold1;
    lowCorInds2 = beamCor < correlationThreshold2;
    badInds = (outlierInds & lowCorInds1) | lowCorInds2;
    
    % Remove bad data (replace with NaN for now)
    beamVel(badInds) = NaN;
    
    % Update data structure
    burstData.(['Burst_VelBeam',num2str(n)]) = beamVel;
end

%% Apply quality metrics to altimeter data
% At this point it is unclear what the threshold values should be in order
% to retain and maximize data for wave processing while still recjecting
% bad data.  For the time being, skip this step.


 
% % Define quality thresholds (determined semi-empirically) 
% minAltimeterCutoff = 20;
% altimQualityThresholdAST = 8000;
% altimQualityThresholdLE = 30000;
% 
% 
% % Apply them:
% averageData.Burst_AltimeterAST(...
%     averageData.Burst_AltimeterAST < minAltimeterCutoff | ...
%     averageData.Burst_AltimeterQualityAST < altimQualityThresholdAST ) = NaN;
% averageData.Burst_AltimeterLE(...
%     averageData.Burst_AltimeterLE < minAltimeterCutoff | ...
%     averageData.Burst_AltimeterQualityLE < altimQualityThresholdLE ) = NaN;


%% Calculate corrected altimeter distances
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
burstData = sodaSoundSpeedCorrection(burstData,mcData,'burst',lon,lat,datelim);

% Calculate altimeter correction factors
soundSpeedCor = burstData.Burst_SpeedOfSound_Corrected./burstData.Burst_SpeedOfSound;
tiltCor = cosd(burstData.Burst_Pitch) .* cosd(burstData.Burst_Roll);
distCor = soundSpeedCor.*tiltCor;

burstData.Burst_AltimeterAST_Corrected = distCor.*burstData.Burst_AltimeterAST;
burstData.Burst_AltimeterLE_Corrected  = distCor.*burstData.Burst_AltimeterLE;

%% Transform BEAM to ENU coordinates

burstData = sigBeamMapping(burstData,Config,'burst','enu');


%% Save output

saveDir = ['../data/SODA_',sodaLabel,'/Processed_Data/'];
fName = ['SODA_',sodaLabel,'_Burst.mat'];

save([saveDir,fName],'burstData','-v7.3');


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