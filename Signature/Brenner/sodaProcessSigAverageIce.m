function sodaProcessSigAverageIce( sodaLabel )
% SODAPROCESSSIGAVERAGEICE Process "AverageIce" data for SODA Signature500
% ADCP's
%
%   sodaProcessSigAverageIce( sodaLabel ), for 'sodaLabel' one of: 'A',
%   'B', or 'C'.
%
%   Note: this function requires a specific relative folder-path setup in
%   order to successfully find, load, and then save relevant data.  The
%   function has no direct output, but saves processed data to be used for
%   further analysis.
%
%   S.D.Brenner, 2019

%% Load data

dataDir = ['../data/SODA_',sodaLabel,'/Concatenated_Matfiles/'];
fNameSig = ['SODA_',sodaLabel,'_AverageIce_raw.mat'];
fNameConfig = ['SODA_',sodaLabel,'_Config.mat'];
% fNameMcat = ['soda_',lower(sodaLabel),'_microcats.mat'];

load('sodaConstants.mat');
load( [dataDir,fNameSig] );
load( [dataDir,fNameConfig] );
% mc = load( ['../data/',fNameMcat] );

iceData = averageiceData; clear averageiceData;

%% For SODA-C: re-create beam velocity record (inverse beam mapping)

if strcmpi( sodaLabel, 'C')
    iceData = sigBeamMappingInv(iceData,Config,'ice');
    iceData.AverageIce_VelBeam1 = iceData.AverageIce_VelBeam1(:);
    iceData.AverageIce_VelBeam2 = iceData.AverageIce_VelBeam2(:);
    iceData.AverageIce_VelBeam3 = iceData.AverageIce_VelBeam3(:);
    iceData.AverageIce_VelBeam4 = iceData.AverageIce_VelBeam4(:);
end

%% Remove "in-air" measurements

% find incides of "in-air" data (based on predefined date limits)
datelim = sodaConsts.(['SODA_',sodaLabel]).datelim;
inAirInd = findbetween( iceData.AverageIce_MatlabTimeStamp, datelim);

% Remove data
iceData = subset(iceData,inAirInd);

%% Apply magnetic offset compass correction

magOffset = sodaConsts.(['SODA_',sodaLabel]).magOffset;
decOffset = sodaConsts.(['SODA_',sodaLabel]).declinationOffset;
decOffsetTime = sodaConsts.(['SODA_',sodaLabel]).declinationOffsetTime;

iceData = sigMagCorrection2D( iceData, 'ice', magOffset);
iceData = sigDeclinationCorrection( iceData, 'ice', decOffset,decOffsetTime);


%% Apply quality metrics to velocity data

% Define quality metrics
fomMax = 1000;
maxVel = 1.5;
      
% Loop through beams
numBeams = 4;
for n = 1:numBeams
    beamFOM = iceData.(['AverageIce_FOMBeam',num2str(n)]);
    beamVel = iceData.(['AverageIce_VelBeam',num2str(n)]);
    
    highFomInds = beamFOM > fomMax;
    highVelInds = abs(beamVel) > maxVel;
    noVelInds =  (beamVel == 0);
    
    badInds = highFomInds | noVelInds | highVelInds ;
    beamVel(badInds) = NaN;
    iceData.(['AverageIce_VelBeam',num2str(n)]) = beamVel;
    
end

%% Perform ensemble averaging

nancount = @(x) sum( ~isnan(x) );

iceDataCnt = sigEnsembleAvg(iceData,'ice',[],nancount);
iceDataStd = sigEnsembleAvg(iceData,'ice',[],@nanstd);
iceData = sigEnsembleAvg(iceData,'ice');

%% Transform BEAM to ENU coordinates

iceDataStd = sigBeamMapping(iceDataStd,Config,'ice','enu');
iceData = sigBeamMapping(iceData,Config,'ice','enu');

%% Make a plot

% plot( iceData.AverageIce_MatlabTimeStamp, iceData.AverageIce_VelNorth,'.');

%% Save output

saveDir  = ['../data/SODA_',sodaLabel,'/Processed_Data/'];
fName    = ['SODA_',sodaLabel,'_Ice.mat'];
fNameStd = ['SODA_',sodaLabel,'_IceStd.mat'];
fNameCnt = ['SODA_',sodaLabel,'_IceCnt.mat'];

save( [saveDir,fName],   'iceData');
save( [saveDir,fNameStd],'iceDataStd');
save( [saveDir,fNameCnt],'iceDataCnt');

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