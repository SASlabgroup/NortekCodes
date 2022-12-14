function sodaIceProduct( sodaLabel )
% SODAICEPRODUCT Process and create ice data product for SODA
% Signature500 ADCP's
%
%   sodaIceProduct( sodaLabel ), for 'sodaLabel' one of: 'A', 'B', or 'C'.
%
%   Note: this function requires a specific relative folder-path setup in
%   order to successfully find, load, and then save relevant data.  The
%   function has no direct output, but saves processed data to be used for
%   further analysis.
%
%   S.D.Brenner, 2019


%% Load data

dataDir = ['../data/SODA_',sodaLabel,'/Processed_Data/'];
% Ocn (ice draft)
fName = ['SODA_',sodaLabel,'_Average.mat'];
fNameStd = ['SODA_',sodaLabel,'_AverageStd.mat'];
fNameCnt = ['SODA_',sodaLabel,'_AverageCnt.mat'];
% Ice (ice velocity)
fNameIce = ['SODA_',sodaLabel,'_Ice.mat'];
fNameIceCnt = ['SODA_',sodaLabel,'_IceCnt.mat'];

load( [dataDir,fName] );
load( [dataDir,fNameStd] );
load( [dataDir,fNameCnt] );
load( [dataDir,fNameIce] );
load( [dataDir,fNameIceCnt] );

dataDir = ['../data/SODA_',sodaLabel,'/Concatenated_Matfiles/'];
fNameConfig = ['SODA_',sodaLabel,'_Config.mat'];
load( [dataDir,fNameConfig] );

load('sodaConstants.mat');

% Load re-analysis data ( ice draft calculations depend on atmospheric
% pressure offsets )
dbar2pa = 1e4;
ERA5 = load('../../Remote_sensing_and_reanalysis/ERA5/mooringsERA5.mat');
atmoMattime = ERA5.(['SODA_',sodaLabel]).mattime;
atmoPressure = ERA5.(['SODA_',sodaLabel]).sp/dbar2pa;


%% Calculate water depth
% Apply atmospheric pressure offset, and make guess to account for density
% variation in hydrostatic pressure integration

% Calculate corrected water pressure
mattime = averageData.Average_MatlabTimeStamp;
pAbsolute = averageData.Average_AltimeterPressure + Config.pressureOffset;
pAtm = interp1( atmoMattime, atmoPressure, mattime);
pWater = pAbsolute - pAtm;

% Estimate bulk density
lon = sodaConsts.(['SODA_',sodaLabel]).lon;
lat = sodaConsts.(['SODA_',sodaLabel]).lat;
SA = gsw_SA_from_SP(averageData.Average_Salinity_Reconstructed,pWater,lon,lat);
CT = gsw_CT_from_t(SA,averageData.Average_WaterTemperature_Corrected,pWater);
rho = gsw_rho(SA,CT,pWater);
waterDepth = dbar2pa * pWater./(rho * g);

%% Calculate ice draft

% Extract data
altimeterDistAst = averageData.Average_AltimeterIceAST_Corrected;
altimeterDistLE  = averageData.Average_AltimeterIceLE_Corrected;

% Calculate
iceDraftAST_raw = waterDepth - altimeterDistAst;
iceDraftLE_raw = waterDepth - altimeterDistLE;

%% Apply quality controls to ice draft

% Define thresholds
minDraft = 0.2;
maxEnsStd = 0.5;
minCnt = 20;
altimeterDiffThreshold = 1;
maxRng = 45; rngDraftLim = 1;
maxFOM = 55e3;

% Define "corrected" variables (in order to not over-write "raw" variables)
iceDraftAST = iceDraftAST_raw;
iceDraftLE = iceDraftLE_raw;

% Remove ice draft that are below a minimum threshold (based on a guess of
% the accuracy of the atmospheric pressure )
iceDraftAST( iceDraftAST < minDraft ) = NaN;
iceDraftLE( iceDraftLE < minDraft ) = NaN;

% Remove ice draft with large standard deviations per ensemble
%   ( these can result from a moving keel, but on a one-minute timescale
%     probably represent noise in the data )
iceDraftAST( averageDataStd.Average_AltimeterIceAST_Corrected > maxEnsStd ) = NaN;
iceDraftLE ( averageDataStd.Average_AltimeterIceLE_Corrected  > maxEnsStd ) = NaN;

% Remove ice draft with low ensemble counts
iceDraftAST( averageDataCnt.Average_AltimeterIceAST_Corrected < minCnt ) = NaN;
iceDraftLE ( averageDataCnt.Average_AltimeterIceLE_Corrected  < minCnt ) = NaN;

% Remove ice draft where the altimeter has maxed out
% ( the range can be higher for a strong reflector than for a weak signal,
%   so the maxRng value given is sort of a lower bound and I'll apply it
%   only when I start getting large draft );
iceDraftAST( iceDraftAST > rngDraftLim & altimeterDistAst > maxRng ) = NaN;
iceDraftLE ( iceDraftLE  > rngDraftLim & altimeterDistLE  > maxRng ) = NaN;
         
% Make a guess that when AST and LE measurements REALLY disagree, they're
% both bad (when probably LE is the bad one in most instances)
iceDraftAST( abs(iceDraftAST_raw-iceDraftLE_raw) > altimeterDiffThreshold ) = NaN;
iceDraftLE ( abs(iceDraftAST_raw-iceDraftLE_raw) > altimeterDiffThreshold ) = NaN;

% Use "AverageIce" Figure-of-Merit (FOM), which indicate data quality for
% ice velocity (from slant-beams).  Open-water pings seems to have very
% hight FOM.  Guess that if all 4 slant-beams see open water (indicated by
% high FOM), then there probably is no ice over the middle either
allFOM(:,:,1) = iceData.AverageIce_FOMBeam1;
allFOM(:,:,2) = iceData.AverageIce_FOMBeam2;
allFOM(:,:,3) = iceData.AverageIce_FOMBeam3;
allFOM(:,:,4) = iceData.AverageIce_FOMBeam4;
allFOMmin = min( allFOM, [], 3);
iceDraftAST( allFOMmin > maxFOM ) = NaN;
iceDraftLE ( allFOMmin > maxFOM ) = NaN;

% Retained points
retndAST = (iceDraftAST_raw == iceDraftAST);
retndLE  = (iceDraftLE_raw  == iceDraftLE );


%% Apply quality controls to ice velocity

% Define quality metric:
minCount = 20;

% Find the number of points remaining from each beam velocity after initial
% QC screening
iceCount = 0;
numBeams = 4;
for n = 1:numBeams
    iceCount = iceCount + iceDataCnt.(['AverageIce_VelBeam',num2str(n)]);
end

% Apply screening
iceData.AverageIce_VelEast ( iceCount < numBeams*minCount ) = NaN;
iceData.AverageIce_VelNorth( iceCount < numBeams*minCount ) = NaN;
iceData.AverageIce_VelUp1  ( iceCount < numBeams*minCount ) = NaN;
iceData.AverageIce_VelUp2  ( iceCount < numBeams*minCount ) = NaN;

%% Create data structure

% Initialize data structure
sigIce.mooringName = ['SODA_',sodaLabel];
sigIce.lat = lat;
sigIce.lon = lon;

% Create universal time vector
%   If ensembles are missing (e.g. SODA-C!) then iceData and averageData
%   may have different timestamps (i.e. differents sets of ensembles might
%   be missing).  Everything should be on a common vector.
datelim = sodaConsts.(['SODA_',sodaLabel]).datelim;
samplingPeriod = 10/1440; %[ hardcoding bad! :(  ]
sigIce.mattime = ( datelim(1) : samplingPeriod : (datelim(end)-samplingPeriod) ).';
numSamples = length(sigIce.mattime);

% Initialize other data structure fields
sigIce.pAbsolute        = NaN( numSamples, 1 );
sigIce.pWater           = NaN( numSamples, 1 );
sigIce.waterDepth       = NaN( numSamples, 1 );
sigIce.altimeterDistAst = NaN( numSamples, 1 );
sigIce.altimeterDistLE  = NaN( numSamples, 1 );
sigIce.iceDraftAST      = NaN( numSamples, 1 );
sigIce.iceDraftLE       = NaN( numSamples, 1 );
sigIce.velEast          = NaN( numSamples, 1 );
sigIce.velNorth         = NaN( numSamples, 1 );
sigIce.velUp1           = NaN( numSamples, 1 );
sigIce.velUp2           = NaN( numSamples, 1 );

% Populate 'Average' data fields
tol = 2/1440;
sigIce_tElapsed = sigIce.mattime - sigIce.mattime(1) ;
avgDat_tElapsed = averageData.Average_MatlabTimeStamp - sigIce.mattime(1) ;
[ avgLia , avgInd] = ismembertol( avgDat_tElapsed, sigIce_tElapsed,  tol, 'DataScale', 1);
avgInd = avgInd(avgLia);

sigIce.pAbsolute(avgInd)        = pAbsolute(avgLia);
sigIce.pWater(avgInd)           = pWater(avgLia);
sigIce.waterDepth(avgInd)       = waterDepth(avgLia);
sigIce.altimeterDistAst(avgInd) = altimeterDistAst(avgLia);
sigIce.altimeterDistLE(avgInd)  = altimeterDistLE(avgLia);
sigIce.iceDraftAST(avgInd)      = iceDraftAST(avgLia);
sigIce.iceDraftLE(avgInd)       = iceDraftLE(avgLia);

% Populate 'AverageIce' data fields
iceDat_tElapsed = iceData.AverageIce_MatlabTimeStamp - sigIce.mattime(1) ;
[ iceLia , iceInd] = ismembertol( iceDat_tElapsed, sigIce_tElapsed, tol, 'DataScale', 1);
iceInd = iceInd(iceLia);

sigIce.velEast(iceInd)  = iceData.AverageIce_VelEast(iceLia);
sigIce.velNorth(iceInd) = iceData.AverageIce_VelNorth(iceLia);
sigIce.velUp1(iceInd)   = iceData.AverageIce_VelUp1(iceLia);
sigIce.velUp2(iceInd)   = iceData.AverageIce_VelUp2(iceLia);

%% Save data product

saveDir = ['../dataProducts/'];
fName = ['SODA_',sodaLabel,'_sigIce'];
save([saveDir,fName],'sigIce');
msgbox('File saved') 


end