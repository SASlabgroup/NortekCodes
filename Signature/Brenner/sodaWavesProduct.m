function sodaWavesProduct( sodaLabel )
% SODAWAVESPROCESS Process and create waves data product for SODA
% Signature500 ADCP's
%
%   sodaWavesProcess( sodaLabel ), for 'sodaLabel' one of: 'A', 'B', or 'C'
%
%   Note: this function requires a specific relative folder-path setup in
%   order to successfully find, load, and then save relevant data.  The
%   function has no direct output, but saves processed data to be used for
%   further analysis.
%
%   S.D.Brenner, 2019

%% Load data
% [ NOTE: loading burst data can take a considerable amount of time ]

dataDir = ['../data/SODA_',sodaLabel,'/Processed_Data/'];
fNameSig = ['SODA_',sodaLabel,'_Burst.mat'];

load('sodaConstants.mat');

% for testing: don't waste time loading burst data variable if it's already
% in the workspace
if ~exist('burstData','var') 
    load( [dataDir,fNameSig] );
end

%% Set constants

lat = sodaConsts.(['SODA_',sodaLabel]).lat;
lon = sodaConsts.(['SODA_',sodaLabel]).lon;
waterdepth = 3800; % [m] approximately; don't need that high of accuracy at this depth.

%% Loop 

% Find indices when the ensemble count stops increasing:
deltaCount = diff(double(burstData.Burst_EnsembleCount));
negInd = find( deltaCount<0 );
negInd = [ 0; negInd; length(burstData.Burst_EnsembleCount) ];
L = length(negInd)-1;


% Loop through ensembles and process
% ( speed improvements could be made by pre-allocating sigWaves, but its
%   not clear to me how to pre-allocate a structure-arry with nested
%   structure fields )
for l = 1:L
    burstInd = (negInd(l)+1) : (negInd(l+1));
    singleBurstData = subsetData( burstData, burstInd);
    sigWaves(l) = sigWavesProcess( singleBurstData ,waterdepth,lon,lat); 
end

%% Apply additional screening

% [ screening criteria to be determined ]


%% Save data product

saveDir = ['../dataProducts/'];
fName = ['SODA_',sodaLabel,'_sigWaves'];
save([saveDir,fName],'sigWaves');
msgbox('File saved') 



end

%% EMBEDDED FUNCTIONS %% ==================================================

function Data = subsetData(Data,ind)
    flds = fields(Data);
    numFlds = length(flds);
    for n = 1:numFlds
        fldName = flds{n};
        [numSamples,~] = size( Data.(fldName) );
        if numSamples > 1
            Data.(fldName) = Data.(fldName)(ind,:);
        end
    end
end