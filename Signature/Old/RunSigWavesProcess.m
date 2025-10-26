% script to call waves processing for Signature500 ADCP's
%
%   S.D.Brenner, 2019
% J. Thomson, Nov 2019... other export format
%       also make more generic, less SODA-specific
%
% could expand to process avg data also
%
clear all
tic

%% Load data
% [ NOTE: loading burst data can take a considerable amount of time ]

dataDir = ['./'];

fNameSigBase = ['CODA_S1A1_Sig_2020_'];

flist = dir([dataDir fNameSigBase '*.mat']);


%% Set constants

lat = 70.48695; lon = -162.28278; waterdepth = 30; % CODA S1-A1

%% Loop through files
counter = 1;

for fi=1:length(flist)

    load([dataDir flist(fi).name])
    % check for ENU velocities 
    if ~isfield(Data,'Burst_VelNorth')
        disp('computing ENU rotations')
        [ Data, Config, T_beam2xyz ] = signatureAD2CP_beam2xyz_enu_altexport( Data, Config, 'burst',0,0);
        save([dataDir flist(fi).name],'Data','Config','Units','Descriptions')
    end
        

%% Loop through bursts

% Find indices when the ensemble count stops increasing:
deltaCount = diff(double(Data.Burst_EnsembleCount));
negInd = find( deltaCount<0 );
negInd = [ negInd; length(Data.Burst_EnsembleCount) ];
L = length(negInd)-1;


% Loop through ensembles and process
% ( speed improvements could be made by pre-allocating sigWaves, but its
%   not clear to me how to pre-allocate a structure-arry with nested
%   structure fields )
for l = 1:L
    burstInd = (negInd(l)+1) : (negInd(l+1));
    singleBurstData = subsetData( Data, burstInd);
    %[sigWaves(counter), astWaves(counter), leWaves(counter), presWaves(counter)] = sigWavesProcess( singleBurstData ,waterdepth,lon,lat); 
    [sigWaves(counter)] = sigWavesProcess( singleBurstData ,waterdepth,lon,lat); 
    counter = counter + 1;
end

end
%% Apply additional screening and sorting

% [ screening criteria to be determined ]
bad = isnan([sigWaves.sigwaveheight]);
sigWaves(bad) = [];
%astWaves(bad) = [];
%leWaves(bad) = [];
%presWaves(bad) = [];

[newtime indices] = sort([sigWaves.time]);
sigWaves = sigWaves(indices);


%% Save data product

saveDir = ['./'];
save([saveDir 'sigWaves_' fNameSigBase ],'sigWaves');
%msgbox('File saved') 

%% plot
plotSWIFT(sigWaves)

toc

%% EMBEDDED FUNCTIONS %% ==================================================

function Data = subsetData(Data,ind)
    flds = fields(Data);
    numFlds = length(flds);
    for n = 1:numFlds
        fldName = flds{n};
        if fldName(1:6) == 'Burst_',
        [numSamples,~] = size( Data.(fldName) );
        if numSamples > 1
            Data.(fldName) = Data.(fldName)(ind,:);
        end
        else
        end
    end
end