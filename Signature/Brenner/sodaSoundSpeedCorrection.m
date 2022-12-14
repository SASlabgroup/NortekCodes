function [Data,linFitAvg,linFitBst] = sodaSoundSpeedCorrection(Data,mcData,mode,lon,lat,datelim,makePlots)
% SODASOUNDSPEEDCORRECTION calculates corrected water temperature and sound
% speed for SODA Signature500 ADCP's
%
%   Data = sodaSoundSpeedCorrection(Data,mcData,lon,lat,datelim)
%   [Data,linFitAvg,linFitBst] = sodaSoundSpeedCorrection(...)

%% Parse inputs

    if isempty(mode); mode = {'avg'}; end
    if nargin < 7 || isempty(makePlots); makePlots = 0; end
    
    
    % Parse 'mode' choice(s)
    mode = lower(mode);
    modeChoices = {'avg','ice','burst'};
    dataWordChoices = {'Average','AverageIce','Burst'};
    [modeLog,modeInd] = ismember( mode , modeChoices );
    if ~modeLog
        error('The input variable ''mode'' must be one of: ''avg'', ''ice'', or ''burst''');
    elseif length(modeLog) > 1
        % If multiple mode words are entered, recursively run this script for
        % each of them individually (this may break something)
        for n = 1:length(lia)
            modeN = modeChoices{locb(n)};
            [Data,linFitAvgN,linFitBstN] = sodaSoundSpeedCorrection(Data,mcData,modeN,lon,lat,datelim);
            linFitAvg(n,:) = linFitAvgN;
            linFitBst(n,:) = linFitBstN;
        end
        return;
    else
        ad2cpstr = '';
        dataModeWord = dataWordChoices{modeInd};
    end

%% Identify concurrent samples indices
% Water temperatures reported by the Signature500 during 'Average'
% ensembles measured during concurruent 'Burst' ensembles have an
% temperature offset compared with those taken when 'Burst' was not being
% measured

    % Defining indices for the 'burst' samples is simple counting IF ensembles
    % are well behaved.  But I know that some SODA-C ensembles aren't.
    % Instead, I'll base the index sorting on sample plan information:
    [~,~,~,H,MI,~] = datevec(Data.([dataModeWord,'_MatlabTimeStamp']) );
    burstInds = find( rem(H,2) == 0 & MI < 18 );

    % Define indices for "non-burst" (avg) samples:
    numSamples = length( Data.([dataModeWord,'_MatlabTimeStamp']) );
    avgInds = 1:numSamples;
    avgInds( burstInds ) = [];
    
    
%% Bin-average data

    % Set limits for bin averaging
    % ( MicroCat data doesn't span the full deployment, so cut off the average
    %   based on how long it lasted )
    binDates = [ datelim(1), mcData.datenum(end) ];
    binSize = 4/24;
    
    % Bin-average temeratures
    [ ~, avgTemps] = bin_avg( Data.([dataModeWord,'_MatlabTimeStamp'])(avgInds),...
                              Data.([dataModeWord,'_WaterTemperature'])(avgInds),...
                              binSize, binDates,@nanmedian);
    [ ~, bstTemps] = bin_avg( Data.([dataModeWord,'_MatlabTimeStamp'])(burstInds),...
                              Data.([dataModeWord,'_WaterTemperature'])(burstInds),...
                              binSize, binDates,@nanmedian);
    [ ~ , mcTemps] = bin_avg( mcData.datenum, mcData.temperature,...
                              binSize, binDates,@nanmedian);
    
    % Bin-average Microcat salinity (for salinity adjustment, below)
    [ mcMattime , mcSP] = bin_avg( mcData.datenum.', real(mcData.salinity).',...
                                   binSize, binDates,@nanmedian);                      
                                   
%% Find regressions and calculated corrected temperature

    % Find regressions
    if ~isempty(avgTemps)
        linFitAvg = polyfit( avgTemps, mcTemps, 1);
    else
        linFitAvg = NaN(1,2);
    end
    if ~isempty(bstTemps)
        linFitBst = polyfit( bstTemps, mcTemps, 1);
    else
        linFitBst = NaN(1,2);
    end
   
    % Calculated corrected temperature
    tempCorrected = NaN( size(Data.([dataModeWord,'_WaterTemperature']) ) );
    tempCorrected(avgInds)  = polyval(linFitAvg, Data.([dataModeWord,'_WaterTemperature'])(avgInds) );
    tempCorrected(burstInds) = polyval(linFitBst, Data.([dataModeWord,'_WaterTemperature'])(burstInds) );
    
    % Updata data structure
    Data.([dataModeWord,'_WaterTemperature_Corrected']) = tempCorrected;
              
%% Plot temperature regressions (if requested)    
    
    if makePlots
        one2one = [-2,8];
        col = get(0,'DefaultAxesColorOrder');
        figure; clf;
        hold on;
        plot( bstTemps, mcTemps,'.','color',col(2,:) );
        plot( avgTemps, mcTemps,'.','color',col(1,:) );
        plot( one2one, polyval(linFitAvg,one2one),....
              'color',col(1,:),'linewidth',1);
        plot( one2one, polyval(linFitBst,one2one),....
              'color',col(2,:),'linewidth',1);  
        plot( one2one, one2one, 'k--','linewidth',2);
        xlabel('Signature-500 measured temperature [{\circ}C]');
        ylabel('SBE-37 measured temperature [{\circ}C]');
        
        set(gca,'xlim',one2one,'ylim',one2one,...
                'xtick',one2one(1):2:one2one(2),... 
                'ytick',one2one(1):2:one2one(2) );
        daspect([1,1,1]);
        grid on;
        
        
    end
    
    
%% Adjust salinity
% Signature500 ADCP's calculate sound speed based on a measured water
% temperature and a fixed salinity.  For SODA moorings, they were
% programmed to use 35 psu as the fixed salinity, when in reality is is
% close to 29.  Sounds speed is much less sensitive to salinity than it is
% to temperature (part of the reason that a fixed salinity is acceptable),
% but the difference between 35 and 29 is high.
% Rather than just use a fixed salinity, I will interpolate Microcat data
% onto the Signature500 timestamp, and use a varying salinity during the
% part of the record that that is available, and a fixed salinity
% there-after.

    
    % Define "default" fixed SP value to use when the microcat didn't
    % collect data
    fixedSP = nanmedian( mcSP );
    % Interpolate onto Signature500 times
    sigSP = interp1( mcMattime, mcSP, Data.([dataModeWord,'_MatlabTimeStamp']),...
                  'linear', fixedSP );
    % Convert practical salinity to absolute salinity
    sigSA = gsw_SA_from_SP(sigSP,Data.([dataModeWord,'_AltimeterPressure']),lon,lat);    
    
    % Updata Data structure (with practical salinity)
    Data.([dataModeWord,'_Salinity_Reconstructed']) = sigSP;
    
%% Calculate corrected sound speed

    Data.([dataModeWord,'_SpeedOfSound_Corrected']) = ...
           gsw_sound_speed_t_exact(sigSA,tempCorrected,Data.([dataModeWord,'_AltimeterPressure']) );

end
