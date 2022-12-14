function sodaConcatRaw( sodaLabel )
% SODACONCATRAW Concatenate raw data for SODA Signature500 ADCP's
%
%   sodaConcatRaw( sodaLabel ), for 'sodaLabel' one of: 'A', 'B', or 'C'.
%
%   Note: this function requires a specific relative folder-path setup in
%   order to successfully find, load, and then save relevant data.  The
%   function has no direct output, but saves concatenated data to be used
%   for further processing.
%
%   S.D.Brenner, 2019

%% Define data and save directories

dataDir = ['../data/SODA_',sodaLabel,'/Converted_MIDAS/'];
saveDir = ['../data/SODA_',sodaLabel,'/Concatenated_Matfiles/'];

%% Generate lists of files and fields

files = dir([dataDir,'*.mat']);
load( [dataDir,files(1).name] );
Data = rmfield(Data,{'Units','Comments'});
flds = fields( Data );


%% Get list of field prefixes ("dataModWords")
% To redude filesizes, each 'prefix' (e.g. "Average_", "Burst_", etc.) will
% be eventually be contained in a separate data structure and saved in a
% different file. Field prefixes are identified by underscores (us); to
% find them, look for the index of the first underscore in each field name.

usInd = strfind(flds,'_');
for n = 1:length(usInd)
    dataModeWord{n} = flds{n}(1:usInd{n}-1);
end


%% Pre-allocate fields sizes for data concatenation

deploymentDays = 380;


numFlds = length(flds);
wbHand = waitbar(0,'Pre-allocating fields...');
for n = 1:numFlds
    waitbar( n/numFlds,wbHand );
    fldName = flds{n};
    % Get field size
    [numSamples,numBins] = size( Data.(fldName) );
    if numSamples > 1 
        % Extrapolate size of numSamples to total deployment length
        % (it's better if this is an overestimate)
        numDays = range(Data.([dataModeWord{n},'_MatlabTimeStamp']) );
        deploymentSamples = ceil( numSamples * deploymentDays/numDays );
        % Pre-allocate field
        allData.(fldName) = NaN(deploymentSamples, numBins);
    else % (static field)
        allData.(fldName) = NaN(1, numBins);
    end
end
waitbar( 1,wbHand,'Pre-allocation complete' );


%% Concatenate ALL raw data

numFiles = length(files);
wbHand = waitbar(0,wbHand,'Concatenating raw data files');
for l = 1:numFiles
    % Output loop progress   
	waitbar(l/numFiles,wbHand,... 
            sprintf('Concatenating raw data files:\n%g of %g',l,numFiles) );    
    % Load data
    fName = files(l).name;
    load( [dataDir,fName] );  
    % Fix burst config data for SODA-C
    if strcmp(sodaLabel,'C')
        Config.burst_nCells = 16;
        Config.burst_nBeams = 4;
        Config.burst_nSamples = 2048;
    end
    % Loop through fields and copy data to concatenated structure
    for n = 1:numFlds
        fldName = flds{n};      
        [numSamples,~] = size( Data.(fldName) );
        if numSamples > 1
            % Define indices to write
            indStart = find(~isnan(allData.(fldName)(:,1)),1,'last') + 1;
            if isempty(indStart)
                indStart = 1;
            end
            indEnd = indStart + numSamples - 1 ;
            dataInd = indStart:indEnd;
            if indEnd > size(allData.(fldName),1)
                disp([fldName,' pre-allocation insufficient'])
            end
            % Add data to concatenated structure in selected index range
            allData.(fldName)(dataInd,:) = Data.(fldName);
        elseif numSamples == 1 && l == 1    
            % If the variable is static (e.g. 'Ranges'), only add it to
            % concatenated structure once
            allData.(fldName) = Data.(fldName);
        end
    end
end
waitbar( 1,wbHand,'Concatenation complete' );

%% Remove extraneous data entries from "static" variables

% Loop through all fields
wbHand = waitbar(0,wbHand,'Checking for static fields');
for n = 1:numFlds
    % Output loop progress   
	waitbar(n/numFlds,wbHand,... 
            sprintf('Checking for static fields:\n%g of %g',n,numFlds) );   
    % Read field, and check if it's static    
    fldName = flds{n};
    if isStatic( allData.(fldName) )
        % if it is static, we only need the first row of values
        allData.(fldName) = allData.(fldName)(1,:);
    end
end

%% Remove trailing NaNs

% Loop through all fields
wbHand = waitbar(0,wbHand,'Removing trailing NaN''s...');
for n = 1:numFlds
    % Output loop progress   
	waitbar(n/numFlds,wbHand,... 
            sprintf('Removing trailing NaN''s:\n%g of %g',n,numFlds) );   
    % Get field name and size
    fldName = flds{n};
    [numSamples,~] = size( allData.(fldName) );    
    % Get indices of NaN's 
    % ( based on the matlab timestamp vectors corresponding to each prefix)
    if n == 1
        % Only get indices if I don't already have them
        nanInd = isnan( allData.([dataModeWord{n},'_MatlabTimeStamp']) );
    elseif n > 1 && ~strcmp(dataModeWord{n},dataModeWord{n-1})
         nanInd = isnan( allData.([dataModeWord{n},'_MatlabTimeStamp']) );
    end
    if numSamples > 1 % Only adjust "non-static" fields
        allData.(fldName)(nanInd,:) = [];
    end
end
waitbar( 1,wbHand,'Removing NaN''s complete' );

%% Sort data chronologically
% File names don't reflect chronological order

wbHand = waitbar(0,wbHand,'Sorting data chronologically...');
for n = 1:numFlds
    % Output loop progress   
	waitbar(n/numFlds,wbHand,... 
            sprintf('Sorting data chronologically:\n%g of %g',n,numFlds) );   
    % Get field name and size
    fldName = flds{n};
    [numSamples,~] = size( allData.(fldName) );      
    % Get sorting indices 
    if n == 1
        % Only get indices if I don't already have them
        [~,sortInd] = sort( allData.([dataModeWord{n},'_MatlabTimeStamp']) );
    elseif n > 1 && ~strcmp(dataModeWord{n},dataModeWord{n-1})
        [~,sortInd] = sort( allData.([dataModeWord{n},'_MatlabTimeStamp']) );
    end
    if numSamples > 1 % Only adjust "non-static" fields
        allData.(fldName) = allData.(fldName)(sortInd,:);
    end
end
waitbar( 1,wbHand,'Sorting complete' );
delete(wbHand);

%% Separate data structure into different structures based on prefixes

prefixes = unique(dataModeWord);
% remove empty "{1×0 char}" result that shows up in prefixes:
% prefixes = prefixes( ~cellfun(@(S) isempty(S), prefixes) );
numPrefixes = length(prefixes);

for m = 1:numPrefixes
    dataFldName = [ lower(prefixes{m}) ,'Data' ];
    sig.(dataFldName) = subsetFields( allData, [prefixes{m},'_'] );
end

%% Save each structure separately
% These files are massive (e.g., SODA-B data builds to approx 25 GB*).

% Include a waitbar that updates based on variable memory sizes (instead of
% numbers)
W = whosFields(sig);
Wconfig = whos('Config');
totalMem = sum([W.bytes])+Wconfig.bytes;

% Loop and save
wbHand = waitbar(0,'Saving files...');
memSaved = 0;
for m = 1:numPrefixes
    dataFldName = [ lower(prefixes{m}) ,'Data' ];
    % output progress
    varSize = byteConv( W(m).bytes );
    waitbar(memSaved/totalMem,wbHand,...
            sprintf('Saving %s...\n( %s variable size )',dataFldName,varSize ) );
    % Save
    saveFileName = ['SODA_',sodaLabel,'_',prefixes{m},'_raw.mat'];
    save([saveDir,saveFileName],'-struct','sig',dataFldName,'-v7.3');
    memSaved = memSaved + W(m).bytes;
    
    waitbar(memSaved/totalMem,wbHand,['Saved ',dataFldName]);
    
end
waitbar(memSaved/totalMem,wbHand,...
        sprintf('Saving %s...\n','Config' ) );
save([saveDir,'SODA_',sodaLabel,'_Config.mat'],'Config');
delete(wbHand);
msgbox('All files saved') 


%%

% %% SPEED TEST!
% 
% % first, I need to go through and mix up all of the time indices so that
% % there will be something to sort!
% % (this might actually be a little tricky...)
% prefixes = unique(dataModeWord);
% M = length( prefixes );
% 
% L = 10;
% T = NaN(2,L*numFlds);
% m = 1;
% for l = 1:L
% 
% for m = 1:M
%     if isfield(Data,[prefixes{m},'_MatlabTimeStamp'])
%         numSamples = length( Data.([prefixes{m},'_MatlabTimeStamp']) );
%         P = randperm( numSamples);
%         Data.([prefixes{m},'_MatlabTimeStamp']) = Data.([prefixes{m},'_MatlabTimeStamp'])(P);
%     end
% end
% 
% 
% for n = 1:numFlds
%     % METHOD 1
%     tic;
%     [~,sortInd] = sort( Data.([dataModeWord{n},'_MatlabTimeStamp']) );
%     T1n = toc;
%     T(1,m) = T1n;
%     
%     % METHOD 2
%     tic;
%     if n>1 && ~strcmp(dataModeWord{n},dataModeWord{n-1})
%         [~,sortInd] = sort( Data.([dataModeWord{n},'_MatlabTimeStamp']) );
%     end
%     T2n = toc;
%     T(2,m) = T2n;
%     m = m + 1;
% end
% end
% 
% minT = min(T(:));
% maxT = max(T(:));
% ts = logspace( floor(log10(minT)), ceil(log10(maxT)), 51 );
% figure(1); clf;
% hold on;
% histogram( T(1,:), ts);
% histogram( T(2,:), ts);
% set(gca,'xscale','log');
% legend('Method 1','Method 2');



end

%% EMBEDDED FUNCTIONS %% ==================================================

function static = isStatic(var)
%     static = ~all( nansum(abs(diff( var ))) );
    ind = ~isnan( var(:,1) );
    static = all( var(ind,:) == var(1,:) ,'all');
end


function [subsetData,subsetFlds] = subsetFields( Data, dataModeWord)
    % list all fields in data structure
    flds = fields(Data);
    % Find and list only fields that start with 'AverageIce_'
    subsetFlds = flds( contains(flds,dataModeWord) );
    % Loop and subset
    for n = 1:length(subsetFlds)
        fldName = subsetFlds{n};
        subsetData.(fldName) = Data.(fldName);
    end
end


function W = whosFields(S)
    flds = fields(S);
    numFlds = length(flds);
    for n = 1:numFlds
        A = getfield( S, flds{n} );
        W(n) = whos('A');
        W(n).name = flds{n};
    end
end


function str = byteConv(n)
    log1024 = @(x) log(x)/log(1024);
    order = round( log1024(n) );
    unitChoices = {'B','KB','MB','GB','TB'};
    order(order>4) = 4;
    order(order<0) = 0;
    str = sprintf('%3.3g %s', n./(1024.^order), unitChoices{[order+1]} );
end