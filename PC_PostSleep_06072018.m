%% Sleep-only post spikes:


%% Select Stage Scored File:
working_dir=pwd;
current_dir='C:\';
cd(current_dir);
scoredCheck = 0;
while isequal(scoredCheck, 0)
    [scoredFile, scoredPath] = uigetfile({'*.xls','Excel 1997-2003 File (*.xls)'},...
        'Select the Sleep Scored File');
    if isequal(scoredFile,0) || isequal(scoredPath,0)
        uiwait(errordlg('You need to select a file. Please try again',...
            'ERROR','modal'));
    else
        cd(working_dir);
        stageScoredFile= fullfile(scoredPath, scoredFile);
        %Load sleep scored file:
        try
            [numData, stringData] = xlsread(stageScoredFile);
            scoredCheck = 1;
        catch %#ok<*CTCH>
            % If file fails to load, it will notify user and prompt to
            % choose another file.
            uiwait(errordlg('Check if the scored file is saved in Microsoft Excel format.',...
             'ERROR','modal'));
         scoredCheck = 0;
        end

    end
end
clear scoredCheck

%% Detect if states are in number or 2-letter format:
if isequal(size(numData,2),3)
    scoredStates = numData(:,2:3);
    clear numData stringData
else
    scoredStates = numData(:,2);
    clear numData
    stringData = stringData(3:end,3);
    [stateNumber] = stateLetter2NumberConverter(stringData);
    scoredStates = [scoredStates stateNumber];
    clear stateNumber stringData
end
epochInSeconds = scoredStates(2,1) - scoredStates(1,1);
% startTime = scoredStates(1,1) * 10^6;
% endTime = (scoredStates(end,1) + epochInSeconds) * 10^6;

%% Reformat spike data into one long array:
numCells = size(PC_ID, 1); %# of units/cells based on # of text files loaded
allSpikeTS = [];
allID = [];
for i = 1:numCells
    spikeTimes = cell2mat(PC_post(i));
    allSpikeTS = [allSpikeTS; spikeTimes];
    allID = [allID; i*ones(size(spikeTimes))];
    clear spikeTimes
end
clear numCells

% Sort all spike data by time stamp:
[allSpikeTS, IX] = sort(allSpikeTS,1);
allID = allID(IX);
clear IX

%% Assign states to each likely word:
numSpikes = size(allSpikeTS,1);
spikeState = zeros(numSpikes,1);
for i = 1:size(scoredStates, 1)
    if isequal(i, size(scoredStates, 1))
        subIntervalIndx = find(allSpikeTS >= scoredStates(i,1) & allSpikeTS < (scoredStates(i,1) + 10));
    else
        subIntervalIndx = find(allSpikeTS >= scoredStates(i,1) & allSpikeTS < scoredStates(i+1,1));
    end
    if ~isempty(subIntervalIndx)
        spikeState(subIntervalIndx) = scoredStates(i,2);
        clear subIntervalIndx
    end
end

targNREM = spikeState== 2 | spikeState== 6;
targREM = spikeState== 3;
targSleep = spikeState== 2 | spikeState== 6 | spikeState== 3;

spikeTS_NREM = allSpikeTS(targNREM);
spikeID_NREM = allID(targNREM);

spikeTS_REM = allSpikeTS(targREM);
spikeID_REM = allID(targREM);

spikeTS_Sleep = allSpikeTS(targSleep);
spikeID_Sleep = allID(targSleep);
clear targNREM targREM targSleep allSpikeTS allID subIntervalIndx spikeState

%% Convert spike data to cell arrays
numPC = size(PC_ID,1);
PC_NREM = cell(numPC,1);
PC_REM = PC_NREM;
PC_Sleep = PC_NREM;
for i = 1:numPC 
    targetPost = ismember(spikeID_NREM, i);
    if ~isempty(targetPost)
        PC_NREM{i,1} = spikeTS_NREM(targetPost);
    end
    
    targetPost = ismember(spikeID_REM, i);
    if ~isempty(targetPost)
        PC_REM{i,1} = spikeTS_REM(targetPost);
    end
    
    targetPost = ismember(spikeID_Sleep, i);
    if ~isempty(targetPost)
        PC_Sleep{i,1} = spikeTS_Sleep(targetPost);
    end
end
clear numPC spikeID_NREM spikeID_REM spikeID_Sleep spikeTS_NREM spikeTS_REM...
    spikeTS_Sleep targetPost i 

