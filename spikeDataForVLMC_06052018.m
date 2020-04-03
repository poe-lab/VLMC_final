function spikeDataForVLMC_06052018(placeCells, PC_Order, sleepBounds)
% Condition data in format for VLMC.
%--INPUT: .MAT file generated by spikeDataCellArray_########.m containing
%cellID and cell array of spike time stamps in usec
%--OUTPUT: .MAT file containing: Tetrode #, Unit # of tetrode, time of each
%spike in seconds

%% Select all of the files to include in RUN part of the spike data:
% These are the files generated in 'MazeTransformRun_06032018' for each set of laps. 
[fileSet,pathName] = uigetfile('*.mat','Select file(s) RUN portion of spikes','MultiSelect', 'on');
numFiles = size(fileSet,2); % # of files selected
% allSets = cell(numFiles, 1);

%% Load RUN and STOP spikes from all sets of laps:
spikesRUN = [];
spikesSTOP = [];
boundsRUN = [];
for i = 1:numFiles
    load(fullfile(pathName, fileSet{1,i}),'-mat', 'BV')
    boundsRUN = [boundsRUN; BV.runBounds];
    spikesRUN = [spikesRUN; BV.RUN_spikes];
    spikesSTOP = [spikesSTOP; BV.STOP_spikes];
    clear BV
end
clear numFiles

%% Convert time from usec to sec:
runSegments = boundsRUN./10^6;
clear boundsRUN
runBounds = [runSegments(1,1) runSegments(end,2)];
spikesRUN(:,1) = spikesRUN(:,1)./10^6;
spikesSTOP(:,1) = spikesSTOP(:,1)./10^6;

%% POST spikes portion:
% Load spike data for getting POST portion:
[fileName,pathName] = uigetfile('*.mat','Select file for Place Field Analyses');
load(fullfile(pathName, fileName),'-mat', 'cellID', 'spikeCellArray');
clear pathName

% Reformat spike data into one long array:
numCells = size(cellID, 1); %# of units/cells based on # of text files loaded
allSpikeTS = [];
allID = [];
for i = 1:numCells
    spikeTimes = cell2mat(spikeCellArray(i));
    allSpikeTS = [allSpikeTS; spikeTimes];
    allID = [allID; i*ones(size(spikeTimes))];
    clear spikeTimes
end
clear numCells spikeCellArray

% Sort all spike data by time stamp:
[allSpikeTS, IX] = sort(allSpikeTS,1);
allSpikeTS = allSpikeTS./10^6; % Convert to sec
allID = allID(IX);
clear IX

%% Post-Run Segments
postIdx = allSpikeTS >= sleepBounds(1,1) & allSpikeTS < sleepBounds(1,2);
postSpikeTS = allSpikeTS(postIdx);
ID_post = allID(postIdx);
postID = cellID(ID_post,:);
clear postIdx allID allSpikeTS ID_post

%% Run Segments
runSpikeTS = spikesRUN(:,1);
runID = cellID(spikesRUN(:,2),:);
clear spikesRUN

%% Stop Segments
stopSpikeTS = spikesSTOP(:,1);
stopID = cellID(spikesSTOP(:,2),:);
clear spikesSTOP cellID

%% Remove cells not in both RUN and Post
% Find units in both Run and Post
uniqueRun = unique(runID, 'rows');
uniquePost = unique(postID, 'rows');
bothRunPost = intersect(uniqueRun,uniquePost,'rows');
clear uniqueRun uniquePost

% Remove units from Run that are not in both Run and Post
isBothRun = ismember(runID, bothRunPost, 'rows');
runSpikeTS = runSpikeTS(isBothRun);
runID = runID(isBothRun,:);
clear isBothRun

% Remove units from Post that are not in both Run and Post
isBothPost = ismember(postID, bothRunPost, 'rows');
postSpikeTS = postSpikeTS(isBothPost);
postID = postID(isBothPost,:);
clear isBothPost bothRunPost

%% Separate place cells from non-place cells in RUN and POST
% Separate for RUN:
isPlaceCell = ismember(runID, placeCells, 'rows');
runSpikeTS_PC = runSpikeTS(isPlaceCell);
runID_PC = runID(isPlaceCell,:);
runSpikeTS_NPC = runSpikeTS(~isPlaceCell);
runID_NPC = runID(~isPlaceCell,:);
clear isPlaceCell runSpikeTS runID

% Separate for POST:
isPlaceCell = ismember(postID, placeCells, 'rows');
postSpikeTS_PC = postSpikeTS(isPlaceCell);
postID_PC = postID(isPlaceCell,:);
postSpikeTS_NPC = postSpikeTS(~isPlaceCell);
postID_NPC = postID(~isPlaceCell,:);
clear isPlaceCell postSpikeTS postID placeCells

%% Convert spike data to cell arrays
PC_ID = unique(postID_PC, 'rows');
numPC = size(PC_ID,1);
PC_run = cell(numPC,1);
PC_post = PC_run;
for i = 1:numPC
    targetRun = ismember(runID_PC, PC_ID(i,:), 'rows');
    PC_run{i,1} = runSpikeTS_PC(targetRun);
    clear targetRun
    
    targetPost = ismember(postID_PC, PC_ID(i,:), 'rows');
    PC_post{i,1} = postSpikeTS_PC(targetPost);
    clear targetPost
end
clear numPC runSpikeTS_PC runID_PC postSpikeTS_PC postID_PC

NPC_ID = unique(postID_NPC, 'rows');
numNPC = size(NPC_ID,1);
NPC_run = cell(numNPC,1);
NPC_post = NPC_run;
for i = 1:numNPC
    targetRun = ismember(runID_NPC, NPC_ID(i,:), 'rows');
    NPC_run{i,1} = runSpikeTS_NPC(targetRun);
    clear targetRun
    
    targetPost = ismember(postID_NPC, NPC_ID(i,:), 'rows');
    NPC_post{i,1} = postSpikeTS_NPC(targetPost);
    clear targetPost
end
clear numNPC runSpikeTS_NPC runID_NPC postSpikeTS_NPC postID_NPC

%% SAVE DATA AND FIGURE
%Request user to name output file:
prompt = {'Enter the filename you want to save it as: (just the name)'};
def = {'Rat#_Day'};
dlgTitle = 'Save .MAT file';
lineNo = 1;
answer = inputdlg(prompt,dlgTitle,lineNo,def);
filename = char(answer(1,:));
clear prompt def dlgTitle lineNo answer
resultsFolder = uigetdir(pwd,'Select Folder to Save Results');
save(fullfile(resultsFolder,['conditionData', filename, '.mat']),...
    'fileSet', 'fileName', 'NPC_ID', 'NPC_run', 'NPC_post', 'PC_ID', 'PC_run',...
    'PC_post', 'runBounds', 'sleepBounds', 'runSegments', 'PC_Order', 'stopSpikeTS', 'stopID'); %, 'state', 'stateTS');
clear runBounds sleepBounds runSegments fileSet fileName NPC_ID NPC_run NPC_post stopSpikeTS stopID i

