function spikeDataForVLMC_02192018
% Condition data in format for VLMC.
%--INPUT: NTT files in the selected folder
%--OUTPUT: .MAT file containing: Tetrode #, Unit # of tetrode, time of each
%spike

%% Load NTT files
% Select folder and get list of NTT files:
fileType = '*.ntt';
[dataFolder, fileList, numberOfDataFiles] = batchLoadFiles(fileType);
clear fileType

%% For each tetrode, load in the data and combine all spikes into one file:
allSpikeTS = [];
allID = [];
for i = 1:numberOfDataFiles
    fileName = strtrim(fileList(i,:)); %Removes any white space at end of file name string.
    tetrodeFile = fullfile(dataFolder,fileName); %Full file path for Neuralynx file to be loaded
    
    %% Import the data:
    % Set up variable for sorted Neuralynx NTT file:
    [spikeTimes, tetrodeNum, cellNumber] = Nlx2MatSpike(tetrodeFile, [1 1 1 0 0], 0, 1, []); % Load only time stamps, cell #, and amplitudes
    
    %% Remove unsorted spikes:
    nonZerosIndex = find(cellNumber);       % Identify spikes with a unit assignment
    cellNumber = cellNumber(nonZerosIndex)'; % Remove cell # unsorted spikes
    spikeTimes = spikeTimes(nonZerosIndex)' ./ 1000000;   % Remove time stamps of unsorted spikes and convert to seconds
    tetrodeNum = tetrodeNum(nonZerosIndex)' + 1; %TT # starts at 0 in the file data
    cellIdentity = [tetrodeNum cellNumber];
    clear nonZerosIndex cellNumber tetrodeNum
    allSpikeTS = [allSpikeTS; spikeTimes];
    allID = [allID; cellIdentity];
    clear spikeTimes cellIdentity fileName
end
clear numberOFDataFiles

%% Sort all spike data by time stamp
[allSpikeTS, IX] = sort(allSpikeTS,1);
allID = allID(IX,:);
clear IX

%% Run Segments
runBounds = [segmentTimeStamps{1,1}(1,1) segmentTimeStamps{1,30}(end,2)];
runIdx = allSpikeTS >= runBounds(1,1) & allSpikeTS < runBounds(1,2);
runSpikeTS = allSpikeTS(runIdx);
runID = allID(runIdx,:);
runSegmentsTS = segmentTimeStamps;
clear runIdx segmentTimeStamps

%% Post-Run Segments
postIdx = allSpikeTS >= sleepBounds(1,1) & allSpikeTS < sleepBounds(1,2);
postSpikeTS = allSpikeTS(postIdx);
postID = allID(postIdx,:);
clear postIdx allID allSpikesTS

%% Remove cells not in both RUN and Post
% Find units in both Run and Post
uniqueRun = unique(runID, 'rows');
uniquePost = unique(postID, 'rows');
bothRunPost = intersect(uniqueRun,uniquePost,'rows');

% Remove units from Run that are not in both Run and Post
isBothRun = ismember(runID, bothRunPost, 'rows');
runSpikeTS = runSpikeTS(isBothRun);
runID = runID(isBothRun,:);
clear isBothPost

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
resultsFolder = 'Z:\Data Analysis\Optogenetic_AnalzedData\VLMC_Analyses';
save(fullfile(resultsFolder,['conditionData', filename, '.mat']),...
    'fileList', 'NPC_ID', 'NPC_run', 'NPC_post', 'PC_ID', 'PC_run',...
    'PC_post', 'runBounds', 'sleepBounds', 'runSegmentsTS'); %, 'state', 'stateTS');


