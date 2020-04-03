function [wordState,scoredFile] = assignStatesToWords_02252018(wordsTimeBins)

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
startTime = scoredStates(1,1) * 10^6;
endTime = (scoredStates(end,1) + epochInSeconds) * 10^6;

%% Assign states to each likely word:
numWords = size(wordsTimeBins,1);
wordState = zeros(numWords,1);
for i = 1:size(scoredStates, 1)
    if isequal(i, size(scoredStates, 1))
        subIntervalIndx = find(wordsTimeBins(:,1) >= scoredStates(i,1) & wordsTimeBins(:,1) < (scoredStates(i,1) + 10));
    else
        subIntervalIndx = find(wordsTimeBins(:,1) >= scoredStates(i,1) & wordsTimeBins(:,1) < scoredStates(i+1,1));
    end
    if ~isempty(subIntervalIndx)
        wordState(subIntervalIndx) = scoredStates(i,2);
        clear subIntervalIndx
    end
end