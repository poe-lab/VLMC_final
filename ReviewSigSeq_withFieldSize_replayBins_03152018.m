function ReviewSigSeq_withFieldSize_replayBins_03152018
% All LFPs must have same sampling rate and time stamps for function to
% work properly.
orderPC = PC_post(PC_Order);
orderPC_ID = PC_ID(PC_Order,:);
PC_cutoff= 0;
numCells = size(orderPC,1)-PC_cutoff;
trackLength = 297; % approximate circumference of the track

numEvents = size(replay.wordsTimeBins,1); %Determine number of detected spindles on the selected LFP/EEG channel
% Determine maximum y-axis values (mV) for the spindle filtered signals plus
% added 10%
maxY_HcSpindle = max(abs(HcSpindleSignal)) + 0.1*max(abs(HcSpindleSignal));
maxY_CtxSpindle = max(abs(smoothdata(ctxSpindleSignal))) + 0.1*max(abs(smoothdata(ctxSpindleSignal)));
timeWindow= 1.5; %3 second window
% maxSigSyllables = max(cellfun('length', replay.syllableIdxAll(logical(replay.inSpindle))));
mazeCodeNew = mazeCode+2;
mazeMarkColor = {'g';'none';'k'};
mazeMarker = {'.';'none';'.'};


%% Calculate and plot results based on each detected spindle centered in a 3 sec window:
for i = 1:numEvents
%     if ~isempty(spindleWordsReplay.idx{i})
    eventStartTS = replay.wordsTimeBins(i,1);
    eventStopTS = replay.wordsTimeBins(i,2);
    addTime = timeWindow - (eventStopTS - eventStartTS);
    if addTime < 0   % If the spindle is the max allowed of 3 sec duration,
        addTime = 0; %do not add any data points to the beginning and end.
    end
    startTime = replay.wordsTimeBins(i,1) - addTime/2;
    stopTime = replay.wordsTimeBins(i,2) + addTime/2;

    spikesInWindow = [];
    for c = 1:numCells                                            % For each place cell
        spikes = orderPC{c};                                      % Keep all its spikes
        spikes = spikes(spikes >= startTime & spikes <= stopTime); % Keep all its spikes that fall within the target interval

        if ~isempty(spikes) % If there are any spikes
            spikesInWindow = [spikesInWindow; [spikes c*ones(size(spikes,1),1)]];
        end   
    end

    ax1 = subplot(4,1,1);
    targetIdx = TimeStamps >= startTime & TimeStamps <= stopTime;
    plot(ax1, TimeStamps(targetIdx), HcSpindleSignal(targetIdx), 'b')
    hold on

    plot(ax1, TimeStamps(targetIdx), HcSignal(targetIdx),...
        'Color', [0.5 0.5 0.5])

    targetIdx = TimeStamps >= eventStartTS &...
        TimeStamps <= eventStopTS;
    plot(ax1, TimeStamps(targetIdx),...
        HcSpindleSignal(targetIdx), 'g')

    targetIdx = TimeStamps >= startTime & TimeStamps <= stopTime;
    plot(ax1, TimeStamps(targetIdx), rippleSignal(targetIdx), 'r')
    hold off

    switch replay.wordState(i)
        case 2
            stage = 'SWS';
        case 3
            stage = 'REM';
        case 4
            stage = 'QW';
        case 6
            stage = 'TR ';
    end
    title(['Event ' num2str(i) ' of ' num2str(numEvents) '          Stage: ' stage])
    xlim([startTime stopTime])
    ylim([-maxY_HcSpindle maxY_HcSpindle])
    ylabel(ax1, 'HC LFP')

    %% Create a plot with raster height on the y-axis representing place field location and width
    ax2 = subplot(4,1,2);
    map=lines(7);

    if isempty(fieldLocations) % field locations are not provided
%             % Create a normal raster plot
        fieldLocations = [PC_Order-0.5, PC_Order+0.5];
    end

    %% Create a field width raster plot
%             colorUnit = map(mod(spikesInSpindle(1,2),7)+1,:);
    for numFields = 1:size(fieldLocations{spikesInWindow(1,2),1},1)
        plot(ax2, [spikesInWindow(1,1) spikesInWindow(1,1)],...
            [fieldLocations{spikesInWindow(1,2),1}(numFields,1)...
            fieldLocations{spikesInWindow(1,2),1}(numFields,2)],...
            'Color', [0.5 0.5 0.5], 'LineWidth', 0.1, 'LineStyle', '-')
%                     'Color', colorUnit, 'LineWidth', 0.25)
    end
    hold on
    for m = 2:size(spikesInWindow,1)
%                 colorUnit = map(mod(spikesInSpindle(m,2),7)+1,:);             
        for numFields = 1:size(fieldLocations{spikesInWindow(m,2),1},1)
            plot(ax2, [spikesInWindow(m,1) spikesInWindow(m,1)],...
                [fieldLocations{spikesInWindow(m,2),1}(numFields,1)...
                fieldLocations{spikesInWindow(m,2),1}(numFields,2)],...
                'Color', [0.5 0.5 0.5], 'LineWidth', 0.1, 'LineStyle', '-')
%                     'Color', colorUnit, 'LineWidth', 0.25)
        end
    end

% Plot dashed borders for each event
    plot(ax2, [eventStartTS eventStartTS], [0 trackLength+1], ':k')
    plot(ax2, [eventStopTS eventStopTS], [0 trackLength+1], ':k')

    %% Get cleaned spikes for  replay bins in the word:
    spikeBurstsInWord = [];
    for k = 1:size(replay.spikeSeq{i,1},1)
        spikeBurstsInWord =[spikeBurstsInWord ; replay.spikeSeq{i,1}{k,1}];
    end

    %% Highlight the burst spikes in the raster:
    for m = 1:size(spikeBurstsInWord,1)
        colorUnit = map(mod(spikeBurstsInWord(m,2),7)+1,:);             
        for numFields = 1:size(fieldLocations{spikeBurstsInWord(m,2),1},1)
            plot(ax2, [spikeBurstsInWord(m,1) spikeBurstsInWord(m,1)],...
                [fieldLocations{spikeBurstsInWord(m,2),1}(numFields,1)...
                fieldLocations{spikeBurstsInWord(m,2),1}(numFields,2)],...
                'Color', colorUnit,'LineWidth', 2,...
                'Marker', char(mazeMarker(mazeCodeNew(spikeBurstsInWord(m,2)))),...
                'MarkerEdgeColor', char(mazeMarkColor(mazeCodeNew(spikeBurstsInWord(m,2)))),...    
                'MarkerSize', 6)
        end
    end    
        

     %% Find any LC spikes in window and plot them:
    numLC = size(lcUnits,1);
    if ~isempty(lcUnits)
        eventLC_spikes = [];
        for iLC = 1:numLC                                                         % For each place cell
            spikesLC = lcUnits{iLC};                                      % Keep all its spikes
            spikesLC = spikesLC(spikesLC >= startTime & spikesLC <= stopTime); % Keep all its spikes that fall within the target interval
            if ~isempty(spikesLC)                                                 % If there are any spikes
                eventLC_spikes = [eventLC_spikes; [spikesLC (trackLength+iLC)*ones(size(spikesLC,1),1)]]; 
            end
        end

        if ~isempty(eventLC_spikes)
            scatter(ax2,eventLC_spikes(:,1),eventLC_spikes(:,2), 'r', '*')
        end   
    end

    hold off
    ylim(ax2,[0 trackLength+numLC])
    xlim(ax2,[startTime stopTime])
    ylabel(ax2, 'Place Fields') 


    %% Plot significant forward and reverse replay
    ax3 = subplot(4,1,3);
    hold on
    maxNumSyll = [];

    %% Get cleaned spikes for  replay bins in the word:
    for k = 1:size(replay.spikeSeq{i,1},1)
        % Get all spike bursts for replay bin:
        tempReplay = replay.spikeSeq{i,1}{k,1};
        % Get start and stop indices of syllables:
        syllIdx = replay.syllableIdxAll{i,1}{k,1};
        tempRho = replay.rho{i,1}{k,1};
        tempPval = replay.pVal{i,1}{k,1};
        tempPval = -tempPval;
        minP = min(tempPval);
        rangeP = max(tempPval) - minP;
        if rangeP == 0
            normRescaleP = ones(length(tempPval),1);
        else
            tempPval = (tempPval - minP) / rangeP;

            % Then scale to [x,y]:
            scaleMaxP = 6;
            scaleMinP = .1; 
            range2 = scaleMaxP - scaleMinP;
            normRescaleP = (tempPval*range2) + scaleMinP;
        end
        numSyll = size(syllIdx,1);
        for n = 1:numSyll
            if tempRho(n) > 0
                replayColor = 'g';
            else
                replayColor = 'r';
            end
            plot(ax3, [tempReplay(syllIdx(n,1),1) tempReplay(syllIdx(n,2),1)],...
            [abs(tempRho(n)) abs(tempRho(n))], replayColor, 'LineWidth', normRescaleP(n))
        end
%         maxNumSyll = [maxNumSyll; numSyll];
    end

    hold off
%     ylim(ax3,[0.9 (max(maxNumSyll)+0.1)])
    ylim(ax3,[0.5 1])
    xlim(ax3,[startTime stopTime])
    ylabel(ax3, 'Replay (Spearman rho)')           


    ax4 = subplot(4,1,4);
    targetIdx = TimeStamps >= startTime & TimeStamps <= stopTime;

    plot(ax4, TimeStamps(targetIdx), ctxSpindleSignal(targetIdx), 'b')
    hold on
    plot(ax4, TimeStamps(targetIdx), ctxSignal(targetIdx), 'Color', [0.5 0.5 0.5])
    targetIdx = TimeStamps >= eventStartTS &...
        TimeStamps <= eventStopTS;

    plot(ax4, TimeStamps(targetIdx),...
        ctxSpindleSignal(targetIdx), 'g')
    hold off
    xlim([startTime stopTime])
    ylim([-maxY_CtxSpindle maxY_CtxSpindle]);
    ylabel(ax4, 'Ctx EEG')
    xlabel(ax4, 'Time(seconds)');
    pause
    
%     end
end