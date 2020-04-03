function ReviewSigSeq_withFieldSize_replayBins_03142018
% All LFPs must have same sampling rate and time stamps for function to
% work properly.
orderPC = PC_post(PC_Order);
orderPC_ID = PC_ID(PC_Order,:);
PC_cutoff= 0;
numCells = size(orderPC,1)-PC_cutoff;
trackLength = 297; % approximate circumference of the track

numSpindles = size(EschenkoSpindle.startIdx,1); %Determine number of detected spindles on the selected LFP/EEG channel
% Determine maximum y-axis values (mV) for the spindle filtered signals plus
% added 10%
maxY_HcSpindle = max(abs(HcSpindleSignal)) + 0.1*max(abs(HcSpindleSignal));
maxY_CtxSpindle = max(abs(smoothdata(ctxSpindleSignal))) + 0.1*max(abs(smoothdata(ctxSpindleSignal)));
timeWindow= 1.5; %3 second window
maxSigSyllables = max(cellfun('length', replay.syllableIdxAll(logical(replay.inSpindle))));
mazeCodeNew = mazeCode+2;
mazeMarkColor = {'g';'none';'k'};
mazeMarker = {'.';'none';'.'};


%% Calculate and plot results based on each detected spindle centered in a 3 sec window:
for i = 1:numSpindles
    if ~isempty(spindleWordsReplay.idx{i})
        addTime = timeWindow - EschenkoSpindle.duration(i);
        if addTime < 0   % If the spindle is the max allowed of 3 sec duration,
            addTime = 0; %do not add any data points to the beginning and end.
        end
        startTime = EschenkoSpindle.timestamp(i) - addTime/2;
        stopTime = EschenkoSpindle.timestamp(i) + EschenkoSpindle.duration(i) + addTime/2;

        spikesInSpindle = [];
        for c = 1:numCells                                            % For each place cell
            spikes = orderPC{c};                                      % Keep all its spikes
            spikes = spikes(spikes >= startTime & spikes <= stopTime); % Keep all its spikes that fall within the target interval

            if ~isempty(spikes) % If there are any spikes
                spikesInSpindle = [spikesInSpindle; [spikes c*ones(size(spikes,1),1)]];
            end   
        end

        ax1 = subplot(4,1,1);
        targetIdx = TimeStamps >= startTime & TimeStamps <= stopTime;
        plot(ax1, TimeStamps(targetIdx), HcSpindleSignal(targetIdx), 'b')
        hold on

        plot(ax1, TimeStamps(targetIdx), HcSignal(targetIdx),...
            'Color', [0.5 0.5 0.5])

        targetIdx = TimeStamps >= EschenkoSpindle.timestamp(i) &...
            TimeStamps <= (EschenkoSpindle.timestamp(i)+EschenkoSpindle.duration(i));
        plot(ax1, TimeStamps(targetIdx),...
            HcSpindleSignal(targetIdx), 'g')

        targetIdx = TimeStamps >= startTime & TimeStamps <= stopTime;
        plot(ax1, TimeStamps(targetIdx), rippleSignal(targetIdx), 'r')
        hold off

        switch EschenkoSpindle.scoring(i)
            case 2
                stage = 'SWS';
            case 3
                stage = 'REM';
            case 4
                stage = 'QW';
            case 6
                stage = 'TR ';
        end
        title(['Spindle ' num2str(i) ' of ' num2str(numSpindles) '          Stage: ' stage])
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
        for numFields = 1:size(fieldLocations{spikesInSpindle(1,2),1},1)
            plot(ax2, [spikesInSpindle(1,1) spikesInSpindle(1,1)],...
                [fieldLocations{spikesInSpindle(1,2),1}(numFields,1)...
                fieldLocations{spikesInSpindle(1,2),1}(numFields,2)],...
                'Color', [0.5 0.5 0.5], 'LineWidth', 0.1, 'LineStyle', '-')
    %                     'Color', colorUnit, 'LineWidth', 0.25)
        end
        hold on
        for m = 2:size(spikesInSpindle,1)
    %                 colorUnit = map(mod(spikesInSpindle(m,2),7)+1,:);             
            for numFields = 1:size(fieldLocations{spikesInSpindle(m,2),1},1)
                plot(ax2, [spikesInSpindle(m,1) spikesInSpindle(m,1)],...
                    [fieldLocations{spikesInSpindle(m,2),1}(numFields,1)...
                    fieldLocations{spikesInSpindle(m,2),1}(numFields,2)],...
                    'Color', [0.5 0.5 0.5], 'LineWidth', 0.1, 'LineStyle', '-')
    %                     'Color', colorUnit, 'LineWidth', 0.25)
            end
        end

        % for each word in a spindle
        for j = 1:length(spindleWordsReplay.idx{i})
            wordIdx = spindleWordsReplay.idx{i}(j);
            plot(ax2, [replay.wordsTimeBins(wordIdx,1)...
                replay.wordsTimeBins(wordIdx,1)],...
                [0 trackLength+1], ':k')

            plot(ax2, [replay.wordsTimeBins(wordIdx,2)...
                replay.wordsTimeBins(wordIdx,2)],...
                [0 trackLength+1], ':k')

            %% Get cleaned spikes for  replay bins in the word:
            spikeBurstsInWord = [];
            for k = 1:size(replay.spikeSeq{wordIdx,1},1)
                spikeBurstsInWord =[spikeBurstsInWord ; replay.spikeSeq{wordIdx,1}{k,1}];
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
        end

         %% Find any LC spikes in window and plot them:
        numLC = size(lcUnits,1);
        if ~isempty(lcUnits)
            spindleLC_spikes = [];
            for iLC = 1:numLC                                                         % For each place cell
                spikesLC = lcUnits{iLC};                                      % Keep all its spikes
                spikesLC = spikesLC(spikesLC >= startTime & spikesLC <= stopTime); % Keep all its spikes that fall within the target interval
                if ~isempty(spikesLC)                                                 % If there are any spikes
                    spindleLC_spikes = [spindleLC_spikes; [spikesLC (trackLength+iLC)*ones(size(spikesLC,1),1)]]; 
                end
            end

            if ~isempty(spindleLC_spikes)
                scatter(ax2,spindleLC_spikes(:,1),spindleLC_spikes(:,2), 'r', '*')
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
        for j = 1:length(spindleWordsReplay.idx{i})
            wordIdx = spindleWordsReplay.idx{i}(j);
            %% Get cleaned spikes for  replay bins in the word:
            for k = 1:size(replay.spikeSeq{wordIdx,1},1)
                % Get all spike bursts for replay bin:
                tempReplay = replay.spikeSeq{wordIdx,1}{k,1};
                % Get start and stop indices of syllables:
                syllIdx = replay.syllableIdxAll{wordIdx,1}{k,1};

                tempRho = replay.rho{wordIdx,1}{k,1};
                tempPval = replay.pVal{wordIdx,1}{k,1};
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
                    [n n], replayColor, 'LineWidth', normRescaleP(n))
                end
                maxNumSyll = [maxNumSyll; numSyll];
            end
               
        end
        hold off
        ylim(ax3,[0.9 (max(maxNumSyll)+0.1)])
        xlim(ax3,[startTime stopTime])
        ylabel(ax3, 'Replay')           


        ax4 = subplot(4,1,4);
        targetIdx = TimeStamps >= startTime & TimeStamps <= stopTime;

        plot(ax4, TimeStamps(targetIdx), ctxSpindleSignal(targetIdx), 'b')
        hold on
        plot(ax4, TimeStamps(targetIdx), ctxSignal(targetIdx), 'Color', [0.5 0.5 0.5])
        targetIdx = TimeStamps >= EschenkoSpindle.timestamp(i) &...
            TimeStamps <= (EschenkoSpindle.timestamp(i)+EschenkoSpindle.duration(i));

        plot(ax4, TimeStamps(targetIdx),...
            ctxSpindleSignal(targetIdx), 'g')
        hold off
        xlim([startTime stopTime])
        ylim([-maxY_CtxSpindle maxY_CtxSpindle]);
        ylabel(ax4, 'Ctx EEG')
        xlabel(ax4, 'Time(seconds)');
        pause
    
    end
end