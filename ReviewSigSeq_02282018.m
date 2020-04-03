function ReviewSigSeq_02282018

orderPC = PC_post(PC_Order);
orderPC_ID = PC_ID(PC_Order,:);
PC_cutoff= 0;
numCells = size(orderPC,1)-PC_cutoff;
numSpindles = size(EschenkoSpindle.startIdx,1); %Determine number of detected spindles on the selected LFP/EEG channel
% Determine maximum y-axis values (mV) for the spindle filtered signals plus
% added 10%
maxY_HcSpindle = max(abs(HcSpindleSignal)) + 0.1*max(abs(HcSpindleSignal));
maxY_CtxSpindle = max(abs(smoothdata(ctxSpindleSignal))) + 0.1*max(abs(smoothdata(ctxSpindleSignal)));
timeWindow= 3; %3 second window
%% Calculate and plot results based on each detected spindle centered in a 3 sec window:
for i = 1:numSpindles
    if ~isempty(spindleWords.idx{i})
    addTime = timeWindow - EschenkoSpindle.duration(i);
    if addTime < 0   % If the spindle is the max allowed of 3 sec duration,
        addTime = 0; %do not add any data points to the beginning and end.
    end
    startTime = EschenkoSpindle.timestamp(i) - addTime/2;
    stopTime = EschenkoSpindle.timestamp(i) + EschenkoSpindle.duration(i) + addTime/2;

    spikesInSpindle = [];
    
    for c = 1:numCells                                                         % For each place cell
        spikes = orderPC{c};                                      % Keep all its spikes
        spikes = spikes(spikes >= startTime & spikes <= stopTime); % Keep all its spikes that fall within the target interval
        if ~isempty(spikes)                                                 % If there are any spikes
            spikesInSpindle = [spikesInSpindle; [spikes c*ones(size(spikes,1),1)]]; 
        end
    end

    if ~isempty(spikesInSpindle)
        lcSpikes = lcSpikeData/1000000;
        lcSpikes = lcSpikes(lcSpikes >= startTime & lcSpikes <= stopTime);
        lcSpikesInSpindle = [];
        if ~isempty(lcSpikes)
            lcSpikesInSpindle = [lcSpikesInSpindle; [lcSpikes (numCells+1)*ones(size(lcSpikes,1),1)]];
        end
        ax1 = subplot(3,1,1);
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
        
        ax2 = subplot(3,1,2);
        map=lines(7);
        colorUnit = map(mod(spikesInSpindle(1,2),7)+1,:);
        plot(ax2, [spikesInSpindle(1,1) spikesInSpindle(1,1)], [spikesInSpindle(1,2)-.5 spikesInSpindle(1,2)+.5], 'Color', colorUnit)
        hold on
        for m = 2:size(spikesInSpindle,1)
            colorUnit = map(mod(spikesInSpindle(m,2),7)+1,:);             %de2bi(mod(spikesInSpindle(m,2),7),3);
            plot(ax2, [spikesInSpindle(m,1) spikesInSpindle(m,1)], [spikesInSpindle(m,2)-.5 spikesInSpindle(m,2)+.5], 'Color', colorUnit)
        end
%         scatter(ax2,spikesInSpindle(:,1),spikesInSpindle(:,2), 'b', 'd')
%         hold on
        for j = 1:length(spindleWords.idx{i})
            plot(ax2, [BH_FDR.wordsTimeBins(spindleWords.idx{i}(j),1)...
                BH_FDR.wordsTimeBins(spindleWords.idx{i}(j),1)],...
                [0 numCells+1], ':g')
            plot(ax2, [BH_FDR.wordsTimeBins(spindleWords.idx{i}(j),2)...
                BH_FDR.wordsTimeBins(spindleWords.idx{i}(j),2)],...
                [0 numCells+1], ':gk')    
        end

        if ~isempty(lcSpikesInSpindle)
            scatter(ax2,lcSpikesInSpindle(:,1),lcSpikesInSpindle(:,2), 'r', 's')
        end
        hold off
        ylim([0 numCells+1])
        xlim([startTime stopTime])
        ylabel(ax2, 'Place Cells')
%         ax2.TickDirMode = 'manual';
%         ax2.TickDir = 'out';
        
        ax3 = subplot(3,1,3);
        targetIdx = TimeStamps >= startTime & TimeStamps <= stopTime;
        
        plot(ax3, TimeStamps(targetIdx), ctxSpindleSignal(targetIdx), 'b')
        hold on
        plot(ax3, TimeStamps(targetIdx), ctxSignal(targetIdx), 'Color', [0.5 0.5 0.5])
        targetIdx = TimeStamps >= EschenkoSpindle.timestamp(i) &...
            TimeStamps <= (EschenkoSpindle.timestamp(i)+EschenkoSpindle.duration(i));
        
        plot(ax3, TimeStamps(targetIdx),...
            ctxSpindleSignal(targetIdx), 'g')
        hold off
        xlim([startTime stopTime])
        ylim([-maxY_CtxSpindle maxY_CtxSpindle]);
        ylabel(ax3, 'Ctx EEG')
        xlabel(ax3, 'Time(seconds)');
        pause
    end
    end
end