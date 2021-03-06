function replayStats = replayStats_11022018(replay)
replayStats.pVal = [];
replayStats.rho = [];
startTS = [];
stopTS = [];
replayStats.replaySeqTS = [];
for i = 1:length(replay.rho)
    for j = 1:length(replay.rho{i,1})
        [minP, idxMinP] = min(replay.pVal{i,1}{j,1});
        replayStats.pVal = [replayStats.pVal; minP];
        replayStats.rho = [replayStats.rho; replay.rho{i,1}{j,1}(idxMinP)];
        startIdx = replay.syllableIdxAll{i,1}{j,1}(idxMinP,1);
        startTS = replay.spikeSeq{i,1}{j,1}(startIdx,1);
        stopIdx = replay.syllableIdxAll{i,1}{j,1}(idxMinP,2);
        stopTS = replay.spikeSeq{i,1}{j,1}(stopIdx,1);
        replayStats.replaySeqTS = [replayStats.replaySeqTS; [startTS stopTS]];
    end    
end

replayStats.revTS = replayStats.replaySeqTS(replayStats.rho<0,:);
replayStats.fwdTS = replayStats.replaySeqTS(replayStats.rho>0,:);
