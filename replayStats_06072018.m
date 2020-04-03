function replayStats_06072018
pVal = [];
rho = [];
startTS = [];
stopTS = [];
replaySeqTS = [];
for i = 1:length(replay.rho)
    for j = 1:length(replay.rho{i,1})
        [minP, idxMinP] = min(replay.pVal{i,1}{j,1});
        pVal = [pVal; minP];
        rho = [rho; replay.rho{i,1}{j,1}(idxMinP)];
        startIdx = replay.syllableIdxAll{i,1}{j,1}(idxMinP,1);
        startTS = replay.spikeSeq{i,1}{j,1}(startIdx,1);
        stopIdx = replay.syllableIdxAll{i,1}{j,1}(idxMinP,2);
        stopTS = replay.spikeSeq{i,1}{j,1}(stopIdx,1);
        replaySeqTS = [replaySeqTS; [startTS stopTS]];
    end    
end

revTS = replaySeqTS(rho<0,:);
fwdTS = replaySeqTS(rho>0,:);
