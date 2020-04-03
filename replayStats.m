% Convert cell array of cell arrays to a vector
function replayStats


pVal = [];
rho = [];
startBinTS = [];
replaySeqTS = [];
for i = 1:length(replay.rho)
    for j = 1:length(replay.rho{i,1})
        [minP, idxMinP] = min(replay.pVal{i,1}{j,1});
        pVal = [pVal; minP];
        rho = [rho; replay.rho{i,1}{j,1}(idxMinP)];
        startBinTS =[startBinTS; replay.spikeSeq{i,1}{j,1}(1)];
        replaySeqTS = [replaySeqTS; [replay.spikeSeq{i,1}{j,1}(1)
    end    
end

revTS = startBinTS(rho<0);
fwdTS = startBinTS(rho>0);
