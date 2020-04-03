% find the closest spindle for each word
closeSpindleIdx = [];
closeSpindleTime = [];

for i = 1:size(BH_FDR.wordsTimeBins,1)
    idx1 = find(EschenkoSpindle.timestamp(:) >=  BH_FDR.wordsTimeBins(i,1), 1);
    idx2 = find(EschenkoSpindle.timestamp(:) <=  BH_FDR.wordsTimeBins(i,1));
    if ~isempty(idx2)
        idx2 = idx2(end);
    end
    idxBoth = [idx1 idx2];
    [spindleTime, idx3] = min(abs(EschenkoSpindle.timestamp(idxBoth) - BH_FDR.wordsTimeBins(i,1)));
    spindleIdx = idxBoth(idx3);
    closeSpindleTime = [closeSpindleTime; spindleTime];
    closeSpindleIdx = [closeSpindleIdx; spindleIdx];
end