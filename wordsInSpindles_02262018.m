function [inSpindle, spindleWords] = wordsInSpindles_02262018(wordsTimeBins, spindle)

%% Determine if words occur during spindles:
numWords = size(wordsTimeBins,1);
inSpindle = zeros(numWords,1);
numSpindles = size(spindle.timestamp, 1);
spindleWords = [];
spindleWords.idx = cell(numSpindles,1);

for i = 1:numSpindles
    idx1 = find(wordsTimeBins(:,1) >= spindle.timestamp(i) & wordsTimeBins(:,1) <= spindle.timestamp(i) + spindle.duration(i));
    idx2 = find(wordsTimeBins(:,2) >= spindle.timestamp(i) & wordsTimeBins(:,2) <= spindle.timestamp(i) + spindle.duration(i));
    idxBoth = [idx1;idx2];
    clear idx1 idx2
    if ~isempty(idxBoth)
        targIdx = sort(unique(idxBoth));
        spindleWords.idx{i} = targIdx';
        inSpindle(targIdx) = 1;
        clear targIdx
    end
    clear idxBoth
end
spindleWords.length = cellfun('length',spindleWords.idx);
spindleWords.ratio = sum(spindleWords.length>0)/numSpindles;