function probWord = sequence_prob_BGmod(P, words)
% prob = sequence_prob(P, seq, P1, P0)
% Compute the probability of a sequence under a 2nd order markov model.

%% Build digrams
usymb = 0:(sqrt(size(P, 1)) - 1);
nsymb = length(usymb);
digrams = list_kgrams((1:nsymb) - 1, 2);

stat_dist = markov_stat_dist(P);

numWords = length(words);
probWord = zeros(numWords,1);

for iWords = 1:numWords
    seq = words{iWords};
    if size(seq, 2) == 1
        seq = seq';
    end
    
    prob = stat_dist(ismember(digrams, seq(1:2), 'rows'));

    for iSeq = 3:length(seq)    
        % Source digram
        src_idx =...
            ismember(digrams, seq((iSeq - 2):(iSeq - 1)), 'rows');

        % Target digram   
        tar_idx =...
            ismember(digrams, seq((iSeq - 1):iSeq), 'rows');

        % Update probability
        prob = prob * P(src_idx, tar_idx);    
    end
    probWord(iWords) = prob;
end


















