%% Condition data for VLMC using 'spikeDataForVLMC_########.m':

% Enter PC IDs:
placeCells = [];

% Enter PC order from place field analysis:
PC_Order = [];

% Enter sleep bound (1st sleep epoch to last sleep epoch in scored file):
sleepBounds = [];

% Run 'spikeDataForVLMC_06052018.m':
spikeDataForVLMC_06052018(placeCells, PC_Order, sleepBounds);

% Load above saved data from spikeDataForVLMC_06052018.m

% Adjust cell order due to PCs that were not in the post-run recording:
% Order vector formed in Excel by labeling numbers 1 to n on left and then
% sorting rows by center of mass (cm).
PC_Order_RunAndPost = PC_Order;

%% Re-order place cells in POST by position on track:
PC_ID=PC_ID(PC_Order_RunAndPost,:);
PC_post=PC_post(PC_Order_RunAndPost);

mazeCode = [];
centerOfMass = [];
fieldSize = [];
%% Run field boundary calculations
fieldLocations = calculateFieldBoundaries_03052018(centerOfMass, fieldSize);

%% Build words from the spike data based on the max ISI setting:
max_isi = 0.05;
[words, W, wordsTimeBins, bin_size] = raster2avlW_02272018(PC_post, max_isi);

%% Find the state scored epoch in which each word begins:
[wordState,scoredFile] = assignStatesToWords_02252018(wordsTimeBins);

%% Generate real Markov chain
d = 2 %0:2;
nfold = 10;
alg = 'PST';
[M, vmm_model] = word2markovchain_02262018(W, nfold, d, alg)

%% Generate 100 null (shuffled) Markov chains
nsamp = 100;
null_P = generate_N_MarkovNulls_02252018(W, nsamp);

%% Get Markov chains
P = M.P;
clear M

% % Below was for original way I saved these variables:
% W = W_Post;
% words = words_Post;
% clear W_Post words_Post

%% Select words of target state(s):
targetLogic = wordState == 2 | wordState == 6;
words = words(targetLogic);
wordsTimeBins = wordsTimeBins(targetLogic,:);
wordState = wordState(targetLogic);
clear targetLogic

%% Remove words less than _ letters
minLetters = 5;
targetLogic = cellfun('length',words)>= minLetters;
words = words(targetLogic);
wordsTimeBins = wordsTimeBins(targetLogic,:);
wordState = wordState(targetLogic);
clear targetLogic

% alph = unique(W);
% digram = list_kgrams(alph, 2);

%% Compute word probabilities
word_prob = sequence_prob_BGmod(P, words);

%% Null probability distributions
null_dist = [];
parfor isamp = 1:length(null_P)    
    % Compute probabilities
    word_null_prob = sequence_prob_BGmod(null_P{isamp}, words);
    null_dist = [null_dist, word_null_prob];
end

%% Z-score
true_z = zeros(size(word_prob));
for i = 1:length(word_prob)
    mn = mean(null_dist(i, :));
    sd = std(null_dist(i, :));
    true_z(i) = (word_prob(i) - mn) / sd;
end
pval = 1 - normcdf(true_z);

%% Compute FDR - Matt's version
null_z = zscore(null_dist, [], 2);    
% figure;
% hist(mean(null_z))
% vline(mean(true_z), 'r');
%
% disp('Mean p-value')
% sum(mean(null_z) > mean(true_z)) / nsamp,
%
% figure;
% for i = 1:size(null_dist, 1)
%
%     subplot(5, 5, i)
%     hist(null_dist(i, :));
%     vline(fo_word_prob(i), 'r');
%
% end

%     figure;
%     hist(null_z(:), 20)
%     vline(true_z, 'r');

[sort_z, idx]  = sort(true_z, 'descend');

ave_fd = zeros(size(sort_z));
for i = 1:length(sort_z)
    ave_fd(i) = mean(sum(null_z >= sort_z(i)));
end

falseDiscovery.fdr = ave_fd ./ (1:length(ave_fd))';

%[sort_z, ave_fd, fdr],
targEnd = sum(falseDiscovery.fdr < 0.1);
falseDiscovery.words = words(idx(1:targEnd));
falseDiscovery.wordsTimeBins = wordsTimeBins(idx(1:targEnd),:);
falseDiscovery.wordState = wordState(idx(1:targEnd));

%% Compute False Discovery Rate - Benjamani-Hochberg Procedure
% The Benjamini–Hochberg procedure (BH step-up procedure) controls the FDR
% at level alpha. It works as follows:
% 1) For a given alpha, find the largest k such that P(k) <= (k/m)*alpha.
% 2) Reject the null hypothesis (i.e., declare discoveries) for all H(i) for i = 1,...,k.
alpha = 0.05;
numTests = length(pval); % m
pThreshold = ((1:1:numTests)/numTests * alpha);
[sort_pVal, idx_pVal]  = sort(pval);
targLogic = sort_pVal <= pThreshold';
targEnd = find(targLogic, 1, 'last');

BH_FDR.words = words(idx_pVal(1:targEnd));
BH_FDR.wordsTimeBins = wordsTimeBins(idx_pVal(1:targEnd),:);
BH_FDR.wordState = wordState(idx_pVal(1:targEnd));

% Sort by time
[~, idx_time]  = sort(BH_FDR.wordsTimeBins(:,1));
BH_FDR.words = BH_FDR.words(idx_time);
BH_FDR.wordsTimeBins = BH_FDR.wordsTimeBins(idx_time,:);
BH_FDR.wordState = BH_FDR.wordState(idx_time);
BH_FDR.wordLength = cellfun('length',BH_FDR.words);







