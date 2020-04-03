function fo_words = firing_order_words(seq)
% fo_words = firing_order_words(seq)
% Compute firing order_words with gaps.

L = length(seq);

fo_words = cell(L, 1);
for i = 1:(2^L - 1)
    curr_subset = dec2bin(i, L) == '1';
    word_size = sum(curr_subset);
    fo_words{word_size} = [fo_words{word_size}; seq(curr_subset)];
end
fo_words = fo_words(2:end);

fo_words = cellfun(...
    @(x) mat2cell(x, ones(size(x, 1), 1), size(x, 2)), fo_words, 'UniformOutput', 0);
temp_fo_words = [];
for i = 1:length(fo_words)
    temp_fo_words = [temp_fo_words; fo_words{i}];
end
fo_words = temp_fo_words;