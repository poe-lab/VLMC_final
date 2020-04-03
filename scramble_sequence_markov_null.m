function [null_P, M] = scramble_sequence_markov_null(W, keep_symb)
% [null_P, M] = scramble_sequence_markov_null(W, keep_symb)
% Build null Markov chain by scrambling sequence.

% zero_idx = find(W == 0);
% 
% null_W = W;
% curr_idx = 1;
% next_idx = 0;
% while next_idx < length(W)
%     
%     %curr_idx / length(W),
%     
%     if W(curr_idx) == 0
%         curr_idx = curr_idx + 1;
%     else
%         next_idx = zero_idx(find(zero_idx > curr_idx, 1, 'first'));
%         stretch = curr_idx:(next_idx - 1);
%         null_W(stretch) = W(stretch(randperm(length(stretch))));
%         curr_idx = next_idx;
%     end
% end

%%
nz_idx = find(W ~= 0);
null_W = W;
null_W(nz_idx) = W(nz_idx(randperm(length(nz_idx))));

if nargin == 2
    null_W = null_W(ismember(null_W, keep_symb));
end

%%
d = 2 %0:2;
nfold = 10;
alg = 'PST';

M = word2markovchain(null_W, nfold, d, alg);

null_P = M.P;

if size(null_P, 1) == length(unique(W))
    null_P = fut_dist2markov(repmat(null_P, [size(null_P, 1), 1]));
end