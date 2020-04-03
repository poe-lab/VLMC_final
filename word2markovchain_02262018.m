function [M, vmm_model] = word2markovchain_02262018(W, nfold, d, alg)
% M = word2markovchain(W, nfold, d, alg)
% Compute Markov chain associated to sequence by cross-validation and
% output its transition matrix and CV data.

if (nargin < 2) || isempty(nfold)
    nfold = 10;
end

%% 

W = W + 1;

% Convert to character sequence
alph = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o',...
    'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z'};
seq = cell2mat(alph(W));

% Compute folds
L = length(seq);
fold_size = floor(L / nfold);
fold = cell(nfold, 1);
for i = 1:nfold
    fold{i} = (i - 1) * fold_size + (1:fold_size);
end

% IID model
cv_iid_entropy = zeros(length(fold), 1);
uword = unique(W);
for i = 1:length(fold)
    curr_W = W(setdiff(1:length(W), fold{i}));
    freq = histc(curr_W, uword) / length(curr_W);
    %cv_iid_entropy(i) = -freq * log2(freq');
    cv_iid_entropy(i) = mean(-log2(freq(...
        integer_label(W(fold{i}), uword))));
end

% Markov models
if (nargin < 3) || isempty(d)
    d = 1:5;
end

if (nargin < 4) || isempty(alg)
    alg = 'PPMC';
end

[vmm_model, cv_bits_per_symb, train_bits_per_symb] = cv_vmm(seq, d, nfold, alg);

%% Convert to transition matrix

ab = alphabet(seq);
params.d = d(find(mean(cv_bits_per_symb, 2,'omitnan') == min(mean(cv_bits_per_symb, 2,'omitnan'))));
params.ab_size = size(ab);
params.pMin = 0.006;
params.alpha= 0;
params.gamma = 1e-3;
params.r = 1.05;

num_state = length(unique(W))^params.d;
if num_state < 1e4
    opt_model = vmm_create(map(ab, seq), alg, params);
    P = vmm_transition_matrix(seq, opt_model, params.d);
else
    disp('State space is too large to fill Markov chain matrix.');
    opt_model = [];
    P = [];
end

%% Output data

M = [];
M.P = P;
M.cv_bits_per_symb = cv_bits_per_symb;
M.train_bits_per_symb = train_bits_per_symb;
 




