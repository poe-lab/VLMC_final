function null_P = generate_N_MarkovNulls_02252018(W, nsamp)

%% Null P
null_P = cell(nsamp,1);
for isamp = 1:nsamp
    isamp
    [null_P{isamp}, ~] = scramble_sequence_markov_null_02262018(W);
end
