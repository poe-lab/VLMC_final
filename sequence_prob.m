function [prob, P0, P1, stat_dist] = sequence_prob(P, seq, P1, P0)
% prob = sequence_prob(P, seq, P1, P0)
% Compute the probability of a sequence under a 2nd order markov model.

%%

if size(seq, 2) == 1
    seq = seq';
end

%% Build digrams
usymb = 0:(sqrt(size(P, 1)) - 1);
nsymb = length(usymb);

% digrams = [];
% for i = 1:nsymb
%     for j = 1:nsymb
%         digrams = [digrams; [usymb(i), usymb(j)]];
%     end
% end

digrams = list_kgrams((1:nsymb) - 1, 2);

%% Construct corresponding first order model

stat_dist = markov_stat_dist(P);

if nargin < 3
    
%     F = [];
%     for i = 1:size(P, 1)
%         F = [F; P(i, P(i, :) ~= 0)];
%     end
    F = markov2fut_dist(P);

    P1 = zeros(nsymb);
    for i = 1:nsymb
        
        idx = find(digrams(:, 2) == usymb(i));
        
        P1(i, :) = sum(diag(stat_dist(idx)) * F(idx, :));
    end
    
end
P1 = diag(1 ./ sum(P1, 2)) * P1;


% Construct corresponding zeroth order model
if nargin < 4
    P0 = markov_stat_dist(P1);
end

%%

% prob = P0(usymb == seq(1));
% 
% prob = prob * P1(usymb == seq(1), usymb == seq(2));

prob = stat_dist(ismember(digrams, seq(1:2), 'rows'));

for i = 3:length(seq)
    
    % Source digram
%     src_idx =...
%         all(digrams == repmat([seq(i - 2), seq(i - 1)],...
%         [size(digrams, 1), 1]), 2);

    src_idx =...
        ismember(digrams, seq((i - 2):(i - 1)), 'rows');
    
    % Target digram
%     tar_idx =...
%         all(digrams == repmat([seq(i - 1), seq(i)],...
%         [size(digrams, 1), 1]), 2);
    
    tar_idx =...
        ismember(digrams, seq((i - 1):i), 'rows');
    
    % Update probability
    prob = prob * P(src_idx, tar_idx);
    
end









