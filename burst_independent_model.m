function [P0, smallP0, Z] = burst_independent_model(W, s1, Z)
% [P0, smallP0, Z] = burst_independent_model(W, s1, Z)
% Compute the "bursty independent model" that preserves burst word and
% silence periods, but has no temporal structure between bursts.
%
% The optional S1 is a HACK to insert an alternative stationary
% distribution than given by W, but retains the 0-0, x-0, 0-x transition
% frequencies.

% Get stationary distribution
alph = unique(W);

if nargin < 2
    s1 = histc(W, unique(W));
    s1 = s1 / sum(s1);
end

% Get transition rate into and out of silence
if nargin < 3
    Z = zeros(2);
    
    last_symb = W(1);
    for i = 2:length(W)
        curr_symb = W(i);
        
        if (last_symb == 0) && (curr_symb == 0)
            Z(1, 1) = Z(1, 1) + 1;
        elseif (last_symb == 0) && (curr_symb > 0)
            Z(1, 2) = Z(1, 2) + 1;
        elseif (last_symb > 0) && (curr_symb == 0)
            Z(2, 1) = Z(2, 1) + 1;
        else
            Z(2, 2) = Z(2, 2) + 1;
        end
        
        last_symb = curr_symb;
    end
    
    Z = diag(1 ./ sum(Z, 2)) * Z;
end

% Fill burst independent model
P0 = zeros(length(alph));

nz_dist = s1(2:end) / sum(s1(2:end));

P0(1, 1) = Z(1, 1);
P0(1, 2:end) = Z(1, 2) * nz_dist;
P0(2:end, 1) = Z(2, 1) * ones(size(P0, 1) - 1, 1);
P0(2:end, 2:end) = Z(2, 2) * repmat(nz_dist, [size(P0, 1) - 1, 1]);

smallP0 = P0;

P0 = fut_dist2markov(repmat(P0, [size(P0, 1), 1]));