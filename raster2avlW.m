function [words, W] = raster2avlW(raster, max_isi)
% [words, W] = raster2avlW(raster, max_isi)
% Convert a raster into an avalanche word.
% Time stamps in raster and max_isi should be in same units.
% Default max_isi = 0.050 (seconds)

all_spike = [];
for i = 1:length(raster)
    all_spike = [all_spike; raster{i}];
end
all_spike = sort(all_spike);
isi = diff(all_spike);

bin_size = mean(isi);

[X, min_time] = bin_raster(raster, bin_size);

avl_time = get_index_blocks(find(sum(X) ~= 0));

offset = min_time - bin_size;

% words = cell(size(avl_time));
% tic;
% h = waitbar(0, 'Computing avalanche words...');
% for i = 1:length(avl_time)
%     waitbar(i / length(avl_time));
%     start_time = bin_size * (avl_time{i}(1) - 1) + offset;
%     end_time = bin_size * avl_time{i}(end) + offset;
%     
%     c_raster = cut_raster(raster, start_time, end_time);
%     
%     words{i} = raster2word_avl(c_raster, max_isi);
%     
% end
% close(h);
% timer = toc,

% Start and end times as cell array
start_time = cellfun(@(x) bin_size * (x(1) - 1) + offset, avl_time,...
    'UniformOutput', false);
end_time = cellfun(@(x) bin_size * x(end) + offset, avl_time,...
    'UniformOutput', false);

% Cut rasters
cut_fun = @(x, y) cut_raster(raster, x, y);

c_raster = cellfun(cut_fun, start_time, end_time,...
    'UniformOutput', false);

% Convert to words
word_fun = @(x) raster2word_avl(x, max_isi);
words = cellfun(word_fun, c_raster,...
    'UniformOutput', false);

%%

if nargout > 1
    
    sil_time = get_index_blocks(find(sum(X) == 0));
    
    start_left = 0;
    if sil_time{1}(1) == 1
        left_buffer = zeros(1, length(sil_time{1}));
        Z = cellfun(@(x) zeros(1, length(x)), sil_time(2:end), 'UniformOutput', 0);
    else
        left_buffer = [];
        Z = cellfun(@(x) zeros(1, length(x)), sil_time(1:end), 'UniformOutput', 0);
    end
    
    W = [words', Z']';
    W = W(:)';
    W = cell2mat([left_buffer, W]);
    
end

end

%% Dependent functions
function W = raster2word_avl(raster, max_isi)
% W = raster2word_avl(raster, max_isi)
% Parse small rasters into words. For use within a spike avalanche.

W = [];
burst_time = [];

for i = 1:length(raster)
    
    spike = raster{i};
    
    if length(spike) >= 1
        isi = diff(spike);
        keep_idx = find([1; isi > max_isi] > 0);
        
        burst_time = [burst_time; spike(keep_idx)];
        
        num_burst = length(keep_idx);
        W = [W, i * ones(1, num_burst)];
    end
end

[~, burst_ord] = sort(burst_time);
W = W(burst_ord);

end







