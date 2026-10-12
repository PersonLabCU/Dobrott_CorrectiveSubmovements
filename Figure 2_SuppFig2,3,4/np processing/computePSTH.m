function [psth_raw, psth_smooth, time_bins] = computePSTH(BinaryMaster, bin_size, smoothing_bins)
    % Inputs:
    % BinaryMaster: E x T binary matrix (E = events/trials, T = time in ms)
    % bin_size: bin size in ms (e.g., 10)
    % smoothing_bins: # of bins for smoothing (e.g., 3)

    [num_events, total_time] = size(BinaryMaster);

    % Number of new bins after downsampling
    num_bins = floor(total_time / bin_size);

    % Preallocate binned BinaryMaster
    BinaryBinned = zeros(num_events, num_bins);

    % Aggregate data into bin_size chunks
    for b = 1:num_bins
        idx_start = (b - 1) * bin_size + 1;
        idx_end = b * bin_size;
        BinaryBinned(:, b) = sum(BinaryMaster(:, idx_start:idx_end), 2);
    end

    % Compute mean firing rate per bin (Hz)
    psth_raw = mean(BinaryBinned, 1) / bin_size * 1000;

    % Smooth using moving average
   % psth_smooth = movmean(psth_raw, smoothing_bins, 'Endpoints', 'discard');
    psth_smooth = movmean(psth_raw, smoothing_bins); % keeps full length 

    % Create time bins centered around 0
    % time_bins = ((1:num_bins) - num_bins / 2) * bin_size;
    time_bins = ((1:num_bins) - num_bins / 2) * bin_size + bin_size/2; % centering the bin 
end
