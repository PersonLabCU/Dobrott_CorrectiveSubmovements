function [mean_psth, ci] = BootstrapPSTH(raster, n_boot)
    n_trials = size(raster, 1);
    n_timebins = size(raster, 2);
    
    bin_size_sec = 0.01;  % 10 ms bins (converted to seconds)
    smooth_win = 5;       % smoothing window in bins (increased window size)

    bootstrapped_psths = zeros(n_boot, n_timebins);

    for i = 1:n_boot
        % Bootstrap resampling
        resample_idx = randi(n_trials, [n_trials, 1]);
        resample_raster = raster(resample_idx, :);
        
        % Compute mean for each bootstrap sample
        bootstrapped_psths(i, :) = mean(resample_raster, 1);
    end

    % Convert bootstrapped PSTHs to firing rate (in Hz)
    bootstrapped_psths = bootstrapped_psths / bin_size_sec;

    % Compute mean and confidence intervals (CI)
    mean_psth = mean(bootstrapped_psths, 1);
    ci = prctile(bootstrapped_psths, [2.5, 97.5]);

    % Smoothing the mean PSTH and CI
    mean_psth = movmean(mean_psth, smooth_win);
    ci(1,:) = movmean(ci(1,:), smooth_win);
    ci(2,:) = movmean(ci(2,:), smooth_win);
end
