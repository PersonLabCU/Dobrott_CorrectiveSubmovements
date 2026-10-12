function [mean_psth, sem_psth, time_bins] = bootstrapPSTH_smooth(BinaryBinned, bin_size, smoothing_bins, nBoot)

    [nTrials, nBins] = size(BinaryBinned);

    % Compute a test smooth to get output length
    test_mean = mean(BinaryBinned, 1) / bin_size * 1000;
    test_smooth = movmean(test_mean, smoothing_bins, 'Endpoints', 'shrink');
    new_nBins = length(test_smooth);

    % Initialize bootstrap matrix with new_nBins
    boot_psths = zeros(nBoot, new_nBins);

    % Bootstrap loop
    for i = 1:nBoot
        % Sample trials with replacement
        boot_idx = randsample(nTrials, nTrials, true);
        boot_trials = BinaryBinned(boot_idx, :);

        % Compute mean firing rate per bin (Hz)
        boot_mean = mean(boot_trials, 1) / bin_size * 1000;

        % Smooth with moving average (shrink)
        boot_smooth = movmean(boot_mean, smoothing_bins, 'Endpoints', 'shrink');

        % Store
        boot_psths(i, :) = boot_smooth;
    end

    % Compute mean and SEM across bootstraps
    mean_psth = mean(boot_psths, 1);
    sem_psth = std(boot_psths, [], 1);

    % Adjust time bins
    time_bins = ((1:new_nBins) - nBins/2) * bin_size; % Approximate alignment

end
