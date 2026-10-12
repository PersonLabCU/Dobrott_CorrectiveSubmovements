    % Calculate confidence intervals for bootstrapped data
function ci = bootstrap_ci(psth_data, n_samples, alpha)
    % Computes confidence intervals from bootstrap resampling
    n_timepoints = length(psth_data);
    boot_means = zeros(n_samples, n_timepoints);

    for i = 1:n_samples
        boot_sample = datasample(psth_data, n_timepoints, 'Replace', true);
        boot_means(i, :) = boot_sample;
    end

    % Compute percentile-based confidence intervals
    ci = prctile(boot_means, [(100 * alpha / 2), (100 * (1 - alpha / 2))], 1);
end