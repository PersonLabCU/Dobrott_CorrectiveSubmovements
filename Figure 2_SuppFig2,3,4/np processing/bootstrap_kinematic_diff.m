function [mean_diff, ci_smooth] = bootstrap_kinematic_diff(data1, data2, sigma, num_boot)
    % Truncate to the smaller number of trials
    min_trials = min(size(data1, 2), size(data2, 2));
    data1 = data1(:, 1:min_trials);
    data2 = data2(:, 1:min_trials);

    % Compute actual mean difference
    mean_diff = mean(data1 - data2, 2);

    % Bootstrap
    n_timepoints = size(data1, 1);
    boot_diffs = zeros(n_timepoints, num_boot);
    for b = 1:num_boot
        idx = randi(min_trials, [1, min_trials]);
        boot_diffs(:, b) = mean(data1(:, idx) - data2(:, idx), 2);
    end

    % Compute 95% confidence intervals
    ci_low = prctile(boot_diffs, 2.5, 2);
    ci_high = prctile(boot_diffs, 97.5, 2);
    ci_smooth = [smoothdata(ci_low, 'gaussian', sigma), smoothdata(ci_high, 'gaussian', sigma)];
end
