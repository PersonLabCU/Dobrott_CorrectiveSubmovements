
%% Bootstrapping Function
function [boot_mean, ci] = bootstrap_psth(psth_raw, n_bootstrap)
    if isempty(psth_raw)
        error("PSTH data is empty.");
    end

    nTrials = size(psth_raw, 1); % Number of trials
    nTime = size(psth_raw, 2);   % Time bins

    boot_means = zeros(n_bootstrap, nTime); % Preallocate

    for i = 1:n_bootstrap
        sample_indices = randi(nTrials, nTrials, 1); % Sample with replacement
        boot_means(i, :) = mean(psth_raw(sample_indices, :), 1); % Mean across trials, keep time intact
    end

    % Final bootstrapped mean
    boot_mean = mean(boot_means, 1);

    % Compute confidence intervals
    ci = prctile(boot_means, [2.5, 97.5], 1);
end