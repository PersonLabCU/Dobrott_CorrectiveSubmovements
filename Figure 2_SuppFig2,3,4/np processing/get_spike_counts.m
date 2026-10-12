function [Trial_spike_times, spike_counts, time_vector] = get_spike_counts(spike_times, event_times, baseline_epoch, event_epoch, bin_size)
    % Convert event times to ms
    event_times = event_times * 1000;
    spike_times = spike_times * 1000;

    % Parameters
    time_vector = -baseline_epoch:bin_size:event_epoch - bin_size;
    num_bins = length(time_vector);
    num_trials = length(event_times);

    % Pre-allocate
    spike_counts = zeros(num_trials, num_bins);
    Trial_spike_times = struct();

    % Loop over trials
    for t = 1:num_trials
        ref_time = event_times(t);
        trial_spikes = spike_times(spike_times >= ref_time - baseline_epoch & spike_times <= ref_time + event_epoch) - ref_time;
        Trial_spike_times.(sprintf('R%d', t)) = trial_spikes;

        for b = 1:num_bins
            bin_start = time_vector(b);
            bin_end = bin_start + bin_size;
            spike_counts(t, b) = sum(trial_spikes >= bin_start & trial_spikes < bin_end);
        end
    end
end