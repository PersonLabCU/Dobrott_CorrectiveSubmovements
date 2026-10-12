function [ReReach, StimRR, MaxVel, ReachesWithoutMaxJerk] = CSMEarlyThresh_jerk_edit7(Reaches, thresh)
    % Function to detect the max jerk before the threshold crossing
    % Saves full metadata into structured ReReach output
    
    ReReach = struct();
    ReachesWithoutMaxJerk = struct();
    StimRR = [];
    MaxVel = [];

    reach_ind = fieldnames(Reaches);  % Get all reach indices

    for reachnum = 1:length(reach_ind)
        reachname = reach_ind{reachnum};
        reach = Reaches.(reachname).kin;

        % 1. Find the threshold crossing (position exceeding the threshold)
        stim_inx = find(reach(:, 2) >= thresh, 1);

        if isempty(stim_inx)
            continue;  % Skip if no threshold crossing is found
        end

        % 2. Define the window around the threshold crossing (0.1s before crossing)
        windowend = reach(stim_inx, 1);  % Time of threshold crossing
        windowstart = windowend - 0.10;

        % Find the closest time point to windowstart
        [~, index] = min(abs(reach(:, 1) - windowstart));

        % Clip the reach data to the time window
        accel_reach = reach(index:stim_inx, :);

        % 3. Extract jerk values and time
        jerk_col = accel_reach(:, 15);
        time_col = accel_reach(:, 1);

        % 4. Smooth jerk using moving average filter (window size 5)
        smoothed_jerk = smooth(jerk_col, 5);

        % 5. Find the maximum jerk within the window
        [max_jerk, max_jerk_idx] = max(smoothed_jerk);

        stim_time = reach(stim_inx, 1);  % Time of threshold crossing

        if isempty(max_jerk) || isnan(max_jerk)
            % No max jerk found, record this reach and append NaN to StimRR
            ReachesWithoutMaxJerk.(reachname).kin = reach;
            ReachesWithoutMaxJerk.(reachname).inx = Reaches.(reachname).inx;
            ReachesWithoutMaxJerk.(reachname).success = Reaches.(reachname).success;
            ReachesWithoutMaxJerk.(reachname).stim = Reaches.(reachname).stim;

            StimRR = [StimRR; NaN];
            MaxVel = [MaxVel; NaN, NaN];
            continue;
        end

        % 6. Get max jerk time and position
        max_jerk_time = accel_reach(max_jerk_idx, 1);
        max_jerk_pos = accel_reach(max_jerk_idx, 2);

        % 7. Save to ReReach struct
        ReReach.(reachname).kin = reach;
        ReReach.(reachname).inx = Reaches.(reachname).inx;
        ReReach.(reachname).success = Reaches.(reachname).success;
        ReReach.(reachname).stim = Reaches.(reachname).stim;
        ReReach.(reachname).jerk_time = max_jerk_time;
        ReReach.(reachname).jerk_pos = max_jerk_pos;

        % 8. Velocity before and after max jerk
        maxvbefore = max(accel_reach(1:max_jerk_idx, 7));
        maxvafter = max(accel_reach(max_jerk_idx:end, 7));

        ReReach.(reachname).maxvbefore = maxvbefore;
        ReReach.(reachname).maxvafter = maxvafter;

        % Store into shared arrays
        MaxVel = [MaxVel; maxvbefore, maxvafter];
        StimRR = [StimRR; max_jerk_time - stim_time];
    end
end
