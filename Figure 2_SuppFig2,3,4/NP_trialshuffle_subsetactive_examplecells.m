% Sampling even and odd trials to look out for spurious tiling

cd("C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\ShamAnimals\ReachandCSMMod_allCSMStogther\CSMAll_wtihMaxV");
MyFolderInfo = dir('*CSMAll.mat');

tic
clearvars -except MyFolderInfo

%% ALL REACHES PSTH + PEAK TIME ANALYSIS (Even/Odd + Early/Late)

NPdata = struct;

PsthMasterCtrl = [];
PsthMasterEven = [];
PsthMasterOdd  = [];

PsthMaster_Before_Ctrl = [];
PsthMaster_After_Ctrl  = [];
PsthMaster_Before_Even = [];
PsthMaster_After_Even  = [];
PsthMaster_Before_Odd  = [];
PsthMaster_After_Odd   = [];

for day = 1:length(MyFolderInfo)
    load(MyFolderInfo(day).name, 'Reaches', 'data', 'filtered_modulated_indices','ExcitRR');

    [Control, Stim] = StimFilter(Reaches);

    data_opto_all = data.units(filtered_modulated_indices);
    if isempty(data_opto_all)
        fprintf('Day %d: No modulated units.\n', day);
        continue
    end

    % Time vector
    timevec = data_opto_all(1).max_vx.time_bins;
    valid_peak_range = find(timevec >= -100 & timevec <= 100);

    % Determine which units have peaks within window
    keep_unit = false(length(data_opto_all), 1);
    peak_times_ctrl_all = NaN(length(data_opto_all), 1);

    for t = 1:length(data_opto_all)
        psth = data_opto_all(t).max_vx.PSTH_smoothbin;
        [~, peak_idx] = max(psth);
        peak_times_ctrl_all(t) = peak_idx;
        if ismember(peak_idx, valid_peak_range)
            keep_unit(t) = true;
        end
    end

    % Save only the filtered units for future plotting
    data_opto = data_opto_all(keep_unit);
    peak_times_window = peak_times_ctrl_all(keep_unit);

    if isempty(data_opto)
        fprintf('Day %d: No units with peaks in window. Skipping.\n', day);
        continue
    end

    % Convert peak indices to actual times (ms)
    peak_times_ms = timevec(peak_times_window);

    % Early vs late relative to 0 ms
    before_peakv_mask = (peak_times_ms < 0);
    after_peakv_mask  = (peak_times_ms >= 0);

    % Preallocate
    nCells = length(data_opto);
    nBins = length(data_opto(1).max_vx.PSTH_smoothbin);

    psth_master_ctrl = NaN(nCells, nBins);
    psth_master_even = NaN(nCells, nBins);
    psth_master_odd  = NaN(nCells, nBins);

    peak_times_ctrl = NaN(nCells,1);
    peak_times_even = NaN(nCells,1);
    peak_times_odd  = NaN(nCells,1);

    % --- Loop through each cell ---
    for t = 1:nCells
        % --- Control PSTH ---
        psth = data_opto(t).max_vx.PSTH_smoothbin;
        psth_master_ctrl(t,:) = psth;
        [~, peak_times_ctrl(t)] = max(psth);

        % --- Even/Odd PSTHs ---
        if isfield(data_opto(t).max_vx, 'Raster') && ~isempty(data_opto(t).max_vx.Raster)
            raster = data_opto(t).max_vx.Raster; % trials x time
            nTrials = size(raster,1);
            even_idx = 2:2:nTrials;
            odd_idx  = 1:2:nTrials;

            raster_even = raster(even_idx,:);
            raster_odd  = raster(odd_idx,:);

            bin_size = 10; % ms
            smoothing_bins = 3;

            [~, psth_smooth_even, ~] = computePSTH(raster_even, bin_size, smoothing_bins);
            [~, psth_smooth_odd,  ~] = computePSTH(raster_odd,  bin_size, smoothing_bins);

            % Normalize
            if max(psth_smooth_even) > 0
               psth_master_even(t,:) = psth_smooth_even / max(psth_smooth_even);
                %psth_master_even(t,:) = psth_smooth_even; % non normlized FR 
            else
               psth_master_even(t,:) = NaN(1,length(psth_smooth_even));
          
            end

            if max(psth_smooth_odd) > 0
                psth_master_odd(t,:)  = psth_smooth_odd / max(psth_smooth_odd);
                %psth_master_odd(t,:)  = psth_smooth_odd;
            else
                psth_master_odd(t,:)  = NaN(1,length(psth_smooth_odd));
            end

            [~, peak_times_even(t)] = max(psth_master_even(t,:));
            [~, peak_times_odd(t)]  = max(psth_master_odd(t,:));
        end
    end

    % --- Store into NPdata ---
    NPdata(day).peak_times_ctrl = peak_times_ctrl;
    NPdata(day).peak_times_even = peak_times_even;
    NPdata(day).peak_times_odd  = peak_times_odd;
    NPdata(day).psth_master_ctrl = psth_master_ctrl;
    NPdata(day).psth_master_even = psth_master_even;
    NPdata(day).psth_master_odd  = psth_master_odd;
    NPdata(day).before_peakv_mask = before_peakv_mask;
    NPdata(day).after_peakv_mask  = after_peakv_mask;
    NPdata(day).kept_indices = filtered_modulated_indices(keep_unit);

    % --- NEW: save full unit struct for later plotting ---
    NPdata(day).kept_units = data_opto;

    % --- Concatenate across days ---
    PsthMasterCtrl = vertcat(PsthMasterCtrl, psth_master_ctrl);
    PsthMasterEven = vertcat(PsthMasterEven, psth_master_even);
    PsthMasterOdd  = vertcat(PsthMasterOdd,  psth_master_odd);

    PsthMaster_Before_Ctrl = vertcat(PsthMaster_Before_Ctrl, psth_master_ctrl(before_peakv_mask,:));
    PsthMaster_After_Ctrl  = vertcat(PsthMaster_After_Ctrl,  psth_master_ctrl(after_peakv_mask,:));

    PsthMaster_Before_Even = vertcat(PsthMaster_Before_Even, psth_master_even(before_peakv_mask,:));
    PsthMaster_After_Even  = vertcat(PsthMaster_After_Even,  psth_master_even(after_peakv_mask,:));

    PsthMaster_Before_Odd  = vertcat(PsthMaster_Before_Odd,  psth_master_odd(before_peakv_mask,:));
    PsthMaster_After_Odd   = vertcat(PsthMaster_After_Odd,   psth_master_odd(after_peakv_mask,:));

    fprintf('Day %d: Added %d units (%d early, %d late)\n', day, nCells, sum(before_peakv_mask), sum(after_peakv_mask));
    clear data Reaches Control filtered_modulated_indices ExcitRR
end


%% NEW EDIT 
% Generate raster + PSTH for each kept unit and save as SVG (white raster background)

output_folder = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\Sham_evenoddtrials\ExampleUnits_max_vx aligned";

oddcolor = rgb('DarkRed');
evencolor = rgb('DarkGreen');

bin_size = 10;         % ms
smoothing_bins = 3;    % for PSTH smoothing
time_range = [-1000, 1000];  % ms, full PSTH range
time_window = [-500, 500];   % ms, plotting window

for day = 1:length(NPdata)
    if ~isfield(NPdata(day), 'kept_units') || isempty(NPdata(day).kept_units)
        continue
    end

    nUnits = length(NPdata(day).kept_units);

    for u = 1:nUnits
        unit = NPdata(day).kept_units(u);
        raster = unit.max_vx.Raster;            % trials x time

        if isempty(raster)
            fprintf('Day %d, Unit %d: empty raster, skipping.\n', day, u);
            continue
        end

        timevec_raster = linspace(unit.max_vx.time_bins(1), unit.max_vx.time_bins(end), size(raster,2));

        nTrials = size(raster, 1);
        even_idx = 2:2:nTrials;
        odd_idx  = 1:2:nTrials;

        raster_even = raster(even_idx, :);
        raster_odd  = raster(odd_idx, :);

        % --- Compute PSTHs (with smoothing) ---
        [~, psth_all, ~]  = computePSTH(raster, bin_size, smoothing_bins);
        [~, psth_even, ~] = computePSTH(raster_even, bin_size, smoothing_bins);
        [~, psth_odd,  ~] = computePSTH(raster_odd,  bin_size, smoothing_bins);

        % --- Construct consistent PSTH time vector ---
        num_time_bins = size(psth_all, 2);
        timevec_psth = linspace(time_range(1), time_range(2), num_time_bins);

        % --- Apply time window mask for plotting ---
        time_mask = timevec_psth >= time_window(1) & timevec_psth <= time_window(2);
        timevec_plot = timevec_psth(time_mask);
        psth_all_plot  = psth_all(time_mask);
        psth_even_plot = psth_even(time_mask);
        psth_odd_plot  = psth_odd(time_mask);

        % --- Prepare spike times for raster safely ---
        [trial_idx, col_idx] = find(raster);  % trial x raster bin index
        spike_times = timevec_raster(col_idx);

        % --- Plot Raster ---
        fig = figure('Visible', 'off', 'Position', [100, 100, 600, 500]);
        subplot(2, 1, 1);
        hold on;
        set(gca, 'Color', 'w'); % white background

        % Plot each spike as a small black dot
        scatter(spike_times, trial_idx, 5, 'k', 'filled'); 
        xlim(time_window);
        ylim([0 nTrials+1]);
        xlabel('Time (ms)');
        ylabel('Trials');
        title(sprintf('Day %d, Unit %d Raster', day, u));
        box off;
        set(gca, 'TickDir', 'out');

        % --- Plot PSTH ---
        subplot(2, 1, 2);
        plot(timevec_plot, psth_all_plot, 'k', 'LineWidth', 1.5); hold on;
        plot(timevec_plot, psth_even_plot, 'Color', evencolor, 'LineWidth', 1.2);
        plot(timevec_plot, psth_odd_plot, 'Color', oddcolor, 'LineWidth', 1.2);
        xlim(time_window);
        xlabel('Time (ms)');
        ylabel('Firing Rate (spikes/s)');
        legend({'All', 'Even', 'Odd'}, 'Location', 'best');
        title(sprintf('Day %d, Unit %d PSTH (Firing Rate)', day, u));
        box off;
        set(gca, 'TickDir', 'out');

        % --- Save figure as SVG ---
        filename = fullfile(output_folder, sprintf('Day%d_Unit%d.svg', day, u));
        saveas(fig, filename);
        close(fig);
    end
end

%% NEw code that says whether it is an erly or late cell
% Generate raster + PSTH for each kept unit and save as SVG (white raster background)
%output_folder = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\Sham_evenoddtrials\ExampleUnits_max_vx aligned\units with labels\tick marks";
%output_folder= "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\Sham_evenoddtrials\Example Units endpoint aligned";

%output_folder="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\Sham_evenoddtrials\pos_thresh aligned";
output_folder="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\Sham_evenoddtrials\ExampleUnits_max_vx aligned\1000ms time window";

oddcolor = rgb('DarkRed');
evencolor = rgb('DarkGreen');

bin_size = 10;         % ms
smoothing_bins = 3;    % for PSTH smoothing
time_range = [-1000, 1000];  % ms, full PSTH range
%time_window = [-500, 500];   % ms, plotting window
time_window = [-1000, 1000];   % ms, plotting window

for day = 1:length(NPdata)
    if ~isfield(NPdata(day), 'kept_units') || isempty(NPdata(day).kept_units)
        continue
    end

    nUnits = length(NPdata(day).kept_units);

    for u = 1:nUnits
        unit = NPdata(day).kept_units(u);
        raster = unit.max_vx.Raster;  % trials x time
        %raster = unit.end_point.Raster;  % trials x time
        %raster = unit.pos_thresh.Raster;  % trials x time

        if isempty(raster)
            fprintf('Day %d, Unit %d: empty raster, skipping.\n', day, u);
            continue
        end

        % Determine early vs late category
        if NPdata(day).before_peakv_mask(u)
            unit_type = 'Early';
        elseif NPdata(day).after_peakv_mask(u)
            unit_type = 'Late';
        else
            unit_type = 'Uncategorized';
        end

        timevec_raster = linspace(unit.max_vx.time_bins(1), unit.max_vx.time_bins(end), size(raster,2));
        %timevec_raster = linspace(unit.end_point.time_bins(1), unit.end_point.time_bins(end), size(raster,2));
        %timevec_raster = linspace(unit.pos_thresh.time_bins(1), unit.pos_thresh.time_bins(end), size(raster,2));

        nTrials = size(raster, 1);
        even_idx = 2:2:nTrials;
        odd_idx  = 1:2:nTrials;

        raster_even = raster(even_idx, :);
        raster_odd  = raster(odd_idx, :);

        % --- Compute PSTHs (with smoothing) ---
        [~, psth_all, ~]  = computePSTH(raster, bin_size, smoothing_bins);
        [~, psth_even, ~] = computePSTH(raster_even, bin_size, smoothing_bins);
        [~, psth_odd,  ~] = computePSTH(raster_odd,  bin_size, smoothing_bins);

        % --- Construct consistent PSTH time vector ---
        num_time_bins = size(psth_all, 2);
        timevec_psth = linspace(time_range(1), time_range(2), num_time_bins);

        % --- Apply time window mask for plotting ---
        time_mask = timevec_psth >= time_window(1) & timevec_psth <= time_window(2);
        timevec_plot = timevec_psth(time_mask);
        psth_all_plot  = psth_all(time_mask);
        psth_even_plot = psth_even(time_mask);
        psth_odd_plot  = psth_odd(time_mask);

        % --- Prepare spike times for raster safely ---
        [trial_idx, col_idx] = find(raster);  % trial x raster bin index
        spike_times = timevec_raster(col_idx);

        % % --- Plot Raster ---
        fig = figure('Visible', 'off', 'Position', [100, 100, 600, 500]);
        subplot(2, 1, 1);
        hold on;
        set(gca, 'Color', 'w'); % white background
       % scatter(spike_times, trial_idx, 5, 'k', 'filled'); 
        scatter(spike_times, trial_idx, 2, 'k', 'filled'); 
        xlim(time_window);
        ylim([0 nTrials+1]);
        xlabel('Time (ms)');
        ylabel('Trials');
        title(sprintf('Day %d, Unit %d Raster (%s)', day, u, unit_type));
        box off;
        set(gca, 'TickDir', 'out');
        % --- Plot Raster ---
        % fig = figure('Visible', 'off', 'Position', [100, 100, 600, 500]);
        % subplot(2, 1, 1);
        % hold on;
        % set(gca, 'Color', 'w'); % white background
        % 
        % % Plot each spike as a short vertical dash
        % spike_height = 0.8;  % length of each dash in trials
        % for k = 1:length(spike_times)
        %     line([spike_times(k) spike_times(k)], [trial_idx(k)-spike_height/2, trial_idx(k)+spike_height/2], ...
        %         'Color', 'k', 'LineWidth', 1);
        % end
        % 
        % xlim(time_window);
        % ylim([0 nTrials+1]);
        % xlabel('Time (ms)');
        % ylabel('Trials');
        % title(sprintf('Day %d, Unit %d Raster (%s)', day, u, unit_type));
        % box off;
        % set(gca, 'TickDir', 'out');


        % --- Plot PSTH ---
        subplot(2, 1, 2);
        plot(timevec_plot, psth_all_plot, 'k', 'LineWidth', 1.5); hold on;
        plot(timevec_plot, psth_even_plot, 'Color', evencolor, 'LineWidth', 1.2);
        plot(timevec_plot, psth_odd_plot, 'Color', oddcolor, 'LineWidth', 1.2);
        xlim(time_window);
        %xlabel('Time (ms)');
        xlabel('Time from endpoint (ms)');
        ylabel('Firing Rate (spikes/s)');
        legend({'All', 'Even', 'Odd'}, 'Location', 'best');
        title(sprintf('Day %d, Unit %d PSTH (%s)', day, u, unit_type));
        box off;
        set(gca, 'TickDir', 'out');

        % --- Save figure as SVG ---
        filename = fullfile(output_folder, sprintf('Day%d_Unit%d_%s.svg', day, u, unit_type));
        saveas(fig, filename);
        close(fig);
    end
end

%% ALIGN TO CSM T0 or T0
%%= NEw code that says whether it is an erly or late cell
% Generate raster + PSTH for each kept unit and save as SVG (white raster background)
%output_folder ="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\Sham_evenoddtrials\Example Units_CSMT0 aligned";
%output_folder="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\Sham_evenoddtrials\Example Units_CSMT2 aligned";
output_folder="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\Sham Neural analysis\Endpoint Aligned";
oddcolor = rgb('Black'); % NO CSM reaches 
evencolor = rgb('DarkMagenta'); % CSM reaches 

bin_size = 10;         % ms
smoothing_bins = 3;    % for PSTH smoothing
time_range = [-1000, 1000];  % ms, full PSTH range
%time_window = [-500, 500];   % ms, plotting window
%time_window = [-200, 200];   % ms, plotting window
time_window = [-1000, 1000];   % ms, plotting window

for day = 1:length(NPdata)
    if ~isfield(NPdata(day), 'kept_units') || isempty(NPdata(day).kept_units)
        continue
    end

    nUnits = length(NPdata(day).kept_units);

    for u = 1:nUnits
        unit = NPdata(day).kept_units(u);
        %raster = unit.max_vx.Raster;  % trials x time

        raster = unit.end_point.Raster;  % trials x time
       % raster = unit.AllCSM_csmT0.Raster;  % trials x time
        %raster = unit.AllCSM_csmT2.Raster;  % trials x time

        if isempty(raster)
            fprintf('Day %d, Unit %d: empty raster, skipping.\n', day, u);
            continue
        end

        % Determine early vs late category
        if NPdata(day).before_peakv_mask(u)
            unit_type = 'Early';
        elseif NPdata(day).after_peakv_mask(u)
            unit_type = 'Late';
        else
            unit_type = 'Uncategorized';
        end

        timevec_raster = linspace(unit.max_vx.time_bins(1), unit.max_vx.time_bins(end), size(raster,2));
 

        nTrials = size(raster, 1);
        even_idx = 2:2:nTrials;
        odd_idx  = 1:2:nTrials;

        raster_even = raster(even_idx, :);
        raster_odd  = raster(odd_idx, :);
        % 
        % raster_even = unit.AllCSM_csmT0.Raster; 
        % raster_odd  = unit.NoCSM_csmT0.Raster; 

        % raster_even = unit.AllCSM_csmT2.Raster; 
        % raster_odd  = unit.NoCSM_csmT2.Raster; 

        % --- Compute PSTHs (with smoothing) ---
        [~, psth_all, ~]  = computePSTH(raster, bin_size, smoothing_bins);
        [~, psth_even, ~] = computePSTH(raster_even, bin_size, smoothing_bins);
        [~, psth_odd,  ~] = computePSTH(raster_odd,  bin_size, smoothing_bins);

        % --- Construct consistent PSTH time vector ---
        num_time_bins = size(psth_all, 2);
        timevec_psth = linspace(time_range(1), time_range(2), num_time_bins);

        % --- Apply time window mask for plotting ---
        time_mask = timevec_psth >= time_window(1) & timevec_psth <= time_window(2);
        timevec_plot = timevec_psth(time_mask);
        psth_all_plot  = psth_all(time_mask);
        psth_even_plot = psth_even(time_mask);
        psth_odd_plot  = psth_odd(time_mask);

        % --- Prepare spike times for raster safely ---
        [trial_idx, col_idx] = find(raster);  % trial x raster bin index
        spike_times = timevec_raster(col_idx);

        % % --- Plot Raster ---
        fig = figure('Visible', 'off', 'Position', [100, 100, 600, 500]);
        subplot(2, 1, 1);
        hold on;
        set(gca, 'Color', 'w'); % white background
        scatter(spike_times, trial_idx, 5, 'k', 'filled'); 
        xlim(time_window);
        ylim([0 nTrials+1]);
        xlabel('Time (ms)');
        ylabel('Trials');
        title(sprintf('Day %d, Unit %d Raster (%s)', day, u, unit_type));
        box off;
        set(gca, 'TickDir', 'out');
        % --- Plot Raster ---
        % fig = figure('Visible', 'off', 'Position', [100, 100, 600, 500]);
        % subplot(2, 1, 1);
        % hold on;
        % set(gca, 'Color', 'w'); % white background
        % 
        % % Plot each spike as a short vertical dash
        % spike_height = 0.8;  % length of each dash in trials
        % for k = 1:length(spike_times)
        %     line([spike_times(k) spike_times(k)], [trial_idx(k)-spike_height/2, trial_idx(k)+spike_height/2], ...
        %         'Color', 'k', 'LineWidth', 1);
        % end
        % 
        % xlim(time_window);
        % ylim([0 nTrials+1]);
        % xlabel('Time (ms)');
        % ylabel('Trials');
        % title(sprintf('Day %d, Unit %d Raster (%s)', day, u, unit_type));
        % box off;
        % set(gca, 'TickDir', 'out');


        % --- Plot PSTH ---
        subplot(2, 1, 2);
        plot(timevec_plot, psth_all_plot, 'k', 'LineWidth', 1.5); hold on;
        plot(timevec_plot, psth_even_plot, 'Color', evencolor, 'LineWidth', 1.2);hold on;
        plot(timevec_plot, psth_odd_plot, 'Color', oddcolor, 'LineWidth', 1.2);
        xlim(time_window);
        %xlabel('Time (ms)');
        xlabel('Time from CSM T0 (ms)');
        %xlabel('Time from CSM T2 (ms)');
        ylabel('Firing Rate (spikes/s)');
        legend({'All', 'CSM', 'NoCSM'}, 'Location', 'best');
        title(sprintf('Day %d, Unit %d PSTH (%s)', day, u, unit_type));
        box off;
        set(gca, 'TickDir', 'out');

        % --- Save figure as SVG ---
        filename = fullfile(output_folder, sprintf('Day%d_Unit%d_%s.svg', day, u, unit_type));
        saveas(fig, filename);
        close(fig);
    end
end
%% ---- STEP 1: Sort each PSTH matrix by its own peak firing time ----

% Compute peak indices for each condition
[~, peak_idx_ctrl] = max(PsthMasterCtrl, [], 2);
[~, peak_idx_even] = max(PsthMasterEven, [], 2);
[~, peak_idx_odd]  = max(PsthMasterOdd, [], 2);

% Sort indices by peak timing
[~, sort_idx_ctrl] = sort(peak_idx_ctrl);
[~, sort_idx_even] = sort(peak_idx_even);
[~, sort_idx_odd]  = sort(peak_idx_odd);

% Sort the PSTH matrices
sorted_ctrl = PsthMasterCtrl(sort_idx_ctrl, :);
sorted_even = PsthMasterEven(sort_idx_even, :);
sorted_odd  = PsthMasterOdd(sort_idx_odd, :);

% Generate time vector for x-axis
num_bins = size(PsthMasterCtrl, 2);
time_vector = linspace(-1000, 1000, num_bins);  % adjust if bins differ
tick_times = -1000:250:1000;
[~, tick_indices] = min(abs(time_vector' - tick_times), [], 1);

% Plot all 3, each sorted by its own peak time
figure('Position', [100 100 1600 600]);

subplot(1,3,1);
imagesc(sorted_ctrl);
colormap(parula);
colorbar;
set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times);
xlabel('Time (ms)');
ylabel('Units (sorted by own peak)');
title('Control (sorted by own peaks)');

subplot(1,3,2);
imagesc(sorted_even);
colormap(parula);
colorbar;
set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times);
xlabel('Time (ms)');
ylabel('Units');
title('Even Trials (sorted by own peaks)');

subplot(1,3,3);
imagesc(sorted_odd);
colormap(parula);
colorbar;
set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times);
xlabel('Time (ms)');
ylabel('Units');
title('Odd Trials (sorted by own peaks)');

sgtitle('Each PSTH Sorted by Its Own Peak Firing Time');


%% ---- STEP 2: Sort all PSTHs according to the Control’s peak order ----

% Apply control’s sorting order to all three
sorted_by_ctrl_ctrl = PsthMasterCtrl(sort_idx_ctrl, :);
sorted_by_ctrl_even = PsthMasterEven(sort_idx_ctrl, :);
sorted_by_ctrl_odd  = PsthMasterOdd(sort_idx_ctrl, :);

% Generate time vector and x-axis tick marks
num_bins = size(PsthMasterCtrl, 2);
time_vector = linspace(-1000, 1000, num_bins);   % adjust range if needed
tick_times = -1000:250:1000;                     % tick positions in ms
[~, tick_indices] = min(abs(time_vector' - tick_times), [], 1);

figure('Position', [100 100 1600 600]);

subplot(1,3,1);
imagesc(sorted_by_ctrl_ctrl);
colormap(parula);
colorbar;
set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times);
xlabel('Time (ms)');
ylabel('Units (sorted by Control peaks)');
title('Control');

subplot(1,3,2);
imagesc(sorted_by_ctrl_even);
colormap(parula);
colorbar;
set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times);
xlabel('Time (ms)');
ylabel('Units');
title('Even Trials (sorted by Control peaks)');

subplot(1,3,3);
imagesc(sorted_by_ctrl_odd);
colormap(parula);
colorbar;
set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times);
xlabel('Time (ms)');
ylabel('Units');
title('Odd Trials (sorted by Control peaks)');

sgtitle('All PSTHs Ordered by Control Peak Firing Time');
%% Plot only a subset of the window 
%% ---- STEP 2: Sort all PSTHs according to the Control’s peak order ----

% Apply control’s sorting order to all three
sorted_by_ctrl_ctrl = PsthMasterCtrl(sort_idx_ctrl, :);
sorted_by_ctrl_even = PsthMasterEven(sort_idx_ctrl, :);
sorted_by_ctrl_odd  = PsthMasterOdd(sort_idx_ctrl, :);

% Generate time vector and define window
num_bins = size(PsthMasterCtrl, 2);
time_vector = linspace(-1000, 1000, num_bins);   % full range
plot_window = [-250 250];                        % window in ms

% Find indices corresponding to the plot window
[~, start_idx] = min(abs(time_vector - plot_window(1)));
[~, end_idx]   = min(abs(time_vector - plot_window(2)));

% Crop PSTHs and time vector
time_vector_crop = time_vector(start_idx:end_idx);
sorted_by_ctrl_ctrl = sorted_by_ctrl_ctrl(:, start_idx:end_idx);
sorted_by_ctrl_even = sorted_by_ctrl_even(:, start_idx:end_idx);
sorted_by_ctrl_odd  = sorted_by_ctrl_odd(:, start_idx:end_idx);

% Define tick marks within the cropped window
tick_times = -250:125:250;
[~, tick_indices] = min(abs(time_vector_crop' - tick_times), [], 1);

figure('Position', [100 100 1600 600]);

subplot(1,3,1);
imagesc(time_vector_crop, 1:size(sorted_by_ctrl_ctrl,1), sorted_by_ctrl_ctrl);
colormap(parula);
colorbar;
set(gca, 'XTick', tick_times);
xlabel('Time (ms)');
ylabel('Units (sorted by Control peaks)');
title('Control');

subplot(1,3,2);
imagesc(time_vector_crop, 1:size(sorted_by_ctrl_even,1), sorted_by_ctrl_even);
colormap(parula);
colorbar;
set(gca, 'XTick', tick_times);
xlabel('Time (ms)');
ylabel('Units');
title('Even Trials (sorted by Control peaks)');

subplot(1,3,3);
imagesc(time_vector_crop, 1:size(sorted_by_ctrl_odd,1), sorted_by_ctrl_odd);
colormap(parula);
colorbar;
set(gca, 'XTick', tick_times);
xlabel('Time (ms)');
ylabel('Units');
title('Odd Trials (sorted by Control peaks)');

sgtitle('All PSTHs Ordered by Control Peak Firing Time (−250 to 250 ms)');

%% Now plot using a similar style as before - Endpoint aligned all trials 

sorted_psth_ctrl = PsthMasterCtrl(sort_idx_ctrl, :);
%sorted_psth_norr = PsthMaster_NoRR(sort_idx_csm, :);

% Generate time vector
num_bins = size(sorted_psth_ctrl, 2);
time_vector = linspace(-1000, 1000, num_bins);  % assuming 10 ms bins over 200 bins
tick_times = -1000:250:1000;  % Set tick marks for x-axis
[~, tick_indices] = min(abs(time_vector' - tick_times), [], 1);  % find closest indices

figure;
imagesc(sorted_psth_ctrl);
colormap(slanCM('GnBu'));
colorbar;
set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times);
xlabel('Time from Endpoint (ms)');
ylabel('Units (sorted by peak firing time)');
title('PSTHs Sorted by Peak Firing Time- All Trials');
%% Window
% Plot Even Trials Sorted by Control Peak Times (−250 to 250 ms)
sorted_psth_ctrl = PsthMasterCtrl(sort_idx_ctrl, :);

% Generate full time vector
num_bins = size(sorted_psth_ctrl, 2);
time_vector = linspace(-1000, 1000, num_bins);  % full 2000 ms window

% Define window and crop indices
plot_window = [-250 250];
[~, start_idx] = min(abs(time_vector - plot_window(1)));
[~, end_idx]   = min(abs(time_vector - plot_window(2)));

% Crop time vector and PSTH
time_vector_crop = time_vector(start_idx:end_idx);
sorted_psth_ctrl_crop = sorted_psth_ctrl(:, start_idx:end_idx);

% Define tick marks within cropped range
tick_times = -250:125:250;
[~, tick_indices] = min(abs(time_vector_crop' - tick_times), [], 1);

% Plot cropped heatmap
figure;
imagesc(time_vector_crop, 1:size(sorted_psth_ctrl_crop,1), sorted_psth_ctrl_crop);
colormap(slanCM('GnBu'));
colorbar;
xline(0)
set(gca, 'XTick', tick_times);
xlabel('Time from Endpoint (ms)');
ylabel('Units (sorted by peak firing time)');
title('PSTHs Sorted by Peak Firing Time - All Trials (−250 to 250 ms)');

%% Now plot even trials according to ctrl peak firing rate idx 

sorted_psth_even = PsthMasterEven(sort_idx_ctrl, :);

% Generate time vector
num_bins = size(sorted_psth_even, 2);
time_vector = linspace(-1000, 1000, num_bins);  % assuming 10 ms bins over 200 bins
tick_times = -1000:250:1000;  % Set tick marks for x-axis
[~, tick_indices] = min(abs(time_vector' - tick_times), [], 1);  % find closest indices


figure;
imagesc(sorted_psth_even);
colormap(slanCM('GnBu'));
colorbar;
set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times);
xlabel('Time from Endpoint (ms)');
ylabel('Units (sorted by peak firing time)');
title('PSTHs Sorted by Peak Firing Time- Even Trials Only');
%% Window
% Plot Even Trials Sorted by Control Peak Times (−250 to 250 ms)
sorted_psth_even = PsthMasterEven(sort_idx_ctrl, :);

% Generate full time vector
num_bins = size(sorted_psth_even, 2);
time_vector = linspace(-1000, 1000, num_bins);  % full 2000 ms window

% Define window and crop indices
plot_window = [-250 250];
[~, start_idx] = min(abs(time_vector - plot_window(1)));
[~, end_idx]   = min(abs(time_vector - plot_window(2)));

% Crop time vector and PSTH
time_vector_crop = time_vector(start_idx:end_idx);
sorted_psth_even_crop = sorted_psth_even(:, start_idx:end_idx);

% Define tick marks within cropped range
tick_times = -250:125:250;
[~, tick_indices] = min(abs(time_vector_crop' - tick_times), [], 1);

% Plot cropped heatmap
figure;
imagesc(time_vector_crop, 1:size(sorted_psth_even_crop,1), sorted_psth_even_crop);
colormap(slanCM('GnBu'));
colorbar;
xline(0)
set(gca, 'XTick', tick_times);
xlabel('Time from Endpoint (ms)');
ylabel('Units (sorted by peak firing time)');
title('PSTHs Sorted by Peak Firing Time - Even Trials Only (−250 to 250 ms)');


%% Now plot even trials according to even peak firing rate idx 

sorted_psth_even = PsthMasterEven(sort_idx_even, :);

% Generate time vector
num_bins = size(sorted_psth_even, 2);
time_vector = linspace(-1000, 1000, num_bins);  % assuming 10 ms bins over 200 bins
tick_times = -1000:250:1000;  % Set tick marks for x-axis
[~, tick_indices] = min(abs(time_vector' - tick_times), [], 1);  % find closest indices


figure;
imagesc(sorted_psth_even);
colormap(slanCM('GnBu'));
colorbar;
set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times);
xlabel('Time from Endpoint (ms)');
ylabel('Units (sorted by peak firing time)');
title('PSTHs Sorted by Peak Firing Time- Even Trials Only');

%% Now plot odd trials according to ctrl peak firing rate idx 

sorted_psth_odd = PsthMasterOdd(sort_idx_ctrl, :);

% Generate time vector
num_bins = size(sorted_psth_odd, 2);
time_vector = linspace(-1000, 1000, num_bins);  % assuming 10 ms bins over 200 bins
tick_times = -1000:250:1000;  % Set tick marks for x-axis
[~, tick_indices] = min(abs(time_vector' - tick_times), [], 1);  % find closest indices


figure;
imagesc(sorted_psth_odd);
colormap(slanCM('GnBu'));
colorbar;
set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times);
xlabel('Time from Endpoint (ms)');
ylabel('Units (sorted by peak firing time)');
title('PSTHs Sorted by Peak Firing Time- Odd Trials Only');


%% Window
% Plot Even Trials Sorted by Control Peak Times (−250 to 250 ms)
sorted_psth_odd = PsthMasterOdd(sort_idx_ctrl, :);

% Generate full time vector
num_bins = size(sorted_psth_odd, 2);
time_vector = linspace(-1000, 1000, num_bins);  % full 2000 ms window

% Define window and crop indices
plot_window = [-250 250];
[~, start_idx] = min(abs(time_vector - plot_window(1)));
[~, end_idx]   = min(abs(time_vector - plot_window(2)));

% Crop time vector and PSTH
time_vector_crop = time_vector(start_idx:end_idx);
sorted_psth_odd_crop = sorted_psth_odd(:, start_idx:end_idx);

% Define tick marks within cropped range
tick_times = -250:125:250;
[~, tick_indices] = min(abs(time_vector_crop' - tick_times), [], 1);

% Plot cropped heatmap
figure;
imagesc(time_vector_crop, 1:size(sorted_psth_odd_crop,1), sorted_psth_odd_crop);
colormap(slanCM('GnBu'));
colorbar;
xline(0)
set(gca, 'XTick', tick_times);
xlabel('Time from Endpoint (ms)');
ylabel('Units (sorted by peak firing time)');
title('PSTHs Sorted by Peak Firing Time - Odd Trials Only (−250 to 250 ms)');

%% Now plot odd trials according to odd peak firing rate idx 

sorted_psth_odd = PsthMasterOdd(sort_idx_odd, :);

% Generate time vector
num_bins = size(sorted_psth_odd, 2);
time_vector = linspace(-1000, 1000, num_bins);  % assuming 10 ms bins over 200 bins
tick_times = -1000:250:1000;  % Set tick marks for x-axis
[~, tick_indices] = min(abs(time_vector' - tick_times), [], 1);  % find closest indices


figure;
imagesc(sorted_psth_odd);
colormap(slanCM('GnBu'));
colorbar;
set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times);
xlabel('Time from Endpoint (ms)');
ylabel('Units (sorted by peak firing time)');
title('PSTHs Sorted by Peak Firing Time- Odd Trials Only');

%% Editing new code

%% ---- STEP 1: Sort each PSTH matrix by its own peak firing time ----
% Compute peak indices for each condition
[~, peak_idx_ctrl] = max(PsthMasterCtrl, [], 2);
[~, peak_idx_even] = max(PsthMasterEven, [], 2);
[~, peak_idx_odd]  = max(PsthMasterOdd, [], 2);

% Sort indices by peak timing
[~, sort_idx_ctrl] = sort(peak_idx_ctrl);
[~, sort_idx_even] = sort(peak_idx_even);
[~, sort_idx_odd]  = sort(peak_idx_odd);

% Separate early/late masks based on 0 ms peak times (apply to ctrl)
peak_times_ctrl_ms = time_vector(peak_idx_ctrl);
before_ctrl_mask = peak_times_ctrl_ms < 0;
after_ctrl_mask  = peak_times_ctrl_ms >= 0;

% Sort the PSTH matrices by their own peak
sorted_ctrl = PsthMasterCtrl(sort_idx_ctrl, :);
sorted_even = PsthMasterEven(sort_idx_even, :);

sorted_odd  = PsthMasterOdd(sort_idx_odd, :);

% Sort early/late separately (own peak sorting)
sorted_ctrl_early = sorted_ctrl(before_ctrl_mask(sort_idx_ctrl), :);
sorted_ctrl_late  = sorted_ctrl(after_ctrl_mask(sort_idx_ctrl), :);

sorted_even_early = sorted_even(before_ctrl_mask(sort_idx_even), :);
sorted_even_late  = sorted_even(after_ctrl_mask(sort_idx_even), :);

sorted_odd_early  = sorted_odd(before_ctrl_mask(sort_idx_odd), :);
sorted_odd_late   = sorted_odd(after_ctrl_mask(sort_idx_odd), :);

% Generate time vector for x-axis
num_bins = size(PsthMasterCtrl, 2);
time_vector = linspace(-1000, 1000, num_bins);
tick_times = -1000:250:1000;
[~, tick_indices] = min(abs(time_vector' - tick_times), [], 1);

% Plot all 3 sorted by own peak
figure('Position',[100 100 1800 600]);
subplot(1,3,1); imagesc(sorted_ctrl); colormap(parula); colorbar;
set(gca,'XTick',tick_indices,'XTickLabel',tick_times);
xlabel('Time (ms)'); ylabel('Units'); title('Control - sorted by own peaks');

subplot(1,3,2); imagesc(sorted_even); colormap(parula); colorbar;
set(gca,'XTick',tick_indices,'XTickLabel',tick_times);
xlabel('Time (ms)'); ylabel('Units'); title('Even - sorted by own peaks');

subplot(1,3,3); imagesc(sorted_odd); colormap(parula); colorbar;
set(gca,'XTick',tick_indices,'XTickLabel',tick_times);
xlabel('Time (ms)'); ylabel('Units'); title('Odd - sorted by own peaks');

sgtitle('Each PSTH Sorted by Its Own Peak Firing Time');

%% ---- STEP 2: Sort all PSTHs by Control’s peak firing order ----
sorted_by_ctrl_ctrl = PsthMasterCtrl(sort_idx_ctrl, :);
sorted_by_ctrl_even = PsthMasterEven(sort_idx_ctrl, :);
sorted_by_ctrl_odd  = PsthMasterOdd(sort_idx_ctrl, :);

% Early/late for Control-sorted
sorted_by_ctrl_early_ctrl = sorted_by_ctrl_ctrl(before_ctrl_mask(sort_idx_ctrl), :);
sorted_by_ctrl_late_ctrl  = sorted_by_ctrl_ctrl(after_ctrl_mask(sort_idx_ctrl), :);

sorted_by_ctrl_early_even = sorted_by_ctrl_even(before_ctrl_mask(sort_idx_ctrl), :);
sorted_by_ctrl_late_even  = sorted_by_ctrl_even(after_ctrl_mask(sort_idx_ctrl), :);

sorted_by_ctrl_early_odd  = sorted_by_ctrl_odd(before_ctrl_mask(sort_idx_ctrl), :);
sorted_by_ctrl_late_odd   = sorted_by_ctrl_odd(after_ctrl_mask(sort_idx_ctrl), :);

% Plot sorted by control peaks
figure('Position',[100 100 1800 600]);
subplot(1,3,1); imagesc(sorted_by_ctrl_ctrl); colormap(parula); colorbar;
set(gca,'XTick',tick_indices,'XTickLabel',tick_times);
xlabel('Time (ms)'); ylabel('Units'); title('Control - sorted by Control peaks');

subplot(1,3,2); imagesc(sorted_by_ctrl_even); colormap(parula); colorbar;
set(gca,'XTick',tick_indices,'XTickLabel',tick_times);
xlabel('Time (ms)'); ylabel('Units'); title('Even - sorted by Control peaks');

subplot(1,3,3); imagesc(sorted_by_ctrl_odd); colormap(parula); colorbar;
set(gca,'XTick',tick_indices,'XTickLabel',tick_times);
xlabel('Time (ms)'); ylabel('Units'); title('Odd - sorted by Control peaks');

sgtitle('All PSTHs Ordered by Control Peak Firing Time');

%% ---- STEP 3: Endpoint-aligned heatmaps (Ctrl / Even / Odd) ----
% Control
figure; imagesc(sorted_by_ctrl_ctrl); colormap(slanCM('GnBu')); colorbar;
set(gca,'XTick',tick_indices,'XTickLabel',tick_times);
xlabel('Time from Endpoint (ms)'); ylabel('Units'); title('PSTHs Sorted by Peak Firing Time- All Trials');

% Even
figure; imagesc(sorted_by_ctrl_even); colormap(slanCM('GnBu')); colorbar;
set(gca,'XTick',tick_indices,'XTickLabel',tick_times);
xlabel('Time from Endpoint (ms)'); ylabel('Units'); title('PSTHs Sorted by Peak Firing Time- Even Trials');

% Odd
figure; imagesc(sorted_by_ctrl_odd); colormap(slanCM('GnBu')); colorbar;
set(gca,'XTick',tick_indices,'XTickLabel',tick_times);
xlabel('Time from Endpoint (ms)'); ylabel('Units'); title('PSTHs Sorted by Peak Firing Time- Odd Trials');


%% Edited code
% ---- STEP 4: Early vs Late PSTHs (Control-sorted) ----
figure('Position',[100 100 1800 600]);

lineStyle = {'--k','LineWidth',1.5};

subplot(2,3,1);
imagesc(time_vector, 1:size(sorted_by_ctrl_early_ctrl,1), sorted_by_ctrl_early_ctrl);
colormap(slanCM('GnBu')); colorbar;
xline(0, lineStyle{:});
title('Ctrl - Early');
set(gca,'XTick',tick_times);

subplot(2,3,4);
imagesc(time_vector, 1:size(sorted_by_ctrl_late_ctrl,1), sorted_by_ctrl_late_ctrl);
colormap(slanCM('GnBu')); colorbar;
xline(0, lineStyle{:});
title('Ctrl - Late');
set(gca,'XTick',tick_times);

subplot(2,3,2);
imagesc(time_vector, 1:size(sorted_by_ctrl_early_even,1), sorted_by_ctrl_early_even);
colormap(slanCM('GnBu')); colorbar;
xline(0, lineStyle{:});
title('Even - Early');
set(gca,'XTick',tick_times);

subplot(2,3,5);
imagesc(time_vector, 1:size(sorted_by_ctrl_late_even,1), sorted_by_ctrl_late_even);
colormap(slanCM('GnBu')); colorbar;
xline(0, lineStyle{:});
title('Even - Late');
set(gca,'XTick',tick_times);

subplot(2,3,3);
imagesc(time_vector, 1:size(sorted_by_ctrl_early_odd,1), sorted_by_ctrl_early_odd);
colormap(slanCM('GnBu')); colorbar;
xline(0, lineStyle{:});
title('Odd - Early');
set(gca,'XTick',tick_times);

subplot(2,3,6);
imagesc(time_vector, 1:size(sorted_by_ctrl_late_odd,1), sorted_by_ctrl_late_odd);
colormap(slanCM('GnBu')); colorbar;
xline(0, lineStyle{:});
title('Odd - Late');
set(gca,'XTick',tick_times);

sgtitle('Early vs Late PSTHs (sorted by Control peak)');
%% -250 to 250 

%% ---- STEP 4: Early vs Late PSTHs (Control-sorted, limited to -250:250 ms) ----
figure('Position',[100 100 1800 600]);

lineStyle = {'--k','LineWidth',1.5};

% Define the time window
time_window = [-250 250];
time_mask = time_vector >= time_window(1) & time_vector <= time_window(2);
time_subset = time_vector(time_mask);

% Subset all PSTHs to that window
early_ctrl_sub = sorted_by_ctrl_early_ctrl(:, time_mask);
late_ctrl_sub  = sorted_by_ctrl_late_ctrl(:, time_mask);
early_even_sub = sorted_by_ctrl_early_even(:, time_mask);
late_even_sub  = sorted_by_ctrl_late_even(:, time_mask);
early_odd_sub  = sorted_by_ctrl_early_odd(:, time_mask);
late_odd_sub   = sorted_by_ctrl_late_odd(:, time_mask);

tick_times = -250:125:250;  % Tick marks for display
[~, tick_indices] = min(abs(time_subset' - tick_times), [], 1);

subplot(2,3,1);
imagesc(time_subset, 1:size(early_ctrl_sub,1), early_ctrl_sub);
colormap(slanCM('GnBu')); colorbar;
xline(0, lineStyle{:});
title('Ctrl - Early');
set(gca, 'XTick', tick_times);

subplot(2,3,4);
imagesc(time_subset, 1:size(late_ctrl_sub,1), late_ctrl_sub);
colormap(slanCM('GnBu')); colorbar;
xline(0, lineStyle{:});
title('Ctrl - Late');
set(gca, 'XTick', tick_times);

subplot(2,3,2);
imagesc(time_subset, 1:size(early_even_sub,1), early_even_sub);
colormap(slanCM('GnBu')); colorbar;
xline(0, lineStyle{:});
title('Even - Early');
set(gca, 'XTick', tick_times);

subplot(2,3,5);
imagesc(time_subset, 1:size(late_even_sub,1), late_even_sub);
colormap(slanCM('GnBu')); colorbar;
xline(0, lineStyle{:});
title('Even - Late');
set(gca, 'XTick', tick_times);

subplot(2,3,3);
imagesc(time_subset, 1:size(early_odd_sub,1), early_odd_sub);
colormap(slanCM('GnBu')); colorbar;
xline(0, lineStyle{:});
title('Odd - Early');
set(gca, 'XTick', tick_times);

subplot(2,3,6);
imagesc(time_subset, 1:size(late_odd_sub,1), late_odd_sub);
colormap(slanCM('GnBu')); colorbar;
xline(0, lineStyle{:});
title('Odd - Late');
set(gca, 'XTick', tick_times);

sgtitle('Early vs Late PSTHs (sorted by Control peak, -250 to 250 ms)');

%% Firing rate consistency 

% ---- CHECK CONSISTENCY OF EARLY/LATE FIRING ACROSS EVEN/ODD TRIALS ----
% Goal: Determine if units classified as 'early' or 'late' in Control also
% fire early/late in Even and Odd trial subsets.

% --- Step 1: Aggregate peak indices and masks across all sessions ---
all_ctrl_peaks = [];
all_even_peaks = [];
all_odd_peaks  = [];
all_early_mask = [];
all_late_mask  = [];

for d = 1:length(NPdata)
    all_ctrl_peaks = [all_ctrl_peaks; NPdata(d).peak_times_ctrl(:)];
    all_even_peaks = [all_even_peaks; NPdata(d).peak_times_even(:)];
    all_odd_peaks  = [all_odd_peaks;  NPdata(d).peak_times_odd(:)];
    all_early_mask = [all_early_mask; NPdata(d).before_peakv_mask(:)];
    all_late_mask  = [all_late_mask;  NPdata(d).after_peakv_mask(:)];
end

% --- Step 2: Convert peak indices to time (ms) if needed ---
if all_ctrl_peaks(1) <= length(timevec)
    ctrl_peak_ms = timevec(all_ctrl_peaks);
    even_peak_ms = timevec(all_even_peaks);
    odd_peak_ms  = timevec(all_odd_peaks);
else
    % Already in ms
    ctrl_peak_ms = all_ctrl_peaks;
    even_peak_ms = all_even_peaks;
    odd_peak_ms  = all_odd_peaks;
end

% --- Step 3: Define early/late based on 0 ms threshold ---
even_early_mask = even_peak_ms < 0;
odd_early_mask  = odd_peak_ms  < 0;
even_late_mask  = even_peak_ms >= 0;
odd_late_mask   = odd_peak_ms  >= 0;

% --- Step 4: Compute consistency percentages ---
% (How many "early" control units stay early in even/odd trials, and same for "late")
early_units_even_consistent = sum(all_early_mask & even_early_mask) / sum(all_early_mask) * 100;
early_units_odd_consistent  = sum(all_early_mask & odd_early_mask)  / sum(all_early_mask) * 100;
late_units_even_consistent  = sum(all_late_mask & even_late_mask)   / sum(all_late_mask) * 100;
late_units_odd_consistent   = sum(all_late_mask & odd_late_mask)    / sum(all_late_mask) * 100;

fprintf('\n---- FIRING CONSISTENCY CHECK ----\n');
fprintf('Early units firing early again: %.1f%% (even), %.1f%% (odd)\n', ...
    early_units_even_consistent, early_units_odd_consistent);
fprintf('Late units firing late again:   %.1f%% (even), %.1f%% (odd)\n\n', ...
    late_units_even_consistent, late_units_odd_consistent);

% --- Step 5: Optional visualization ---
figure('Position', [100 100 1200 500]);

subplot(1,2,1);
scatter(ctrl_peak_ms, even_peak_ms, 40, all_early_mask, 'filled');
xlabel('Control Peak Time (ms)');
ylabel('Even Peak Time (ms)');
title('Even vs Control Peak Times');
xline(0,'--k'); yline(0,'--k');
axis equal; grid on;
xlim([-250 250]); ylim([-250 250]);
colormap(slanCM('GnBu'));

subplot(1,2,2);
scatter(ctrl_peak_ms, odd_peak_ms, 40, all_early_mask, 'filled');
xlabel('Control Peak Time (ms)');
ylabel('Odd Peak Time (ms)');
title('Odd vs Control Peak Times');
xline(0,'--k'); yline(0,'--k');
axis equal; grid on;
xlim([-250 250]); ylim([-250 250]);
colormap(slanCM('GnBu'));

sgtitle('Peak Timing Consistency: Control vs Even/Odd Trials');



% --- Step 5: Visualization by Early/Late Units ---
figure('Position', [100 100 1200 500]);

% Colors for early vs late
c_early = [0 0.4470 0.7410];  % blue
c_late  = [0.8500 0.3250 0.0980]; % red

% --- Even vs Control ---
subplot(1,2,1); hold on;
scatter(ctrl_peak_ms(all_early_mask), even_peak_ms(all_early_mask), 50, c_early, 'filled', 'o');
scatter(ctrl_peak_ms(all_late_mask),  even_peak_ms(all_late_mask), 50, c_late, 'filled', 's');
xlabel('Control Peak Time (ms)');
ylabel('Even Peak Time (ms)');
title('Even vs Control Peak Times');
xline(0,'--k'); yline(0,'--k');
axis equal; grid on;
xlim([-250 250]); ylim([-250 250]);
legend({'Early Units','Late Units'}, 'Location', 'best');

% --- Odd vs Control ---
subplot(1,2,2); hold on;
scatter(ctrl_peak_ms(all_early_mask), odd_peak_ms(all_early_mask), 50, c_early, 'filled', 'o');
scatter(ctrl_peak_ms(all_late_mask),  odd_peak_ms(all_late_mask), 50, c_late, 'filled', 's');
xlabel('Control Peak Time (ms)');
ylabel('Odd Peak Time (ms)');
title('Odd vs Control Peak Times');
xline(0,'--k'); yline(0,'--k');
axis equal; grid on;
xlim([-250 250]); ylim([-250 250]);
legend({'Early Units','Late Units'}, 'Location', 'best');

sgtitle('Peak Timing Consistency: Control vs Even/Odd Trials');

%% Second attempt with edis 
%% Firing rate consistency 

% ---- CHECK CONSISTENCY OF EARLY/LATE FIRING ACROSS EVEN/ODD TRIALS ----
% Goal: Determine if units classified as 'early' or 'late' in Control also
% fire early/late in Even and Odd trial subsets.

% --- Step 1: Aggregate peak indices and masks across all sessions ---
all_ctrl_peaks = [];
all_even_peaks = [];
all_odd_peaks  = [];
all_early_mask = [];
all_late_mask  = [];

for d = 1:length(NPdata)
    all_ctrl_peaks = [all_ctrl_peaks; NPdata(d).peak_times_ctrl(:)];
    all_even_peaks = [all_even_peaks; NPdata(d).peak_times_even(:)];
    all_odd_peaks  = [all_odd_peaks;  NPdata(d).peak_times_odd(:)];
    all_early_mask = [all_early_mask; NPdata(d).before_peakv_mask(:)];
    all_late_mask  = [all_late_mask;  NPdata(d).after_peakv_mask(:)];
end

% Ensure masks are logical
all_early_mask = logical(all_early_mask);
all_late_mask  = logical(all_late_mask);

% --- Step 2: Convert peak indices to time (ms) if needed ---
if all_ctrl_peaks(1) <= length(timevec)
    ctrl_peak_ms = timevec(all_ctrl_peaks);
    even_peak_ms = timevec(all_even_peaks);
    odd_peak_ms  = timevec(all_odd_peaks);
else
    % Already in ms
    ctrl_peak_ms = all_ctrl_peaks;
    even_peak_ms = all_even_peaks;
    odd_peak_ms  = all_odd_peaks;
end

% --- Step 3: Define early/late based on 0 ms threshold ---
even_early_mask = logical(even_peak_ms < 0);
odd_early_mask  = logical(odd_peak_ms  < 0);
even_late_mask  = logical(even_peak_ms >= 0);
odd_late_mask   = logical(odd_peak_ms  >= 0);

% --- Step 4: Compute consistency percentages ---
early_units_even_consistent = sum(all_early_mask & even_early_mask) / sum(all_early_mask) * 100;
early_units_odd_consistent  = sum(all_early_mask & odd_early_mask)  / sum(all_early_mask) * 100;
late_units_even_consistent  = sum(all_late_mask & even_late_mask)   / sum(all_late_mask) * 100;
late_units_odd_consistent   = sum(all_late_mask & odd_late_mask)    / sum(all_late_mask) * 100;

fprintf('\n---- FIRING CONSISTENCY CHECK ----\n');
fprintf('Early units firing early again: %.1f%% (even), %.1f%% (odd)\n', ...
    early_units_even_consistent, early_units_odd_consistent);
fprintf('Late units firing late again:   %.1f%% (even), %.1f%% (odd)\n\n', ...
    late_units_even_consistent, late_units_odd_consistent);

% --- Step 5: Visualization by Early/Late Units ---
figure('Position', [100 100 1200 500]);

% Colors for early vs late
c_early = rgb('LightSteelBlue');  % blue
c_late  =rgb('DarkGreen') % red

% --- Even vs Control ---
subplot(1,2,1); hold on;
scatter(ctrl_peak_ms(all_early_mask), even_peak_ms(all_early_mask), 50, c_early, 'filled', 'o');
scatter(ctrl_peak_ms(all_late_mask),  even_peak_ms(all_late_mask), 50, c_late, 'filled', 'o');
xlabel('Control Peak Time (ms)');
ylabel('Even Peak Time (ms)');
title('Even vs Control Peak Times');
xline(0,'--k'); yline(0,'--k');
axis equal; grid on;
xlim([-250 250]); ylim([-250 250]);
legend({'Early Units','Late Units'}, 'Location', 'best');

% --- Odd vs Control ---
subplot(1,2,2); hold on;
scatter(ctrl_peak_ms(all_early_mask), odd_peak_ms(all_early_mask), 50, c_early, 'filled', 'o');
scatter(ctrl_peak_ms(all_late_mask),  odd_peak_ms(all_late_mask), 50, c_late, 'filled', 'o');
xlabel('Control Peak Time (ms)');
ylabel('Odd Peak Time (ms)');
title('Odd vs Control Peak Times');
xline(0,'--k'); yline(0,'--k');
axis equal; grid on;
xlim([-250 250]); ylim([-250 250]);
legend({'Early Units','Late Units'}, 'Location', 'best');

sgtitle('Peak Timing Consistency: Control vs Even/Odd Trials');

%%
%% Firing rate consistency 

% ---- CHECK CONSISTENCY OF EARLY/LATE FIRING ACROSS EVEN/ODD TRIALS ----
% Goal: Compare peak times in Even vs Odd trials for units classified as early/late in Control.

% --- Step 1: Aggregate peak indices and masks across all sessions ---
all_even_peaks = [];
all_odd_peaks  = [];
all_early_mask = [];
all_late_mask  = [];

for d = 1:length(NPdata)
    all_even_peaks = [all_even_peaks; NPdata(d).peak_times_even(:)];
    all_odd_peaks  = [all_odd_peaks;  NPdata(d).peak_times_odd(:)];
    all_early_mask = [all_early_mask; NPdata(d).before_peakv_mask(:)];
    all_late_mask  = [all_late_mask;  NPdata(d).after_peakv_mask(:)];
end

% Ensure masks are logical
all_early_mask = logical(all_early_mask);
all_late_mask  = logical(all_late_mask);

% --- Step 2: Convert peak indices to time (ms) if needed ---
if all_even_peaks(1) <= length(timevec)
    even_peak_ms = timevec(all_even_peaks);
    odd_peak_ms  = timevec(all_odd_peaks);
else
    even_peak_ms = all_even_peaks;
    odd_peak_ms  = all_odd_peaks;
end

% --- Step 3: Define early/late based on 0 ms threshold (for reference) ---
even_early_mask = logical(even_peak_ms < 0);
odd_early_mask  = logical(odd_peak_ms  < 0);
even_late_mask  = logical(even_peak_ms >= 0);
odd_late_mask   = logical(odd_peak_ms  >= 0);

% --- Step 4: Optional: Compute consistency percentages ---
early_units_consistent = sum(all_early_mask & even_early_mask & odd_early_mask) / sum(all_early_mask) * 100;
late_units_consistent  = sum(all_late_mask & even_late_mask & odd_late_mask) / sum(all_late_mask) * 100;

fprintf('\n---- FIRING CONSISTENCY CHECK ----\n');
fprintf('Early units firing early in both Even/Odd: %.1f%%\n', early_units_consistent);
fprintf('Late units firing late in both Even/Odd:   %.1f%%\n\n', late_units_consistent);

% --- Step 5: Visualization by Early/Late Units ---
figure('Position', [100 100 600 500]);

% Colors for early vs late
c_early = rgb('LightSteelBlue');  
c_late  = rgb('DarkGreen');

hold on;
scatter(even_peak_ms(all_early_mask), odd_peak_ms(all_early_mask), 50, c_early, 'filled', 'o');
scatter(even_peak_ms(all_late_mask),  odd_peak_ms(all_late_mask), 50, c_late, 'filled', 'o');

xlabel('Even Peak Time (ms)');
ylabel('Odd Peak Time (ms)');
title('Even vs Odd Peak Times (Early/Late Units)');
xline(0,'--k'); yline(0,'--k');
axis equal; grid on;
xlim([-250 250]); ylim([-250 250]);
legend({'Early Units','Late Units'}, 'Location', 'best');
