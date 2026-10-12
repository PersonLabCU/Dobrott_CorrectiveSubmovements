%% ================== SHAM ANALYSIS: NoCSM vs AllCSM (Even/Odd + Early/Late) ==================
% Filters units by NoCSM peak within valid window, stores kept_units for plotting

% Sampling even and odd trials to look out for spurious tiling

cd("C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\ShamAnimals\ReachandCSMMod_allCSMStogther\CSMAll_wtihMaxV")
MyFolderInfo = dir('*CSMAll.mat');

tic
clearvars -except MyFolderInfo

%% ================== SHAM: ALL REACHES PSTH + PEAK TIME ANALYSIS (Even/Odd + Early/Late) ==================

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

% Get all unit indices
All_filtered_modulated_indices = [];

% ---------- params ----------
bin_size       = 10;   % ms
smoothing_bins = 3;    % for computePSTH
valid_peak_range_ms = [-250 250];   % your choice

for day = 1:length(MyFolderInfo)

    fprintf('\n========== SESSION %d / %d ==========\n', day, length(MyFolderInfo));
    load(MyFolderInfo(day).name, 'Reaches','data','filtered_modulated_indices','ExcitRR');

    if ~exist('filtered_modulated_indices','var') || isempty(filtered_modulated_indices)
        fprintf('Day %d: filtered_modulated_indices empty. Skipping...\n', day);
        continue;
    end

    All_filtered_modulated_indices = [All_filtered_modulated_indices; filtered_modulated_indices(:)];

    % ---- pull units ----
    data_units_all = data.units(filtered_modulated_indices);

    % ---- sanity: required field exists? ----
    if ~isfield(data_units_all(1), 'AllCSM_maxvx') || ~isfield(data_units_all(1).max_vx, 'time_bins')
        fprintf('Day %d: Missing AllCSM_maxvx.time_bins. Skipping...\n', day);
        continue;
    end

    % Time vector
    timevec = data_units_all(1).max_vx.time_bins(:);
    valid_peak_range = find(timevec >= valid_peak_range_ms(1) & timevec <= valid_peak_range_ms(2));

    if isempty(valid_peak_range)
        fprintf('Day %d: valid_peak_range empty. Check timevec scale.\n', day);
        continue;
    end

    % ---- Determine which units have peaks within window (based on AllCSM PSTH) ----
    keep_unit = false(length(data_units_all), 1);
    peak_idx_all = NaN(length(data_units_all), 1);

    for t = 1:length(data_units_all)

        if ~isfield(data_units_all(t).AllCSM_maxvx,'PSTH_smoothbin') || isempty(data_units_all(t).max_vx.PSTH_smoothbin)
            continue;
        end

        psth = data_units_all(t).max_vx.PSTH_smoothbin(:)';
        [~, peak_idx] = max(psth);
        peak_idx_all(t) = peak_idx;

        if ismember(peak_idx, valid_peak_range)
            keep_unit(t) = true;
        end
    end

    % Save only the filtered units for future plotting
    data_kept = data_units_all(keep_unit);
    peak_idx_kept = peak_idx_all(keep_unit);

    if isempty(data_kept)
        fprintf('Day %d: No units with peaks in window. Skipping.\n', day);
        continue
    end

    % Convert peak indices to actual times (ms)
    peak_times_ms = timevec(peak_idx_kept);

    % Early vs late relative to 0 ms
    before_peakv_mask = (peak_times_ms < 0);
    after_peakv_mask  = (peak_times_ms >= 0);

    % Preallocate
    nCells = length(data_kept);
    nBins  = length(data_kept(1).max_vx.PSTH_smoothbin);

    psth_master_ctrl = NaN(nCells, nBins);
    psth_master_even = NaN(nCells, nBins);
    psth_master_odd  = NaN(nCells, nBins);

    peak_times_ctrl = NaN(nCells,1);
    peak_times_even = NaN(nCells,1);
    peak_times_odd  = NaN(nCells,1);

    % --- Loop through each cell ---
    for t = 1:nCells

        % --- "Control" PSTH (here = AllCSM, since Sham has AllCSM/NoCSM elsewhere) ---
        psth = data_kept(t).max_vx.PSTH_smoothbin(:)';
        psth_master_ctrl(t,:) = psth;
        [~, peak_times_ctrl(t)] = max(psth);

        % --- Even/Odd PSTHs from AllCSM raster ---
        if isfield(data_kept(t).max_vx, 'Raster') && ~isempty(data_kept(t).max_vx.Raster)

            raster = data_kept(t).max_vx.Raster;  % trials x time
            nTrials = size(raster,1);

            if nTrials >= 2
                even_idx = 2:2:nTrials;
                odd_idx  = 1:2:nTrials;

                raster_even = raster(even_idx,:);
                raster_odd  = raster(odd_idx,:);

                [~, psth_smooth_even, ~] = computePSTH(raster_even, bin_size, smoothing_bins);
                [~, psth_smooth_odd,  ~] = computePSTH(raster_odd,  bin_size, smoothing_bins);

                % keep raw FR (non-normalized), like you prefer
                if max(psth_smooth_even) > 0
                    psth_master_even(t,:) = psth_smooth_even;
                else
                    psth_master_even(t,:) = NaN(1, length(psth_smooth_even));
                end

                if max(psth_smooth_odd) > 0
                    psth_master_odd(t,:)  = psth_smooth_odd;
                else
                    psth_master_odd(t,:)  = NaN(1, length(psth_smooth_odd));
                end

                [~, peak_times_even(t)] = max(psth_master_even(t,:));
                [~, peak_times_odd(t)]  = max(psth_master_odd(t,:));
            end
        end
    end

    % --- Store into NPdata ---
    NPdata(day).peak_times_ctrl   = peak_times_ctrl;
    NPdata(day).peak_times_even   = peak_times_even;
    NPdata(day).peak_times_odd    = peak_times_odd;

    NPdata(day).psth_master_ctrl  = psth_master_ctrl;
    NPdata(day).psth_master_even  = psth_master_even;
    NPdata(day).psth_master_odd   = psth_master_odd;

    NPdata(day).before_peakv_mask = before_peakv_mask;
    NPdata(day).after_peakv_mask  = after_peakv_mask;

    % kept indices (relative to data.units)
    NPdata(day).kept_indices = filtered_modulated_indices(keep_unit);

    % save full unit struct for later plotting
    NPdata(day).kept_units = data_kept;

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

    clear Reaches data filtered_modulated_indices ExcitRR
end


%% SHAM: Plot NoCSM_csmT2 vs AllCSM_csmT2 (Raster + PSTH) and save SVGs
% Generate raster + PSTH for each kept unit and save as SVG (white raster background)

% ---------- OUTPUT FOLDER ----------
output_folder = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\Sham Neural analysis\Modulation index Analysis\CSMT2aligned units";
%output_folder= "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\Sham Neural analysis\Modulation index Analysis\CSMT0aligned units";

if ~exist(output_folder,'dir'); mkdir(output_folder); end

% ---------- PSTH PARAMS ----------
bin_size       = 10;         % ms
smoothing_bins = 3;          % for PSTH smoothing
time_range     = [-1000, 1000];  % ms, full PSTH range (for PSTH timevec)
time_window    = [-1000, 1000];    % ms, plotting window

for day = 1:length(NPdata)

    if ~isfield(NPdata(day), 'kept_units') || isempty(NPdata(day).kept_units)
        continue
    end
    if ~isfield(NPdata(day), 'kept_indices') || isempty(NPdata(day).kept_indices)
        continue
    end

    nUnits = length(NPdata(day).kept_units);

    for u = 1:nUnits

        unit   = NPdata(day).kept_units(u);
        unitID = NPdata(day).kept_indices(u);   % <-- ACTUAL unit index in data.units

        % % % ---------- Pull rasters ---------- T2
        % if ~isfield(unit,'NoCSM_csmT2') || ~isfield(unit.NoCSM_csmT2,'Raster') || isempty(unit.NoCSM_csmT2.Raster)
        %     continue
        % end
        % if ~isfield(unit,'AllCSM_csmT2') || ~isfield(unit.AllCSM_csmT2,'Raster') || isempty(unit.AllCSM_csmT2.Raster)
        %     continue
        % end

        % raster_NoCSM  = unit.NoCSM_csmT2.Raster;    % trials x time
        % raster_AllCSM = unit.AllCSM_csmT2.Raster;   % trials x time

        raster_NoCSM  = unit.NoCSM_pos_thresh.Raster;    % trials x time
        raster_AllCSM = unit.AllCSM_pos_thresh.Raster;   % trials x time


        % ---------- Time vector for raster ----------
        % if isfield(unit.NoCSM_csmT2,'time_bins') && ~isempty(unit.NoCSM_csmT2.time_bins)
        %     timevec_raster = linspace(unit.NoCSM_csmT2.time_bins(1), unit.NoCSM_csmT2.time_bins(end), size(raster_NoCSM,2));
        % else
        %     % fallback if time_bins missing
        %     timevec_raster = linspace(time_range(1), time_range(2), size(raster_NoCSM,2));
        % end

            % ---------- Time vector for raster ----------
        if isfield(unit.NoCSM_pos_thresh,'time_bins') && ~isempty(unit.NoCSM_pos_thresh.time_bins)
            timevec_raster = linspace(unit.NoCSM_pos_thresh.time_bins(1), unit.NoCSM_pos_thresh.time_bins(end), size(raster_NoCSM,2));
        else
            % fallback if time_bins missing
            timevec_raster = linspace(time_range(1), time_range(2), size(raster_NoCSM,2));
        end



        %    % ---------- Pull rasters ----------
        % if ~isfield(unit,'NoCSM_csmT0') || ~isfield(unit.NoCSM_csmT0,'Raster') || isempty(unit.NoCSM_csmT0.Raster)
        %     continue
        % end
        % if ~isfield(unit,'AllCSM_csmT0') || ~isfield(unit.AllCSM_csmT0,'Raster') || isempty(unit.AllCSM_csmT0.Raster)
        %     continue
        % end
        % 
        % raster_NoCSM  = unit.NoCSM_csmT0.Raster;    % trials x time
        % raster_AllCSM = unit.AllCSM_csmT0.Raster;   % trials x time
        % 
        % % ---------- Time vector for raster ----------
        % if isfield(unit.NoCSM_csmT0,'time_bins') && ~isempty(unit.NoCSM_csmT0.time_bins)
        %     timevec_raster = linspace(unit.NoCSM_csmT0.time_bins(1), unit.NoCSM_csmT0.time_bins(end), size(raster_NoCSM,2));
        % else
        %     % fallback if time_bins missing
        %     timevec_raster = linspace(time_range(1), time_range(2), size(raster_NoCSM,2));
        % end

        % ---------- Compute PSTHs (with smoothing) ----------
        [~, psth_nocsm, ~]  = computePSTH(raster_NoCSM,  bin_size, smoothing_bins);
        [~, psth_allcsm, ~] = computePSTH(raster_AllCSM, bin_size, smoothing_bins);

        % ---------- Construct consistent PSTH time vector ----------
        num_time_bins = size(psth_nocsm, 2);
        timevec_psth  = linspace(time_range(1), time_range(2), num_time_bins);

        % ---------- Apply time window mask for plotting ----------
        time_mask        = timevec_psth >= time_window(1) & timevec_psth <= time_window(2);
        timevec_plot     = timevec_psth(time_mask);
        psth_nocsm_plot  = psth_nocsm(time_mask);
        psth_allcsm_plot = psth_allcsm(time_mask);

        % ---------- Raster points: NoCSM ----------
        [trial_idx_NoCSM, col_idx_NoCSM] = find(raster_NoCSM);
        spike_times_NoCSM = timevec_raster(col_idx_NoCSM);
        nTrials_NoCSM     = size(raster_NoCSM, 1);

        % ---------- Raster points: AllCSM ----------
        [trial_idx_AllCSM, col_idx_AllCSM] = find(raster_AllCSM);
        spike_times_AllCSM = timevec_raster(col_idx_AllCSM);
        nTrials_AllCSM     = size(raster_AllCSM, 1);

        % ---------- Plot ----------
        fig = figure('Visible', 'off', 'Position', [100, 100, 500, 500]);

        % --- NoCSM Raster ---
        subplot(3, 1, 1); hold on;
        set(gca, 'Color', 'w');
        scatter(spike_times_NoCSM, trial_idx_NoCSM, 1.5, 'k', 'filled');

        
        % xlim(time_window);
        % ylim([0 nTrials_NoCSM+1]);

        xlim(time_window);
        xticks(time_window(1):200:time_window(2));   % <-- 200 ms ticks
        ylim([0 nTrials_total+1]);
        xline(0,'--'); xline(50,'--');

        xlabel('Time from CSM T2 (ms)');
        %xlabel('Time from CSM T0 (ms)');
        ylabel('Trials');
        title(sprintf('Day %d | UnitID %d (u=%d) | NoCSM Raster', day, unitID, u));
        box off; set(gca, 'TickDir', 'out');

        % --- AllCSM Raster ---
        subplot(3, 1, 2); hold on;
        set(gca, 'Color', 'w');
        scatter(spike_times_AllCSM, trial_idx_AllCSM, 1.5, 'r', 'filled');
        % xlim(time_window);
        % ylim([0 nTrials_AllCSM+1]);
        xlim(time_window);
        xticks(time_window(1):200:time_window(2));   % <-- 200 ms ticks
        ylim([0 nTrials_total+1]);
        xline(0,'--'); xline(50,'--');
        xlabel('Time from CSM T2 (ms)');
        %xlabel('Time from CSM T0 (ms)');
        ylabel('Trials');
        title(sprintf('Day %d | UnitID %d (u=%d) | AllCSM Raster', day, unitID, u));
        box off; set(gca, 'TickDir', 'out');

        % --- PSTH ---
        subplot(3, 1, 3); hold on;
        plot(timevec_plot, psth_nocsm_plot,  'k', 'LineWidth', 0.5);
        plot(timevec_plot, psth_allcsm_plot, 'r', 'LineWidth', 0.5);
        xlim(time_window);
        xticks(time_window(1):200:time_window(2));   % <-- 200 ms ticks
        xline(0,'--');
        xline(50,'--');
        % xlim(time_window);
        % xline(0,'--'); xline(50,'--');
        xlabel('Time from threshold crossing(ms)');
        %xlabel('Time from CSM T2 (ms)');
        %xlabel('Time from CSM T0 (ms)');
        ylabel('Firing Rate (spikes/s)');
        title(sprintf('Day %d | UnitID %d (u=%d) | PSTH (NoCSM vs AllCSM)', day, unitID, u));
        box off; set(gca, 'TickDir', 'out');

        % --- Save figure as SVG (include UnitID) ---
        filename = fullfile(output_folder, sprintf('Day%d_UnitID%d_u%d_NoCSM_vs_AllCSM.svg', day, unitID, u));
        saveas(fig, filename);
        close(fig);

    end
end


%% SHAM: Plot NoCSM vs AllCSM (SINGLE raster + PSTH) and save SVGs
% Raster: NoCSM trials first (black), then AllCSM trials (Dark Magenta), all on SAME axis.

% ---------- OUTPUT FOLDER ----------
%output_folder = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\Sham Neural analysis\Heatmaps and Ex units\PosThreshAligned";
output_folder="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\Sham Neural analysis\Heatmaps and Ex units\PosthreshAligned_1000";

if ~exist(output_folder,'dir'); mkdir(output_folder); end

% ---------- PSTH PARAMS ----------
bin_size       = 10;            % ms
smoothing_bins = 3;             % for PSTH smoothing
time_range     = [-1000, 1000]; % ms, full PSTH range (for PSTH timevec)
%time_window    = [-500, 500];   % ms, plotting window
time_window    = [-1000, 1000];   % ms, plotting window

% ---------- COLORS ----------
cNo   = [0 0 0];        % black
cCSM  = [0.55 0 0.55];  % Dark Magenta (RGB-ish). Tweak if you want darker/lighter.

for day = 1:length(NPdata)

    if ~isfield(NPdata(day), 'kept_units') || isempty(NPdata(day).kept_units)
        continue
    end
    if ~isfield(NPdata(day), 'kept_indices') || isempty(NPdata(day).kept_indices)
        continue
    end

    nUnits = length(NPdata(day).kept_units);

    for u = 1:nUnits

        unit   = NPdata(day).kept_units(u);
        unitID = NPdata(day).kept_indices(u);   % <-- ACTUAL unit index in data.units

        % ---------- Pull rasters ----------
        if ~isfield(unit,'NoCSM_pos_thresh') || ~isfield(unit.NoCSM_pos_thresh,'Raster') || isempty(unit.NoCSM_pos_thresh.Raster)
            continue
        end
        if ~isfield(unit,'AllCSM_pos_thresh') || ~isfield(unit.AllCSM_pos_thresh,'Raster') || isempty(unit.AllCSM_pos_thresh.Raster)
            continue
        end

        raster_NoCSM  = unit.NoCSM_pos_thresh.Raster;    % trials x time
        raster_AllCSM = unit.AllCSM_pos_thresh.Raster;   % trials x time

        % ---------- Time vector for raster ----------
        if isfield(unit.NoCSM_pos_thresh,'time_bins') && ~isempty(unit.NoCSM_pos_thresh.time_bins)
            timevec_raster = linspace(unit.NoCSM_pos_thresh.time_bins(1), unit.NoCSM_pos_thresh.time_bins(end), size(raster_NoCSM,2));
        else
            timevec_raster = linspace(time_range(1), time_range(2), size(raster_NoCSM,2));
        end

        % ---------- Compute PSTHs (with smoothing) ----------
        [~, psth_nocsm, ~]  = computePSTH(raster_NoCSM,  bin_size, smoothing_bins);
        [~, psth_allcsm, ~] = computePSTH(raster_AllCSM, bin_size, smoothing_bins);

        % ---------- Construct consistent PSTH time vector ----------
        num_time_bins = size(psth_nocsm, 2);
        timevec_psth  = linspace(time_range(1), time_range(2), num_time_bins);

        % ---------- Apply time window mask for plotting ----------
        time_mask        = timevec_psth >= time_window(1) & timevec_psth <= time_window(2);
        timevec_plot     = timevec_psth(time_mask);
        psth_nocsm_plot  = psth_nocsm(time_mask);
        psth_allcsm_plot = psth_allcsm(time_mask);

        % ---------- Raster points ----------
        [trial_idx_NoCSM,  col_idx_NoCSM]  = find(raster_NoCSM);
        spike_times_NoCSM = timevec_raster(col_idx_NoCSM);
        nTrials_NoCSM     = size(raster_NoCSM, 1);

        [trial_idx_AllCSM, col_idx_AllCSM] = find(raster_AllCSM);
        spike_times_AllCSM = timevec_raster(col_idx_AllCSM);
        nTrials_AllCSM     = size(raster_AllCSM, 1);

        % Shift AllCSM trial indices so they appear AFTER NoCSM on the same raster
        trial_idx_AllCSM_shift = trial_idx_AllCSM + nTrials_NoCSM;

        nTrials_total = nTrials_NoCSM + nTrials_AllCSM;

        % ---------- Plot ----------
        fig = figure('Visible', 'off', 'Position', [100, 100, 500, 500]);

        % --- SINGLE Raster (NoCSM then AllCSM) ---
        subplot(3, 1, [1 2]); hold on;
        set(gca, 'Color', 'w');

        scatter(spike_times_NoCSM,  trial_idx_NoCSM,        1.5, cNo,  'filled');
        scatter(spike_times_AllCSM, trial_idx_AllCSM_shift, 1.5, cCSM, 'filled');

        % xlim(time_window);
        % ylim([0 nTrials_total+1]);
        % xline(0,'--'); xline(50,'--');

        xlim(time_window);
        xticks(time_window(1):250:time_window(2));   % <-- 200 ms ticks
        ylim([0 nTrials_total+1]);
        xline(0,'--'); xline(50,'--');
        xlabel('Time from threshold crossing (ms)');
        %xlabel('Time from CSM T0 (ms)');
        ylabel('Trials');

        % Optional divider line between NoCSM and CSM blocks
        yline(nTrials_NoCSM + 0.5, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.75);

        title(sprintf('Day %d | UnitID %d (u=%d) | Raster (NoCSM black, AllCSM magenta)', day, unitID, u));
        box off; set(gca, 'TickDir', 'out');

        % --- PSTH ---
        subplot(3, 1, 3); hold on;
        plot(timevec_plot, psth_nocsm_plot,  'Color', cNo,  'LineWidth', 0.5);
        plot(timevec_plot, psth_allcsm_plot, 'Color', cCSM, 'LineWidth', 0.5);
        % xlim(time_window);
        % xline(0,'--'); xline(50,'--');

        xlim(time_window);
        xticks(time_window(1):250:time_window(2));   % <-- 200 ms ticks
        xline(0,'--');
        xline(50,'--');

        xlabel('Time from threshold crossing(ms)');
        %xlabel('Time from CSM T0 (ms)');
        ylabel('Firing Rate (spikes/s)');
        title(sprintf('Day %d | UnitID %d (u=%d) | PSTH (NoCSM vs AllCSM)', day, unitID, u));
        box off; set(gca, 'TickDir', 'out');

        % --- Save figure as SVG (include UnitID) ---
        filename = fullfile(output_folder, sprintf('Day%d_UnitID%d_u%d_NoCSM_then_AllCSM_SingleRaster.svg', day, unitID, u));
        saveas(fig, filename);
        close(fig);

    end
end