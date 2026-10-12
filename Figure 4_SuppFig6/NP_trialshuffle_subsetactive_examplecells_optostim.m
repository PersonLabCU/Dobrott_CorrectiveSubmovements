% Sampling even and odd trials to look out for spurious tiling

cd("C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\StimAnimals\BlockSeparated_CSMvsNoCSM_inducedalignment")
MyFolderInfo = dir('*Block_CSMvsNoCSM.mat'); 


tic
clearvars -except MyFolderInfo

%% ALL REACHES PSTH + PEAK TIME ANALYSIS (Even/Odd + Early/Late)

NPdata = struct;

PsthMasterCtrl = [];
PsthMasterEven = [];
PsthMasterOdd  = [];
PsthMasterOpto  = [];

PsthMaster_Before_Ctrl = [];
PsthMaster_After_Ctrl  = [];
PsthMaster_Before_Even = [];
PsthMaster_After_Even  = [];
PsthMaster_Before_Odd  = [];
PsthMaster_After_Odd   = [];

PsthMaster_Before_Opto = [];
PsthMaster_After_Opto= [];

% Get all unit indices 
All_Salt_index = [];
All_filtered_modulated_indices = [];
All_opto_index = [];

for day = 1:length(MyFolderInfo)
    

     fprintf('\n========== SESSION %d / %d ==========\n', day, length(MyFolderInfo));
     load(MyFolderInfo(day).name, 'Reaches','data','filtered_modulated_indices','ExcitRR_BL','ExcitRR_AD','ExcitRR_WOE','ExcitRR_WO','Salt_index','opto_index');

     if isempty(Salt_index)
         fprintf('Day %d: Salt_index is empty. Skipping...\n', day);
         continue;
     end
     % 
     % if isempty(opto_index)
     %     fprintf('Day %d: Opto_index is empty. Skipping...\n', day);
     %     continue;
     % end


     % --- concatenate per-session index values ---
     if exist('Salt_index','var') && ~isempty(Salt_index)
         All_Salt_index = [All_Salt_index; Salt_index(:)];
     end

     if exist('filtered_modulated_indices','var') && ~isempty(filtered_modulated_indices)
         All_filtered_modulated_indices = [All_filtered_modulated_indices; filtered_modulated_indices(:)];
     end

     if exist('opto_index','var') && ~isempty(opto_index)
         All_opto_index = [All_opto_index; opto_index(:)];
     end



    data_opto_all = data.units(Salt_index);
    %data_opto_all = data.units(filtered_modulated_indices);
    %data_opto_all = data.units(opto_index);

    % Time vector
    timevec = data_opto_all(1).BLCSM_max_vx.time_bins;
    valid_peak_range = find(timevec >= -250 & timevec <= 250);
    

    % Determine which units have peaks within window
    keep_unit = false(length(data_opto_all), 1);
    peak_times_ctrl_all = NaN(length(data_opto_all), 1);

    for t = 1:length(data_opto_all)
        psth = data_opto_all(t).BLCSM_max_vx.PSTH_smoothbin;
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
    nBins = length(data_opto(1).BLCSM_max_vx.PSTH_smoothbin);

    psth_master_ctrl = NaN(nCells, nBins);
    psth_master_even = NaN(nCells, nBins);
    psth_master_odd  = NaN(nCells, nBins);
    psth_opto=NaN(nCells, nBins);

    peak_times_ctrl = NaN(nCells,1);
    peak_times_even = NaN(nCells,1);
    peak_times_odd  = NaN(nCells,1);
    peak_times_oopto  = NaN(nCells,1);

    % --- Loop through each cell ---
    for t = 1:nCells
        % --- Control PSTH ---
        psth = data_opto(t).BLCSM_max_vx.PSTH_smoothbin;
        psth_master_ctrl(t,:) = psth;
        [~, peak_times_ctrl(t)] = max(psth);
        
        %OptoPSTH
        raster_opto=data_opto(t).eventopto.Raster;
        bin_size = 10; % ms
        smoothing_bins = 3;

        [~, psth_smooth_opto, time_bins_opto] = computePSTH(raster_opto, bin_size, smoothing_bins);
        psth_master_opto(t,:)  = psth_smooth_opto;
        [~, peak_times_opto(t)] = max(psth_smooth_opto);

        % --- Even/Odd PSTHs ---
        if isfield(data_opto(t).BLCSM_max_vx, 'Raster') && ~isempty(data_opto(t).BLCSM_max_vx.Raster)
            raster = data_opto(t).BLCSM_max_vx.Raster; % trials x time
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
               % psth_master_even(t,:) = psth_smooth_even / max(psth_smooth_even);
                psth_master_even(t,:) = psth_smooth_even; % non normlized FR 
            else
               psth_master_even(t,:) = NaN(1,length(psth_smooth_even));
          
            end

            if max(psth_smooth_odd) > 0
                %psth_master_odd(t,:)  = psth_smooth_odd / max(psth_smooth_odd);
                psth_master_odd(t,:)  = psth_smooth_odd;
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
    NPdata(day).peak_times_opto = peak_times_opto;
    NPdata(day).psth_master_ctrl = psth_master_ctrl;
    NPdata(day).psth_master_even = psth_master_even;
    NPdata(day).psth_master_odd  = psth_master_odd;
    NPdata(day).psth_master_opto  = psth_master_opto;
    NPdata(day).before_peakv_mask = before_peakv_mask;
    NPdata(day).after_peakv_mask  = after_peakv_mask;
    NPdata(day).kept_indices = Salt_index(keep_unit);
    %NPdata(day).kept_indices = filtered_modulated_indices(keep_unit);
    %NPdata(day).kept_indices = opto_index(keep_unit);

    % --- NEW: save full unit struct for later plotting ---
    NPdata(day).kept_units = data_opto;
   
    NPdata(day).time_bins_opto = time_bins_opto;

    % --- Concatenate across days ---
    PsthMasterCtrl = vertcat(PsthMasterCtrl, psth_master_ctrl);
    PsthMasterEven = vertcat(PsthMasterEven, psth_master_even);
    PsthMasterOdd  = vertcat(PsthMasterOdd,  psth_master_odd);
    PsthMasterOpto = vertcat(PsthMasterOpto, psth_master_opto);

    PsthMaster_Before_Ctrl = vertcat(PsthMaster_Before_Ctrl, psth_master_ctrl(before_peakv_mask,:));
    PsthMaster_After_Ctrl  = vertcat(PsthMaster_After_Ctrl,  psth_master_ctrl(after_peakv_mask,:));

    PsthMaster_Before_Even = vertcat(PsthMaster_Before_Even, psth_master_even(before_peakv_mask,:));
    PsthMaster_After_Even  = vertcat(PsthMaster_After_Even,  psth_master_even(after_peakv_mask,:));

    PsthMaster_Before_Odd  = vertcat(PsthMaster_Before_Odd,  psth_master_odd(before_peakv_mask,:));
    PsthMaster_After_Odd   = vertcat(PsthMaster_After_Odd,   psth_master_odd(after_peakv_mask,:));

    PsthMaster_Before_Opto  = vertcat(PsthMaster_Before_Opto,  psth_master_opto(before_peakv_mask,:));
    PsthMaster_After_Opto   = vertcat(PsthMaster_After_Opto,   psth_master_opto(after_peakv_mask,:));

    fprintf('Day %d: Added %d units (%d early, %d late)\n', day, nCells, sum(before_peakv_mask), sum(after_peakv_mask));
    clear Reaches data filtered_modulated_indices ExcitRR_BL ExcitRR_AD ExcitRR_WOE ExcitRR_WO Salt_index
end


%% NEW EDIT BEHAVIOR AS WELL 
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
%output_folder = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\ExampleUnits_optotagged";
%output_folder="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\ExampleUnits_filteredmodulated";
output_folder="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\ExampleUnits_Opto_index";

oddcolor = rgb('DarkRed');
evencolor = rgb('DarkGreen');

bin_size = 10;         % ms
smoothing_bins = 3;    % for PSTH smoothing
%time_range = [-1000, 1000];  % ms, full PSTH range
%time_window = [-500, 500];   % ms, plotting window
time_range = [-500, 500];  % ms, full PSTH range for optotagged units 
time_window = [-250, 250];   % ms, plotting window

for day = 1:length(NPdata)
    if ~isfield(NPdata(day), 'kept_units') || isempty(NPdata(day).kept_units)
        continue
    end

    nUnits = length(NPdata(day).kept_units);

    for u = 1:nUnits
        unit = NPdata(day).kept_units(u);
        raster = unit.eventopto.Raster;  % trials x time

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

        timevec_raster = linspace(NPdata(1).time_bins_opto(1), NPdata(1).time_bins_opto(end), size(raster,2));

        nTrials = size(raster, 1);
        even_idx = 2:2:nTrials;
        odd_idx  = 1:2:nTrials;

        raster_even = raster(even_idx, :);
        raster_odd  = raster(odd_idx, :);

        % --- Compute PSTHs (with smoothing) ---
        [~, psth_all, ~]  = computePSTH(raster, bin_size, smoothing_bins);
        [~, psth_even, ~] = computePSTH(raster_even, bin_size, smoothing_bins);
        [~, psth_odd,  ~] = computePSTH(raster_odd,  bin_size, smoothing_bins);
        % 
        % psth_raw,psth_all = Psth1opto(raster,2); % was 5 before 
        % psth_raw,psth_even = Psth1opto(raster_even,2); % was 5 before
        % psth_raw,psth_odd = Psth1opto(raster_odd,2); % was 5 before
        % 

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
        scatter(spike_times, trial_idx, 1, 'k', 'filled'); 
        xlim(time_window);
        ylim([0 nTrials+1]);
        xlabel('Time from stimulation (ms)');
        ylabel('Trials');
        title(sprintf('Day %d, Unit %d Raster (%s)', day, u, unit_type));
        box off;
        set(gca, 'TickDir', 'out');
        % % --- Plot Raster ---
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
        xlabel('Time from stimulation (ms)');
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

%% Correct PSTH alignment using smoothed PSTH output (scatter raster preserved)

output_folder = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\ExampleUnits_optotagged\1ms time bins";

oddcolor = rgb('DarkRed'); 
evencolor = rgb('DarkGreen');

time_window = [-250 250];   % plotting window

for day = 1:length(NPdata)

    if ~isfield(NPdata(day),'kept_units') || isempty(NPdata(day).kept_units)
        continue
    end

    nUnits = length(NPdata(day).kept_units);

    for u = 1:nUnits

        unit = NPdata(day).kept_units(u);
        raster = unit.eventopto.Raster;

        if isempty(raster)
            continue
        end

        % assign type
        if NPdata(day).before_peakv_mask(u)
            unit_type = 'Early';
        elseif NPdata(day).after_peakv_mask(u)
            unit_type = 'Late';
        else
            unit_type = 'Uncategorized';
        end

        % --- RASTER INFO ---
        nTrials = size(raster,1);

        % odd/even split
        odd_idx  = 1:2:nTrials;
        even_idx = 2:2:nTrials;

        raster_odd  = raster(odd_idx,:);
        raster_even = raster(even_idx,:);

        % compute PSTH (smoothed) for all / even / odd trials
        [psth_raw, psth_all_smooth]  = Psth1opto(raster,5);   % smoothed PSTH
        [psth_raw, psth_even_smooth] = Psth1opto(raster_even,5);
        [psth_raw, psth_odd_smooth]  = Psth1opto(raster_odd,5);

        % --- Construct time vector from smoothed PSTH output ---
        nBins_all  = length(psth_all_smooth);
        nBins_even = length(psth_even_smooth);
        nBins_odd  = length(psth_odd_smooth);

        timevec_all  = linspace(-500, 499, nBins_all);
        timevec_even = linspace(-500, 499, nBins_even);
        timevec_odd  = linspace(-500, 499, nBins_odd);

        % --- Apply time window mask for plotting ---
        mask_all  = timevec_all  >= time_window(1) & timevec_all  <= time_window(2);
        mask_even = timevec_even >= time_window(1) & timevec_even <= time_window(2);
        mask_odd  = timevec_odd  >= time_window(1) & timevec_odd  <= time_window(2);

        % cropped PSTHs and time vectors
        tplot_all  = timevec_all(mask_all);
        tplot_even = timevec_even(mask_even);
        tplot_odd  = timevec_odd(mask_odd);

        psth_all_p  = psth_all_smooth(mask_all);
        psth_even_p = psth_even_smooth(mask_even);
        psth_odd_p  = psth_odd_smooth(mask_odd);

        % --- Extract spikes for scatter raster (use original full raster) ---
        [trial_idx, col_idx] = find(raster);
        spike_times = linspace(-500, 499, size(raster,2));
        spike_times = spike_times(col_idx);

        % --- FIGURE ---
        fig = figure('Visible','off','Position',[200 200 600 500]);

        % -----------------------
        % RASTER (SCATTER)
        % -----------------------
        subplot(2,1,1); hold on
        scatter(spike_times, trial_idx, 1, 'k', 'filled');
        xlim(time_window);
        ylim([0 nTrials+1]);
        xlabel('Time (ms)');
        ylabel('Trials');
        title(sprintf('Day %d Unit %d Raster (%s)', day, u, unit_type));
        box off; set(gca,'TickDir','out')

        % -----------------------
        % PSTH PLOTS
        % -----------------------
        subplot(2,1,2); hold on
        plot(tplot_all,  psth_all_p,  'k',          'LineWidth',1.6);
        plot(tplot_even, psth_even_p, 'Color', evencolor, 'LineWidth',1.3);
        plot(tplot_odd,  psth_odd_p,  'Color', oddcolor,  'LineWidth',1.3);
        xlim(time_window);
        xlabel('Time (ms)');
        ylabel('Firing Rate (Hz)');
        legend({'All','Even','Odd'},'Location','best');
        title(sprintf('Day %d Unit %d PSTH (%s)', day, u, unit_type));
        box off; set(gca,'TickDir','out')

        % save
        filename = fullfile(output_folder, sprintf('Day%d_Unit%d_%s.svg', day, u, unit_type));
        saveas(fig, filename);
        close(fig);

    end
end

%% PLoting PST hsimilar to behfore 
output_folder="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\ExampleUnits_filteredmodulated";

%output_folder = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\ExampleUnits_optotagged\1ms time bins\Less smooth";

oddcolor = rgb('DarkRed'); 
evencolor = rgb('DarkGreen');

time_window = [-250 250];   % plotting window

for day = 1:length(NPdata)

    if ~isfield(NPdata(day),'kept_units') || isempty(NPdata(day).kept_units)
        continue
    end

    nUnits = length(NPdata(day).kept_units);

    for u = 1:nUnits

        unit = NPdata(day).kept_units(u);
        raster = unit.eventopto.Raster;

        if isempty(raster)
            continue
        end

        % assign type
        if NPdata(day).before_peakv_mask(u)
            unit_type = 'Early';
        elseif NPdata(day).after_peakv_mask(u)
            unit_type = 'Late';
        else
            unit_type = 'Uncategorized';
        end

        % --- RASTER INFO ---
        nTrials = size(raster,1);

        % odd/even split
        odd_idx  = 1:2:nTrials;
        even_idx = 2:2:nTrials;

        raster_odd  = raster(odd_idx,:);
        raster_even = raster(even_idx,:);

        % compute PSTH (smoothed) for all / even / odd trials
        % [psth_raw, psth_all_smooth]  = Psth1opto(raster,5);   % smoothed PSTH
        % [psth_raw, psth_even_smooth] = Psth1opto(raster_even,5);
        % [psth_raw, psth_odd_smooth]  = Psth1opto(raster_odd,5);

        [psth_raw, psth_all_smooth]  = Psth1opto(raster,2);   % smoothed PSTH
        [psth_raw, psth_even_smooth] = Psth1opto(raster_even,2);
        [psth_raw, psth_odd_smooth]  = Psth1opto(raster_odd,2);

        % --- Construct time vector from smoothed PSTH output ---
        nBins_all  = length(psth_all_smooth);
        nBins_even = length(psth_even_smooth);
        nBins_odd  = length(psth_odd_smooth);

        timevec_all  = linspace(-500, 499, nBins_all);
        timevec_even = linspace(-500, 499, nBins_even);
        timevec_odd  = linspace(-500, 499, nBins_odd);

        % --- Apply time window mask for plotting ---
        mask_all  = timevec_all  >= time_window(1) & timevec_all  <= time_window(2);
        mask_even = timevec_even >= time_window(1) & timevec_even <= time_window(2);
        mask_odd  = timevec_odd  >= time_window(1) & timevec_odd  <= time_window(2);

        % cropped PSTHs and time vectors
        tplot_all  = timevec_all(mask_all);
        tplot_even = timevec_even(mask_even);
        tplot_odd  = timevec_odd(mask_odd);

        psth_all_p  = psth_all_smooth(mask_all);
        psth_even_p = psth_even_smooth(mask_even);
        psth_odd_p  = psth_odd_smooth(mask_odd);

        % --- Extract spikes for scatter raster (use original full raster) ---
        [trial_idx, col_idx] = find(raster);
        spike_times = linspace(-500, 499, size(raster,2));
        spike_times = spike_times(col_idx);

        % --- FIGURE ---
        fig = figure('Visible','off','Position',[200 200 600 500]);

        % -----------------------
        % RASTER (SCATTER)
        % -----------------------
        subplot(2,1,1); hold on
        scatter(spike_times, trial_idx, 1, 'k', 'filled');
        xlim(time_window);
        ylim([0 nTrials+1]);
        xlabel('Time from stimulation (ms)');
        ylabel('Trials');
        title(sprintf('Day %d Unit %d Raster (%s)', day, u, unit_type));
        box off; set(gca,'TickDir','out')
        % 
        % % -----------------------
        % % PSTH PLOTS
        % % -----------------------
        subplot(2,1,2); hold on
        plot(tplot_all,  psth_all_p,  'k',          'LineWidth',1.6);
        plot(tplot_even, psth_even_p, 'Color', evencolor, 'LineWidth',1.3);
        plot(tplot_odd,  psth_odd_p,  'Color', oddcolor,  'LineWidth',1.3);
        xlim(time_window);
        xlabel('Time from stimulation (ms)');
        ylabel('Firing Rate (Hz)');
        legend({'All','Even','Odd'},'Location','best');
        title(sprintf('Day %d Unit %d PSTH (%s)', day, u, unit_type));
        box off; set(gca,'TickDir','out')

        % subplot(2,1,2); hold on    
        % epo=250;
        % td1 = length(psth_raw)-length(psth_all_smooth)+1;
        % td2 = td1/2;
        % %plot((-epo:epo),psth_smooth(((1000-epo)-td2):((2000-(1000-epo))-td2)),'Color',[.2 .4 .8],'LineWidth',3)
        % plot((-epo:epo),psth_all_smooth(((500-epo)-td2):((1000-(500-epo))-td2)),'Color','k','LineWidth',3)
        % plot((-epo:epo),psth_even_smooth(((500-epo)-td2):((1000-(500-epo))-td2)),'Color',evencolor,'LineWidth',3)
        % plot((-epo:epo),psth_odd_smooth(((500-epo)-td2):((1000-(500-epo))-td2)),'Color',oddcolor,'LineWidth',3)
        % xlim([-epo,epo])
        % xlabel('Time from Stim')
        % ylabel('FR (Hz)')
        % title(sprintf('Day %d Unit %d PSTH (%s)', day, u, unit_type));
        % box off; set(gca,'TickDir','out')

        % save
        filename = fullfile(output_folder, sprintf('Day%d_Unit%d_%s.svg', day, u, unit_type));
        saveas(fig, filename);
        close(fig);

    end
end



%% Behavior aligned

%output_folder = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\Behvaior Aligned\Posthresh";
%output_folder="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\NoStim_CSMModulation\IndividualUnits_Csmt2";
%output_folder="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\CSMvsNoCSM_BLandWO\Baseline Trials";
output_folder="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\CSMvsNoCSM_BLandWO\Washout";

oddcolor = rgb('DarkRed'); 
evencolor = rgb('DarkGreen');

time_window = [-500 500];   % plotting window

for day = 1:length(NPdata)

    if ~isfield(NPdata(day),'kept_units') || isempty(NPdata(day).kept_units)
        continue
    end

    nUnits = length(NPdata(day).kept_units);

    for u = 1:nUnits

        unit = NPdata(day).kept_units(u);
        % raster_BL = unit.BLCSM_pos_thresh.Raster;
        % raster_AD = unit.ADCSM_pos_thresh.Raster;
        % 
        % raster_BL = unit.StimoffNoCSM_csmT2.Raster;
        % raster_AD = unit.StimoffCSM_csmT2.Raster;

        % raster_BL = unit.BLNoCSM_csmT2.Raster;
        % raster_AD = unit.BLCSM_csmT2.Raster;

        raster_BL = unit.WONoCSM_csmT2.Raster;
        raster_AD = unit.WOCSM_csmT2.Raster;

        if isempty(raster)
            continue
        end

        % assign type
        if NPdata(day).before_peakv_mask(u)
            unit_type = 'Early';
        elseif NPdata(day).after_peakv_mask(u)
            unit_type = 'Late';
        else
            unit_type = 'Uncategorized';
        end

       

        % compute PSTH (smoothed) for all / even / odd trials
        [psth_raw, psth_bl_smooth]  = Psth1(raster_BL,5);   % smoothed PSTH
        [psth_raw, psth_ad_smooth] = Psth1(raster_AD,5);
        %[psth_raw, psth_ad_smooth] = Psth1opto(raster_AD,5);
  

        % --- Construct time vector from smoothed PSTH output ---
        nBins_all  = length(psth_bl_smooth);
        nBins_even = length(psth_ad_smooth);
 

        timevec_all  = linspace(-1000, 999, nBins_all);
        timevec_even = linspace(-1000, 999, nBins_even);
  

        % --- Apply time window mask for plotting ---
        mask_all  = timevec_all  >= time_window(1) & timevec_all  <= time_window(2);
        mask_even = timevec_even >= time_window(1) & timevec_even <= time_window(2);
       

        % cropped PSTHs and time vectors
        tplot_all  = timevec_all(mask_all);
        tplot_even = timevec_even(mask_even);
       
        psth_all_p  = psth_bl_smooth(mask_all);
        psth_even_p = psth_ad_smooth(mask_even);
     

        % --- Extract spikes for scatter raster (use original full raster) ---
        [trial_idx, col_idx] = find(raster_BL);
        spike_times = linspace(-1000, 999, size(raster_BL,2));
        spike_times = spike_times(col_idx);

        % --- FIGURE ---
        fig = figure('Visible','off','Position',[200 200 600 500]);

        % -----------------------
        % RASTER (SCATTER)
        % -----------------------
        subplot(3,1,1); hold on % baseline trials 
        scatter(spike_times, trial_idx, 1, 'k', 'filled');
        xlim(time_window);
        ylim([0 nTrials+1]);
        %xlabel('Time from threshold crossing (ms)');
        xlabel('Time from CSMT2 (ms)');
        ylabel('Trials');
        title(sprintf('Day %d Unit %d Raster (%s) NoCSM', day, u, unit_type));
        box off; set(gca,'TickDir','out')
        
        % Adapt raster 
        [trial_idx_ad, col_idx_ad] = find(raster_AD);
        spike_times_ad = linspace(-1000, 999, size(raster_AD,2));
        spike_times_ad = spike_times_ad(col_idx_ad);
        
    
        subplot(3,1,2); hold on
        scatter(spike_times_ad, trial_idx_ad, 1, 'm', 'filled');
        xlim(time_window);
        ylim([0 nTrials+1]);
        %xlabel('Time from threshold crossing  (ms)');
         xlabel('Time from CSMT2 (ms)');
        ylabel('Trials');
        title(sprintf('Day %d Unit %d Raster (%s) CSM', day, u, unit_type));
        box off; set(gca,'TickDir','out')
        % -----------------------
        % PSTH PLOTS
        % -----------------------
        Stimcolor=rgb('DodgerBlue');
        subplot(3,1,3); hold on
        plot(tplot_all,  psth_all_p,  'k',          'LineWidth',1.6);
        plot(tplot_even, psth_even_p, 'Color', 'm', 'LineWidth',1.3);
       
        xlim(time_window);
        %xlabel('Time from threshold crossing (ms)');
         xlabel('Time from CSMT2 (ms)');
        ylabel('Firing Rate (Hz)');
        legend({'NoCSM','CSM'},'Location','best');
        title(sprintf('Day %d Unit %d PSTH (%s)', day, u, unit_type));
        box off; set(gca,'TickDir','out')

        % save
        filename = fullfile(output_folder, sprintf('Day%d_Unit%d_%s.svg', day, u, unit_type));
        saveas(fig, filename);
        close(fig);

    end
end

%% Behavior aligned — BL + AD on SAME raster axis (BL first, then AD), plus PSTHs overlaid
% Rasters on the same plot 

%output_folder = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\Behvaior Aligned\Posthresh";
%output_folder="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\NoStim_CSMModulation\IndividualUnits_Csmt2";
%output_folder="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\CSMvsNoCSM_BLandWO\Baseline Trials";
output_folder="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\ExampleUnits_optotagged\OptoTagged_SameRasterPlot";
if ~exist(output_folder,'dir'); mkdir(output_folder); end

BLcolor = rgb('Black');   % BL trials
ADcolor = rgb('DodgerBlue');     % AD trials

time_window = [-500 500];   % plotting window (ms)

for day = 1:length(NPdata)

    if ~isfield(NPdata(day),'kept_units') || isempty(NPdata(day).kept_units)
        continue
    end

    nUnits = length(NPdata(day).kept_units);

    for u = 1:nUnits

        unit = NPdata(day).kept_units(u);

        % ---------------------------
        % SELECT RASTERS (EDIT THESE)
        % ---------------------------
        % Example (your current lines):
        raster_BL = unit.BLCSM_pos_thresh.Raster;
        raster_AD = unit.ADCSM_pos_thresh.Raster;

        % ---- basic checks ----
        if isempty(raster_BL) || isempty(raster_AD)
            continue
        end

        % use nTrials safely (they can differ)
        nTrialsBL = size(raster_BL,1);
        nTrialsAD = size(raster_AD,1);

        % ---------------------------
        % assign type
        % ---------------------------
        if isfield(NPdata(day),'before_peakv_mask') && NPdata(day).before_peakv_mask(u)
            unit_type = 'Early';
        elseif isfield(NPdata(day),'after_peakv_mask') && NPdata(day).after_peakv_mask(u)
            unit_type = 'Late';
        else
            unit_type = 'Uncategorized';
        end

        % ---------------------------
        % PSTHs (smoothed)
        % ---------------------------
        [~, psth_bl_smooth] = Psth1(raster_BL,8);
        [~, psth_ad_smooth] = Psth1(raster_AD,8);

        timevec_bl = linspace(-1000, 999, numel(psth_bl_smooth));
        timevec_ad = linspace(-1000, 999, numel(psth_ad_smooth));

        mask_bl = timevec_bl >= time_window(1) & timevec_bl <= time_window(2);
        mask_ad = timevec_ad >= time_window(1) & timevec_ad <= time_window(2);

        tplot_bl = timevec_bl(mask_bl);
        tplot_ad = timevec_ad(mask_ad);

        psth_bl_p = psth_bl_smooth(mask_bl);
        psth_ad_p = psth_ad_smooth(mask_ad);

        % ---------------------------
        % SPIKES for raster scatter
        % ---------------------------
        [trial_idx_bl, col_idx_bl] = find(raster_BL);
        spike_times_axis_bl = linspace(-1000, 999, size(raster_BL,2));
        spike_times_bl = spike_times_axis_bl(col_idx_bl);

        [trial_idx_ad, col_idx_ad] = find(raster_AD);
        spike_times_axis_ad = linspace(-1000, 999, size(raster_AD,2));
        spike_times_ad = spike_times_axis_ad(col_idx_ad);

        % Offset AD trial indices so BL is on top and AD is below (single axis)
        trial_idx_ad_plot = trial_idx_ad + nTrialsBL;

        % ---------------------------
        % FIGURE (single raster axis)
        % ---------------------------
        fig = figure('Visible','off','Position',[200 200 650 520]);

        % ---- RASTER: BL then AD on SAME subplot ----
        subplot(2,1,1); hold on
        scatter(spike_times_bl, trial_idx_bl,      1, 'MarkerEdgeColor', BLcolor, 'MarkerFaceColor', BLcolor);
        scatter(spike_times_ad, trial_idx_ad_plot, 1, 'MarkerEdgeColor', ADcolor, 'MarkerFaceColor', ADcolor);

        xlim(time_window);
        ylim([0 (nTrialsBL+nTrialsAD)+1]);
        xline(0,'--')
        xline(50,'--')

        xlabel('Time from threshold crossing (ms)');
        ylabel('Trials');

        % divider line between BL and AD blocks
        yline(nTrialsBL+0.5,'k-','LineWidth',0.75);

        title(sprintf('Day %d Unit %d Raster (%s): BL (top) + AD (bottom)', day, u, unit_type));
        box off; set(gca,'TickDir','out')

        % ---- PSTH: BL + AD overlaid ----
        subplot(2,1,2); hold on
        plot(tplot_bl, psth_bl_p, 'Color', BLcolor, 'LineWidth',1.6);
        plot(tplot_ad, psth_ad_p, 'Color', ADcolor, 'LineWidth',1.3);

        xlim(time_window);
        xline(0,'--')
        xline(50,'--')
        xlabel('Time from threshold crossing (ms)');
        ylabel('Firing Rate (Hz)');
        legend({'BL','AD'},'Location','best');
        title(sprintf('Day %d Unit %d PSTH (%s)', day, u, unit_type));
        box off; set(gca,'TickDir','out')

        % save
        filename = fullfile(output_folder, sprintf('Day%d_Unit%d_%s_BLvsAD.svg', day, u, unit_type));
        saveas(fig, filename);
        close(fig);

    end
end

%% ===================== Behavior aligned: 4 rasters + 1 PSTH overlay =====================
% DOES NOT SKIP if any single raster is empty:
% - each raster subplot prints "NO TRIALS" if empty
% - PSTH overlays only the conditions that exist

%output_folder="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\CSMvsNoCSM_BLandWO\BL and WO";
output_folder="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\CSMvsNoCSM_BLandWO\BL and WO with UNit IDS";

time_window = [-500 500];   % plotting window
full_window = [-1000 999];  % assumed timebase for rasters/Psth1
smooth_bins = 5;

% colors
cBL_No = rgb('DimGray');
cWO_No = rgb('Plum');
cBL_C  = rgb('Black');
cWO_C  = rgb('Magenta');

for day = 1:length(NPdata)

    if ~isfield(NPdata(day),'kept_units') || isempty(NPdata(day).kept_units)
        continue
    end

    nUnits = length(NPdata(day).kept_units);

    for u = 1:nUnits

        unit = NPdata(day).kept_units(u);

        % -----------------------
        % Pull rasters safely (may be empty)
        % -----------------------
        raster_BL_No = [];
        raster_WO_No = [];
        raster_BL_C  = [];
        raster_WO_C  = [];

        if isfield(unit,'BLNoCSM_csmT2') && isfield(unit.BLNoCSM_csmT2,'Raster')
            raster_BL_No = unit.BLNoCSM_csmT2.Raster;
        end
        if isfield(unit,'WONoCSM_csmT2') && isfield(unit.WONoCSM_csmT2,'Raster')
            raster_WO_No = unit.WONoCSM_csmT2.Raster;
        end
        if isfield(unit,'BLCSM_csmT2') && isfield(unit.BLCSM_csmT2,'Raster')
            raster_BL_C = unit.BLCSM_csmT2.Raster;
        end
        if isfield(unit,'WOCSM_csmT2') && isfield(unit.WOCSM_csmT2,'Raster')
            raster_WO_C = unit.WOCSM_csmT2.Raster;
        end

        % If literally everything is missing, still plot placeholders (no skip)
        % (You can skip here if you want; leaving as "plot anyway".)

        % -----------------------
        % assign type
        % -----------------------
        if isfield(NPdata(day),'before_peakv_mask') && numel(NPdata(day).before_peakv_mask) >= u && NPdata(day).before_peakv_mask(u)
            unit_type = 'Early';
        elseif isfield(NPdata(day),'after_peakv_mask') && numel(NPdata(day).after_peakv_mask) >= u && NPdata(day).after_peakv_mask(u)
            unit_type = 'Late';
        else
            unit_type = 'Uncategorized';
        end

        % -----------------------
        % PSTHs (compute only if raster exists)
        % -----------------------
        psth_BL_No = []; tvec_BL_No = [];
        psth_WO_No = []; tvec_WO_No = [];
        psth_BL_C  = []; tvec_BL_C  = [];
        psth_WO_C  = []; tvec_WO_C  = [];

        if ~isempty(raster_BL_No)
            [~, psth_BL_No] = Psth1(raster_BL_No, smooth_bins);
            tvec_BL_No = linspace(full_window(1), full_window(2), numel(psth_BL_No));
        end
        if ~isempty(raster_WO_No)
            [~, psth_WO_No] = Psth1(raster_WO_No, smooth_bins);
            tvec_WO_No = linspace(full_window(1), full_window(2), numel(psth_WO_No));
        end
        if ~isempty(raster_BL_C)
            [~, psth_BL_C]  = Psth1(raster_BL_C, smooth_bins);
            tvec_BL_C = linspace(full_window(1), full_window(2), numel(psth_BL_C));
        end
        if ~isempty(raster_WO_C)
            [~, psth_WO_C]  = Psth1(raster_WO_C, smooth_bins);
            tvec_WO_C = linspace(full_window(1), full_window(2), numel(psth_WO_C));
        end

        % -----------------------
        % FIGURE
        % -----------------------
        fig = figure('Visible','off','Position',[200 200 750 900]);

        % ===== 1) BL NoCSM raster =====
        subplot(5,1,1); hold on
        if isempty(raster_BL_No)
            xlim(time_window); ylim([0 1]);
            text(mean(time_window), 0.5, 'NO TRIALS', 'HorizontalAlignment','center', 'FontWeight','bold');
        else
            [tri, ci] = find(raster_BL_No);
            spk_t = linspace(full_window(1), full_window(2), size(raster_BL_No,2));
            spk_t = spk_t(ci);
            scatter(spk_t, tri, 1, 'MarkerFaceColor', cBL_No, 'MarkerEdgeColor', cBL_No);
            ylim([0 size(raster_BL_No,1)+1]);
            xlim(time_window);
        end
        ylabel('Trials');
        title(sprintf('Day %d Unit %d (%s) BL NoCSM', day, u, unit_type));
        box off; set(gca,'TickDir','out')

        % ===== 2) WO NoCSM raster =====
        subplot(5,1,2); hold on
        if isempty(raster_WO_No)
            xlim(time_window); ylim([0 1]);
            text(mean(time_window), 0.5, 'NO TRIALS', 'HorizontalAlignment','center', 'FontWeight','bold');
        else
            [tri, ci] = find(raster_WO_No);
            spk_t = linspace(full_window(1), full_window(2), size(raster_WO_No,2));
            spk_t = spk_t(ci);
            scatter(spk_t, tri, 1, 'MarkerFaceColor', cWO_No, 'MarkerEdgeColor', cWO_No);
            ylim([0 size(raster_WO_No,1)+1]);
            xlim(time_window);
        end
        ylabel('Trials');
        title(sprintf('Day %d Unit %d (%s) WO NoCSM', day, u, unit_type));
        box off; set(gca,'TickDir','out')

        % ===== 3) BL CSM raster =====
        subplot(5,1,3); hold on
        if isempty(raster_BL_C)
            xlim(time_window); ylim([0 1]);
            text(mean(time_window), 0.5, 'NO TRIALS', 'HorizontalAlignment','center', 'FontWeight','bold');
        else
            [tri, ci] = find(raster_BL_C);
            spk_t = linspace(full_window(1), full_window(2), size(raster_BL_C,2));
            spk_t = spk_t(ci);
            scatter(spk_t, tri, 1, 'MarkerFaceColor', cBL_C, 'MarkerEdgeColor', cBL_C);
            ylim([0 size(raster_BL_C,1)+1]);
            xlim(time_window);
        end
        ylabel('Trials');
        title(sprintf('Day %d Unit %d (%s) BL CSM', day, u, unit_type));
        box off; set(gca,'TickDir','out')

        % ===== 4) WO CSM raster =====
        subplot(5,1,4); hold on
        if isempty(raster_WO_C)
            xlim(time_window); ylim([0 1]);
            text(mean(time_window), 0.5, 'NO TRIALS', 'HorizontalAlignment','center', 'FontWeight','bold');
        else
            [tri, ci] = find(raster_WO_C);
            spk_t = linspace(full_window(1), full_window(2), size(raster_WO_C,2));
            spk_t = spk_t(ci);
            scatter(spk_t, tri, 1, 'MarkerFaceColor', cWO_C, 'MarkerEdgeColor', cWO_C);
            ylim([0 size(raster_WO_C,1)+1]);
            xlim(time_window);
        end
        xlabel('Time from CSMT2 (ms)');
        ylabel('Trials');
        title(sprintf('Day %d Unit %d (%s) WO CSM', day, u, unit_type));
        box off; set(gca,'TickDir','out')

        % ===== 5) PSTH overlay =====
        subplot(5,1,5); hold on
        leg = {};

        if ~isempty(psth_BL_No)
            m = tvec_BL_No >= time_window(1) & tvec_BL_No <= time_window(2);
            plot(tvec_BL_No(m), psth_BL_No(m), 'LineWidth', 0.5, 'Color', cBL_No);
            leg{end+1} = 'BL NoCSM';
        end
        if ~isempty(psth_WO_No)
            m = tvec_WO_No >= time_window(1) & tvec_WO_No <= time_window(2);
            plot(tvec_WO_No(m), psth_WO_No(m), 'LineWidth', 0.5, 'Color', cWO_No);
            leg{end+1} = 'WO NoCSM';
        end
        if ~isempty(psth_BL_C)
            m = tvec_BL_C >= time_window(1) & tvec_BL_C <= time_window(2);
            plot(tvec_BL_C(m), psth_BL_C(m), 'LineWidth', 0.5, 'Color', cBL_C);
            leg{end+1} = 'BL CSM';
        end
        if ~isempty(psth_WO_C)
            m = tvec_WO_C >= time_window(1) & tvec_WO_C <= time_window(2);
            plot(tvec_WO_C(m), psth_WO_C(m), 'LineWidth', 0.5, 'Color', cWO_C);
            leg{end+1} = 'WO CSM';
        end

        xlim(time_window);
        xlabel('Time from CSMT2 (ms)');
        ylabel('Firing Rate (Hz)');
        if ~isempty(leg)
            legend(leg,'Location','best');
        else
            text(mean(time_window), 0.5, 'NO PSTH (no trials in any condition)', ...
                'HorizontalAlignment','center','FontWeight','bold');
        end
        title(sprintf('Day %d Unit %d PSTH (%s)', day, u, unit_type));
        box off; set(gca,'TickDir','out')
        % 
        % % save
        % filename = fullfile(output_folder, sprintf('Day%d_Unit%d_%s_4cond.svg', day, u, unit_type));
        % saveas(fig, filename);
        % close(fig);

        % -----------------------
        % save (include actual unit index if available)
        % -----------------------
        if isfield(unit,'unitIdx')
            unit_id = unit.unitIdx;   % actual Neuropixels unit ID
        else
            unit_id = u;              % fallback to loop index
        end

        filename = fullfile(output_folder, ...
            sprintf('Day%d_Unit%d_UID%d_%s_4cond.svg', day, u, unit_id, unit_type));

        saveas(fig, filename);
        close(fig);


    end
end

%% NEW EDIT BEHAVIOR AS WELL - with WO raster
% Generate raster + PSTH for each kept unit and save as SVG (white raster background)

%output_folder = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\Behvaior Aligned\Posthresh";

%output_folder="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\ExampleUnits_filteredmodulated\pos_threshaligned";
%output_folder="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\Behvaior Aligned\Pos_thresh_wWO";
%output_folder="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\Behvaior Aligned\Pos_thresh_wWO\With WO raster";
% output_folder="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\Behvaior Aligned\Pos_thresh_wWO\With 5 WO trials";
% output_folder="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\Behvaior Aligned\CSMT2aligned";

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
        % raster_BL = unit.BLCSM_pos_thresh.Raster;            % trials x time
        % raster_AD = unit.ADCSM_pos_thresh.Raster;
        % 
        % raster_WO_full = unit.WOCSM_pos_thresh.Raster;
        % 
        raster_BL = unit.BLCSM_csmT2.Raster;            % trials x time
        raster_AD = unit.ADCSM_csmT2.Raster;

        raster_WO_full = unit.WOCSM_csmT2.Raster;

        % --- Keep only first n WO trials ---
        nWO_plot = 10;   % number of WO trials to include
        %nWO_plot = 5;   % number of WO trials to include

        nWO = min(nWO_plot, size(raster_WO_full,1));
        raster_WO = raster_WO_full(1:nWO, :);

     

        %timevec_raster = linspace(unit.BLCSM_pos_thresh.time_bins(1), unit.BLCSM_pos_thresh.time_bins(end), size(raster_BL,2));
        timevec_raster = linspace(unit.BLCSM_csmT2.time_bins(1), unit.BLCSM_csmT2.time_bins(end), size(raster_BL,2));

        nTrials = size(raster_BL, 1);
       
        % --- Compute PSTHs (with smoothing) ---
        [~, psth_bl, ~]  = computePSTH(raster_BL, bin_size, smoothing_bins);
        [~, psth_ad, ~] = computePSTH(raster_AD, bin_size, smoothing_bins);
        [~, psth_wo, ~] = computePSTH(raster_WO, bin_size, smoothing_bins);


        % --- Construct consistent PSTH time vector ---
        num_time_bins = size(psth_bl, 2);
        timevec_psth = linspace(time_range(1), time_range(2), num_time_bins);

        % --- Apply time window mask for plotting ---
        time_mask = timevec_psth >= time_window(1) & timevec_psth <= time_window(2);
        timevec_plot = timevec_psth(time_mask);
        psth_bl_plot  = psth_bl(time_mask);
        psth_ad_plot = psth_ad(time_mask);
        psth_wo_plot = psth_wo(time_mask);
 

        %Baseline trials --
        [trial_idx_BL, col_idx_BL] = find(raster_BL);  % trial x raster bin index
        spike_times_BL = timevec_raster(col_idx_BL);

        % --- Plot Raster ---
        %fig = figure('Visible', 'off', 'Position', [100, 100, 600, 500]);
        fig = figure('Visible', 'off', 'Position', [100, 100, 500, 550]);
        subplot(4, 1, 1);
        hold on;
        set(gca, 'Color', 'w'); % white background

        % Plot each spike as a small black dot
        scatter(spike_times_BL, trial_idx_BL, 1.5, 'k', 'filled'); 
        xlim(time_window);
        ylim([0 nTrials+1]);
        xline(0,'--')
        xline(50,'--')
        %xlabel('Time from threshold crossing  (ms)');
        xlabel('Time from CSM T2 (ms)');
        ylabel('Trials');
        title(sprintf('Day %d, Unit %d Raster', day, u));
        box off;
        set(gca, 'TickDir', 'out');


        % Stimulated trials 
        [trial_idx_AD, col_idx_AD] = find(raster_AD);  % trial x raster bin index
        spike_times_AD = timevec_raster(col_idx_AD);

        % --- Plot Raster ---
        nTrials = size(raster_AD, 1);
        subplot(4, 1, 2);
        hold on;
        set(gca, 'Color', 'w'); % white background

        % Plot each spike as a small black dot
        scatter(spike_times_AD, trial_idx_AD, 1.5, 'b', 'filled'); 
        xlim(time_window);
        ylim([0 nTrials+1]);
        xline(0,'--')
        xline(50,'--')
        %xlabel('Time from threshold crossing  (ms)');
        xlabel('Time from CSM T2 (ms)');
        ylabel('Trials');
        title(sprintf('Day %d, Unit %d Raster', day, u));
        box off;
        set(gca, 'TickDir', 'out');

          % Washout  trials 
        [trial_idx_WO, col_idx_WO] = find(raster_WO);  % trial x raster bin index
        spike_times_WO = timevec_raster(col_idx_WO);

        % --- Plot Raster ---
        nTrials = size(raster_WO, 1);
        subplot(4, 1, 3);
        hold on;
        set(gca, 'Color', 'w'); % white background

        % Plot each spike as a small black dot
        scatter(spike_times_WO, trial_idx_WO, 1.5, 'm', 'filled'); 
        xlim(time_window);
        ylim([0 nTrials+1]);
        xline(0,'--')
        xline(50,'--')
        %xlabel('Time from threshold crossing  (ms)');
        xlabel('Time from CSM T2 (ms)');
        ylabel('Trials');
        title(sprintf('Day %d, Unit %d Raster', day, u));
        box off;
        set(gca, 'TickDir', 'out');

        % --- Plot PSTH ---
        Stimcolor=rgb('DodgerBlue');
        Washcolor=rgb('Magenta');
        subplot(4, 1, 4);
        plot(timevec_plot, psth_bl_plot, 'k', 'LineWidth', 1.5); hold on;
        plot(timevec_plot, psth_ad_plot, 'Color', Stimcolor, 'LineWidth', 1.2);hold on;
        plot(timevec_plot, psth_wo_plot, 'Color', Washcolor, 'LineWidth', 1.2);
      
        xlim(time_window);
        xline(0,'--')
        xline(50,'--')
        %xlabel('Time from threshold crossing  (ms)');
        xlabel('Time from CSM T2 (ms)');
        ylabel('Firing Rate (spikes/s)');
     
        title(sprintf('Day %d, Unit %d PSTH (Firing Rate)', day, u));
        box off;
        set(gca, 'TickDir', 'out');

        % --- Save figure as SVG ---
        filename = fullfile(output_folder, sprintf('Day%d_Unit%d.svg', day, u));
        saveas(fig, filename);
        close(fig);
    end
end

%% Edited code to plot all rasters on the same plot 

%% ALL 3 RASTERS ON THE SAME PLOT + PSTH BELOW
% Generate combined raster + PSTH for each kept unit and save as SVG

% output_folder = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\Behvaior Aligned\Posthresh";
% output_folder = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\ExampleUnits_filteredmodulated\pos_threshaligned";
% output_folder = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\Behvaior Aligned\Pos_thresh_wWO";
% output_folder = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\Behvaior Aligned\Pos_thresh_wWO\With WO raster";
% output_folder = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\Behvaior Aligned\Pos_thresh_wWO\With 5 WO trials";
% output_folder = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\Behvaior Aligned\CSMT2aligned";
%output_folder="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\Behvaior Aligned\CSMT2AlignedEdit";
output_folder="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\Behvaior Aligned\PosThreshEdit";

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

bin_size       = 10;               % ms
smoothing_bins = 5;                % PSTH smoothing
time_range     = [-1000, 1000];    % full PSTH range (ms)
%time_window    = [-1000, 1000]; 
time_window    = [-500, 500];      % plotting window (ms)
%time_window    = [-250, 250];   

Stimcolor = rgb('DodgerBlue');
Washcolor = rgb('Magenta');

for day = 1:length(NPdata)

    if ~isfield(NPdata(day), 'kept_units') || isempty(NPdata(day).kept_units)
        continue
    end

    nUnits = length(NPdata(day).kept_units);

    for u = 1:nUnits

        unit = NPdata(day).kept_units(u);

        % ---------------- LOAD RASTERS ----------------
         raster_BL      = unit.BLCSM_pos_thresh.Raster;
        raster_AD      = unit.ADCSM_pos_thresh.Raster;
        raster_WO_full = unit.WOCSM_pos_thresh.Raster;
        % raster_BL      = unit.BLCSM_csmT2.Raster;
        % raster_AD      = unit.ADCSM_csmT2.Raster;
        % raster_WO_full = unit.WOCSM_csmT2.Raster;

        % ---------------- LIMIT WO TRIALS ----------------
        nWO_plot = 10;   % set to 5 if desired
        nWO      = min(nWO_plot, size(raster_WO_full,1));
        raster_WO = raster_WO_full(1:nWO, :);

        % ---------------- TIME VECTOR ----------------
        timevec_raster = linspace(unit.BLCSM_csmT2.time_bins(1), ...
                                  unit.BLCSM_csmT2.time_bins(end), ...
                                  size(raster_BL,2));

        % ---------------- COMPUTE PSTHs ----------------
          
        [~, psth_bl, ~] = computePSTH(raster_BL, bin_size, smoothing_bins);
        [~, psth_ad, ~] = computePSTH(raster_AD, bin_size, smoothing_bins);
        [~, psth_wo, ~] = computePSTH(raster_WO, 10, 6);

        num_time_bins = size(psth_bl, 2);
        timevec_psth  = linspace(time_range(1), time_range(2), num_time_bins);

        time_mask     = timevec_psth >= time_window(1) & timevec_psth <= time_window(2);
        timevec_plot  = timevec_psth(time_mask);
        psth_bl_plot  = psth_bl(time_mask);
        psth_ad_plot  = psth_ad(time_mask);
        psth_wo_plot  = psth_wo(time_mask);

        % ---------------- FIND SPIKES FOR EACH CONDITION ----------------
        [trial_idx_BL, col_idx_BL] = find(raster_BL);
        spike_times_BL = timevec_raster(col_idx_BL);

        [trial_idx_AD, col_idx_AD] = find(raster_AD);
        spike_times_AD = timevec_raster(col_idx_AD);

        [trial_idx_WO, col_idx_WO] = find(raster_WO);
        spike_times_WO = timevec_raster(col_idx_WO);

        % ---------------- STACK TRIALS VERTICALLY ----------------
        nBL = size(raster_BL,1);
        nAD = size(raster_AD,1);
        nWO = size(raster_WO,1);

        gap = 2;  % blank rows between groups

        % BL on top, then AD, then WO
        trial_idx_BL_plot = trial_idx_BL;
        trial_idx_AD_plot = trial_idx_AD + nBL + gap;
        trial_idx_WO_plot = trial_idx_WO + nBL + gap + nAD + gap;

        total_trials = nBL + nAD + nWO + 2*gap;

        % centers for labels
        yBL = (1 + nBL) / 2;
        yAD = nBL + gap + (1 + nAD) / 2;
        yWO = nBL + gap + nAD + gap + (1 + nWO) / 2;

        % separator lines
        sep1 = nBL + gap/2;
        sep2 = nBL + gap + nAD + gap/2;

        % ---------------- MAKE FIGURE ----------------
        fig = figure('Visible','off','Color','w','Position',[100,100,650,520]);

        % =========================================================
        % TOP: COMBINED RASTER
        % =========================================================
        ax1 = subplot(2,1,1);
        hold(ax1,'on');
        set(ax1,'Color','w');

        scatter(ax1, spike_times_BL, trial_idx_BL_plot, 1.5, 'k', 'filled');
        scatter(ax1, spike_times_AD, trial_idx_AD_plot, 1.5, Stimcolor, 'filled');
        scatter(ax1, spike_times_WO, trial_idx_WO_plot, 1.5, Washcolor, 'filled');


        xline(ax1, 0,  '--k');
        xline(ax1, -81.5, '--k');
        xline(ax1, -136.5, '--k');

        yline(ax1, sep1, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
        yline(ax1, sep2, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);

        xlim(ax1, time_window);
        ylim(ax1, [0 total_trials+1]);

        xlabel(ax1, 'Time from CSM T2 (ms)');
        ylabel(ax1, 'Trials');
        title(ax1, sprintf('Day %d, Unit %d Raster', day, u));

        set(ax1, 'YTick', [yBL yAD yWO], ...
                 'YTickLabel', {sprintf('BL (n=%d)',nBL), ...
                                sprintf('AD (n=%d)',nAD), ...
                                sprintf('WO (n=%d)',nWO)});
        box(ax1,'off');
        set(ax1,'TickDir','out');

        % =========================================================
        % BOTTOM: PSTH
        % =========================================================
        ax2 = subplot(2,1,2);
        hold(ax2,'on');

        plot(ax2, timevec_plot, psth_bl_plot, 'k', 'LineWidth', 1.5);
        plot(ax2, timevec_plot, psth_ad_plot, 'Color', Stimcolor, 'LineWidth', 1.2);
        plot(ax2, timevec_plot, psth_wo_plot, 'Color', Washcolor, 'LineWidth', 1.2);

        xline(ax2, 0,  '--k');
        
        xline(ax2, 0,  '--k');
        xline(ax2, -81.5, '--k');
        xline(ax2, -136.5, '--k');

        xlim(ax2, time_window);
        xlabel(ax2, 'Time from CSM T2 (ms)');
        ylabel(ax2, 'Firing Rate (spikes/s)');
        title(ax2, sprintf('Day %d, Unit %d PSTH', day, u));

        legend(ax2, {'BL','AD','WO'}, 'Location','best', 'Box','off');
        box(ax2,'off');
        set(ax2,'TickDir','out');

        % ---------------- SAVE ----------------
        filename = fullfile(output_folder, sprintf('Day%d_Unit%d.svg', day, u));
        saveas(fig, filename);
        close(fig);

    end
end
%% Smooth PSTH 
%% ALL 3 RASTERS ON THE SAME PLOT + PSTH BELOW
% Generate combined raster + PSTH for each kept unit and save as SVG

% output_folder = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\Behvaior Aligned\Posthresh";
% output_folder = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\ExampleUnits_filteredmodulated\pos_threshaligned";
% output_folder = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\Behvaior Aligned\Pos_thresh_wWO";
% output_folder = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\Behvaior Aligned\Pos_thresh_wWO\With WO raster";
% output_folder = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\Behvaior Aligned\Pos_thresh_wWO\With 5 WO trials";
% output_folder = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\Behvaior Aligned\CSMT2aligned";
%output_folder="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\Behvaior Aligned\CSMT2AlignedEdit";
output_folder="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\Behvaior Aligned\PosThreshEdit";

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

bin_size       = 10;               % ms
smoothing_bins = 5;                % initial PSTH smoothing inside computePSTH
extra_smooth   = 2;                % additional smoothing after PSTH extraction
time_range     = [-1000, 1000];    % full PSTH range (ms)
%time_window    = [-1000, 1000];
time_window    = [-500, 500];      % plotting window (ms)
%time_window    = [-250, 250];

Stimcolor = rgb('DodgerBlue');
Washcolor = rgb('Magenta');

for day = 1:length(NPdata)

    if ~isfield(NPdata(day), 'kept_units') || isempty(NPdata(day).kept_units)
        continue
    end

    nUnits = length(NPdata(day).kept_units);

    for u = 1:nUnits

        unit = NPdata(day).kept_units(u);

        % ---------------- LOAD RASTERS ----------------
        % raster_BL      = unit.BLCSM_csmT2.Raster;
        % raster_AD      = unit.ADCSM_csmT2.Raster;
        % raster_WO_full = unit.WOCSM_csmT2.Raster;

         raster_BL      = unit.BLCSM_pos_thresh.Raster;
        raster_AD      = unit.ADCSM_pos_thresh.Raster;
        raster_WO_full = unit.WOCSM_pos_thresh.Raster;

        % ---------------- LIMIT WO TRIALS ----------------
        nWO_plot = 10;   % set to 5 if desired
        nWO      = min(nWO_plot, size(raster_WO_full,1));
        raster_WO = raster_WO_full(1:nWO, :);

        % ---------------- TIME VECTOR ----------------
        timevec_raster = linspace(unit.BLCSM_csmT2.time_bins(1), ...
                                  unit.BLCSM_csmT2.time_bins(end), ...
                                  size(raster_BL,2));

        % ---------------- COMPUTE PSTHs ----------------
        [~, psth_bl, ~] = computePSTH(raster_BL, bin_size, smoothing_bins);
        [~, psth_ad, ~] = computePSTH(raster_AD, bin_size, smoothing_bins);
        [~, psth_wo, ~] = computePSTH(raster_WO, bin_size, smoothing_bins);

        num_time_bins = size(psth_bl, 2);
        timevec_psth  = linspace(time_range(1), time_range(2), num_time_bins);

        time_mask     = timevec_psth >= time_window(1) & timevec_psth <= time_window(2);
        timevec_plot  = timevec_psth(time_mask);
        psth_bl_plot  = psth_bl(time_mask);
        psth_ad_plot  = psth_ad(time_mask);
        psth_wo_plot  = psth_wo(time_mask);

        % ---------------- EXTRA SMOOTHING ----------------
        psth_bl_plot = smoothdata(psth_bl_plot, 'gaussian', extra_smooth);
        psth_ad_plot = smoothdata(psth_ad_plot, 'gaussian', extra_smooth);
        psth_wo_plot = smoothdata(psth_wo_plot, 'gaussian', extra_smooth);

        % ---------------- FIND SPIKES FOR EACH CONDITION ----------------
        [trial_idx_BL, col_idx_BL] = find(raster_BL);
        spike_times_BL = timevec_raster(col_idx_BL);

        [trial_idx_AD, col_idx_AD] = find(raster_AD);
        spike_times_AD = timevec_raster(col_idx_AD);

        [trial_idx_WO, col_idx_WO] = find(raster_WO);
        spike_times_WO = timevec_raster(col_idx_WO);

        % ---------------- STACK TRIALS VERTICALLY ----------------
        nBL = size(raster_BL,1);
        nAD = size(raster_AD,1);
        nWO = size(raster_WO,1);

        gap = 1;  % blank rows between groups

        % BL on top, then AD, then WO
        % trial_idx_BL_plot = trial_idx_BL;
        % trial_idx_AD_plot = trial_idx_AD + nBL + gap;
        % trial_idx_WO_plot = trial_idx_WO + nBL + gap + nAD + gap;
        % 
        % BL on top, then AD, then WO
        trial_idx_WO_plot = trial_idx_WO;
        trial_idx_AD_plot = trial_idx_AD + nWO + gap;
        trial_idx_BL_plot = trial_idx_BL + nWO + gap + nAD + gap;

        total_trials = nBL + nAD + nWO + 2*gap;

        % % centers for labels
        % yBL = (1 + nBL) / 2;
        % yAD = nBL + gap + (1 + nAD) / 2;
        % yWO = nBL + gap + nAD + gap + (1 + nWO) / 2;

          % centers for labels
        yWO = (1 + nBL) / 2;
        yAD = nBL + gap + (1 + nAD) / 2;
        yBL = nBL + gap + nAD + gap + (1 + nWO) / 2;

        % separator lines
        sep1 = nBL + gap/2;
        sep2 = nBL + gap + nAD + gap/2;

        % ---------------- MAKE FIGURE ----------------
        fig = figure('Visible','off','Color','w','Position',[100,100,650,520]);

        % =========================================================
        % TOP: COMBINED RASTER
        % =========================================================
        ax1 = subplot(2,1,1);
        hold(ax1,'on');
        set(ax1,'Color','w');

        scatter(ax1, spike_times_BL, trial_idx_BL_plot, 1.5, 'k', 'filled');
        scatter(ax1, spike_times_AD, trial_idx_AD_plot, 1.5, Stimcolor, 'filled');
        scatter(ax1, spike_times_WO, trial_idx_WO_plot, 1.5, Washcolor, 'filled');

        xline(ax1, 0, '--k');
        xline(ax1, 81.5, '--k');
        xline(ax1, 50 ,'--k');

        yline(ax1, sep1, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
        yline(ax1, sep2, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);

        xlim(ax1, time_window);
        ylim(ax1, [0 total_trials+1]);

        xlabel(ax1, 'Time from threshold crossing (ms)');
        ylabel(ax1, 'Trials');
        title(ax1, sprintf('Day %d, Unit %d Raster', day, u));

        % set(ax1, 'YTick', [yBL yAD yWO], ...
        %          'YTickLabel', {sprintf('BL (n=%d)',nBL), ...
        %                         sprintf('AD (n=%d)',nAD), ...
        %                         sprintf('WO (n=%d)',nWO)});
        box(ax1,'off');
        set(ax1,'TickDir','out');

        % =========================================================
        % BOTTOM: PSTH
        % =========================================================
        ax2 = subplot(2,1,2);
        hold(ax2,'on');

        plot(ax2, timevec_plot, psth_bl_plot, 'k', 'LineWidth', 1.5);
        plot(ax2, timevec_plot, psth_ad_plot, 'Color', Stimcolor, 'LineWidth', 1.2);
        plot(ax2, timevec_plot, psth_wo_plot, 'Color', Washcolor, 'LineWidth', 1.2);

        xline(ax2, 0, '--k');
        xline(ax2, 81.5, '--k');
        xline(ax2, 50 ,'--k');

        xlim(ax2, time_window);
        xlabel(ax2, 'Time from threshold corsing (ms)');
        ylabel(ax2, 'Firing Rate (spikes/s)');
        title(ax2, sprintf('Day %d, Unit %d PSTH', day, u));

        legend(ax2, {'BL','AD','WO'}, 'Location','best', 'Box','off');
        box(ax2,'off');
        set(ax2,'TickDir','out');

        % ---------------- SAVE ----------------
        filename = fullfile(output_folder, sprintf('Day%d_Unit%d.svg', day, u));
        saveas(fig, filename);
        close(fig);

    end
end
%% Different version of behavior aligned 
% Generate raster + PSTH for each kept unit and save as SVG (white raster background)
output_folder = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\ExampleUnits_optotagged";

oddcolor = rgb('DarkRed');
evencolor = rgb('DarkGreen');

bin_size = 10;         % ms
smoothing_bins = 3;    % for PSTH smoothing
%time_range = [-1000, 1000];  % ms, full PSTH range
%time_window = [-500, 500];   % ms, plotting window
time_range = [-500, 500];  % ms, full PSTH range for optotagged units 
time_window = [-250, 250];   % ms, plotting window

for day = 1:length(NPdata)
    if ~isfield(NPdata(day), 'kept_units') || isempty(NPdata(day).kept_units)
        continue
    end

    nUnits = length(NPdata(day).kept_units);

    for u = 1:nUnits
        unit = NPdata(day).kept_units(u);
        raster = unit.eventopto.Raster;  % trials x time

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

        timevec_raster = linspace(NPdata(1).time_bins_opto(1), NPdata(1).time_bins_opto(end), size(raster,2));

        nTrials = size(raster, 1);
        even_idx = 2:2:nTrials;
        odd_idx  = 1:2:nTrials;

        raster_even = raster(even_idx, :);
        raster_odd  = raster(odd_idx, :);

        % --- Compute PSTHs (with smoothing) ---
        [~, psth_all, ~]  = computePSTH(raster, bin_size, smoothing_bins);
        [~, psth_even, ~] = computePSTH(raster_even, bin_size, smoothing_bins);
        [~, psth_odd,  ~] = computePSTH(raster_odd,  bin_size, smoothing_bins);
        % 
        % psth_raw,psth_all = Psth1opto(raster,2); % was 5 before 
        % psth_raw,psth_even = Psth1opto(raster_even,2); % was 5 before
        % psth_raw,psth_odd = Psth1opto(raster_odd,2); % was 5 before
        % 

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
        scatter(spike_times, trial_idx, 1, 'k', 'filled'); 
        xlim(time_window);
        ylim([0 nTrials+1]);
        xlabel('Time (ms)');
        ylabel('Trials');
        title(sprintf('Day %d, Unit %d Raster (%s)', day, u, unit_type));
        box off;
        set(gca, 'TickDir', 'out');
        % % --- Plot Raster ---
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
        xlabel('Time (ms)');
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
output_folder ="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\Sham_evenoddtrials\Example Units_CSMT0 aligned";
%output_folder="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\Sham_evenoddtrials\Example Units_CSMT2 aligned";
oddcolor = rgb('Black'); % NO CSM reaches 
evencolor = rgb('DarkMagenta'); % CSM reaches 

bin_size = 10;         % ms
smoothing_bins = 3;    % for PSTH smoothing
time_range = [-1000, 1000];  % ms, full PSTH range
%time_window = [-500, 500];   % ms, plotting window
time_window = [-200, 200];   % ms, plotting window

for day = 1:length(NPdata)
    if ~isfield(NPdata(day), 'kept_units') || isempty(NPdata(day).kept_units)
        continue
    end

    nUnits = length(NPdata(day).kept_units);

    for u = 1:nUnits
        unit = NPdata(day).kept_units(u);
        %raster = unit.max_vx.Raster;  % trials x time

        raster = unit.AllCSM_csmT0.Raster;  % trials x time
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

        % raster_even = raster(even_idx, :);
        % raster_odd  = raster(odd_idx, :);
        % 
        raster_even = unit.AllCSM_csmT0.Raster; 
        raster_odd  = unit.NoCSM_csmT0.Raster; 

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
        %plot(timevec_plot, psth_all_plot, 'k', 'LineWidth', 1.5); hold on;
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
