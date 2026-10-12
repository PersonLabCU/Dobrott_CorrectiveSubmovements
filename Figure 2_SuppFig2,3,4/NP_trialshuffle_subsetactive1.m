% Sampling even and odd trials to look out for spurious tiling

cd("C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\ShamAnimals\ReachandCSMMod_allCSMStogther\CSMAll_wtihMaxV");
MyFolderInfo = dir('*CSMAll.mat');

tic
clearvars -except MyFolderInfo

%% ALL REACHES PSTH + PEAK TIME ANALYSIS (Even/Odd + Early/Late)
NPdata = struct;

PsthMasterCtrl = [];
PsthMasterNoCSM=[];
PsthMasterCSM=[];
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
    %valid_peak_range = find(timevec >= -100 & timevec <= 100);
    valid_peak_range = find(timevec >= -250 & timevec <= 250);
   %valid_peak_range = find(timevec >= -1000 & timevec <= 1000);

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
    psth_master_csm = NaN(nCells, nBins);
    psth_master_nocsm = NaN(nCells, nBins);
    psth_master_even = NaN(nCells, nBins);
    psth_master_odd  = NaN(nCells, nBins);

    peak_times_ctrl = NaN(nCells,1);
    peak_times_csm = NaN(nCells,1);
    peak_times_nocsm = NaN(nCells,1);
    peak_times_even = NaN(nCells,1);
    peak_times_odd  = NaN(nCells,1);

    % --- Loop through each cell ---
    for t = 1:nCells
        % --- Control PSTH ---
        %psth = data_opto(t).max_vx.PSTH_smoothbin;
        psth_raw = data_opto(t).pos_thresh.PSTH_smoothbin;
        psth = psth_raw / max(psth_raw);              % normalize to max
        psth_master_ctrl(t,:) = psth;
        [~, peak_times_ctrl(t)] = max(psth);

        % CSM vs NOCSM triasl 
        psth_csm_raw = data_opto(t).AllCSM_pos_thresh.PSTH_smoothbin;
        %psth_csm = data_opto(t).AllCSM_maxvx.PSTH_smoothbin;
        psth_csm = psth_csm_raw / max(psth_csm_raw);              % normalize to max
        %psth_csm=psth_csm_raw;
        %psth_csm = psth_csm / max(psth_raw);   
        psth_master_csm(t,:) = psth_csm;
        [~, peak_times_csm(t)] = max(psth_csm);
        
        %NoCSM trials
        psth_nocsm = data_opto(t).NoCSM_pos_thresh.PSTH_smoothbin;
        %psth_nocsm = data_opto(t).NoCSM_maxvx.PSTH_smoothbin;
        psth_nocsm = psth_nocsm / max(psth_nocsm);              % normalize to max
        %psth_nocsm = psth_nocsm / max(psth_csm_raw); 
        %psth_nocsm=psth_nocsm;
        psth_master_nocsm(t,:) = psth_nocsm;
        [~, peak_times_nocsm(t)] = max(psth_nocsm);


        % --- Even/Odd PSTHs ---
        % if isfield(data_opto(t).max_vx, 'Raster') && ~isempty(data_opto(t).max_vx.Raster)
        %     raster = data_opto(t).max_vx.Raster; % trials x time

        if isfield(data_opto(t).max_vx, 'Raster') && ~isempty(data_opto(t).max_vx.Raster)
            raster = data_opto(t).max_vx.Raster; % trials x time
            nTrials = size(raster,1);
            even_idx = 2:2:nTrials;
            odd_idx  = 1:2:nTrials;

            raster_even = raster(even_idx,:);
            raster_odd  = raster(odd_idx,:);

            bin_size = 10; % ms
            smoothing_bins = 3;

            [~, psth_smooth_even, ~] = computePSTH_old(raster_even, bin_size, smoothing_bins);
            [~, psth_smooth_odd,  ~] = computePSTH_old(raster_odd,  bin_size, smoothing_bins);

            % Normalize
            if max(psth_smooth_even) > 0
                psth_master_even(t,:) = psth_smooth_even / max(psth_smooth_even);
            else
                psth_master_even(t,:) = NaN(1,length(psth_smooth_even));
            end

            if max(psth_smooth_odd) > 0
                psth_master_odd(t,:)  = psth_smooth_odd / max(psth_smooth_odd);
            else
                psth_master_odd(t,:)  = NaN(1,length(psth_smooth_odd));
            end

            [~, peak_times_even(t)] = max(psth_master_even(t,:));
            [~, peak_times_odd(t)]  = max(psth_master_odd(t,:));
        end
    end

    % --- Store into NPdata ---
    NPdata(day).peak_times_ctrl = peak_times_ctrl;
    NPdata(day).peak_times_csm = peak_times_csm;
    NPdata(day).peak_times_nocsm = peak_times_nocsm;

    NPdata(day).peak_times_even = peak_times_even;
    NPdata(day).peak_times_odd  = peak_times_odd;
    NPdata(day).psth_master_ctrl = psth_master_ctrl;
    NPdata(day).psth_master_even = psth_master_even;
    NPdata(day).psth_master_odd  = psth_master_odd;

    % --- Concatenate across days ---
    PsthMasterCtrl = vertcat(PsthMasterCtrl, psth_master_ctrl);
    PsthMasterCSM = vertcat(PsthMasterCSM, psth_master_csm);
    PsthMasterNoCSM = vertcat(PsthMasterNoCSM, psth_master_nocsm);

    PsthMasterEven = vertcat(PsthMasterEven, psth_master_even);
    PsthMasterOdd  = vertcat(PsthMasterOdd,  psth_master_odd);

    % --- Separate early/late for all three types ---
    PsthMaster_Before_Ctrl = vertcat(PsthMaster_Before_Ctrl, psth_master_ctrl(before_peakv_mask,:));
    PsthMaster_After_Ctrl  = vertcat(PsthMaster_After_Ctrl,  psth_master_ctrl(after_peakv_mask,:));

    PsthMaster_Before_Even = vertcat(PsthMaster_Before_Even, psth_master_even(before_peakv_mask,:));
    PsthMaster_After_Even  = vertcat(PsthMaster_After_Even,  psth_master_even(after_peakv_mask,:));

    PsthMaster_Before_Odd  = vertcat(PsthMaster_Before_Odd,  psth_master_odd(before_peakv_mask,:));
    PsthMaster_After_Odd   = vertcat(PsthMaster_After_Odd,   psth_master_odd(after_peakv_mask,:));

    % Store additional info
    NPdata(day).before_peakv_mask = before_peakv_mask;
    NPdata(day).after_peakv_mask  = after_peakv_mask;
    NPdata(day).kept_indices = filtered_modulated_indices(keep_unit);

    fprintf('Day %d: Added %d units (%d early, %d late)\n', day, nCells, sum(before_peakv_mask), sum(after_peakv_mask));
    clear data Reaches Control filtered_modulated_indices ExcitRR
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

%% ---- CSM vs NoCSM
% Compute peak indices for each condition
[~, peak_idx_ctrl] = max(PsthMasterCtrl, [], 2);
[~, peak_idx_even] = max(PsthMasterCSM, [], 2);
[~, peak_idx_odd]  = max(PsthMasterNoCSM, [], 2);

% Sort indices by peak timing
[~, sort_idx_ctrl] = sort(peak_idx_ctrl);
[~, sort_idx_even] = sort(peak_idx_even);
[~, sort_idx_odd]  = sort(peak_idx_odd);

% % % % Sort the PSTH matrices
% sorted_ctrl = PsthMasterCtrl(sort_idx_ctrl, :);
% sorted_even = PsthMasterCSM(sort_idx_even, :);
% sorted_odd  = PsthMasterNoCSM(sort_idx_odd, :);

% Sort the PSTH matrices
% sorted_ctrl = PsthMasterCtrl(sort_idx_ctrl, :);
% sorted_even = PsthMasterCSM(sort_idx_ctrl, :);
% sorted_odd  = PsthMasterNoCSM(sort_idx_ctrl, :);

% Sort the PSTH matrices
sorted_ctrl = PsthMasterCtrl(sort_idx_ctrl, :);
sorted_even = PsthMasterCSM(sort_idx_even, :);
sorted_odd  = PsthMasterNoCSM(sort_idx_even, :);

% Generate time vector for x-axis
num_bins = size(PsthMasterCtrl, 2);
time_vector = linspace(-1000, 1000, num_bins);  % adjust if bins differ
tick_times = -1000:250:1000;
[~, tick_indices] = min(abs(time_vector' - tick_times), [], 1);

% Plot all 3, each sorted by its own peak time
figure('Position', [100 100 1600 600]);

subplot(1,3,1);
imagesc(sorted_ctrl);
%colormap(parula);
colormap(slanCM('GnBu'));
colorbar;
set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times);
xlabel('Time from threshold crossing (ms)');
%ylabel('Units (sorted by own peak)');
ylabel('Units (sorted by Ctrl All peak)');
title('Control (sorted by own peaks)');

% subplot(1,3,2);
% imagesc(sorted_even);
% %colormap(parula);
% colormap(slanCM('GnBu'));
% colorbar;
% set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times);
% xlabel('Time from threshold crossing (ms)');
% %ylabel('Units');
% ylabel('Units (sorted by Ctrl All peak)');
% title('CSM Trials');
% 
% subplot(1,3,3);
% imagesc(sorted_odd);
% %colormap(parula);
% colormap(slanCM('GnBu'));
% colorbar;
% set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times);
% xlabel('Time from threshold crossing (ms)');
% %ylabel('Units');
% ylabel('Units (sorted by Ctrl All peak)');
% title('NoCSM Trials ');

% --- get global color limits ---
cmin = min([sorted_even(:); sorted_odd(:)]);
cmax = max([sorted_even(:); sorted_odd(:)]);

% --- CSM ---
subplot(1,3,2);
imagesc(sorted_even);
colormap(slanCM('GnBu'));
caxis([cmin cmax]);   % <-- KEY LINE
colorbar;
set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times);
xlabel('Time from threshold crossing (ms)');
ylabel('Units (sorted by Ctrl All peak)');
title('CSM Trials');

% --- NoCSM ---
subplot(1,3,3);
imagesc(sorted_odd);
colormap(slanCM('GnBu'));
caxis([cmin cmax]);   % <-- SAME LIMITS
colorbar;
set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times);
xlabel('Time from threshold crossing (ms)');
ylabel('Units (sorted by Ctrl All peak)');
title('NoCSM Trials');

sgtitle('Each PSTH Sorted by Its Own Peak Firing Time');

%% ---- DELTA HEATMAP (CSM - NoCSM), sorted by CTRL peak

% --- sort using CTRL peak timing ---
[~, peak_idx_ctrl] = max(PsthMasterCtrl, [], 2);
[~, sort_idx_ctrl] = sort(peak_idx_ctrl);

% --- apply sorting ---
sorted_csm   = PsthMasterCSM(sort_idx_ctrl, :);
sorted_nocsm = PsthMasterNoCSM(sort_idx_ctrl, :);

% sorted_csm   = PsthMasterCSM(sort_idx_even, :); % even is CSM and odd is noCSM peak timing
% sorted_nocsm = PsthMasterNoCSM(sort_idx_even, :);

% --- delta ---
sorted_delta = sorted_csm - sorted_nocsm;


figure('Position',[200 200 700 600]);

imagesc(sorted_delta);
colormap(slanCM('RdBu'));   % diverging colormap
colorbar;

set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times);
xlabel('Time from threshold crossing (ms)');
ylabel('Units (sorted by Ctrl peak)');
title('\Delta PSTH (CSM - NoCSM)');
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
%xlabel('Time from Endpoint (ms)');
xlabel('Time from Max Outward Vel (ms)');
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
