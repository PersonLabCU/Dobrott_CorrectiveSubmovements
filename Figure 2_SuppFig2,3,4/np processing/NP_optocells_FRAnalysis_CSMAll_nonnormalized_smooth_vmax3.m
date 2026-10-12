%% Get data - COde for CSM FR analysis - bursters and decreases 


cd("C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\ShamAnimals\ReachandCSMMod_allCSMStogther\CSMAll_wtihMaxV")
MyFolderInfo = dir('*CSMAll.mat');

tic
clearvars -except MyFolderInfo


%% ALL REACHES PSTH + PEAK TIME ANALYSIS
NPdata = struct;

PsthMasterCtrl = [];
PsthMaster_CSM = [];
PsthMaster_NoRR = [];

PsthMaster_Raw_Ctrl = [];
PsthMaster_Raw_CSM = [];
PsthMaster_Raw_NoRR = [];

% NEW: Separate PSTH collections for before/after peak velocity
PsthMaster_Before_Ctrl = [];
PsthMaster_After_Ctrl = [];

% NEW: Separate for CSM and NoRR
PsthMaster_Before_CSM = [];
PsthMaster_After_CSM = [];
PsthMaster_Before_NoRR = [];
PsthMaster_After_NoRR = [];

for day = 1:length(MyFolderInfo)
    load(MyFolderInfo(day).name, 'Reaches', 'data', 'filtered_modulated_indices','neg_corr_cells','pos_corr_cells','FRBins','FRCorrStats_edit','ExcitRR')

    [Control, Stim] = StimFilter(Reaches);

    % Get modulated data
    data_opto_all = data.units(filtered_modulated_indices);

    % === NEW: Determine valid peak window between -100 and 1000 ms ===
    timevec = data.units(filtered_modulated_indices(1)).max_vx.time_bins;
   %timevec = data.units(filtered_modulated_indices(1)).end_point.time_bins;
   valid_peak_range = find(timevec >= -100 & timevec <= 100);
    %valid_peak_range = find(timevec >= -200 & timevec <= 200);
   % valid_peak_range = find(timevec >= -500 & timevec <= 500);

    % === Filter based on whether peak is within -100 to 100 ms ===
    keep_unit = false(length(data_opto_all), 1);
    peak_times_ctrl_all = NaN(length(data_opto_all), 1);

    for t = 1:length(data_opto_all)
        psth = data_opto_all(t).max_vx.PSTH_smoothbin;
        %psth = data_opto_all(t).end_point.PSTH_smoothbin;
        norm_psth = psth;
        [~, peak_idx] = max(norm_psth);
        peak_times_ctrl_all(t) = peak_idx;

        if ismember(peak_idx, valid_peak_range)
            keep_unit(t) = true;
        end
    end

    % === Keep only those units with peak in window ===
    data_opto = data_opto_all(keep_unit);
    peak_times_ctrl = peak_times_ctrl_all(keep_unit);

    if isempty(data_opto)
        fprintf('Day %d: No modulated units with peak between -100 and 100 ms. Skipping...\n', day);
        continue
    end

    % Convert peak indices to actual time values (for kept units)
    peak_times_ms = timevec(peak_times_ctrl);  % in ms
    before_peakv_mask = (peak_times_ms < 0);
    after_peakv_mask = (peak_times_ms >= 0);

    % Global indices
    kept_indices_global = filtered_modulated_indices(keep_unit);
    before_peakv_indices = kept_indices_global(before_peakv_mask);
    after_peakv_indices = kept_indices_global(after_peakv_mask);

    % Preallocate
    psth_master_ctrl = NaN(length(data_opto), length(data_opto(1).max_vx.PSTH_smoothbin));
    psth_master_raw_ctrl = NaN(length(data_opto), length(data_opto(1).max_vx.PSTH_rawbin));

    for t = 1:length(data_opto)
        %psth = data_opto(t).max_vx.PSTH_smoothbin;
        %psth = data_opto(t).end_point.PSTH_smoothbin;
        psth = data_opto(t).pos_thresh.PSTH_smoothbin;
        psth=psth/max(psth);
        psth_raw = data_opto(t).max_vx.PSTH_rawbin;
        psth_master_ctrl(t,:) = psth;
        psth_master_raw_ctrl(t,:) = psth_raw / max(psth_raw);
    end

    PsthMaster_Before_Ctrl = vertcat(PsthMaster_Before_Ctrl, psth_master_ctrl(before_peakv_mask,:));
    PsthMaster_After_Ctrl = vertcat(PsthMaster_After_Ctrl, psth_master_ctrl(after_peakv_mask,:));

    % NoCSM PSTHs
    psth_master_norr = NaN(length(data_opto), length(data_opto(1).NoCSM_maxvx.PSTH_smoothbin));
    psth_master_raw_norr = NaN(length(data_opto), length(data_opto(1).NoCSM_maxvx.PSTH_rawbin));
    peak_times_norr = NaN(length(data_opto), 1);

    for g = 1:length(data_opto)
       %psthnorr = data_opto(g).NoCSM_maxvx.PSTH_smoothbin;
       psthnorr = data_opto(g).NoCSM_pos_thresh.PSTH_smoothbin;
       %psthnorr = data_opto(g).NoCSM_csmT2.PSTH_smoothbin;
       %psthnorr = data_opto(g).NoCSM_csmT0.PSTH_smoothbin;
       % psthnorr = psthnorr/max(psthnorr);
        psth_rawnorr = data_opto(g).NoCSM_maxvx.PSTH_rawbin;

        psth_master_norr(g,:) = psthnorr;
        psth_master_raw_norr(g,:) = psth_rawnorr / max(psth_rawnorr);
        [~, peak_idxnorr] = max(psthnorr);
        peak_times_norr(g) = peak_idxnorr;
    end

    % NEW: Separate NoRR by before/after peak
    PsthMaster_Before_NoRR = vertcat(PsthMaster_Before_NoRR, psth_master_norr(before_peakv_mask,:));
    PsthMaster_After_NoRR = vertcat(PsthMaster_After_NoRR, psth_master_norr(after_peakv_mask,:));

    % CSM PSTHs
    psth_master_csm = NaN(length(data_opto), length(data_opto(1).AllCSM_maxvx.PSTH_smoothbin));
    psth_master_raw_csm = NaN(length(data_opto), length(data_opto(1).AllCSM_maxvx.PSTH_rawbin));
    peak_times_csm = NaN(length(data_opto), 1);

    for b = 1:length(data_opto)
      % psthcsm = data_opto(b).AllCSM_maxvx.PSTH_smoothbin;
       psthcsm = data_opto(b).AllCSM_pos_thresh.PSTH_smoothbin;
        %psthcsm = data_opto(b).AllCSM_csmT2.PSTH_smoothbin;
        %psthcsm = data_opto(b).AllCSM_csmT0.PSTH_smoothbin;
       % psthcsm=psthcsm/max(psthcsm);
        psth_rawcsm = data_opto(b).AllCSM_maxvx.PSTH_rawbin;

        psth_master_csm(b,:) = psthcsm;
        psth_master_raw_csm(b,:) = psth_rawcsm / max(psth_rawcsm);
        [~, peak_idxcsm] = max(psthcsm);
        peak_times_csm(b) = peak_idxcsm;
    end

    % NEW: Separate CSM by before/after peak
    PsthMaster_Before_CSM = vertcat(PsthMaster_Before_CSM, psth_master_csm(before_peakv_mask,:));
    PsthMaster_After_CSM = vertcat(PsthMaster_After_CSM, psth_master_csm(after_peakv_mask,:));

    % Store in structure
    NPdata(day).peak_times_ctrl = peak_times_ctrl;
    NPdata(day).psth_master_ctrl = psth_master_ctrl;
    PsthMasterCtrl = vertcat(PsthMasterCtrl, psth_master_ctrl);
    PsthMaster_Raw_Ctrl = vertcat(PsthMaster_Raw_Ctrl, psth_master_raw_ctrl);

    NPdata(day).peak_times_csm = peak_times_csm;
    NPdata(day).psth_master_csm = psth_master_csm;
    PsthMaster_CSM = vertcat(PsthMaster_CSM, psth_master_csm);
    PsthMaster_Raw_CSM = vertcat(PsthMaster_Raw_CSM, psth_master_raw_csm);

    NPdata(day).peak_times_norr = peak_times_norr;
    NPdata(day).psth_master_norr = psth_master_norr;
    PsthMaster_NoRR = vertcat(PsthMaster_NoRR, psth_master_norr);
    PsthMaster_Raw_NoRR = vertcat(PsthMaster_Raw_NoRR, psth_master_raw_norr);

    NPdata(day).cells = data;
    NPdata(day).Reaches = Reaches;
    NPdata(day).ControlR = Control;

    NPdata(day).NoCSM = ExcitRR.NoCSM;
    NPdata(day).ReReaches = ExcitRR.ReReaches;
    NPdata(day).ExcitRR = ExcitRR;

    NPdata(day).Reachmodidx = filtered_modulated_indices;
    NPdata(day).kept_indices = kept_indices_global;
    NPdata(day).before_peakv_indices = before_peakv_indices;
    NPdata(day).after_peakv_indices = after_peakv_indices;

    fprintf('Day %d: Kept %d/%d modulated units with peak between -100 and 100 ms\n', ...
            day, sum(keep_unit), length(keep_unit));

    clear data Reaches Control filtered_modulated_indices neg_corr_cells pos_corr_cells FRBins FRCorrStats_edit ExcitRR
end
%% Plotting individual cells - before Peak V
% === PLOT INDIVIDUAL PSTHs ===
folder_path="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\ShamAnimals\Figures\CSMaligned_072925\T2CSMAligned_indicidualcells_before";
% Assuming timevec is the same as used before

num_bins = size(PsthMaster_Before_CSM, 2);
timevec = linspace(-1000, 1000, num_bins);

nCells = size(PsthMaster_Before_CSM, 1);
csmcolor=rgb('SlateBlue');
nocsmcolor=rgb('Black');

smooth_window = 3;

for i = 1:nCells
%for i=1;
    figure;
    hold on;

    psth_csm = movmean(PsthMaster_Before_CSM(i,:), smooth_window);
    psth_norr = movmean(PsthMaster_Before_NoRR(i,:), smooth_window);
    
    % Plot CSM PSTH
    plot(timevec, psth_csm, 'Color', csmcolor, 'LineWidth', 1.5);
    
    % Plot NoRR PSTH
    plot(timevec, psth_norr, 'Color', nocsmcolor, 'LineWidth', 1.5);
    
    
    % % Plot CSM PSTH
    % plot(timevec, PsthMaster_Before_CSM(i,:),'Color',csmcolor, 'LineWidth', 1.5); 
    % % Plot NoRR PSTH
    % plot(timevec, PsthMaster_Before_NoRR(i,:),'Color', nocsmcolor, 'LineWidth', 1.5);
    
    
    % Labels and titles
    xlim([-250 250])
    xlabel('Time (ms)');
    ylabel('Firing Rate (Hz)');
    title(sprintf('Cell %d PSTH (Before Peak)', i));
    legend({'Before CSM', 'Before NoRR'}, 'Location', 'Best');
    
    grid off;
    hold off;
    % 
     filename=sprintf('BeforeV_CSMvsNoCSM%d.svg',i);

    full_path=fullfile(folder_path,filename);
    saveas(gcf,full_path);

    close;
end
%% Plotting individual cells - After Peak V 
% === PLOT INDIVIDUAL PSTHs ===
folder_path="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\ShamAnimals\Figures\CSMaligned_072925\T2CSMAligned_indicidualcells_after";
% Assuming timevec is the same as used before

num_bins = size(PsthMaster_After_CSM, 2);
timevec = linspace(-1000, 1000, num_bins);

nCells = size(PsthMaster_After_CSM, 1);
csmcolor=rgb('SlateBlue');
nocsmcolor=rgb('Black');

smooth_window = 3;

for i = 1:nCells
    figure;
    hold on;

    psth_csm = movmean(PsthMaster_After_CSM(i,:), smooth_window);
    psth_norr = movmean(PsthMaster_After_NoRR(i,:), smooth_window);
    
    % Plot CSM PSTH
    plot(timevec, psth_csm, 'Color', csmcolor, 'LineWidth', 1.5);
    
    % Plot NoRR PSTH
    plot(timevec, psth_norr, 'Color', nocsmcolor, 'LineWidth', 1.5);
    
    
    % % Plot CSM PSTH
    % plot(timevec, PsthMaster_Before_CSM(i,:),'Color',csmcolor, 'LineWidth', 1.5); 
    % % Plot NoRR PSTH
    % plot(timevec, PsthMaster_Before_NoRR(i,:),'Color', nocsmcolor, 'LineWidth', 1.5);
    
    
    % Labels and titles
    xlim([-250 250])
    xlabel('Time (ms)');
    ylabel('Firing Rate (Hz)');
    title(sprintf('Cell %d PSTH (After Peak)', i));
    legend({'Before CSM', 'Before NoRR'}, 'Location', 'Best');
    
    grid off;
    hold off;

     filename=sprintf('AfterV_CSMvsNoCSM%d.svg',i);
 
    full_path=fullfile(folder_path,filename);
    saveas(gcf,full_path);

    close;
end
%% Calculate the SEM for all the PSTHS for individual cells 
nBoot = 1000;
bin_size = 10; % 10ms bin size
smoothing_bins = 5; % Moving average over 3 bins

%BinaryMaster = NPdata(1).cells.units(1).NoCSM_csmT2.Raster; %before V
%BinaryMaster = NPdata(3).cells.units(46).NoCSM_csmT2.Raster; %before V
%BinaryMaster = NPdata(3).cells.units(46).NoCSM_csmT0.Raster; %before V
%BinaryMaster = NPdata(3).cells.units(46).NoCSM_pos_thresh.Raster; %before V
%BinaryMaster = NPdata(14).cells.units(7).NoCSM_pos_thresh.Raster; %before V
%BinaryMaster = NPdata(4).cells.units(35).NoCSM_pos_thresh.Raster; %before V
BinaryMaster = NPdata(8).cells.units(1).NoCSM_pos_thresh.Raster; %before V

%BinaryMaster = NPdata(1).cells.units(6).NoCSM_csmT2.Raster; %after V
%BinaryMaster = NPdata(1).cells.units(6).NoCSM_pos_thresh.Raster; %after V
%BinaryMaster = NPdata(8).cells.units(11).NoCSM_pos_thresh.Raster; %after V
%BinaryMaster = NPdata(8).cells.units(9).NoCSM_pos_thresh.Raster; %after V

[num_events, total_time] = size(BinaryMaster);

% Number of new bins after downsampling
num_bins = floor(total_time / bin_size);

% Preallocate binned BinaryMaster
BinaryBinned = zeros(num_events, num_bins);

% Aggregate data into bin_size chunks
for b = 1:num_bins
    idx_start = (b - 1) * bin_size + 1;
    idx_end = b * bin_size;
    BinaryBinned(:, b) = sum(BinaryMaster(:, idx_start:idx_end), 2);
end



% Also compute bootstrapped SEM
% Let's use 'same' to avoid shrinking
[mean_psth_boot, sem_psth_boot, time_bins] = bootstrapPSTH_smooth(BinaryBinned, bin_size, smoothing_bins, nBoot);



%% Calculate the SEM for all the PSTHS for individual cells 
nBoot = 1000;
bin_size = 10; % 10ms bin size
smoothing_bins = 5; % Moving average over 3 bins

%BinaryMaster = NPdata(1).cells.units(1).AllCSM_csmT2.Raster; %before
%BinaryMaster = NPdata(1).cells.units(6).AllCSM_csmT2.Raster;%after
%BinaryMaster = NPdata(3).cells.units(46).AllCSM_csmT2.Raster; %before
%BinaryMaster = NPdata(3).cells.units(46).AllCSM_pos_thresh.Raster; %before
%BinaryMaster = NPdata(3).cells.units(46).AllCSM_pos_thresh.Raster; %before
%BinaryMaster = NPdata(14).cells.units(7).AllCSM_pos_thresh.Raster; %before
%BinaryMaster = NPdata(4).cells.units(34).AllCSM_pos_thresh.Raster; %before
BinaryMaster = NPdata(8).cells.units(1).AllCSM_pos_thresh.Raster; %before

%BinaryMaster = NPdata(1).cells.units(6).AllCSM_pos_thresh.Raster; %before
%BinaryMaster = NPdata(8).cells.units(11).AllCSM_pos_thresh.Raster; %after
%BinaryMaster = NPdata(8).cells.units(9).AllCSM_pos_thresh.Raster; %after

[num_events, total_time] = size(BinaryMaster);

% Number of new bins after downsampling
num_bins = floor(total_time / bin_size);

% Preallocate binned BinaryMaster
BinaryBinned = zeros(num_events, num_bins);

% Aggregate data into bin_size chunks
for b = 1:num_bins
    idx_start = (b - 1) * bin_size + 1;
    idx_end = b * bin_size;
    BinaryBinned(:, b) = sum(BinaryMaster(:, idx_start:idx_end), 2);
end


% Also compute bootstrapped SEM
% Let's use 'same' to avoid shrinking
[mean_psth_boot_CSM, sem_psth_boot_CSM, time_bins_CSM] = bootstrapPSTH_smooth(BinaryBinned, bin_size, smoothing_bins, nBoot);




%%
%figure;hold on;
figure('Position', [100, 100, 600, 300]);hold on

% Define SlateBlue color
%slateBlue = [106, 90, 205] / 255;
slateBlue=rgb('DarkMagenta');

% Plot NoCSM mean PSTH with shaded SEM (black)
shadedErrorBar(time_bins, mean_psth_boot, sem_psth_boot, ...
    'lineprops', {'-k', 'LineWidth', 1.5}, 'transparent', true);

% Plot CSM mean PSTH with shaded SEM (SlateBlue)
shadedErrorBar(time_bins_CSM, mean_psth_boot_CSM, sem_psth_boot_CSM, ...
    'lineprops', {'-', 'Color', slateBlue, 'LineWidth', 1.5}, 'transparent', true);

% Formatting
%xlabel('Time from CSM (ms)');
xlabel('Time from threshold crossing (ms)');
ylabel('Firing Rate (Hz)');
title('Bootstrapped Smoothed PSTH');
xline(0, '--r', 'Event Onset'); % Mark event time (0 ms)
xline(-55)
xline(87)

set(gca, 'FontSize', 12, 'LineWidth', 1);

% Set axis limits: -250 to 250 ms
xlim([-250, 250]);
%ylim([0 110])
%ylim([0, max([mean_psth_boot + sem_psth_boot, mean_psth_boot_CSM + sem_psth_boot_CSM]) * 1.2]);

box off;

%% Calculate the SEM for all the PSTHS for individual cells - AFTER 
nBoot = 1000;
bin_size = 10; % 10ms bin size
smoothing_bins = 3; % Moving average over 3 bins

%BinaryMaster = NPdata(8).cells.units(11).NoCSM_csmT2.Raster;
BinaryMaster = NPdata(8).cells.units(11).NoCSM_pos_thresh.Raster;
%BinaryMaster = NPdata(8).cells.units(8).NoCSM_csmT2.Raster;
%BinaryMaster = NPdata(4).cells.units(2).NoCSM_csmT2.Raster;
%BinaryMaster = NPdata(4).cells.units(2).NoCSM_pos_thresh.Raster;

[num_events, total_time] = size(BinaryMaster);

% Number of new bins after downsampling
num_bins = floor(total_time / bin_size);

% Preallocate binned BinaryMaster
BinaryBinned = zeros(num_events, num_bins);

% Aggregate data into bin_size chunks
for b = 1:num_bins
    idx_start = (b - 1) * bin_size + 1;
    idx_end = b * bin_size;
    BinaryBinned(:, b) = sum(BinaryMaster(:, idx_start:idx_end), 2);
end

% Also compute bootstrapped SEM
% Let's use 'same' to avoid shrinking
[mean_psth_boot, sem_psth_boot, time_bins] = bootstrapPSTH_smooth(BinaryBinned, bin_size, smoothing_bins, nBoot);

%% Calculate the SEM for all the PSTHS for individual cells 
nBoot = 1000;
bin_size = 10; % 10ms bin size
smoothing_bins = 3; % Moving average over 3 bins

%BinaryMaster = NPdata(8).cells.units(11).AllCSM_csmT2.Raster;
BinaryMaster = NPdata(8).cells.units(11).AllCSM_pos_thresh.Raster;
%BinaryMaster = NPdata(8).cells.units(8).AllCSM_csmT2.Raster;
%BinaryMaster = NPdata(4).cells.units(2).AllCSM_csmT2.Raster;
%BinaryMaster = NPdata(4).cells.units(2).AllCSM_pos_thresh.Raster;

[num_events, total_time] = size(BinaryMaster);

% Number of new bins after downsampling
num_bins = floor(total_time / bin_size);

% Preallocate binned BinaryMaster
BinaryBinned = zeros(num_events, num_bins);

% Aggregate data into bin_size chunks
for b = 1:num_bins
    idx_start = (b - 1) * bin_size + 1;
    idx_end = b * bin_size;
    BinaryBinned(:, b) = sum(BinaryMaster(:, idx_start:idx_end), 2);
end

% Also compute bootstrapped SEM
% Let's use 'same' to avoid shrinking
[mean_psth_boot_CSM, sem_psth_boot_CSM, time_bins_CSM] = bootstrapPSTH_smooth(BinaryBinned, bin_size, smoothing_bins, nBoot);


%%
slateBlue=rgb('DarkMagenta');
%figure;hold on;
figure('Position', [100, 100, 600, 300]);hold on

% Define SlateBlue color
slateBlue = [106, 90, 205] / 255;

% Plot NoCSM mean PSTH with shaded SEM (black)
shadedErrorBar(time_bins, mean_psth_boot, sem_psth_boot, ...
    'lineprops', {'-k', 'LineWidth', 1.5}, 'transparent', true);

% Plot CSM mean PSTH with shaded SEM (SlateBlue)
shadedErrorBar(time_bins_CSM, mean_psth_boot_CSM, sem_psth_boot_CSM, ...
    'lineprops', {'-', 'Color', slateBlue, 'LineWidth', 1.5}, 'transparent', true);

% Formatting
%xlabel('Time from CSM (ms)');
xlabel('Time from threshold crossing (ms)');
ylabel('Firing Rate (Hz)');
title('Bootstrapped Smoothed PSTH');
xline(0, '--r', 'Event Onset'); % Mark event time (0 ms)
xline(-55)
xline(87)
set(gca, 'FontSize', 12, 'LineWidth', 1);

% Set axis limits: -250 to 250 ms
xlim([-250, 250]);
ylim([0, max([mean_psth_boot + sem_psth_boot, mean_psth_boot_CSM + sem_psth_boot_CSM]) * 1.2]);

box off;

%%
% === Peak Time Histogram Plotting ===

% Initialize arrays to collect all peak times
all_peak_times_ctrl = [];
all_peak_times_csm = [];
all_peak_times_norr = [];

% Loop over NPdata to collect all peak time indices and timevec
for day = 1:length(NPdata)
    if isfield(NPdata(day), 'peak_times_ctrl')
        all_peak_times_ctrl = [all_peak_times_ctrl; NPdata(day).peak_times_ctrl];
        all_peak_times_csm = [all_peak_times_csm; NPdata(day).peak_times_csm];
        all_peak_times_norr = [all_peak_times_norr; NPdata(day).peak_times_norr];
    end
end

% Get timevec (assuming it's the same across all days)

% Example: timevec from -500 to 1000 ms, 10 ms bins
timevec = -1000:10:1000;


% Convert peak time indices to actual time (ms)
peak_times_ctrl_ms = timevec(all_peak_times_ctrl);
peak_times_csm_ms = timevec(all_peak_times_csm);
peak_times_norr_ms = timevec(all_peak_times_norr);

% === Plotting ===

figure;  hold on;
histogram(peak_times_ctrl_ms(peak_times_ctrl_ms < 0), 'Normalization', 'probability', 'FaceColor', [0.2 0.6 1], 'FaceAlpha', 0.5);
histogram(peak_times_ctrl_ms(peak_times_ctrl_ms >= 0), 'Normalization', 'probability', 'FaceColor', [1 0.4 0.4], 'FaceAlpha', 0.5);
xline(0, '--k', 'Peak Velocity');
xlim([-250 250])
xlabel('Peak Time (ms)'); ylabel('Probability');
title('CTRL Units: Peak Time Distribution');
legend({'Before PeakV', 'After PeakV'});


%%

% Control
plot_sorted_psth(PsthMasterCtrl, 'Control PSTHs Sorted by Peak Firing Time');

plot_sorted_psth_window(PsthMasterCtrl, 'PSTH Sorted by Peak-allreaches', [-250, 250]);

% NoCSM 
plot_sorted_psth(PsthMaster_NoRR, 'NoCSM PSTHs Sorted by Peak Firing Time');

% NoCSM - active before peak velcotiy 
plot_sorted_psth(PsthMaster_Before_NoRR, 'NoCSM PSTHs Sorted by Peak Firing Time-active BEFORE peak');

% NoCSM - active before peak velcotiy 
plot_sorted_psth(PsthMaster_After_NoRR, 'NoCSM PSTHs Sorted by Peak Firing Time-active AFTER peak');

% CSM
plot_sorted_psth(PsthMaster_CSM, 'CSM PSTHs Sorted by Peak Firing Time');

% CSM - cells active BEFORE peak velocity 
plot_sorted_psth(PsthMaster_Before_CSM, 'CSM PSTHs Sorted by Peak Firing Time- active BEFORE peak');

% CSM - cells active AFTER peak velocity 
plot_sorted_psth(PsthMaster_After_CSM, 'CSM PSTHs Sorted by Peak Firing Time- active AFTER peak');

% Ctrl - peak activity before 
plot_sorted_psth(PsthMaster_Before_Ctrl, 'Ctrl PSTHs-Peak act. before PeakV');

% Ctrl - peak activity after 
plot_sorted_psth(PsthMaster_After_Ctrl, 'Ctrl PSTHs-Peak act. after PeakV');

%%
%%%% STEP 1: Calculate Peak Firing Time Indices for PsthMaster_NoRRT2
all_peak_times_norr = [];

% Loop through each day in NPdata
for i = 1:length(NPdata)
    if isfield(NPdata(i), 'peak_times_norr') && ~isempty(NPdata(i).peak_times_norr)
        all_peak_times_norr = [all_peak_times_norr; NPdata(i).peak_times_norr(:)];
    end
end

% Ensure number of peak times matches the number of units in PsthMaster_NoRRT2
if length(all_peak_times_norr) ~= size(PsthMaster_NoRR, 1)
    error('Mismatch between number of units and number of peak times for NoRRT2.');
end

% Sort the peak times and get the sorting index
[~, sort_idx_norr] = sort(all_peak_times_norr);



%% STEP 1: Calculate Peak Firing Time Indices for PsthMaster_CSM
all_peak_times_csm = [];

% Loop through each day in NPdata
for i = 1:length(NPdata)
    if isfield(NPdata(i), 'peak_times_norr') && ~isempty(NPdata(i).peak_times_csm)
        all_peak_times_csm = [all_peak_times_csm; NPdata(i).peak_times_csm(:)];
    end
end

% Ensure number of peak times matches the number of units in PsthMaster_NoRRT2
if length(all_peak_times_csm) ~= size(PsthMaster_CSM, 1)
    error('Mismatch between number of units and number of peak times for NoRRT2.');
end

% Sort the peak times and get the sorting index
[~, sort_idx_csm] = sort(all_peak_times_csm);


%% STEP 2: Sort PsthMaster_NoRRT2 Using the Sorting Index
sorted_psth_norr = PsthMaster_NoRR(sort_idx_norr, :);
%sorted_psth_norr = PsthMaster_NoRR(sort_idx_csm, :);

% Generate time vector
num_bins = size(sorted_psth_norr, 2);
time_vector = linspace(-1000, 1000, num_bins);  % assuming 10 ms bins over 200 bins
tick_times = -1000:250:1000;  % Set tick marks for x-axis
[~, tick_indices] = min(abs(time_vector' - tick_times), [], 1);  % find closest indices

% Plot the sorted NoRRT2 PSTH heatmap
figure;
imagesc(sorted_psth_norr);
colormap(slanCM('GnBu'));
colorbar;
set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times);
xlabel('Time (ms)');
ylabel('Units (sorted by peak firing time)');
title('NoRR PSTHs Sorted by Peak Firing Time');



%% STEP 3: Sort PsthMaster_T2CSM Using the Same Sorting Index
sorted_psth_csm = PsthMaster_CSM(sort_idx_norr, :);
%sorted_psth_csm = PsthMaster_CSM(sort_idx_csm, :);

% Plot the sorted T2 CSM PSTH heatmap
figure;
imagesc(sorted_psth_csm);
colormap(slanCM('GnBu'));
colorbar;
set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times);
xlabel('Time (ms)');
ylabel('Units (sorted by peak firing time)');
title('CSM PSTHs Sorted by Peak Firing Time- No CSM IDX');


%% STEP 4: Compute Difference and Sort
% Subtract T2CSM from NoRRT2
psth_diff = sorted_psth_csm - sorted_psth_norr;


% Create a time vector matching psth_diff
num_time_bins = size(psth_diff, 2);
time_vector = linspace(-500, 500, num_time_bins);  % Adjust range as needed

% % Now make the time mask
time_window = [-250, 250];
time_mask = time_vector >= time_window(1) & time_vector <= time_window(2);

% Subset all relevant arrays
psth_diff_window = psth_diff(:, time_mask);
tick_times_window = time_vector(time_mask);  % use same time vector for ticks

figure;
imagesc(psth_diff_window);
%imagesc(psth_diff;
colormap(flipud(slanCM('RdBu')));
colorbar;
caxis([-max(abs(psth_diff_window(:))), max(abs(psth_diff_window(:)))]);
set(gca, 'XTick', linspace(1, size(psth_diff_window, 2), 5), ...
         'XTickLabel', round(linspace(time_window(1), time_window(2), 5)));  % Optional custom ticks
xlabel('Time (ms)');
ylabel('Units (sorted by NoRRT2 peak firing time)');
title('Difference: CSM - NoCSM (Sorted by NoRRT2 Peaks)');

%%
CSMcolor = rgb('DarkOrchid');

% Compute average PSTHs (mean across units)
avg_psth_csm = mean(PsthMaster_CSM, 1);  % 1 = average across rows (units)
avg_psth_norr = mean(PsthMaster_NoRR, 1);

% Optional: Compute standard error of the mean (SEM)
sem_psth_csm = std(PsthMaster_CSM, 0, 1) / sqrt(size(PsthMaster_CSM, 1));
sem_psth_norr = std(PsthMaster_NoRR, 0, 1) / sqrt(size(PsthMaster_NoRR, 1));

% Generate time vector (assuming same binning and length)
num_bins = size(PsthMaster_CSM, 2);
time_vector = linspace(-1000, 1000, num_bins);  % Update if your window or binning is different

% Statistical testing (paired test across units)
p_vals = nan(1, num_bins);
for t = 1:num_bins
    [p_vals(t), ~] = signrank(PsthMaster_CSM(:, t), PsthMaster_NoRR(:, t));
end

% FDR correction (requires fdr_bh.m)
[h_fdr, ~, ~] = fdr_bh(p_vals);

% Plot average PSTHs
% figure;
% hold on;

figure('Position', [100, 100, 600, 300]); 
hold on;

% Plot mean lines
plot(time_vector, avg_psth_csm, 'Color', CSMcolor, 'LineWidth', 2);
plot(time_vector, avg_psth_norr, 'k', 'LineWidth', 2);

% Shaded error bars
fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_csm + sem_psth_csm, fliplr(avg_psth_csm - sem_psth_csm)], ...
     CSMcolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_norr + sem_psth_norr, fliplr(avg_psth_norr - sem_psth_norr)], ...
     'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Plot significance markers
sig_bins = find(p_vals<=0.05);  % Significant time bins
plot(time_vector(sig_bins), avg_psth_csm(sig_bins), '.', 'Color', 'k', 'MarkerSize', 10);

% Labels and limits
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
legend('CSM', 'NoCSM');
title('Average PSTH across Units with Significance');

% Optional y-limits
% ylim([0 1]);
xlim([-250 250]);
%% Before and after 

Beforecolor = rgb('MediumOrchid');
Aftercolor = rgb('Indigo');

% Compute average PSTHs (mean across units)
avg_psth_before = mean(PsthMaster_Before_Ctrl, 1);  % 1 = average across rows (units)
avg_psth_after = mean(PsthMaster_After_Ctrl, 1);

% Optional: Compute standard error of the mean (SEM)
sem_psth_before = std(PsthMaster_Before_Ctrl, 0, 1) / sqrt(size(PsthMaster_Before_Ctrl, 1));
sem_psth_after= std(PsthMaster_After_Ctrl, 0, 1) / sqrt(size(PsthMaster_After_Ctrl, 1));

% Generate time vector (assuming same binning and length)
num_bins = size(PsthMaster_After_Ctrl, 2);
time_vector = linspace(-1000, 1000, num_bins);  % Update if your window or binning is different

% Plot average PSTHs
figure;
%figure('Position', [100, 100, 600, 300]);  
hold on;

% Plot mean lines
plot(time_vector, avg_psth_before, 'Color', Beforecolor, 'LineWidth', 2);
plot(time_vector, avg_psth_after,'Color',Aftercolor, 'LineWidth', 2);

% Shaded error bars
fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_before + sem_psth_before, fliplr(avg_psth_before - sem_psth_before)], ...
     Beforecolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_after + sem_psth_after, fliplr(avg_psth_after - sem_psth_after)], ...
     Aftercolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');



% Labels and limits
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
legend('Before', 'After');
title('Average PSTH across Units with Significance');

% Optional y-limits
% ylim([0 1]);
%xlim([-250 250]);
%% Before and after - NOCSM

Beforecolornorr = rgb('CornflowerBlue');
Aftercolornorr = rgb('DarkBlue');
Beforecolorcsm = rgb('MediumOrchid');
Aftercolorcsm = rgb('Indigo');

% Compute average PSTHs (mean across units)
avg_psth_before_norr = mean(PsthMaster_Before_NoRR, 1);  % 1 = average across rows (units)
avg_psth_after_norr = mean(PsthMaster_After_NoRR, 1);

avg_psth_before_csm = mean(PsthMaster_Before_CSM, 1);  % 1 = average across rows (units)
avg_psth_after_csm = mean(PsthMaster_After_CSM, 1);

% Optional: Compute standard error of the mean (SEM)
sem_psth_before_norr = std(PsthMaster_Before_NoRR, 0, 1) / sqrt(size(PsthMaster_Before_NoRR, 1));
sem_psth_after_norr= std(PsthMaster_After_NoRR, 0, 1) / sqrt(size(PsthMaster_After_NoRR, 1));

sem_psth_before_csm = std(PsthMaster_Before_CSM, 0, 1) / sqrt(size(PsthMaster_Before_CSM, 1));
sem_psth_after_csm= std(PsthMaster_After_CSM, 0, 1) / sqrt(size(PsthMaster_After_CSM, 1));

% Generate time vector (assuming same binning and length)
num_bins = size(PsthMaster_After_NoRR, 2);
time_vector = linspace(-1000, 1000, num_bins);  % Update if your window or binning is different

% Plot average PSTHs
% figure;
% hold on;
figure('Position', [100, 100, 600, 300]); 
hold on;

% Plot mean lines
plot(time_vector, avg_psth_before_norr, 'Color', Beforecolornorr, 'LineWidth', 2);
plot(time_vector, avg_psth_after_norr,'Color',Aftercolornorr, 'LineWidth', 2);
plot(time_vector, avg_psth_before_csm, 'Color', Beforecolorcsm, 'LineWidth', 2);
plot(time_vector, avg_psth_after_csm,'Color',Aftercolorcsm, 'LineWidth', 2);

% Shaded error bars
fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_before_norr + sem_psth_before_norr, fliplr(avg_psth_before_norr - sem_psth_before_norr)], ...
     Beforecolornorr, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_after_norr + sem_psth_after_norr, fliplr(avg_psth_after_norr - sem_psth_after_norr)], ...
     Aftercolornorr, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_before_csm + sem_psth_before_csm, fliplr(avg_psth_before_csm - sem_psth_before_csm)], ...
     Beforecolorcsm, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_after_csm + sem_psth_after_csm, fliplr(avg_psth_after_csm - sem_psth_after_csm)], ...
     Aftercolorcsm, 'FaceAlpha', 0.2, 'EdgeColor', 'none');



% Labels and limits
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
legend('Before-NoCSM', 'After-NoCSM','Before-CSM', 'After-CSM');
title('Average PSTH across Units with Significance');

% Optional y-limits
% ylim([0 1]);
%xlim([-500 500]);
xlim([-250 250]);
%%
Beforecolornorr = rgb('CornflowerBlue');
Aftercolornorr = rgb('DarkBlue');
Beforecolorcsm = rgb('MediumOrchid');
Aftercolorcsm = rgb('Indigo');

% Compute average PSTHs (mean across units)
avg_psth_before_norr = mean(PsthMaster_Before_NoRR, 1);
avg_psth_after_norr = mean(PsthMaster_After_NoRR, 1);

avg_psth_before_csm = mean(PsthMaster_Before_CSM, 1);
avg_psth_after_csm = mean(PsthMaster_After_CSM, 1);

% Compute standard error of the mean (SEM)
sem_psth_before_norr = std(PsthMaster_Before_NoRR, 0, 1) / sqrt(size(PsthMaster_Before_NoRR, 1));
sem_psth_after_norr  = std(PsthMaster_After_NoRR, 0, 1)  / sqrt(size(PsthMaster_After_NoRR, 1));

sem_psth_before_csm = std(PsthMaster_Before_CSM, 0, 1) / sqrt(size(PsthMaster_Before_CSM, 1));
sem_psth_after_csm  = std(PsthMaster_After_CSM, 0, 1)  / sqrt(size(PsthMaster_After_CSM, 1));

% Generate time vector
num_bins = size(PsthMaster_After_NoRR, 2);
time_vector = linspace(-1000, 1000, num_bins);

% Paired tests: Before CSM vs. Before NoCSM; After CSM vs. After NoCSM
p_vals_before = nan(1, num_bins);
p_vals_after  = nan(1, num_bins);

for t = 1:num_bins
    [p_vals_before(t), ~] = signrank(PsthMaster_Before_CSM(:, t), PsthMaster_Before_NoRR(:, t));
    [p_vals_after(t), ~]  = signrank(PsthMaster_After_CSM(:, t), PsthMaster_After_NoRR(:, t));
end

% FDR correction
[h_fdr_before, ~, ~] = fdr_bh(p_vals_before);
[h_fdr_after, ~, ~]  = fdr_bh(p_vals_after);

% Plot average PSTHs
figure('Position', [100, 100, 600, 300]); 
hold on;

% Plot mean lines
plot(time_vector, avg_psth_before_norr, 'Color', Beforecolornorr, 'LineWidth', 2);
plot(time_vector, avg_psth_after_norr,  'Color', Aftercolornorr,  'LineWidth', 2);
plot(time_vector, avg_psth_before_csm, 'Color', Beforecolorcsm, 'LineWidth', 2);
plot(time_vector, avg_psth_after_csm,  'Color', Aftercolorcsm,  'LineWidth', 2);

% Shaded error bars
fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_before_norr + sem_psth_before_norr, fliplr(avg_psth_before_norr - sem_psth_before_norr)], ...
     Beforecolornorr, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_after_norr + sem_psth_after_norr, fliplr(avg_psth_after_norr - sem_psth_after_norr)], ...
     Aftercolornorr, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_before_csm + sem_psth_before_csm, fliplr(avg_psth_before_csm - sem_psth_before_csm)], ...
     Beforecolorcsm, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_after_csm + sem_psth_after_csm, fliplr(avg_psth_after_csm - sem_psth_after_csm)], ...
     Aftercolorcsm, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% --- Plot significance markers ---
sig_bins_before = find(p_vals_before <= 0.05);
sig_bins_after  = find(p_vals_after <= 0.05);

% Set y-position for significance markers
max_val = max([avg_psth_after_csm + sem_psth_after_csm, avg_psth_after_norr + sem_psth_after_norr]);
y_offset_before = max_val + 0.1 * max_val;
y_offset_after  = max_val + 0.2 * max_val;

plot(time_vector(sig_bins_before), y_offset_before * ones(size(sig_bins_before)), '.', 'Color', Beforecolorcsm, 'MarkerSize', 12);
plot(time_vector(sig_bins_after),  y_offset_after  * ones(size(sig_bins_after)),  '.', 'Color', Aftercolorcsm,  'MarkerSize', 12);

% Labels and limits
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
legend('Before-NoCSM', 'After-NoCSM', 'Before-CSM', 'After-CSM');
title('Average PSTH across Units with Significance');

xlim([-250 250]);
%%
Beforecolornorr = rgb('CornflowerBlue');
Aftercolornorr = rgb('DarkBlue');
Beforecolorcsm = rgb('MediumOrchid');
Aftercolorcsm = rgb('Indigo');

% Compute average PSTHs (mean across units)
avg_psth_before_norr = mean(PsthMaster_Before_NoRR, 1);
avg_psth_after_norr  = mean(PsthMaster_After_NoRR, 1);
avg_psth_before_csm  = mean(PsthMaster_Before_CSM, 1);
avg_psth_after_csm   = mean(PsthMaster_After_CSM, 1);

% Compute SEM
sem_psth_before_norr = std(PsthMaster_Before_NoRR, 0, 1) / sqrt(size(PsthMaster_Before_NoRR, 1));
sem_psth_after_norr  = std(PsthMaster_After_NoRR, 0, 1)  / sqrt(size(PsthMaster_After_NoRR, 1));
sem_psth_before_csm  = std(PsthMaster_Before_CSM, 0, 1)  / sqrt(size(PsthMaster_Before_CSM, 1));
sem_psth_after_csm   = std(PsthMaster_After_CSM, 0, 1)   / sqrt(size(PsthMaster_After_CSM, 1));

% Time vector
num_bins = size(PsthMaster_After_NoRR, 2);
time_vector = linspace(-1000, 1000, num_bins);

% Paired tests: Before CSM vs. Before NoCSM; After CSM vs. After NoCSM
p_vals_before = nan(1, num_bins);
p_vals_after  = nan(1, num_bins);

for t = 1:num_bins
    [p_vals_before(t), ~] = signrank(PsthMaster_Before_CSM(:, t), PsthMaster_Before_NoRR(:, t));
    [p_vals_after(t), ~]  = signrank(PsthMaster_After_CSM(:, t), PsthMaster_After_NoRR(:, t));
end

% Plot PSTHs
figure('Position', [100, 100, 600, 300]); 
hold on;

% Plot mean lines
plot(time_vector, avg_psth_before_norr, 'Color', Beforecolornorr, 'LineWidth', 2);
plot(time_vector, avg_psth_after_norr,  'Color', Aftercolornorr,  'LineWidth', 2);
plot(time_vector, avg_psth_before_csm, 'Color', Beforecolorcsm, 'LineWidth', 2);
plot(time_vector, avg_psth_after_csm,  'Color', Aftercolorcsm,  'LineWidth', 2);

% Shaded error bars
fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_before_norr + sem_psth_before_norr, fliplr(avg_psth_before_norr - sem_psth_before_norr)], ...
     Beforecolornorr, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_after_norr + sem_psth_after_norr, fliplr(avg_psth_after_norr - sem_psth_after_norr)], ...
     Aftercolornorr, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_before_csm + sem_psth_before_csm, fliplr(avg_psth_before_csm - sem_psth_before_csm)], ...
     Beforecolorcsm, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_after_csm + sem_psth_after_csm, fliplr(avg_psth_after_csm - sem_psth_after_csm)], ...
     Aftercolorcsm, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Plot significance: dots, lines, and p-values
sig_bins_before = find(p_vals_before <= 0.05);
sig_bins_after  = find(p_vals_after <= 0.05);

for i = 1:length(sig_bins_before)
    t_idx = sig_bins_before(i);
    x_val = time_vector(t_idx);
    y_val = avg_psth_before_csm(t_idx); % Dot at CSM PSTH
    % Dot
    plot(x_val, y_val, '.', 'Color', Beforecolorcsm, 'MarkerSize', 12);
    % Line down
    plot([x_val x_val], [0 y_val], 'Color', Beforecolorcsm, 'LineWidth', 1);
    % Text
    text(x_val, y_val + 0.05*max(y_val, 1), sprintf('p=%.3f', p_vals_before(t_idx)), ...
        'HorizontalAlignment', 'center', 'FontSize', 8);
end

for i = 1:length(sig_bins_after)
    t_idx = sig_bins_after(i);
    x_val = time_vector(t_idx);
    y_val = avg_psth_after_csm(t_idx); % Dot at CSM PSTH
    % Dot
    plot(x_val, y_val, '.', 'Color', Aftercolorcsm, 'MarkerSize', 12);
    % Line down
    plot([x_val x_val], [0 y_val], 'Color', Aftercolorcsm, 'LineWidth', 1);
    % Text
    text(x_val, y_val + 0.05*max(y_val, 1), sprintf('p=%.3f', p_vals_after(t_idx)), ...
        'HorizontalAlignment', 'center', 'FontSize', 8);
end

% Labels
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
legend('Before-NoCSM', 'After-NoCSM', 'Before-CSM', 'After-CSM');
title('Average PSTH across Units with p-values');

xlim([-250 250]);
ylim([40 130])

%% Diff version
Beforecolornorr = rgb('CornflowerBlue');
Aftercolornorr = rgb('DarkBlue');
Beforecolorcsm = rgb('MediumOrchid');
Aftercolorcsm = rgb('Indigo');

% Compute average PSTHs (mean across units)
avg_psth_before_norr = mean(PsthMaster_Before_NoRR, 1);
avg_psth_after_norr = mean(PsthMaster_After_NoRR, 1);
avg_psth_before_csm = mean(PsthMaster_Before_CSM, 1);
avg_psth_after_csm = mean(PsthMaster_After_CSM, 1);

% Compute SEM
sem_psth_before_norr = std(PsthMaster_Before_NoRR, 0, 1) / sqrt(size(PsthMaster_Before_NoRR, 1));
sem_psth_after_norr  = std(PsthMaster_After_NoRR, 0, 1)  / sqrt(size(PsthMaster_After_NoRR, 1));
sem_psth_before_csm = std(PsthMaster_Before_CSM, 0, 1) / sqrt(size(PsthMaster_Before_CSM, 1));
sem_psth_after_csm  = std(PsthMaster_After_CSM, 0, 1)  / sqrt(size(PsthMaster_After_CSM, 1));

% Generate time vector
num_bins = size(PsthMaster_After_NoRR, 2);
time_vector = linspace(-1000, 1000, num_bins);

% Paired tests: Before CSM vs. Before NoCSM; After CSM vs. After NoCSM
p_vals_before = nan(1, num_bins);
p_vals_after  = nan(1, num_bins);
for t = 1:num_bins
    [p_vals_before(t), ~] = signrank(PsthMaster_Before_CSM(:, t), PsthMaster_Before_NoRR(:, t));
    [p_vals_after(t), ~]  = signrank(PsthMaster_After_CSM(:, t), PsthMaster_After_NoRR(:, t));
end

% Create figure
figure('Position', [100, 100, 900, 500]);
hold on;

% Plot PSTH mean lines
plot(time_vector, avg_psth_before_norr, 'Color', Beforecolornorr, 'LineWidth', 2);
plot(time_vector, avg_psth_after_norr,  'Color', Aftercolornorr,  'LineWidth', 2);
plot(time_vector, avg_psth_before_csm, 'Color', Beforecolorcsm, 'LineWidth', 2);
plot(time_vector, avg_psth_after_csm,  'Color', Aftercolorcsm,  'LineWidth', 2);

% Plot shaded SEM areas
fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_before_norr + sem_psth_before_norr, fliplr(avg_psth_before_norr - sem_psth_before_norr)], ...
     Beforecolornorr, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_after_norr + sem_psth_after_norr, fliplr(avg_psth_after_norr - sem_psth_after_norr)], ...
     Aftercolornorr, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_before_csm + sem_psth_before_csm, fliplr(avg_psth_before_csm - sem_psth_before_csm)], ...
     Beforecolorcsm, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_after_csm + sem_psth_after_csm, fliplr(avg_psth_after_csm - sem_psth_after_csm)], ...
     Aftercolorcsm, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Move heatmap above traces
yl = ylim;
heatmap_y_start = yl(2) + 0.02*diff(yl);
heatmap_height = 0.05 * diff(yl);
bin_width = time_vector(2) - time_vector(1);

% Get colormaps - use slanCM 'OrRd' for Before and 'PuRd' for After
before_cmap = slanCM('Greens', 256);
after_cmap = slanCM('Blues', 256);

% Plot Before p-value heatmap
for t = 1:num_bins
    p_val = min(p_vals_before(t), 0.05); % clip at 0.05
    color_idx = max(1, round((1 - p_val/0.05) * 255) + 1);
    colorVal = before_cmap(color_idx, :);
    x_patch = [time_vector(t)-bin_width/2, time_vector(t)+bin_width/2, time_vector(t)+bin_width/2, time_vector(t)-bin_width/2];
    y_patch = [heatmap_y_start, heatmap_y_start, heatmap_y_start + heatmap_height, heatmap_y_start + heatmap_height];
    fill(x_patch, y_patch, colorVal, 'EdgeColor', 'none');
end

% Plot After p-value heatmap above that
heatmap_y_start2 = heatmap_y_start + 1.2*heatmap_height;
for t = 1:num_bins
    p_val = min(p_vals_after(t), 0.05); % clip at 0.05
    color_idx = max(1, round((1 - p_val/0.05) * 255) + 1);
    colorVal = after_cmap(color_idx, :);
    x_patch = [time_vector(t)-bin_width/2, time_vector(t)+bin_width/2, time_vector(t)+bin_width/2, time_vector(t)-bin_width/2];
    y_patch = [heatmap_y_start2, heatmap_y_start2, heatmap_y_start2 + heatmap_height, heatmap_y_start2 + heatmap_height];
    fill(x_patch, y_patch, colorVal, 'EdgeColor', 'none');
end

% Labels for heatmaps
text(min(time_vector), heatmap_y_start + heatmap_height/2, 'p-value Before', ...
    'FontSize', 9, 'VerticalAlignment', 'middle', 'Color', Beforecolorcsm);
text(min(time_vector), heatmap_y_start2 + heatmap_height/2, 'p-value After', ...
    'FontSize', 9, 'VerticalAlignment', 'middle', 'Color', Aftercolorcsm);

% Legend
%legend('Before-NoCSM', 'After-NoCSM', 'Before-CSM', 'After-CSM', 'Location', 'northwest');
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
title('Average PSTH across Units with p-value heatmap');
xlim([-250 250]);


% Separate colorbars for Before and After heatmaps

% Before p-value colorbar (red/orange)
cb_ax_red = axes('Position', [0.92 0.55 0.02 0.35]);
imagesc(cb_ax_red, [0 1], [0 0.05], linspace(0.05, 0, 256)');
set(cb_ax_red, 'YDir', 'normal', 'XTick', [], 'YTick', [0 0.025 0.05]);
colormap(cb_ax_red, before_cmap);
ylabel(cb_ax_red, 'p-value Before');
caxis(cb_ax_red, [0 0.05]);
box(cb_ax_red, 'off');

% After p-value colorbar (purple/blue)
cb_ax_blue = axes('Position', [0.92 0.1 0.02 0.35]);
imagesc(cb_ax_blue, [0 1], [0 0.05], linspace(0.05, 0, 256)');
set(cb_ax_blue, 'YDir', 'normal', 'XTick', [], 'YTick', [0 0.025 0.05]);
colormap(cb_ax_blue, after_cmap);
ylabel(cb_ax_blue, 'p-value After');
caxis(cb_ax_blue, [0 0.05]);
box(cb_ax_blue, 'off');

% Return focus to main axes
axes(gca);


%% Diff version
Beforecolornorr = rgb('CornflowerBlue');
Aftercolornorr = rgb('DarkBlue');
Beforecolorcsm = rgb('MediumOrchid');
Aftercolorcsm = rgb('Indigo');

% Compute average PSTHs (mean across units)
avg_psth_before_norr = mean(PsthMaster_Before_NoRR, 1);
avg_psth_after_norr = mean(PsthMaster_After_NoRR, 1);
avg_psth_before_csm = mean(PsthMaster_Before_CSM, 1);
avg_psth_after_csm = mean(PsthMaster_After_CSM, 1);

% Compute SEM
sem_psth_before_norr = std(PsthMaster_Before_NoRR, 0, 1) / sqrt(size(PsthMaster_Before_NoRR, 1));
sem_psth_after_norr  = std(PsthMaster_After_NoRR, 0, 1)  / sqrt(size(PsthMaster_After_NoRR, 1));
sem_psth_before_csm = std(PsthMaster_Before_CSM, 0, 1) / sqrt(size(PsthMaster_Before_CSM, 1));
sem_psth_after_csm  = std(PsthMaster_After_CSM, 0, 1)  / sqrt(size(PsthMaster_After_CSM, 1));

% Generate time vector
num_bins = size(PsthMaster_After_NoRR, 2);
time_vector = linspace(-1000, 1000, num_bins);

% Paired tests: Before CSM vs. Before NoCSM; After CSM vs. After NoCSM
p_vals_before = nan(1, num_bins);
p_vals_after  = nan(1, num_bins);
for t = 1:num_bins
    [p_vals_before(t), ~] = signrank(PsthMaster_Before_CSM(:, t), PsthMaster_Before_NoRR(:, t));
    [p_vals_after(t), ~]  = signrank(PsthMaster_After_CSM(:, t), PsthMaster_After_NoRR(:, t));
end

% Create figure
%figure('Position', [100, 100, 900, 500]);
figure;
hold on;

% Plot PSTH mean lines
plot(time_vector, avg_psth_before_norr, 'Color', Beforecolornorr, 'LineWidth', 2);
plot(time_vector, avg_psth_after_norr,  'Color', Aftercolornorr,  'LineWidth', 2);
plot(time_vector, avg_psth_before_csm, 'Color', Beforecolorcsm, 'LineWidth', 2);
plot(time_vector, avg_psth_after_csm,  'Color', Aftercolorcsm,  'LineWidth', 2);

% Plot shaded SEM areas
fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_before_norr + sem_psth_before_norr, fliplr(avg_psth_before_norr - sem_psth_before_norr)], ...
     Beforecolornorr, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_after_norr + sem_psth_after_norr, fliplr(avg_psth_after_norr - sem_psth_after_norr)], ...
     Aftercolornorr, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_before_csm + sem_psth_before_csm, fliplr(avg_psth_before_csm - sem_psth_before_csm)], ...
     Beforecolorcsm, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_after_csm + sem_psth_after_csm, fliplr(avg_psth_after_csm - sem_psth_after_csm)], ...
     Aftercolorcsm, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Set y limits as requested
ylim([40 140]);
yl = ylim;

% Move heatmaps down below y-axis lower limit (40)
heatmap_height = 0.05 * diff(yl);
heatmap_y_start = yl(1) - 1.5 * heatmap_height;
heatmap_y_start2 = heatmap_y_start - 1.2 * heatmap_height;

bin_width = time_vector(2) - time_vector(1);

% Get colormaps - use slanCM 'OrRd' for Before and 'PuRd' for After
before_cmap = slanCM('Blues', 256);
after_cmap = slanCM('Greens', 256);

% Plot Before p-value heatmap
for t = 1:num_bins
    p_val = min(p_vals_before(t), 0.05); % clip at 0.05
    color_idx = max(1, round((1 - p_val/0.05) * 255) + 1);
    colorVal = before_cmap(color_idx, :);
    x_patch = [time_vector(t)-bin_width/2, time_vector(t)+bin_width/2, time_vector(t)+bin_width/2, time_vector(t)-bin_width/2];
    y_patch = [heatmap_y_start, heatmap_y_start, heatmap_y_start + heatmap_height, heatmap_y_start + heatmap_height];
    fill(x_patch, y_patch, colorVal, 'EdgeColor', 'none');
end

% Plot After p-value heatmap below that
for t = 1:num_bins
    p_val = min(p_vals_after(t), 0.05); % clip at 0.05
    color_idx = max(1, round((1 - p_val/0.05) * 255) + 1);
    colorVal = after_cmap(color_idx, :);
    x_patch = [time_vector(t)-bin_width/2, time_vector(t)+bin_width/2, time_vector(t)+bin_width/2, time_vector(t)-bin_width/2];
    y_patch = [heatmap_y_start2, heatmap_y_start2, heatmap_y_start2 + heatmap_height, heatmap_y_start2 + heatmap_height];
    fill(x_patch, y_patch, colorVal, 'EdgeColor', 'none');
end

% Labels for heatmaps
text(min(time_vector), heatmap_y_start + heatmap_height/2, 'p-value Before', ...
    'FontSize', 9, 'VerticalAlignment', 'middle', 'Color', Beforecolorcsm);
text(min(time_vector), heatmap_y_start2 + heatmap_height/2, 'p-value After', ...
    'FontSize', 9, 'VerticalAlignment', 'middle', 'Color', Aftercolorcsm);

% Legend
%legend('Before-NoCSM', 'After-NoCSM', 'Before-CSM', 'After-CSM', 'Location', 'northwest');
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
title('Average PSTH across Units with p-value heatmap');
xlim([-250 250]);

% Expand ylim downward so heatmaps are visible
ylim([heatmap_y_start2 - heatmap_height*0.5, yl(2)]);

% Separate colorbars for Before and After heatmaps

% Before p-value colorbar (green)
cb_ax_red = axes('Position', [0.92 0.55 0.02 0.35]);
imagesc(cb_ax_red, [0 1], [0 0.05], linspace(0.05, 0, 256)');
set(cb_ax_red, 'YDir', 'normal', 'XTick', [], 'YTick', [0 0.025 0.05]);
colormap(cb_ax_red, before_cmap);
ylabel(cb_ax_red, 'p-value Before');
caxis(cb_ax_red, [0 0.05]);
box(cb_ax_red, 'off');

% After p-value colorbar (blue)
cb_ax_blue = axes('Position', [0.92 0.1 0.02 0.35]);
imagesc(cb_ax_blue, [0 1], [0 0.05], linspace(0.05, 0, 256)');
set(cb_ax_blue, 'YDir', 'normal', 'XTick', [], 'YTick', [0 0.025 0.05]);
colormap(cb_ax_blue, after_cmap);
ylabel(cb_ax_blue, 'p-value After');
caxis(cb_ax_blue, [0 0.05]);
box(cb_ax_blue, 'off');

% Return focus to main axes
axes(gca);


%% Second method for plotting the mean PSTH and the SEM 
CSMcolor = rgb('DarkOrchid');

% Parameters
num_bootstrap = 1000;
smooth_window = 5;  % Number of bins for moving average

% Smooth each unit's PSTH before averaging
smoothed_CSM = movmean(PsthMaster_CSM, 3, 2);  % per-unit smoothing
smoothed_NoRR = movmean(PsthMaster_NoRR, 3, 2);

% Average PSTH across units
avg_psth_csm = mean(smoothed_CSM, 1);
avg_psth_norr = mean(smoothed_NoRR, 1);

% Smooth the average PSTH
smooth_avg_csm = movmean(avg_psth_csm, smooth_window);
smooth_avg_norr = movmean(avg_psth_norr, smooth_window);

% Bootstrap SEM
boot_sem = @(data) std(bootstrp(num_bootstrap, @mean, data));  % function to get SEM via bootstrap

sem_psth_csm = zeros(1, size(smoothed_CSM, 2));
sem_psth_norr = zeros(1, size(smoothed_NoRR, 2));

for t = 1:size(smoothed_CSM, 2)
    sem_psth_csm(t) = boot_sem(smoothed_CSM(:, t));
    sem_psth_norr(t) = boot_sem(smoothed_NoRR(:, t));
end

% Smooth the bootstrapped SEMs to match the smoothed average
smooth_sem_csm = movmean(sem_psth_csm, smooth_window);
smooth_sem_norr = movmean(sem_psth_norr, smooth_window);

% Time vector (adjust to match binning)
num_bins = size(PsthMaster_CSM, 2);
time_vector = linspace(-1000, 1000, num_bins);  % adjust if needed

% Statistical testing (paired signrank test)
p_vals = nan(1, num_bins);
for t = 1:num_bins
    [p_vals(t), ~] = signrank(PsthMaster_CSM(:, t), PsthMaster_NoRR(:, t));
end

% FDR correction
[h_fdr, ~, ~] = fdr_bh(p_vals);

% Plot
figure;
hold on;

% Mean lines
plot(time_vector, smooth_avg_csm, 'Color', CSMcolor, 'LineWidth', 2);
plot(time_vector, smooth_avg_norr, 'k', 'LineWidth', 2);

% Shaded error bars using smoothed SEM
fill([time_vector, fliplr(time_vector)], ...
     [smooth_avg_csm + smooth_sem_csm, fliplr(smooth_avg_csm - smooth_sem_csm)], ...
     CSMcolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([time_vector, fliplr(time_vector)], ...
     [smooth_avg_norr + smooth_sem_norr, fliplr(smooth_avg_norr - smooth_sem_norr)], ...
     'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Significance markers
sig_bins = find(h_fdr);
plot(time_vector(sig_bins), smooth_avg_csm(sig_bins), '.', 'Color', 'k', 'MarkerSize', 10);

% Labels and formatting
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
legend('CSM', 'NoCSM');
title('Smoothed Average PSTH with Bootstrapped SEM and Significance');
box off;

%% Second method for plotting the mean PSTH and the SEM - Before and After peak v 

% Parameters
num_bootstrap = 1000;
smooth_window = 5;  % Number of bins for moving average

% Smooth each unit's PSTH before averaging
smoothed_before = movmean(PsthMaster_Before_Ctrl, 3, 2);  % per-unit smoothing
smoothed_after = movmean(PsthMaster_After_Ctrl, 3, 2);

% Average PSTH across units
avg_psth_before = mean(smoothed_before, 1);
avg_psth_after = mean(smoothed_after, 1);

% Smooth the average PSTH
smooth_avg_before = movmean(avg_psth_before, smooth_window);
smooth_avg_after = movmean(avg_psth_after, smooth_window);

% Bootstrap SEM
boot_sem = @(data) std(bootstrp(num_bootstrap, @mean, data));  % function to get SEM via bootstrap

sem_psth_before = zeros(1, size(smoothed_before, 2));
sem_psth_after = zeros(1, size(smoothed_after, 2));

for t = 1:size(smoothed_before, 2)
    sem_psth_before(t) = boot_sem(smoothed_before(:, t));
    sem_psth_after(t) = boot_sem(smoothed_after(:, t));
end

% Smooth the bootstrapped SEMs to match the smoothed average
smooth_sem_before = movmean(sem_psth_before, smooth_window);
smooth_sem_after = movmean(sem_psth_after, smooth_window);

% Time vector (adjust to match binning)
num_bins = size(PsthMaster_Before_Ctrl, 2);
time_vector = linspace(-1000, 1000, num_bins);  % adjust if needed

% Plot
%figure;
figure('Position', [100, 100, 600, 300]); 
hold on;

% Mean lines
plot(time_vector, smooth_avg_before, 'Color', Beforecolor, 'LineWidth', 2);
plot(time_vector, smooth_avg_after,'Color',Aftercolor, 'LineWidth', 2);

% Shaded error bars using smoothed SEM
fill([time_vector, fliplr(time_vector)], ...
     [smooth_avg_before + smooth_sem_before, fliplr(smooth_avg_before - smooth_sem_before)], ...
     Beforecolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([time_vector, fliplr(time_vector)], ...
     [smooth_avg_after + smooth_sem_after, fliplr(smooth_avg_after - smooth_sem_after)], ...
     Aftercolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Labels and formatting
xlim([-250 250])
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
legend('Before', 'After');
title('Smoothed Average PSTH with Bootstrapped SEM and Significance');
box off;


%% 
Color_BL_CSM      = rgb('DarkMagenta');
Color_BL_NoCSM    = rgb('Black');
% Remove from both conditions to keep them matched
Psth_CSM_before_cleaned = PsthMaster_Before_CSM;
Psth_NoCSM_before_cleaned = PsthMaster_Before_NoRR;
bad_units=[18; 25]; %Badunis identified from pos thresh crossing
%bad_units=[7;18;25];
% bad_units=[3];
%bad_units=[2;7;9;10;12;13;14;16;17;18;21;22;24;25];
Psth_CSM_before_cleaned(bad_units, :) = [];
Psth_NoCSM_before_cleaned(bad_units, :) = [];

% Recompute mean and SEM
avg_before_BL_CSM = nanmean(Psth_CSM_before_cleaned, 1);
avg_before_BL_NoCSM = nanmean(Psth_NoCSM_before_cleaned, 1);

sem_before_BL_CSM = nanstd(Psth_CSM_before_cleaned, 0, 1) / sqrt(size(Psth_CSM_before_cleaned, 1));
sem_before_BL_NoCSM = nanstd(Psth_NoCSM_before_cleaned, 0, 1) / sqrt(size(Psth_NoCSM_before_cleaned, 1));

% Compute average PSTHs and SEMs — AFTER PEAK
avg_after_BL_CSM    = nanmean(Psth_CSM_before_cleaned, 1);
avg_after_BL_NoCSM  = nanmean(Psth_NoCSM_before_cleaned, 1);


sem_after_BL_CSM    = nanstd(Psth_CSM_before_cleaned, 0, 1) / sqrt(size(Psth_CSM_before_cleaned, 1));
sem_after_BL_NoCSM  = nanstd(Psth_NoCSM_before_cleaned, 0, 1) / sqrt(size(Psth_NoCSM_before_cleaned, 1));

% Time vector
num_bins = size(Psth_CSM_before_cleaned, 2);
time_vector = linspace(-1000, 1000, num_bins);

figure('Position', [100, 100, 600, 300]);hold on
% === PLOT: BEFORE PEAK UNITS ===
%figure; hold on;

% Lines
plot(time_vector, avg_before_BL_CSM,  'Color', Color_BL_CSM, 'LineWidth', 2);

plot(time_vector, avg_before_BL_NoCSM,  'Color', Color_BL_NoCSM,  'LineWidth', 2);

% Shading
fill([time_vector, fliplr(time_vector)], ...
     [avg_before_BL_CSM + sem_before_BL_CSM, fliplr(avg_before_BL_CSM - sem_before_BL_CSM)], ...
     Color_BL_CSM, 'FaceAlpha', 0.2, 'EdgeColor', 'none');


fill([time_vector, fliplr(time_vector)], ...
     [avg_before_BL_NoCSM + sem_before_BL_NoCSM, fliplr(avg_before_BL_NoCSM - sem_before_BL_NoCSM)], ...
     Color_BL_NoCSM, 'FaceAlpha', 0.2, 'EdgeColor', 'none');


%xline(csmmean,'--')
xline(0,'--')
title('PSTH – Before Peak Firing');

xlabel('Time from Pos Thresh (stim)(ms)');
ylabel('Firing rate (zscore)');
xlim([-250 250])
%ylim([0 10])
box off;




%% 
Color_BL_CSM      = rgb('DarkMagenta');
Color_BL_NoCSM    = rgb('Black');

Psth_CSM_after=PsthMaster_After_CSM;
Psth_NoCSM_after=PsthMaster_After_NoRR;

% Remove from both conditions to keep them matched
Psth_CSM_after_cleaned = Psth_CSM_after;
Psth_NoCSM_after_cleaned = Psth_NoCSM_after;

bad_units=[3,11,13,18,20,21];
Psth_CSM_after_cleaned(bad_units, :) = [];
Psth_NoCSM_after_cleaned(bad_units, :) = [];

% Recompute mean and SEM
avg_after_BL_CSM = nanmean(Psth_CSM_after_cleaned, 1);
avg_after_BL_NoCSM = nanmean(Psth_NoCSM_after_cleaned, 1);

sem_after_BL_CSM = nanstd(Psth_CSM_after_cleaned, 0, 1) / sqrt(size(Psth_CSM_after_cleaned, 1));
sem_after_BL_NoCSM = nanstd(Psth_NoCSM_after_cleaned, 0, 1) / sqrt(size(Psth_NoCSM_after_cleaned, 1));

% Compute average PSTHs and SEMs — AFTER PEAK
avg_after_BL_CSM    = nanmean(Psth_CSM_after_cleaned, 1);
avg_after_BL_NoCSM  = nanmean(Psth_NoCSM_after_cleaned, 1);

sem_after_BL_CSM    = nanstd(Psth_CSM_after_cleaned, 0, 1) / sqrt(size(Psth_CSM_after_cleaned, 1));
sem_after_BL_NoCSM  = nanstd(Psth_NoCSM_after_cleaned, 0, 1) / sqrt(size(Psth_NoCSM_after_cleaned, 1));

% Time vector
num_bins = size(Psth_CSM_after, 2);
time_vector = linspace(-1000, 1000, num_bins);


% === PLOT: AFTER PEAK UNITS ===
figure('Position', [100, 100, 600, 300]);hold on

% Lines
plot(time_vector, avg_after_BL_CSM,  'Color', Color_BL_CSM, 'LineWidth', 2);
plot(time_vector, avg_after_BL_NoCSM,  'Color', Color_BL_NoCSM,  'LineWidth', 2);

% Shading
fill([time_vector, fliplr(time_vector)], ...
     [avg_after_BL_CSM + sem_after_BL_CSM, fliplr(avg_after_BL_CSM - sem_after_BL_CSM)], ...
     Color_BL_CSM, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

fill([time_vector, fliplr(time_vector)], ...
     [avg_after_BL_NoCSM + sem_after_BL_NoCSM, fliplr(avg_after_BL_NoCSM - sem_after_BL_NoCSM)], ...
     Color_BL_NoCSM, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

xline(0,'--')
title('PSTH – After Peak Firing');
xlabel('Time from Pos Thresh (stim) (ms)');
ylabel('Firing rate (zscore)');
xlim([-250 250])
%ylim([0 18])
box off;
%%
% Define the time window of interest (0 to 100 ms)
bin_start = 0;
bin_end = 300;

% Generate time vector
num_bins = size(PsthMaster_CSM, 2);
time_vector = linspace(-1000, 1000, num_bins);

% Find indices corresponding to 0–100 ms
time_indices = find(time_vector >= bin_start & time_vector <= bin_end);

% Average firing rate in this window for each unit
mean_FR_CSM = mean(PsthMaster_CSM(:, time_indices), 2);
mean_FR_NoRR = mean(PsthMaster_NoRR(:, time_indices), 2);

% Perform paired test (use Wilcoxon if non-normal, or t-test)
[p, h, stats] = signrank(mean_FR_CSM, mean_FR_NoRR);

% Display result
fprintf('Signrank test between CSM and NoRR (0–100ms): p = %.4f\n', p);

%% Diff
diffcolor=rgb('LightSeaGreen');
% Compute average PSTHs (mean across units)
psth_diff=PsthMaster_CSM-PsthMaster_NoRR;
avg_psth_diff = mean(psth_diff, 1);  % 1 = average across rows (units)

% Optional: Compute standard error of the mean (SEM)
sem_psth_diff = std(psth_diff, 0, 1) / sqrt(size(psth_diff, 1));

% Generate time vector (assuming same binning and length)
num_bins = size(psth_diff, 2);
time_vector = linspace(-1000, 1000, num_bins);  % Update if your window or binning is different

% Plot average PSTHs
figure;
hold on;
plot(time_vector, avg_psth_diff,'Color',diffcolor, 'LineWidth', 2);

% Optional: Add shaded error bars for SEM
fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_diff + sem_psth_diff, fliplr(avg_psth_diff - sem_psth_diff)], ...
     diffcolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

%ylim([0 1]);
yline(0,'--')
xlim([-250 250])
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
title('Av. Diff in FR CSM-NoCSM');

%% Diff for before and after csm vs nocsm
beforecolor=rgb('DeepSkyBlue');
aftercolor=rgb('DarkMagenta');

% Compute average PSTHs (mean across units)
diff_before=PsthMaster_Before_CSM-PsthMaster_Before_NoRR;
diff_after=PsthMaster_After_CSM-PsthMaster_After_NoRR;

avg_psth_diff_before = mean(diff_before, 1);  % 1 = average across rows (units)
avg_psth_diff_after = mean(diff_after, 1);

% Optional: Compute standard error of the mean (SEM)
sem_psth_diff_before = std(diff_before, 0, 1) / sqrt(size(diff_before, 1));
sem_psth_diff_after = std(diff_after, 0, 1) / sqrt(size(diff_after, 1));

% Generate time vector (assuming same binning and length)
num_bins = size(diff_before, 2);
time_vector = linspace(-1000, 1000, num_bins);  % Update if your window or binning is different

% Plot average PSTHs
%figure;
figure('Position', [100, 100, 600, 300]); 
hold on;
plot(time_vector, avg_psth_diff_before,'Color',beforecolor, 'LineWidth', 2);
plot(time_vector, avg_psth_diff_after,'Color',aftercolor, 'LineWidth', 2);

% Optional: Add shaded error bars for SEM
fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_diff_before + sem_psth_diff_before, fliplr(avg_psth_diff_before - sem_psth_diff_before)], ...
     beforecolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_diff_after + sem_psth_diff_after, fliplr(avg_psth_diff_after - sem_psth_diff_after)], ...
     aftercolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

ylim([-15 20]);
yline(0,'--')
%xlim([-250 250])
xlim([-200 150])
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
title('Av. Diff in FR CSM-NoCSM');
%% Cleaned up edits: 
%% Diff for before and after csm vs nocsm
beforecolor=rgb('OrangeRed');
aftercolor=rgb('DarkBlue');


% Compute average PSTHs (mean across units)
% diff_before=PsthMaster_WOE_CSM-PsthMaster_BL_CSM;
% diff_after=PsthMaster_WOE_NoCSM-PsthMaster_BL_NoCSM;
% 
%diff_before=PsthMaster_Before_CSM-PsthMaster_Before_NoRR;
diff_before=Psth_CSM_before_cleaned-Psth_NoCSM_before_cleaned; % version without the 2 units 
%diff_after=Psth_CSM_after-Psth_NoCSM_after;
diff_after=Psth_CSM_after_cleaned-Psth_NoCSM_after_cleaned;

% diff_before=Psth_AD_CSM_before-Psth_BL_CSM_before;
% diff_after=Psth_AD_CSM_after-Psth_BL_CSM_after;

avg_psth_diff_before = nanmean(diff_before, 1);  % 1 = average across rows (units)
avg_psth_diff_after = nanmean(diff_after, 1);

% Optional: Compute standard error of the mean (SEM)
sem_psth_diff_before = nanstd(diff_before, 0, 1) / sqrt(size(diff_before, 1));
sem_psth_diff_after = nanstd(diff_after, 0, 1) / sqrt(size(diff_after, 1));

% Generate time vector (assuming same binning and length)
num_bins = size(diff_before, 2);
time_vector = linspace(-1000, 1000, num_bins);  % Update if your window or binning is different
x1=0 ;x2=50;
y_upper=4; y_lower=-6;

% Find the index closest to x = 106
target_x = 106; % from peak and valley data T2 aligned
%target_x = 25; % from peak and valley data T2 aligned
[~, idx_closest] = min(abs(time_vector - target_x));

y_after = avg_psth_diff_after(idx_closest);
%y_after = avg_psth_diff_before(idx_closest);

% Plot average PSTHs
%figure;
figure('Position', [100, 100, 600, 300]); 
hold on;
plot(time_vector, avg_psth_diff_before,'Color',beforecolor, 'LineWidth', 2);
plot(time_vector, avg_psth_diff_after,'Color',aftercolor, 'LineWidth', 2);

% Optional: Add shaded error bars for SEM
fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_diff_before + sem_psth_diff_before, fliplr(avg_psth_diff_before - sem_psth_diff_before)], ...
     beforecolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

 fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_diff_after + sem_psth_diff_after, fliplr(avg_psth_diff_after - sem_psth_diff_after)], ...
     aftercolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
 hold on
plot(target_x, y_after, 'o', 'Color', aftercolor, 'MarkerFaceColor', aftercolor);

%hold on 
% Fill the area between the vertical lines
%fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
%ylim([0 1]);
yline(0,'--')
yline(-2,':')
yline(2,':')
%xline(csmmean,'--')
xline(0,'--')
%xlim([-250 250])
%xlim([-100 250])
xlim([-200 150])
%ylim([-3 5])
xline(106)
%xline(25)
%ylim([-25 25])
text(target_x + 5, y_after, sprintf('%.2f', y_after), 'Color', aftercolor, 'FontSize', 10);
xlabel('Time from CSM(ms)');
ylabel('Firing rate (Hz)');
title('Av. Diff in FR WOE-BL CSM for before and after');
%title('Av. Diff in FR AD-BL CSM for before and after');
%% calculating inflection points 
% SETTINGS
plot_window = [-250, 250];  % display window
time_window = [-200, 200];  % peak detection window

% Colors
beforecolor = rgb('SteelBlue');
aftercolor = rgb('DarkGreen');

% Compute averages and SEMs
avg_psth_diff_before = mean(diff_before, 1);
avg_psth_diff_after  = mean(diff_after, 1);

sem_psth_diff_before = std(diff_before, 0, 1) / sqrt(size(diff_before, 1));
sem_psth_diff_after  = std(diff_after, 0, 1) / sqrt(size(diff_after, 1));

% Generate time vector
num_bins = size(diff_before, 2);
time_vector = linspace(-1000, 1000, num_bins);

% Logical and index masks
plot_idx = time_vector >= plot_window(1) & time_vector <= plot_window(2);
time_win_idx = find(time_vector >= time_window(1) & time_vector <= time_window(2));
time_win_vec = time_vector(time_win_idx);

% Initialize inflections struct
inflections = struct();

% ========== BEFORE PLOT ==========
figure('Position', [100, 100, 600, 300]); hold on;

% Extract curves
curve_before = avg_psth_diff_before(plot_idx);
sem_before = sem_psth_diff_before(plot_idx);
times_plot = time_vector(plot_idx);

% Shaded error
fill([times_plot, fliplr(times_plot)], ...
     [curve_before + sem_before, fliplr(curve_before - sem_before)], ...
     beforecolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Plot line
plot(times_plot, curve_before, 'Color', beforecolor, 'LineWidth', 2);

% --- Find inflections ---
curve_win = avg_psth_diff_before(time_win_idx);
sem_win = sem_psth_diff_before(time_win_idx);
[pks, locs_pks] = findpeaks(curve_win);
[vals_neg, locs_vals] = findpeaks(-curve_win); vals = -vals_neg;

inflections.before.peaks = [];
inflections.before.valleys = [];

for i = 1:length(pks)
    t = time_win_vec(locs_pks(i)); y = pks(i); s = sem_win(locs_pks(i));
    plot(t, y, 'k^', 'MarkerFaceColor', beforecolor, 'MarkerSize', 6);
    text(t+5, y, sprintf('%.2f ± %.2f', y, s), 'Color', beforecolor, 'FontSize', 9);
    inflections.before.peaks(end+1,:) = [t, y, s];
end
for i = 1:length(vals)
    t = time_win_vec(locs_vals(i)); y = vals(i); s = sem_win(locs_vals(i));
    plot(t, y, 'kv', 'MarkerFaceColor', beforecolor, 'MarkerSize', 6);
    text(t+5, y, sprintf('%.2f ± %.2f', y, s), 'Color', beforecolor, 'FontSize', 9);
    inflections.before.valleys(end+1,:) = [t, y, s];
end
% 
% % Format
% yline(0, '--', 'Color', [0.4 0.4 0.4]);
% xline(0, '--');
% xlabel('Time from Stim (ms)');
% ylabel('Δ Firing Rate (Hz)');
% title('Avg ΔFR: CSM vs NoRR (Before)');
% xlim(plot_window);
% ylim([min(curve_before - sem_before) - 5, max(curve_before + sem_before) + 5]);
% legend({'Before SEM', 'Before Mean'}, 'Location', 'best');
% hold off;

% ========== AFTER PLOT ==========
%figure('Position', [100, 100, 600, 300]); hold on;

curve_after = avg_psth_diff_after(plot_idx);
sem_after = sem_psth_diff_after(plot_idx);

% Shaded error
fill([times_plot, fliplr(times_plot)], ...
     [curve_after + sem_after, fliplr(curve_after - sem_after)], ...
     aftercolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Plot line
plot(times_plot, curve_after, 'Color', aftercolor, 'LineWidth', 2);

% --- Find inflections ---
curve_win = avg_psth_diff_after(time_win_idx);
sem_win = sem_psth_diff_after(time_win_idx);
[pks, locs_pks] = findpeaks(curve_win);
[vals_neg, locs_vals] = findpeaks(-curve_win); vals = -vals_neg;

inflections.after.peaks = [];
inflections.after.valleys = [];

for i = 1:length(pks)
    t = time_win_vec(locs_pks(i)); y = pks(i); s = sem_win(locs_pks(i));
    plot(t, y, 'k^', 'MarkerFaceColor', aftercolor, 'MarkerSize', 6);
    text(t+5, y, sprintf('%.2f ± %.2f', y, s), 'Color', aftercolor, 'FontSize', 9);
    inflections.after.peaks(end+1,:) = [t, y, s];
end
for i = 1:length(vals)
    t = time_win_vec(locs_vals(i)); y = vals(i); s = sem_win(locs_vals(i));
    plot(t, y, 'kv', 'MarkerFaceColor', aftercolor, 'MarkerSize', 6);
    text(t+5, y, sprintf('%.2f ± %.2f', y, s), 'Color', aftercolor, 'FontSize', 9);
    inflections.after.valleys(end+1,:) = [t, y, s];
end

% Format
%ylim([-15 20])
yline(0, '--', 'Color', [0.4 0.4 0.4]);
xline(0, '--');
xlabel('Time from Stim (ms)');
ylabel('Δ Firing Rate (Hz)');
title('Avg ΔFR: CSM vs NoRR (After)');
xlim(plot_window);
%ylim([min(curve_after - sem_after) - 5, max(curve_after + sem_after) + 5]);
legend({'After SEM', 'After Mean'}, 'Location', 'best');
hold off;

%% Different way to calculate change in firing rate 
% Parameters
num_bootstrap = 1000;
smooth_window = 5;

% Difference PSTHs (after smoothing each unit)
unit_diff = movmean(PsthMaster_CSM, 3, 2) - movmean(PsthMaster_NoRR, 3, 2);  % smooth then subtract

% Average difference
avg_diff = mean(unit_diff, 1);

% Smooth the average difference
smooth_avg_diff = movmean(avg_diff, smooth_window);

% Bootstrapped SEM of the difference
boot_sem = @(data) std(bootstrp(num_bootstrap, @mean, data));
sem_diff = zeros(1, size(unit_diff, 2));
for t = 1:size(unit_diff, 2)
    sem_diff(t) = boot_sem(unit_diff(:, t));
end

% Smooth the SEM
smooth_sem_diff = movmean(sem_diff, smooth_window);

% Time vector
num_bins = size(PsthMaster_CSM, 2);
time_vector = linspace(-1000, 1000, num_bins);  % adjust if needed

% Statistical testing on the difference
p_vals_diff = nan(1, num_bins);
for t = 1:num_bins
    [p_vals_diff(t), ~] = signrank(PsthMaster_CSM(:, t), PsthMaster_NoRR(:, t));
end
[h_fdr_diff, ~, ~] = fdr_bh(p_vals_diff);

% Plot
figure;
hold on;

% Mean difference line
plot(time_vector, smooth_avg_diff, 'k', 'LineWidth', 2);

% Shaded bootstrapped SEM around the difference
fill([time_vector, fliplr(time_vector)], ...
     [smooth_avg_diff + smooth_sem_diff, fliplr(smooth_avg_diff - smooth_sem_diff)], ...
     'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Plot significance markers
sig_bins_diff = find(h_fdr_diff);
plot(time_vector(sig_bins_diff), smooth_avg_diff(sig_bins_diff), '.', 'Color', 'r', 'MarkerSize', 10);

% Labels and formatting
xlabel('Time (ms)');
ylabel('Δ Firing rate (Hz)');
title('Smoothed PSTH Difference with Bootstrapped SEM and Significance');
line([min(time_vector), max(time_vector)], [0, 0], 'Color', [0.5 0.5 0.5], 'LineStyle', '--');  % baseline
box off;
%xlim([-250 250])

%% Now look at kinematics
% Collect all NoCSM reaches and CSM
NoCSMAll= struct; 
CSMAll=struct;
AllReaches=struct;
% Combine all reachse across animals 
for i = 1:length(NPdata)
 
    Reach=NPdata(i).NoCSM;
    NoCSMAll(i).NoCSM=Reach; 

    ReachCSM=NPdata(i).ReReaches;
    CSMAll(i).CSM=ReachCSM; 

    ReachAll=NPdata(i).Reaches;
    AllReaches(i).All=ReachAll; 
   
    clearvars Reach ReachCSM ReachAll
end

%% collect reaches and get average 
NoCSM_stimvel=struct;
for i=1:length(NoCSMAll)
   
    NoCSMR=NoCSMAll(i).NoCSM;

    if isempty(NoCSMR)
        C=[];
    else
        dim=7; %max outward velocity
        dir=1;
        thresh= 0.1;
        C= AlignPeakV(NoCSMR,dim,1000);
        %C= AlignPeakV(NoCSMR,dim,250);
        %C= AlignThreshCSM_T2HF(NoCSMR,thresh,dir,250);
        %C= AlignThreshCSM_T0HF(NoCSMR,thresh,dir,250);
        %C= AlignThresh(NoCSMR,thresh,dir,500);
        C=stimvel_accel(C);
        C=stimvel_accel_jerk(C);
    end
    NoCSM_stimvel(i).stimvel=C;
end
%%
CSM_stimvel=struct;
for i=1:length(CSMAll)
   
   CSMR=CSMAll(i).CSM;
   if isempty(CSMR)
       RR=[];
   else
      RR= AlignPeakV(CSMR,dim,1000);
      % RR= AlignPeakV(CSMR,dim,250);
      %RR= AlignThreshCSM_T2HF(CSMR,thresh,dir,250);
       %RR= AlignThreshCSM_T0HF(CSMR,thresh,dir,250);
      % RR= AlignThresh(CSMR,thresh,dir,500);
       RR=stimvel_accel(RR);
       RR=stimvel_accel_jerk(RR);
   end
    CSM_stimvel(i).stimvel=RR;
end 

%%
dim=7;
All_stimvel=struct;
for i=1:length(AllReaches)
   
   AllR=AllReaches(i).All;
   if isempty(AllR)
       RR=[];
   else
      RR= AlignPeakV(AllR,dim,1000);
      % RR= AlignPeakV(CSMR,dim,250);
      %RR= AlignThreshCSM_T2HF(CSMR,thresh,dir,250);
       %RR= AlignThreshCSM_T0HF(CSMR,thresh,dir,250);
      % RR= AlignThresh(CSMR,thresh,dir,500);
       RR=stimvel_accel(RR);
       RR=stimvel_accel_jerk(RR);
   end
    All_stimvel(i).stimvel=RR;
end 
%% REmove entries with no cells active 
CSM_stimvel(9) = [];
NoCSM_stimvel(9) = [];
All_stimvel(9) = [];
%%
NRR_big=struct;

for a = 1:length(NoCSM_stimvel) 
    NRR_big.stimvel(:,:,a) = NoCSM_stimvel(a).stimvel.Average;
end
%NoReReach
NRR_big.Average = nanmean(NRR_big.stimvel,3);
NRR_big.Std = nanstd(NRR_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(NRR_big.Average),1); 
for nnum = 1:length(NRR_big.Average)
    samp = squeeze(NRR_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
NRR_big.N = N;
%%
RR_big=struct;

for a = 1:length(CSM_stimvel) 
    RR_big.stimvel(:,:,a) = CSM_stimvel(a).stimvel.Average;
end
%NoReReach
RR_big.Average = nanmean(RR_big.stimvel,3);
RR_big.Std = nanstd(RR_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(RR_big.Average),1); 
for nnum = 1:length(RR_big.Average)
    samp = squeeze(RR_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
RR_big.N = N;
%%
All_big=struct;

for a = 1:length(All_stimvel) 
    All_big.stimvel(:,:,a) = All_stimvel(a).stimvel.Average;
end
%NoReReach
All_big.Average = nanmean(All_big.stimvel,3);
All_big.Std = nanstd(All_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(All_big.Average),1); 
for nnum = 1:length(All_big.Average)
    samp = squeeze(All_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
All_big.N = N;
%% Stimvel Plot 
CSMColor=rgb('DarkMagenta');
dim=5;
 StimvelPlot(All_big,dim,'k',1);
% StimvelPlot(NRR_big,dim,'k',1);
% StimvelPlot(RR_big,dim,CSMColor,0);
yline(0,'--')
%ylim('auto')
ylim([-1.5 2.5])
%xlim([-0.25 0.25])
%xlim('auto')
xlim([-1 1])
ylabel('Outward Velocity (cm/s)')
xlabel('Time from Peak Velocity (s)')
f = gcf; % Get current figure handle
f.Units = 'inches';
pos = f.Position;
f.Position = [pos(1:2), pos(3), pos(4)/2]; % Halve the height
%% With special ticks 
%% Stimvel Plot 
CSMColor = rgb('DarkMagenta');
dim = 5;

StimvelPlot(All_big, dim, 'k', 1);
% StimvelPlot(NRR_big, dim, 'k', 1);
% StimvelPlot(RR_big, dim, CSMColor, 0);

yline(0,'--')

ylim([-1.5 2.5])
xlim([-1 1])

% ---- ADD THIS LINE ----
xticks(-1:0.25:1)

ylabel('Outward Velocity (cm/s)')
xlabel('Time from Peak Velocity (s)')

f = gcf; % Get current figure handle
f.Units = 'inches';
pos = f.Position;
f.Position = [pos(1:2), pos(3), pos(4)/2]; % Halve the height



%%
figure;
subplot(3,1,1); hold on;
StimvelPlot(NRR_big,8,[0 0 0],0);
StimvelPlot(RR_big,8,[0 0 1],0);

ylabel('Out. Velocity (cm/s)')
yline(0,'--k')
xlabel('Time from threshold crossing (ms)')
xlim([-1 1])
ylim('auto')

subplot(3,1,2); hold on
imagesc(sorted_psth_norr);
colormap(slanCM('GnBu'));
colorbar;
% Adjust x-axis to reflect -1000 to 1000 ms
num_bins = size(sorted_psth_norr, 2);
bin_size = 10;  % ms per bin
time_vector = linspace(-1000, 1000, num_bins);  % create time axis

% Set xticks at reasonable intervals (e.g., every 250 ms)
tick_times = -1000:250:1000;
[~, tick_indices] = min(abs(time_vector' - tick_times), [], 1);  % find closest indices

set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times);

xlabel('Time (ms)');
ylabel('Units (sorted by peak firing time)');
title('Control PSTHs Sorted by Peak Firing Time');

subplot(3,1,3); hold on
imagesc(sorted_psth_csm);
colormap(slanCM('GnBu'));
colorbar;
% Adjust x-axis to reflect -1000 to 1000 ms
num_bins = size(sorted_psth_csm, 2);
bin_size = 10;  % ms per bin
time_vector = linspace(-1000, 1000, num_bins);  % create time axis

[~, tick_indices] = min(abs(time_vector' - tick_times), [], 1);  % find closest indices

set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times);

xlabel('Time (ms)');
ylabel('Units (sorted by peak firing time)');
title('Control PSTHs Sorted by Peak Firing Time');
%%
figure('Units', 'normalized', 'Position', [0.3 0.1 0.25 0.8]);  % Narrow, tall figure

% Layout parameters
ax_width = 0.65;             % narrow width
heatmap_height = 0.35;       % taller height for square-ish heatmaps
gap = 0.05;                  % vertical gap between plots
cb_width = 0.02;

% Velocity plot (top) — 1/3 the heatmap height
top_height = heatmap_height / 3;
top_bottom = 1 - top_height - 0.05;
ax1 = axes('Position', [0.15, top_bottom, ax_width, top_height]); hold on;
StimvelPlot(NRR_big, 8, [0 0 0], 0);
StimvelPlot(RR_big, 8, [0 0 1], 0);
ylabel('Out. Acceleration (cm/s)');
yline(0, '--k');
xlabel('Time from threshold crossing (ms)');
xlim([-0.25 0.25]);
%xlim([-1 1]);
% ylim('auto');
ylim([-40 40])
title('Acce;');

% Middle heatmap (Control)
mid_bottom = top_bottom - gap - heatmap_height;
ax2 = axes('Position', [0.15, mid_bottom, ax_width, heatmap_height]); hold on;
imagesc(sorted_psth_norr);
colormap(ax2, slanCM('GnBu'));
axis tight;
set(ax2, 'YDir', 'reverse');  % Flip Y-axis so first unit is at top

% Limit to -250 to 250 ms
num_bins = size(sorted_psth_norr, 2);
time_vector = linspace(-1000, 1000, num_bins);
plot_range = [-250 250];
%plot_range = [-1000 1000];
[~, x_start] = min(abs(time_vector - plot_range(1)));
[~, x_end] = min(abs(time_vector - plot_range(2)));

set(ax2, 'XLim', [x_start x_end], ...
         'XTick', linspace(x_start, x_end, 5), ...
         'XTickLabel', linspace(plot_range(1), plot_range(2), 5));
xlabel('Time (ms)');
ylabel('Units');
title('Control PSTHs');
colorbar(ax2, 'Position', [0.15 + ax_width + 0.01, mid_bottom, cb_width, heatmap_height]);

% Bottom heatmap (CSM)
bot_bottom = mid_bottom - gap - heatmap_height;
ax3 = axes('Position', [0.15, bot_bottom, ax_width, heatmap_height]); hold on;
imagesc(sorted_psth_csm);
colormap(ax3, slanCM('GnBu'));
axis tight;
set(ax3, 'YDir', 'reverse'); 

num_bins = size(sorted_psth_csm, 2);
time_vector = linspace(-1000, 1000, num_bins);
[~, x_start] = min(abs(time_vector - plot_range(1)));
[~, x_end] = min(abs(time_vector - plot_range(2)));

set(ax3, 'XLim', [x_start x_end], ...
         'XTick', linspace(x_start, x_end, 5), ...
         'XTickLabel', linspace(plot_range(1), plot_range(2), 5));
xlabel('Time (ms)');
ylabel('Units');
title('CSM PSTHs');
colorbar(ax3, 'Position', [0.15 + ax_width + 0.01, bot_bottom, cb_width, heatmap_height]);
%% 
diff_avg = RR_big.Average(:,8) - NRR_big.Average(:,8);
diff_avg=diff_avg*10;
time=RR_big.Average(:,1);

% Standard error of the difference (element-wise)
SE_diff = sqrt((RR_big.Std.^2 ./ RR_big.N) + (NRR_big.Std.^2 ./ NRR_big.N));
%%
figure; hold on
plot(time,diff_avg)
xlim([-0.25 0.25])
yline(0,'--')
xlabel('Time from postional threshold (ms)')
ylabel('Change in Outward Acceleration')

%% Plot difference between CSM and NoCSM with shaded error bars

% 1. Compute average difference
diff_avg = RR_big.Average - NRR_big.Average;

% 2. Compute standard error of the difference
SE_diff = sqrt((RR_big.Std.^2 ./ RR_big.N) + (NRR_big.Std.^2 ./ NRR_big.N));

% 3. Define time axis (adjust based on your binning)
time = (1:length(diff_avg))';  % or multiply by bin size if you have it in ms

% 4. Plot with shadedErrorBar
figure; hold on;
shadedErrorBar(time, diff_avg, SE_diff, 'lineProps', {'k', 'LineWidth', 2});
xlabel('Time (ms)');
ylabel('Velocity Difference (CSM - NoCSM)');
title('Kinematic Difference Between CSM and NoCSM');
yline(0, '--', 'Color', [0.5 0.5 0.5]); % optional: zero line
%%
% For example, just plot the y-dimension (2nd column)
diff_avg_y = RR_big.Average(:,8) - NRR_big.Average(:,8);
diff_avg_y=diff_avg_y*10;
SE_diff_y = sqrt((RR_big.Std(:,8).^2 ./ RR_big.N) + (NRR_big.Std(:,8).^2 ./ NRR_big.N));

% diff_avg_y = RR_big.Average(:,5) - NRR_big.Average(:,5);
% diff_avg_y=diff_avg_y*10;
% SE_diff_y = sqrt((RR_big.Std(:,5).^2 ./ RR_big.N) + (NRR_big.Std(:,5).^2 ./ NRR_big.N));

time = (RR_big.Average(:,1))'; % ensure it's a column vector

% Now plot
figure('Position', [100, 100, 600, 300]);  
shadedErrorBar(time, diff_avg_y, SE_diff_y, 'lineProps', {'k', 'LineWidth', 2});
xlabel('Time (ms)');
%xlim([-0.2 0.15])
xlim('auto')
ylabel('Y-Velocity Difference (CSM - NoCSM)');
title('CSM vs NoCSM Kinematic Difference (Y-Direction)');
yline(0, '--');
%% Generate stimvel zscore
% Initialize arrays for p-values and Z-scores
parray0 = NaN(size(RR_big.stimvel, 1), 1);  % Store p-values
zarray0 =   NaN(size(RR_big.stimvel, 1), 1)   ;    % Store Z-scores (same shape as stimvel)
Zscores0= NaN(size(RR_big.stimvel, 1), 5);

for pnum = 1:length(parray0)
    % Extract trial data for C_big (control) and T0_big (experimental)
    a4b = squeeze(RR_big.stimvel(pnum, 5, :));  % Experimental data (T0_big)
    b4b = squeeze(NRR_big.stimvel(pnum,5, :));   % Control data (C_big), matching trial

    % Calculate p-value using Wilcoxon signed-rank test (trial-wise)
    [parray0(pnum, 1), ~, pstat] = signrank(a4b, b4b, 'method', 'approximate');
    parray0(pnum, 2) = pstat.zval;  % Store Z-value from signrank

    % Calculate mean and standard deviation for the control trial (b4b)
    mu_b4b = mean(b4b);             % Mean of the current trial in C_big
    sigma_b4b = std(b4b);           % Standard deviation of the current trial in C_big

    % Calculate Z-scores for each time point in a4b relative to b4b stats
    zscores = (a4b - mu_b4b) ./ sigma_b4b;  % Point-wise Z-score
    Zscores(pnum,:)=zscores;
    mean_zscore=mean(zscores);
    zarray0(pnum,:)=mean_zscore;
end

%%

CSMcolor=rgb('Orchid');
%zscore plot 
% Adjust the figure to provide more vertical space
set(gcf, 'Position', [100, 100, 800, 600]); % [left, bottom, width, height]

dim = 5;
StimvelPlot(RR_big, dim,CSMcolor, 1);
StimvelPlot(NRR_big, dim,'k', 0);

% Plot reference line
yline(0, '--k');

% Axis limits and labels
ylim('auto')
%ylim([-20 20])
xlabel('Time from MaxV (s)');

ylabel('Outward Acceleration');


% Adjust the main plot's position to create space below
ax1 = gca;
ax1.Position(2) = ax1.Position(2) + 0.1;
ax1.Position(4) = ax1.Position(4) - 0.1;

% Add Z-score color bar below
ax2 = axes('Position', [ax1.Position(1), ax1.Position(2) - 0.2, ax1.Position(3), 0.05]);

% Plot Z-scores as a heatmap
imagesc(zarray0', 'Parent', ax2);  % Transpose to match time x trial layout
colormap(ax2, brewermap([], 'PrGn'));  % Diverging colormap
caxis([-2, 2]);  % Adjust the Z-score color limits
colorbar(ax2, 'Position', [0.92, ax2.Position(2), 0.02, ax2.Position(4)]);

% Hide ticks and labels for the Z-score axis
set(ax2, 'XTick', [], 'YTick', [], 'YColor', 'none', 'XColor', 'none');


%%
CSMcolor = rgb('Orchid');

% zscore plot 
set(gcf, 'Position', [100, 100, 800, 600]); % [left, bottom, width, height]

dim = 5;

% Plot both datasets on the same axes
StimvelPlot(RR_big, dim, CSMcolor, 1);
hold on;
StimvelPlot(NRR_big, dim, 'k', 0);
hold off;

% Plot reference line
yline(0, '--k');

% Axis limits and labels
ylim('auto')
xlabel('Time from MaxV (s)');
ylabel('Outward Acceleration');

% Adjust the main plot's position to make it shorter
ax1 = gca;
ax1.Position(2) = ax1.Position(2) + 0.15;  % Raise the plot up
ax1.Position(4) = ax1.Position(4) * 0.6;   % Shrink the height (e.g., 60% of original)

% Add Z-score color bar below
ax2 = axes('Position', [ax1.Position(1), ax1.Position(2) - 0.18, ax1.Position(3), 0.05]);

% Plot Z-scores as a heatmap
imagesc(zarray0', 'Parent', ax2);  % Transpose to match time x trial layout
colormap(ax2, brewermap([], 'PrGn'));  % Diverging colormap
caxis([-2, 2]);  % Adjust the Z-score color limits
colorbar(ax2, 'Position', [0.92, ax2.Position(2), 0.02, ax2.Position(4)]);

% Hide ticks and labels for the Z-score axis
set(ax2, 'XTick', [], 'YTick', [], 'YColor', 'none', 'XColor', 'none');
