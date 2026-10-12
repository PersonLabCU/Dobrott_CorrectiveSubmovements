%% Get data - COde for CSM FR analysis - bursters and decreases 


% %using opto index and reaches are separated and csm alignment is present
cd("C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\StimAnimals\BlockSeparated_CSMvsNoCSM_inducedalignment")
MyFolderInfo = dir('*Block_CSMvsNoCSM.mat'); 

tic
clearvars -except MyFolderInfo
%%
% === INITIALIZE MASTER MATRICES ===================================== %%
NPdata = struct;

% -- BL
PsthMaster_BL_CSM     = [];
PsthMaster_BL_NoCSM   = [];
% -- AD
PsthMaster_AD_CSM     = [];
PsthMaster_AD_NoCSM   = [];
% -- WOE
PsthMaster_WOE_CSM    = [];
PsthMaster_WOE_NoCSM  = [];

% --- New: CSM/NoCSM split by BEFORE/AFTER peak time ----------------- %%
% BL
Psth_BL_CSM_before = []; Psth_BL_CSM_after = [];
Psth_BL_NoCSM_before = []; Psth_BL_NoCSM_after = [];
% AD
Psth_AD_CSM_before = []; Psth_AD_CSM_after = [];
Psth_AD_NoCSM_before = []; Psth_AD_NoCSM_after = [];
% WOE
Psth_WOE_CSM_before = []; Psth_WOE_CSM_after = [];
Psth_WOE_NoCSM_before = []; Psth_WOE_NoCSM_after = [];

% === MAIN LOOP ======================================================= %%
for day = 1:length(MyFolderInfo)
    load(MyFolderInfo(day).name, ...
         'Reaches', 'data', 'filtered_modulated_indices', ...
         'Salt_index','neg_corr_cells','pos_corr_cells','FRBins', ...
         'FRCorrStats_edit','optotagged_modulated_indices', ...
         'ExcitRR_BL','ExcitRR_AD','ExcitRR_WO','ExcitRR_WOE');

    if isempty(Salt_index)
        fprintf('Day %d: Salt_index is empty. Skipping...\n', day);
        continue;
    end

    data_opto_all = data.units(Salt_index);
    timevec = data.units(Salt_index(1)).BLCSM_max_vx.time_bins;

    [~, zero_idx] = min(abs(timevec));  % index closest to 0
    valid_peak_range = find(timevec >= -250 & timevec <= 250);

    keep_unit          = false(numel(data_opto_all),1);
    peak_times_ctrl_all = NaN(numel(data_opto_all),1);

    for t = 1:numel(data_opto_all)
        psth               = data_opto_all(t).BLCSM_max_vx.PSTH_smoothbin;
        [~, peak_idx]      = max(psth);
        peak_times_ctrl_all(t) = peak_idx;
        keep_unit(t)       = ismember(peak_idx, valid_peak_range);
    end

    data_opto         = data_opto_all(keep_unit);
    peak_times_ctrl   = peak_times_ctrl_all(keep_unit);

    is_before = peak_times_ctrl < zero_idx;
    is_after  = peak_times_ctrl >= zero_idx;

    if isempty(data_opto)
        fprintf('Day %d: No modulated units with peak between -250 and 250 ms. Skipping...\n', day);
        continue
    end

    % =================== BL CONDITION =============================== %%
    psth_bl_len  = numel(data_opto(1).BLCSM_csmT2.PSTH_smoothbin);
    bl_csm_mat   = NaN(numel(data_opto), psth_bl_len);
    bl_nocsm_mat = NaN(numel(data_opto), psth_bl_len);

    for t = 1:numel(data_opto)
        % bl_csm_mat(t,:)   = data_opto(t).BLCSM_csmT2.PSTH_smoothbin;
        % bl_nocsm_mat(t,:) = data_opto(t).BLNoCSM_csmT2.PSTH_smoothbin;
        % 
        bl_csm_mat(t,:)   = data_opto(t).BLCSM_pos_thresh.PSTH_smoothbin;
        bl_nocsm_mat(t,:) = data_opto(t).BLNoCSM_pos_thresh.PSTH_smoothbin;
        % 
        % bl_csm   = data_opto(t).BLCSM_pos_thresh.PSTH_smoothbin;
        % bl_nocsm = data_opto(t).BLNoCSM_pos_thresh.PSTH_smoothbin;
        % 
        % bl_csm_mat(t,:)   = bl_csm   ./ max(bl_csm);
        % bl_nocsm_mat(t,:)   = bl_nocsm   ./ max(bl_nocsm);
    end

    PsthMaster_BL_CSM   = vertcat(PsthMaster_BL_CSM,   bl_csm_mat);
    PsthMaster_BL_NoCSM = vertcat(PsthMaster_BL_NoCSM, bl_nocsm_mat);

    % Append by before/after
    Psth_BL_CSM_before     = vertcat(Psth_BL_CSM_before,     bl_csm_mat(is_before,:));
    Psth_BL_CSM_after      = vertcat(Psth_BL_CSM_after,      bl_csm_mat(is_after,:));
    Psth_BL_NoCSM_before   = vertcat(Psth_BL_NoCSM_before,   bl_nocsm_mat(is_before,:));
    Psth_BL_NoCSM_after    = vertcat(Psth_BL_NoCSM_after,    bl_nocsm_mat(is_after,:));

    % =================== AD CONDITION =============================== %%
    psth_ad_len  = numel(data_opto(1).ADCSM_csmT2.PSTH_smoothbin);
    ad_csm_mat   = NaN(numel(data_opto), psth_ad_len);
    ad_nocsm_mat = NaN(numel(data_opto), psth_ad_len);

    for t = 1:numel(data_opto)
        % ad_csm_mat(t,:)   = data_opto(t).ADCSM_csmT2.PSTH_smoothbin;
        % ad_nocsm_mat(t,:) = data_opto(t).ADNoCSM_csmT2.PSTH_smoothbin;

        ad_csm_mat(t,:)   = data_opto(t).ADCSM_pos_thresh.PSTH_smoothbin;
        ad_nocsm_mat(t,:) = data_opto(t).ADNoCSM_pos_thresh.PSTH_smoothbin;

        % ad_csm  = data_opto(t).ADCSM_pos_thresh.PSTH_smoothbin;
        % ad_nocsm= data_opto(t).ADNoCSM_pos_thresh.PSTH_smoothbin;
        % 
        % ad_csm_mat(t,:)   = ad_csm   ./ max(ad_csm);
        % ad_nocsm_mat(t,:)   = ad_nocsm   ./ max(ad_nocsm);
    end

    PsthMaster_AD_CSM   = vertcat(PsthMaster_AD_CSM,   ad_csm_mat);
    PsthMaster_AD_NoCSM = vertcat(PsthMaster_AD_NoCSM, ad_nocsm_mat);

    Psth_AD_CSM_before     = vertcat(Psth_AD_CSM_before,     ad_csm_mat(is_before,:));
    Psth_AD_CSM_after      = vertcat(Psth_AD_CSM_after,      ad_csm_mat(is_after,:));
    Psth_AD_NoCSM_before   = vertcat(Psth_AD_NoCSM_before,   ad_nocsm_mat(is_before,:));
    Psth_AD_NoCSM_after    = vertcat(Psth_AD_NoCSM_after,    ad_nocsm_mat(is_after,:));

    % =================== WOE CONDITION ============================== %%
    psth_woe_len  = numel(data_opto(1).WOECSM_csmT2.PSTH_smoothbin);
    woe_csm_mat   = NaN(numel(data_opto), psth_woe_len);
    woe_nocsm_mat = NaN(numel(data_opto), psth_woe_len);

    for t = 1:numel(data_opto)
        % woe_csm_mat(t,:)   = data_opto(t).WOECSM_csmT2.PSTH_smoothbin;
        % woe_nocsm_mat(t,:) = data_opto(t).WOENoCSM_csmT2.PSTH_smoothbin;

        woe_csm_mat(t,:)   = data_opto(t).WOECSM_pos_thresh.PSTH_smoothbin;
        woe_nocsm_mat(t,:) = data_opto(t).WOENoCSM_pos_thresh.PSTH_smoothbin;

        % woe_csm  = data_opto(t).WOECSM_pos_thresh.PSTH_smoothbin;
        % woe_nocsm= data_opto(t).WOENoCSM_pos_thresh.PSTH_smoothbin;
        % 
        % woe_csm_mat(t,:)   = woe_csm   ./ max(woe_csm);
        % woe_nocsm_mat(t,:)   = woe_nocsm   ./ max(woe_nocsm);
    end

    PsthMaster_WOE_CSM   = vertcat(PsthMaster_WOE_CSM,   woe_csm_mat);
    PsthMaster_WOE_NoCSM = vertcat(PsthMaster_WOE_NoCSM, woe_nocsm_mat);

    Psth_WOE_CSM_before     = vertcat(Psth_WOE_CSM_before,     woe_csm_mat(is_before,:));
    Psth_WOE_CSM_after      = vertcat(Psth_WOE_CSM_after,      woe_csm_mat(is_after,:));
    Psth_WOE_NoCSM_before   = vertcat(Psth_WOE_NoCSM_before,   woe_nocsm_mat(is_before,:));
    Psth_WOE_NoCSM_after    = vertcat(Psth_WOE_NoCSM_after,    woe_nocsm_mat(is_after,:));

    % ============= METADATA STORAGE ================================ %%
    NPdata(day).peak_times_ctrl = peak_times_ctrl;
    NPdata(day).kept_indices    = Salt_index(keep_unit);
    NPdata(day).ExcitRR_BL      = ExcitRR_BL;
    NPdata(day).ExcitRR_AD      = ExcitRR_AD;
    NPdata(day).ExcitRR_WOE     = ExcitRR_WOE;
    NPdata(day).cells           = data;
    NPdata(day).is_before       = is_before;
    NPdata(day).is_after        = is_after;

    fprintf('Day %d: Kept %d/%d modulated units with peak between -250 and 250 ms\n', ...
        day, sum(keep_unit), numel(keep_unit));
end



%% Plotting individual cells - After Peak V 
% === PLOT INDIVIDUAL PSTHs ===
folder_path="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\StimAnimals\Figures\posthreshaligned_individualunits\BeforeCells";
% Assuming timevec is the same as used before

num_bins = size(Psth_BL_CSM_before, 2);
timevec = linspace(-1000, 1000, num_bins);

nCells = size(Psth_BL_CSM_before, 1);
blcolor=rgb('Black');
stimcolor=rgb('DodgerBlue');
washcolor=rgb('Magenta');

smooth_window = 3;

for i = 1:nCells
    figure;
    hold on;

    psth_csm = movmean(Psth_BL_CSM_before(i,:), smooth_window);
    psth_norr = movmean(Psth_AD_CSM_before(i,:), smooth_window);
    psth_wo = movmean(Psth_WOE_CSM_before(i,:), smooth_window);
    
    % Plot CSM PSTH
    plot(timevec, psth_csm, 'Color', blcolor, 'LineWidth', 1.5);
    
    % Plot NoRR PSTH
    plot(timevec, psth_norr, 'Color', stimcolor, 'LineWidth', 1.5);

    plot(timevec, psth_wo, 'Color', washcolor, 'LineWidth', 1.5);
    
    
    % % Plot CSM PSTH
    % plot(timevec, PsthMaster_Before_CSM(i,:),'Color',csmcolor, 'LineWidth', 1.5); 
    % % Plot NoRR PSTH
    % plot(timevec, PsthMaster_Before_NoRR(i,:),'Color', nocsmcolor, 'LineWidth', 1.5);
    
    
    % Labels and titles
    xlabel('Time frim Stim (ms)');
    xline(0,'--')
    xline(50,'--')
    ylabel('Firing Rate (Hz)');
    title(sprintf('Cell %d PSTH (Before Peak)', i));
    legend({'Baseline', 'Stim','Wash'}, 'Location', 'Best');
    
    grid off;
    hold off;

     filename=sprintf('Before_BaselinevsStim%d.svg',i);
 
    full_path=fullfile(folder_path,filename);
    saveas(gcf,full_path);

    close;
end


%% Plotting individual cells - After Peak V 
% === PLOT INDIVIDUAL PSTHs ===
folder_path="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\StimAnimals\Figures\posthreshaligned_individualunits\AfterCells";
% Assuming timevec is the same as used before

num_bins = size(Psth_BL_CSM_after, 2);
timevec = linspace(-1000, 1000, num_bins);

nCells = size(Psth_BL_CSM_after, 1);
blcolor=rgb('Black');
stimcolor=rgb('DodgerBlue');
washcolor=rgb('Magenta');

smooth_window = 3;

for i = 1:nCells
    figure;
    hold on;

    psth_csm = movmean(Psth_BL_CSM_after(i,:), smooth_window);
    psth_norr = movmean(Psth_AD_CSM_after(i,:), smooth_window);
    psth_wo = movmean(Psth_WOE_CSM_after(i,:), smooth_window);
    
    % Plot CSM PSTH
    plot(timevec, psth_csm, 'Color', blcolor, 'LineWidth', 1.5);
    
    % Plot NoRR PSTH
    plot(timevec, psth_norr, 'Color', stimcolor, 'LineWidth', 1.5);

    plot(timevec, psth_wo, 'Color', washcolor, 'LineWidth', 1.5);
    
    
    % % Plot CSM PSTH
    % plot(timevec, PsthMaster_Before_CSM(i,:),'Color',csmcolor, 'LineWidth', 1.5); 
    % % Plot NoRR PSTH
    % plot(timevec, PsthMaster_Before_NoRR(i,:),'Color', nocsmcolor, 'LineWidth', 1.5);
    
    
    % Labels and titles
    xlabel('Time frim Stim (ms)');
    xline(0,'--')
    xline(50,'--')
    ylabel('Firing Rate (Hz)');
    title(sprintf('Cell %d PSTH (After Peak)', i));
    legend({'Baseline', 'Stim','Wash'}, 'Location', 'Best');
    
    grid off;
    hold off;

     filename=sprintf('After_BaselinevsStim%d.svg',i);
 
    full_path=fullfile(folder_path,filename);
    saveas(gcf,full_path);

    close;
end

%% Calculate the SEM for all the PSTHS for individual cells - AFTER 
nBoot = 1000;
bin_size = 10; % 10ms bin size
smoothing_bins = 3; % Moving average over 3 bins

%BinaryMaster = NPdata(10).cells.units(28).BLCSM_pos_thresh.Raster; % Good after cells 
BinaryMaster = NPdata(10).cells.units(33).BLCSM_pos_thresh.Raster; % Good after cells 
%BinaryMaster = NPdata(12).cells.units(43).BLCSM_pos_thresh.Raster; % Good after cells 

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


%BinaryMaster = NPdata(10).cells.units(28).ADCSM_pos_thresh.Raster; % good after cell 
BinaryMaster = NPdata(10).cells.units(33).ADCSM_pos_thresh.Raster; % goodbefore cell 
%BinaryMaster = NPdata(12).cells.units(43).ADCSM_pos_thresh.Raster; % goodbefore cell 

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

%% Calculate the SEM for all the PSTHS for individual cells 
nBoot = 1000;
bin_size = 10; % 10ms bin size
smoothing_bins = 3; % Moving average over 3 bins

%BinaryMaster =NPdata(10).cells.units(28).WOCSM_pos_thresh.Raster;%aftercell
BinaryMaster = NPdata(10).cells.units(33).WOCSM_pos_thresh.Raster; %before
%BinaryMaster = NPdata(12).cells.units(43).WOCSM_pos_thresh.Raster; %before

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
[mean_psth_boot_CSMWO, sem_psth_boot_CSMWO, time_bins_CSMWO] = bootstrapPSTH_smooth(BinaryBinned, bin_size, smoothing_bins, nBoot);


%%
x1=0 ;x2=50;
%y_upper=250; y_lower=0;
y_upper=160; y_lower=0;
%figure;hold on;
figure('Position', [100, 100, 600, 300]);hold on

% Define SlateBlue color
slateBlue = rgb('DodgerBlue');
WOcolor=rgb('Violet');

% Plot NoCSM mean PSTH with shaded SEM (black)
shadedErrorBar(time_bins, mean_psth_boot, sem_psth_boot, ...
    'lineprops', {'-k', 'LineWidth', 1.5}, 'transparent', true);

% Plot CSM mean PSTH with shaded SEM (SlateBlue)
shadedErrorBar(time_bins_CSM, mean_psth_boot_CSM, sem_psth_boot_CSM, ...
    'lineprops', {'-', 'Color', slateBlue, 'LineWidth', 1.5}, 'transparent', true);
% 
% % Plot CSM mean PSTH with shaded SEM (SlateBlue)
shadedErrorBar(time_bins_CSMWO, mean_psth_boot_CSMWO, sem_psth_boot_CSMWO, ...
    'lineprops', {'-', 'Color', WOcolor, 'LineWidth', 1.5}, 'transparent', true);

hold on 
% Fill the area between the vertical lines
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');

% Formatting
xlabel('Time from stim (ms)');
ylabel('Firing Rate (Hz)');
xline(81.5)
title('Bootstrapped Smoothed PSTH');


set(gca, 'FontSize', 12, 'LineWidth', 1);

% Set axis limits: -250 to 250 ms
%ylim([0 250])
ylim([0 160])
xlim([-250, 250]);
%ylim([0, max([mean_psth_boot + sem_psth_boot, mean_psth_boot_CSM + sem_psth_boot_CSM]) * 1.2]);

box off;




%% Define colors
csmmean=81.5;
Color_BL_CSM      = rgb('Black');
Color_AD_CSM      = rgb('DodgerBlue');
Color_WOE_CSM     = rgb('Violet');

Color_BL_NoCSM    = rgb('Gray');
Color_AD_NoCSM    = rgb('LightSteelBlue');
Color_WOE_NoCSM   = rgb('Magenta');

% Compute average PSTHs and SEMs — BEFORE PEAK
avg_before_BL_CSM    = nanmean(Psth_BL_CSM_before, 1);
avg_before_AD_CSM    = nanmean(Psth_AD_CSM_before, 1);
avg_before_WOE_CSM   = nanmean(Psth_WOE_CSM_before, 1);

avg_before_BL_NoCSM  = nanmean(Psth_BL_NoCSM_before, 1);
avg_before_AD_NoCSM  = nanmean(Psth_AD_NoCSM_before, 1);
avg_before_WOE_NoCSM = nanmean(Psth_WOE_NoCSM_before, 1);

sem_before_BL_CSM    = nanstd(Psth_BL_CSM_before, 0, 1) / sqrt(size(Psth_BL_CSM_before, 1));
sem_before_AD_CSM    = nanstd(Psth_AD_CSM_before, 0, 1) / sqrt(size(Psth_AD_CSM_before, 1));
sem_before_WOE_CSM   = nanstd(Psth_WOE_CSM_before, 0, 1) / sqrt(size(Psth_WOE_CSM_before, 1));

sem_before_BL_NoCSM  = nanstd(Psth_BL_NoCSM_before, 0, 1) / sqrt(size(Psth_BL_NoCSM_before, 1));
sem_before_AD_NoCSM  = nanstd(Psth_AD_NoCSM_before, 0, 1) / sqrt(size(Psth_AD_NoCSM_before, 1));
sem_before_WOE_NoCSM = nanstd(Psth_WOE_NoCSM_before, 0, 1) / sqrt(size(Psth_WOE_NoCSM_before, 1));

% Compute average PSTHs and SEMs — AFTER PEAK
avg_after_BL_CSM    = nanmean(Psth_BL_CSM_after, 1);
avg_after_AD_CSM    = nanmean(Psth_AD_CSM_after, 1);
avg_after_WOE_CSM   = nanmean(Psth_WOE_CSM_after, 1);

avg_after_BL_NoCSM  = nanmean(Psth_BL_NoCSM_after, 1);
avg_after_AD_NoCSM  = nanmean(Psth_AD_NoCSM_after, 1);
avg_after_WOE_NoCSM = nanmean(Psth_WOE_NoCSM_after, 1);

sem_after_BL_CSM    = nanstd(Psth_BL_CSM_after, 0, 1) / sqrt(size(Psth_BL_CSM_after, 1));
sem_after_AD_CSM    = nanstd(Psth_AD_CSM_after, 0, 1) / sqrt(size(Psth_AD_CSM_after, 1));
sem_after_WOE_CSM   = nanstd(Psth_WOE_CSM_after, 0, 1) / sqrt(size(Psth_WOE_CSM_after, 1));

sem_after_BL_NoCSM  = nanstd(Psth_BL_NoCSM_after, 0, 1) / sqrt(size(Psth_BL_NoCSM_after, 1));
sem_after_AD_NoCSM  = nanstd(Psth_AD_NoCSM_after, 0, 1) / sqrt(size(Psth_AD_NoCSM_after, 1));
sem_after_WOE_NoCSM = nanstd(Psth_WOE_NoCSM_after, 0, 1) / sqrt(size(Psth_WOE_NoCSM_after, 1));

% Time vector
num_bins = size(Psth_BL_CSM_before, 2);
time_vector = linspace(-1000, 1000, num_bins);

x1=-50 ;x2=50;
y_upper=120; y_lower=20;

figure('Position', [100, 100, 600, 300]);hold on
% === PLOT: BEFORE PEAK UNITS ===
%figure; hold on;

% Lines
plot(time_vector, avg_before_BL_CSM,  'Color', Color_BL_CSM, 'LineWidth', 2);
plot(time_vector, avg_before_AD_CSM,  'Color', Color_AD_CSM, 'LineWidth', 2);
plot(time_vector, avg_before_WOE_CSM, 'Color', Color_WOE_CSM,'LineWidth', 2);

% plot(time_vector, avg_before_BL_NoCSM,  'Color', Color_BL_NoCSM,  'LineWidth', 2);
% plot(time_vector, avg_before_AD_NoCSM,  'Color', Color_AD_NoCSM,  'LineWidth', 2);
% plot(time_vector, avg_before_WOE_NoCSM, 'Color', Color_WOE_NoCSM, 'LineWidth', 2);

% Shading
fill([time_vector, fliplr(time_vector)], ...
     [avg_before_BL_CSM + sem_before_BL_CSM, fliplr(avg_before_BL_CSM - sem_before_BL_CSM)], ...
     Color_BL_CSM, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([time_vector, fliplr(time_vector)], ...
     [avg_before_AD_CSM + sem_before_AD_CSM, fliplr(avg_before_AD_CSM - sem_before_AD_CSM)], ...
     Color_AD_CSM, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([time_vector, fliplr(time_vector)], ...
     [avg_before_WOE_CSM + sem_before_WOE_CSM, fliplr(avg_before_WOE_CSM - sem_before_WOE_CSM)], ...
     Color_WOE_CSM, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% fill([time_vector, fliplr(time_vector)], ...
%      [avg_before_BL_NoCSM + sem_before_BL_NoCSM, fliplr(avg_before_BL_NoCSM - sem_before_BL_NoCSM)], ...
%      Color_BL_NoCSM, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% fill([time_vector, fliplr(time_vector)], ...
%      [avg_before_AD_NoCSM + sem_before_AD_NoCSM, fliplr(avg_before_AD_NoCSM - sem_before_AD_NoCSM)], ...
%      Color_AD_NoCSM, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% fill([time_vector, fliplr(time_vector)], ...
%      [avg_before_WOE_NoCSM + sem_before_WOE_NoCSM, fliplr(avg_before_WOE_NoCSM - sem_before_WOE_NoCSM)], ...
%      Color_WOE_NoCSM, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on 
% Fill the area between the vertical lines
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
%ylim([0 14])
% xline(csmmean,'--')
% xline(0,'--')
title('PSTH – Before Peak Firing');
xlabel('Time from CSM (ms)');
%xlabel('Time from Stim (ms)');
%xlabel('Time from CSM(ms)');
ylabel('Firing rate (Hz)');
xlim([-250 250])
box off;

% === PLOT: AFTER PEAK UNITS ===
%figure; hold on;

figure('Position', [100, 100, 600, 300]);hold on
y_upper=160; y_lower=40;

% Lines
plot(time_vector, avg_after_BL_CSM,  'Color', Color_BL_CSM, 'LineWidth', 2);
plot(time_vector, avg_after_AD_CSM,  'Color', Color_AD_CSM, 'LineWidth', 2);
plot(time_vector, avg_after_WOE_CSM, 'Color', Color_WOE_CSM,'LineWidth', 2);

% plot(time_vector, avg_after_BL_NoCSM,  'Color', Color_BL_NoCSM,  'LineWidth', 2);
% plot(time_vector, avg_after_AD_NoCSM,  'Color', Color_AD_NoCSM,  'LineWidth', 2);
% plot(time_vector, avg_after_WOE_NoCSM, 'Color', Color_WOE_NoCSM, 'LineWidth', 2);

% Shading
fill([time_vector, fliplr(time_vector)], ...
     [avg_after_BL_CSM + sem_after_BL_CSM, fliplr(avg_after_BL_CSM - sem_after_BL_CSM)], ...
     Color_BL_CSM, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([time_vector, fliplr(time_vector)], ...
     [avg_after_AD_CSM + sem_after_AD_CSM, fliplr(avg_after_AD_CSM - sem_after_AD_CSM)], ...
     Color_AD_CSM, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([time_vector, fliplr(time_vector)], ...
     [avg_after_WOE_CSM + sem_after_WOE_CSM, fliplr(avg_after_WOE_CSM - sem_after_WOE_CSM)], ...
     Color_WOE_CSM, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% fill([time_vector, fliplr(time_vector)], ...
%      [avg_after_BL_NoCSM + sem_after_BL_NoCSM, fliplr(avg_after_BL_NoCSM - sem_after_BL_NoCSM)], ...
%      Color_BL_NoCSM, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% fill([time_vector, fliplr(time_vector)], ...
%      [avg_after_AD_NoCSM + sem_after_AD_NoCSM, fliplr(avg_after_AD_NoCSM - sem_after_AD_NoCSM)], ...
%      Color_AD_NoCSM, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% fill([time_vector, fliplr(time_vector)], ...
%      [avg_after_WOE_NoCSM + sem_after_WOE_NoCSM, fliplr(avg_after_WOE_NoCSM - sem_after_WOE_NoCSM)], ...
%      Color_WOE_NoCSM, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on 
% Fill the area between the vertical lines
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
%ylim([0 14])
%xline(csmmean,'--')
xline(-50,'--')
xline(50,'--')
title('PSTH – After Peak Firing');
%xlabel('Time from Stim (ms)');
xlabel('Time from CSM (ms)');
ylabel('Firing rate (Hz)');
xlim([-250 250])
box off;
%% 
% === Heatmap Panel for Normalized PSTHs (CSM Only) ===================== %%
% Crop time vector to -250 to 250 ms
time_window_idx = find(time_vector >= -250 & time_vector <= 250);
%time_window_idx = find(time_vector >= -1000 & time_vector <= 1000);
time_crop = time_vector(time_window_idx);

% % Normalize color scale
% vmin = 0;
% vmax = 1;

% --- Sort by time of peak within cropped window
[~, peak_sort_BL_b]  = max(Psth_BL_CSM_before(:, time_window_idx), [], 2);
[~, peak_sort_AD_b]  = max(Psth_AD_CSM_before(:, time_window_idx), [], 2);
[~, peak_sort_WOE_b] = max(Psth_WOE_CSM_before(:, time_window_idx), [], 2);
[~, idx_sort_BL_b]  = sort(peak_sort_BL_b);
[~, idx_sort_AD_b]  = sort(peak_sort_AD_b);
[~, idx_sort_WOE_b] = sort(peak_sort_WOE_b);

[~, peak_sort_BL_a]  = max(Psth_BL_CSM_after(:, time_window_idx), [], 2);
[~, peak_sort_AD_a]  = max(Psth_AD_CSM_after(:, time_window_idx), [], 2);
[~, peak_sort_WOE_a] = max(Psth_WOE_CSM_after(:, time_window_idx), [], 2);
[~, idx_sort_BL_a]  = sort(peak_sort_BL_a);
[~, idx_sort_AD_a]  = sort(peak_sort_AD_a);
[~, idx_sort_WOE_a] = sort(peak_sort_WOE_a);

% --- Initialize figure
figure('Position', [100, 100, 1000, 600]);

% === TOP ROW: BEFORE PEAK UNITS ========================================
subplot(2,3,1);
imagesc(time_crop, 1:size(Psth_BL_CSM_before,1), ...
        Psth_BL_CSM_before(idx_sort_BL_b, time_window_idx));
title('BL – Before Peak');
xlabel('Time from CSM (ms)');
ylabel('Unit #');
xline(0, '--k'); set(gca, 'YDir', 'normal');

subplot(2,3,2);
imagesc(time_crop, 1:size(Psth_AD_CSM_before,1), ...
        Psth_AD_CSM_before(idx_sort_AD_b, time_window_idx));
title('AD – Before Peak');
xlabel('Time from CSM (ms)');
xline(0, '--k'); set(gca, 'YDir', 'normal'); yticklabels([]);

subplot(2,3,3);
imagesc(time_crop, 1:size(Psth_WOE_CSM_before,1), ...
        Psth_WOE_CSM_before(idx_sort_WOE_b, time_window_idx));
title('WOE – Before Peak');
xlabel('Time from CSM (ms)');
xline(0, '--k'); set(gca, 'YDir', 'normal'); yticklabels([]);

% === BOTTOM ROW: AFTER PEAK UNITS ======================================
subplot(2,3,4);
imagesc(time_crop, 1:size(Psth_BL_CSM_after,1), ...
        Psth_BL_CSM_after(idx_sort_BL_a, time_window_idx));
title('BL – After Peak');
xlabel('Time from CSM (ms)');
ylabel('Unit #');
xline(0, '--k'); set(gca, 'YDir', 'normal');

subplot(2,3,5);
imagesc(time_crop, 1:size(Psth_AD_CSM_after,1), ...
        Psth_AD_CSM_after(idx_sort_AD_a, time_window_idx));
title('AD – After Peak');
xlabel('Time from CSM (ms)');
xline(0, '--k'); set(gca, 'YDir', 'normal'); yticklabels([]);

subplot(2,3,6);
imagesc(time_crop, 1:size(Psth_WOE_CSM_after,1), ...
        Psth_WOE_CSM_after(idx_sort_WOE_a, time_window_idx));
title('WOE – After Peak');
xlabel('Time from CSM (ms)');
xline(0, '--k'); set(gca, 'YDir', 'normal'); yticklabels([]);

% Shared colorbar
cb = colorbar('Position', [0.92, 0.11, 0.02, 0.77]);
ylabel(cb, 'Normalized Firing Rate (max = 1)');

sgtitle('Normalized Unit Activity Aligned to CSM (–250 to 250 ms)');

%%
% --- Clip time range to -250 to 250 ms ---
time_window = [-250, 250];
time_mask = (timevec >= time_window(1)) & (timevec <= time_window(2));
xvec = timevec(time_mask);

% --- Clip data to time window ---
Psth_BL_CSM_before1 = Psth_BL_CSM_before(:, time_mask);
Psth_AD_CSM_before1 = Psth_AD_CSM_before(:, time_mask);
Psth_WOE_CSM_before1 = Psth_WOE_CSM_before(:, time_mask);

Psth_BL_CSM_after1 = Psth_BL_CSM_after(:, time_mask);
Psth_AD_CSM_after1 = Psth_AD_CSM_after(:, time_mask);
Psth_WOE_CSM_after1 = Psth_WOE_CSM_after(:, time_mask);

% % --- Combine data and compute global limits ---
% all_data_combined = [Psth_BL_CSM_before1(:); Psth_AD_CSM_before1(:); Psth_WOE_CSM_before1(:);
%                      Psth_BL_CSM_after1(:);  Psth_AD_CSM_after1(:);  Psth_WOE_CSM_after1(:)];
% vmin = min(all_data_combined);
% vmax = max(all_data_combined);

% --- Sorting ---
[~, idx_sort_BL_before]  = sort(max(Psth_BL_CSM_before1, [], 2));
[~, idx_sort_AD_before]  = sort(max(Psth_AD_CSM_before1, [], 2));
[~, idx_sort_WOE_before] = sort(max(Psth_WOE_CSM_before1, [], 2));
[~, idx_sort_BL_after]   = sort(max(Psth_BL_CSM_after1, [], 2));
[~, idx_sort_AD_after]   = sort(max(Psth_AD_CSM_after1, [], 2));
[~, idx_sort_WOE_after]  = sort(max(Psth_WOE_CSM_after1, [], 2));

% --- Plot ---
figure('Position', [100, 100, 1200, 600]);

% Top Row: BEFORE
subplot(2,3,1); imagesc(xvec, 1:size(Psth_BL_CSM_before1,1), Psth_BL_CSM_before1(idx_sort_BL_before,:));
title('BL - Before'); xlabel('Time (ms)'); ylabel('Units'); formatPlot();
% 
% subplot(2,3,2); imagesc(xvec, 1:size(Psth_AD_CSM_before,1), Psth_AD_CSM_before(idx_sort_AD_before,:), [vmin vmax]);
% title('AD - Before'); xlabel('Time (ms)'); ylabel('Units'); formatPlot();
subplot(2,3,2); imagesc(xvec, 1:size(Psth_AD_CSM_before1,1), Psth_AD_CSM_before1(idx_sort_BL_before,:));
title('AD - Before'); xlabel('Time (ms)'); ylabel('Units'); formatPlot();

% subplot(2,3,3); imagesc(xvec, 1:size(Psth_WOE_CSM_before,1), Psth_WOE_CSM_before(idx_sort_WOE_before,:), [vmin vmax]);
% title('WOE - Before'); xlabel('Time (ms)'); ylabel('Units'); formatPlot();
subplot(2,3,3); imagesc(xvec, 1:size(Psth_WOE_CSM_before1,1), Psth_WOE_CSM_before1(idx_sort_BL_before,:));
title('WOE - Before'); xlabel('Time (ms)'); ylabel('Units'); formatPlot();

% Bottom Row: AFTER
subplot(2,3,4); imagesc(xvec, 1:size(Psth_BL_CSM_after1,1), Psth_BL_CSM_after1(idx_sort_BL_after,:));
title('BL - After'); xlabel('Time (ms)'); ylabel('Units'); formatPlot();

% subplot(2,3,5); imagesc(xvec, 1:size(Psth_AD_CSM_after,1), Psth_AD_CSM_after(idx_sort_AD_after,:), [vmin vmax]);
% title('AD - After'); xlabel('Time (ms)'); ylabel('Units'); formatPlot();
subplot(2,3,5); imagesc(xvec, 1:size(Psth_AD_CSM_after1,1), Psth_AD_CSM_after1(idx_sort_BL_after,:));
title('AD - After'); xlabel('Time (ms)'); ylabel('Units'); formatPlot();
% 
% subplot(2,3,6); imagesc(xvec, 1:size(Psth_WOE_CSM_after,1), Psth_WOE_CSM_after(idx_sort_WOE_after,:), [vmin vmax]);
% title('WOE - After'); xlabel('Time (ms)'); ylabel('Units'); formatPlot();
subplot(2,3,6); imagesc(xvec, 1:size(Psth_WOE_CSM_after1,1), Psth_WOE_CSM_after1(idx_sort_BL_after,:));
title('WOE - After'); xlabel('Time (ms)'); ylabel('Units'); formatPlot();

% Colorbar
cb = colorbar('Position', [0.93, 0.11, 0.02, 0.75]);
ylabel(cb, 'Normalized firing rate');
sgtitle('Normalized Firing Rates by Group and Peak Timing');


%% === Compute Max Firing Per Unit: BEFORE PEAK ===
max_FR_before_BL  = max(Psth_BL_CSM_before, [], 2, 'omitnan');
max_FR_before_AD  = max(Psth_AD_CSM_before, [], 2, 'omitnan');
max_FR_before_WOE = max(Psth_WOE_CSM_before, [], 2, 'omitnan');

% === Compute Max Firing Per Unit: AFTER PEAK ===
max_FR_after_BL  = max(Psth_BL_CSM_after, [], 2, 'omitnan');
max_FR_after_AD  = max(Psth_AD_CSM_after, [], 2, 'omitnan');
max_FR_after_WOE = max(Psth_WOE_CSM_after, [], 2, 'omitnan');

% === Plot Max Firing Rate — BEFORE PEAK ===
figure('Position', [100, 100, 400, 300]); hold on;

% X positions for each group
x_BL = 1; x_AD = 2; x_WOE = 3;

% Plot individual points
scatter(repelem(x_BL, length(max_FR_before_BL)), max_FR_before_BL, 20, [0.5 0.5 0.5], 'filled');
scatter(repelem(x_AD, length(max_FR_before_AD)), max_FR_before_AD, 20, [0.5 0.5 0.5], 'filled');
scatter(repelem(x_WOE, length(max_FR_before_WOE)), max_FR_before_WOE, 20, [0.5 0.5 0.5], 'filled');

% Plot group means with SEM
errorbar(x_BL, mean(max_FR_before_BL, 'omitnan'), sem(max_FR_before_BL), 'o', ...
         'Color', 'k', 'MarkerFaceColor', Color_BL_CSM, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerSize', 10);
errorbar(x_AD, mean(max_FR_before_AD, 'omitnan'), sem(max_FR_before_AD), 'o', ...
         'Color', 'k', 'MarkerFaceColor', Color_AD_CSM, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerSize', 10);
errorbar(x_WOE, mean(max_FR_before_WOE, 'omitnan'), sem(max_FR_before_WOE), 'o', ...
         'Color', 'k', 'MarkerFaceColor', Color_WOE_CSM, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerSize', 10);

xlim([0.5 3.5])
xticks([1 2 3])
xticklabels({'BL', 'AD', 'WOE'})
ylabel('Max firing rate (Hz)')
title('Max FR Before Peak')
box off

% === Plot Max Firing Rate — AFTER PEAK ===
figure('Position', [100, 100, 400, 300]); hold on;

% Plot individual points
scatter(repelem(x_BL, length(max_FR_after_BL)), max_FR_after_BL, 20, [0.5 0.5 0.5], 'filled');
scatter(repelem(x_AD, length(max_FR_after_AD)), max_FR_after_AD, 20, [0.5 0.5 0.5], 'filled');
scatter(repelem(x_WOE, length(max_FR_after_WOE)), max_FR_after_WOE, 20, [0.5 0.5 0.5], 'filled');

% Plot group means with SEM
errorbar(x_BL, mean(max_FR_after_BL, 'omitnan'), sem(max_FR_after_BL), 'o', ...
         'Color', 'k', 'MarkerFaceColor', Color_BL_CSM, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerSize', 10);
errorbar(x_AD, mean(max_FR_after_AD, 'omitnan'), sem(max_FR_after_AD), 'o', ...
         'Color', 'k', 'MarkerFaceColor', Color_AD_CSM, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerSize', 10);
errorbar(x_WOE, mean(max_FR_after_WOE, 'omitnan'), sem(max_FR_after_WOE), 'o', ...
         'Color', 'k', 'MarkerFaceColor', Color_WOE_CSM, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerSize', 10);

xlim([0.5 3.5])
xticks([1 2 3])
xticklabels({'BL', 'AD', 'WOE'})
ylabel('Max firing rate (Hz)')
title('Max FR After Peak')
box off


%% === Compute Max Firing Per Unit: BEFORE PEAK ===
max_FR_before_BL  = max(Psth_BL_CSM_before, [], 2, 'omitnan');
max_FR_before_AD  = max(Psth_AD_CSM_before, [], 2, 'omitnan');
max_FR_before_WOE = max(Psth_WOE_CSM_before, [], 2, 'omitnan');

WOE_relative_before=max_FR_before_WOE-max_FR_before_BL;
AD_relative_before=max_FR_before_AD-max_FR_before_BL;

figure; hold on
scatter(AD_relative_before,WOE_relative_before)


% === Compute Max Firing Per Unit: AFTER PEAK ===
max_FR_after_BL  = max(Psth_BL_CSM_after, [], 2, 'omitnan');
max_FR_after_AD  = max(Psth_AD_CSM_after, [], 2, 'omitnan');
max_FR_after_WOE = max(Psth_WOE_CSM_after, [], 2, 'omitnan');

WOE_relative_after=max_FR_after_WOE-max_FR_after_BL;
AD_relative_after=max_FR_after_AD-max_FR_after_BL;

figure; hold on
scatter(AD_relative_after,WOE_relative_after)



%% === Define Time Window ===
% Adjust to match your actual time vector variable name
window_start = 0;  % in ms
window_end   = 50; % in ms

% Assume `time_vector` corresponds to columns in Psth matrices
window_idx = (time_vector >= window_start & time_vector <= window_end);

% === Compute Mean FR in Time Window ===
mean_FR_before_BL  = mean(Psth_BL_CSM_before(:, window_idx), 2, 'omitnan');
mean_FR_before_AD  = mean(Psth_AD_CSM_before(:, window_idx), 2, 'omitnan');
mean_FR_before_WOE = mean(Psth_WOE_CSM_before(:, window_idx), 2, 'omitnan');

mean_FR_after_BL   = mean(Psth_BL_CSM_after(:, window_idx), 2, 'omitnan');
mean_FR_after_AD   = mean(Psth_AD_CSM_after(:, window_idx), 2, 'omitnan');
mean_FR_after_WOE  = mean(Psth_WOE_CSM_after(:, window_idx), 2, 'omitnan');

% === Compute Modulation Index in Window ===
MI_AD_before  = (mean_FR_before_AD  - mean_FR_before_BL)  ./ (mean_FR_before_AD  + mean_FR_before_BL);
MI_WOE_before = (mean_FR_before_WOE - mean_FR_before_BL) ./ (mean_FR_before_WOE + mean_FR_before_BL);

MI_AD_after   = (mean_FR_after_AD  - mean_FR_after_BL)  ./ (mean_FR_after_AD  + mean_FR_after_BL);
MI_WOE_after  = (mean_FR_after_WOE - mean_FR_after_BL) ./ (mean_FR_after_WOE + mean_FR_after_BL);

% === Prepare Data for Box Plot ===
MI_values = [MI_AD_before; MI_WOE_before; MI_AD_after; MI_WOE_after];
group_labels = [...
    repmat({'AD Before'},  length(MI_AD_before),  1);
    repmat({'WOE Before'}, length(MI_WOE_before), 1);
    repmat({'AD After'},   length(MI_AD_after),   1);
    repmat({'WOE After'},  length(MI_WOE_after),  1)];

% === Plot Box Plot ===
figure;
boxplot(MI_values, group_labels, 'Colors', 'k', 'Symbol', 'o', 'Whisker', 1.5);
ylabel('Modulation Index (Mean FR in Window)');
title(sprintf('Modulation Index in %d–%d ms Window', window_start, window_end));
ylim([-1 1]);
grid on;
set(gca, 'FontSize', 12);



%% Diff for before and after csm vs nocsm
beforecolor=rgb('OrangeRed');
aftercolor=rgb('DarkBlue');

% % 
% diff_before=Psth_WOE_CSM_before-Psth_BL_CSM_before;
% diff_after=Psth_WOE_CSM_after-Psth_BL_CSM_after;
% 
diff_before=Psth_AD_CSM_before-Psth_BL_CSM_before;
diff_after=Psth_AD_CSM_after-Psth_BL_CSM_after;

avg_psth_diff_before = nanmean(diff_before, 1);  % 1 = average across rows (units)
avg_psth_diff_after = nanmean(diff_after, 1);

% Optional: Compute standard error of the mean (SEM)
sem_psth_diff_before = nanstd(diff_before, 0, 1) / sqrt(size(diff_before, 1));
sem_psth_diff_after = nanstd(diff_after, 0, 1) / sqrt(size(diff_after, 1));

% Generate time vector (assuming same binning and length)
num_bins = size(diff_before, 2);
time_vector = linspace(-1000, 1000, num_bins);  % Update if your window or binning is different
x1=-50 ;x2=50;
y_upper=20; y_lower=-50;

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
% Fill the area between the vertical lines
%fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
%ylim([0 1]);
yline(0,'--')
% yline(-2,':')
% yline(2,':')
% xline(csmmean,'--')
% xline(0,'--')
xlim([-250 250])
%ylim([-6 4])
%ylim([-25 25])
xlabel('Time from CSM(ms)');
%xlabel('Time From Stim (ms)')
ylabel('Firing rate (Hz)');
%title('Av. Diff in FR WOE-BL CSM for before and after');
title('Av. Diff in FR AD-BL CSM for before and after');
%% Diff for before and after csm vs nocsm


stimcolor=rgb('DodgerBlue');
washcolor=rgb('Violet');


% % 
diff_before=Psth_WOE_CSM_before-Psth_BL_CSM_before;
diff_after=Psth_WOE_CSM_after-Psth_BL_CSM_after;
% 
diff_before_AD=Psth_AD_CSM_before-Psth_BL_CSM_before;
diff_after_AD=Psth_AD_CSM_after-Psth_BL_CSM_after;

avg_psth_diff_before = nanmean(diff_before, 1);  % 1 = average across rows (units)
avg_psth_diff_after = nanmean(diff_after, 1);

avg_psth_diff_before_AD = nanmean(diff_before_AD, 1);  % 1 = average across rows (units)
avg_psth_diff_after_AD = nanmean(diff_after_AD, 1);

% Optional: Compute standard error of the mean (SEM)
sem_psth_diff_before = nanstd(diff_before, 0, 1) / sqrt(size(diff_before, 1));
sem_psth_diff_after = nanstd(diff_after, 0, 1) / sqrt(size(diff_after, 1));


sem_psth_diff_before_AD = nanstd(diff_before_AD, 0, 1) / sqrt(size(diff_before_AD, 1));
sem_psth_diff_after_AD = nanstd(diff_after_AD, 0, 1) / sqrt(size(diff_after_AD, 1));

% Generate time vector (assuming same binning and length)
num_bins = size(diff_before, 2);
time_vector = linspace(-1000, 1000, num_bins);  % Update if your window or binning is different
x1=-50 ;x2=50;
y_upper=40; y_lower=-50;

% Plot average PSTHs
%figure;
figure('Position', [100, 100, 600, 300]); 
hold on;
plot(time_vector, avg_psth_diff_before,'Color',washcolor, 'LineWidth', 2);
plot(time_vector, avg_psth_diff_before_AD,'Color',stimcolor, 'LineWidth', 2);

% Optional: Add shaded error bars for SEM
fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_diff_before + sem_psth_diff_before, fliplr(avg_psth_diff_before - sem_psth_diff_before)], ...
     washcolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_diff_before_AD + sem_psth_diff_before_AD, fliplr(avg_psth_diff_before_AD - sem_psth_diff_before_AD)], ...
     stimcolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on 
% Fill the area between the vertical lines
%fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
%fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.5    0.5    0.5],'FaceAlpha',0.2,'EdgeColor','none');
%ylim([0 1]);
yline(0,'--')
% yline(-2,':')
% yline(2,':')
%xline(csmmean,'--')
% xline(0,'--')
xlim([-250 250])
%ylim([-6 4])
%ylim([-50 40])
xlabel('Time from CSM(ms)');
%xlabel('Time From Stim (ms)')
ylabel('Firing rate (Hz)');
title('Av. Diff for Before cells' );
%title('Av. Diff in FR AD-BL CSM for before and after');
%%
figure('Position', [100, 100, 600, 300]); 
hold on;
plot(time_vector, avg_psth_diff_after,'Color',washcolor, 'LineWidth', 2);
plot(time_vector, avg_psth_diff_after_AD,'Color',stimcolor, 'LineWidth', 2);

% Optional: Add shaded error bars for SEM
fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_diff_after + sem_psth_diff_after, fliplr(avg_psth_diff_after - sem_psth_diff_after)], ...
     washcolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

fill([time_vector, fliplr(time_vector)], ...
     [avg_psth_diff_after_AD + sem_psth_diff_after_AD, fliplr(avg_psth_diff_after_AD - sem_psth_diff_after_AD)], ...
     stimcolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on 
% Fill the area between the vertical lines
%fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.5    0.5    0.5],'FaceAlpha',0.2,'EdgeColor','none');
%ylim([0 1]);
yline(0,'--')
% yline(-2,':')
% yline(2,':')
%xline(csmmean,'--')
% xline(0,'--')
xlim([-250 250])
%ylim([-6 4])
ylim([-50 40])
xlabel('Time from CSM(ms)');
%xlabel('Time From Stim (ms)')
ylabel('Firing rate (Hz)');
title('Av. Diff for after cells' );
%title('Av. Diff in FR AD-BL CSM for before and after');

%%
% === Set up ===
% Define time window and time vector
time_zoom_idx = time_vector >= -250 & time_vector <= 250;
time_zoom = time_vector(time_zoom_idx);

% Choose diverging colormap
cmap = colormap(flipud(slanCM('RdBu')));

% === Compute differences ===
diff_WOE_before = Psth_WOE_CSM_before - Psth_BL_CSM_before;
diff_WOE_after  = Psth_WOE_CSM_after  - Psth_BL_CSM_after;
diff_AD_before  = Psth_AD_CSM_before  - Psth_BL_CSM_before;
diff_AD_after   = Psth_AD_CSM_after   - Psth_BL_CSM_after;

% === Plotting ===
figure('Position', [100, 100, 1200, 800]);

% --- Subplot 1: WOE-BL Before ---
subplot(2,2,1);
data = diff_WOE_before(:, time_zoom_idx);
[~, sortIdx] = sort(max(abs(data), [], 2), 'descend');
imagesc(time_zoom, 1:size(data,1), data(sortIdx,:));
set(gca, 'YDir', 'normal');
colormap(cmap);
colorbar;
%caxis([-8 8]);
hold on;
xline(0,'--k','LineWidth',1.5);
yl = ylim;
rectangle('Position',[50, yl(1), 50, yl(2)-yl(1)], ...
    'EdgeColor','k','LineWidth',1.2,'LineStyle','--', ...
    'FaceColor',[0.8 0.8 1], 'FaceAlpha', 0.4);
xlabel('Time from CSM (ms)');
ylabel('Units (sorted)');
title('WOE - BL Before CSM');

% --- Subplot 2: WOE-BL After ---
subplot(2,2,2);
data = diff_WOE_after(:, time_zoom_idx);
[~, sortIdx] = sort(max(abs(data), [], 2), 'descend');
imagesc(time_zoom, 1:size(data,1), data(sortIdx,:));
set(gca, 'YDir', 'normal');
colormap(cmap);
colorbar;
%caxis([-8 8]);
hold on;
xline(0,'--k','LineWidth',1.5);
yl = ylim;
rectangle('Position',[50, yl(1), 50, yl(2)-yl(1)], ...
    'EdgeColor','k','LineWidth',1.2,'LineStyle','--', ...
    'FaceColor',[0.8 0.8 1], 'FaceAlpha', 0.4);
xlabel('Time from CSM (ms)');
ylabel('Units (sorted)');
title('WOE - BL After CSM');

% --- Subplot 3: AD-BL Before ---
subplot(2,2,3);
data = diff_AD_before(:, time_zoom_idx);
[~, sortIdx] = sort(max(abs(data), [], 2), 'descend');
imagesc(time_zoom, 1:size(data,1), data(sortIdx,:));
set(gca, 'YDir', 'normal');
colormap(cmap);
colorbar;
%caxis([-8 8]);
hold on;
xline(0,'--k','LineWidth',1.5);
yl = ylim;
rectangle('Position',[50, yl(1), 50, yl(2)-yl(1)], ...
    'EdgeColor','k','LineWidth',1.2,'LineStyle','--', ...
    'FaceColor',[0.8 0.8 1], 'FaceAlpha', 0.4);
xlabel('Time from CSM (ms)');
ylabel('Units (sorted)');
title('AD - BL Before CSM');

% --- Subplot 4: AD-BL After ---
subplot(2,2,4);
data = diff_AD_after(:, time_zoom_idx);
[~, sortIdx] = sort(max(abs(data), [], 2), 'descend');
imagesc(time_zoom, 1:size(data,1), data(sortIdx,:));
set(gca, 'YDir', 'normal');
colormap(cmap);
colorbar;
%caxis([-8 8]);
hold on;
xline(0,'--k','LineWidth',1.5);
yl = ylim;
rectangle('Position',[50, yl(1), 50, yl(2)-yl(1)], ...
    'EdgeColor','k','LineWidth',1.2,'LineStyle','--', ...
    'FaceColor',[0.8 0.8 1], 'FaceAlpha', 0.4);
xlabel('Time from CSM (ms)');
ylabel('Units (sorted)');
title('AD - BL After CSM');
%%
% Define plot window
plot_window = [-250, 250];
plot_idx = time_vector >= plot_window(1) & time_vector <= plot_window(2);

% Also find peak window (unchanged, -50 to 50 ms)
time_window = [-100, -50];
%time_window = [0, 50];
window_idx = time_vector >= time_window(1) & time_vector <= time_window(2);

% ====== BEFORE Figure =======
figure('Position', [100, 100, 600, 300]); 
hold on;

% Curves and SEM (subset for plotting)
curve_AD_before = avg_psth_diff_before_AD(plot_idx);
sem_AD_before = sem_psth_diff_before_AD(plot_idx);
times_plot = time_vector(plot_idx);

curve_WOE_before = avg_psth_diff_before(plot_idx);
sem_WOE_before = sem_psth_diff_before(plot_idx);

% Plot shaded error bars
fill([times_plot, fliplr(times_plot)], ...
     [curve_AD_before + sem_AD_before, fliplr(curve_AD_before - sem_AD_before)], ...
     stimcolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([times_plot, fliplr(times_plot)], ...
     [curve_WOE_before + sem_WOE_before, fliplr(curve_WOE_before - sem_WOE_before)], ...
     washcolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Plot mean curves
plot(times_plot, curve_AD_before, 'Color', stimcolor, 'LineWidth', 2);
plot(times_plot, curve_WOE_before, 'Color', washcolor, 'LineWidth', 2);

% Find peaks in peak window (-50 to 50 ms) **index relative to full time_vector**
peak_win_idx = find(window_idx);  % indices of peak window relative to full vector
peak_win_idx_before = peak_win_idx(peak_win_idx >= find(plot_idx, 1, 'first') & peak_win_idx <= find(plot_idx, 1, 'last'));

% AD Before: max in window
%[peak_AD_before_val, peak_AD_before_idx_rel] = max(avg_psth_diff_before_AD(window_idx));
[peak_AD_before_val, peak_AD_before_idx_rel] = min(avg_psth_diff_before_AD(window_idx));
peak_AD_before_idx = find(window_idx, 1, 'first') - 1 + peak_AD_before_idx_rel;
peak_AD_before_time = time_vector(peak_AD_before_idx);
peak_AD_before_sem = sem_psth_diff_before_AD(peak_AD_before_idx);

% WOE Before: min in window
[peak_WOE_before_val, peak_WOE_before_idx_rel] = min(avg_psth_diff_before(window_idx));
peak_WOE_before_idx = find(window_idx, 1, 'first') - 1 + peak_WOE_before_idx_rel;
peak_WOE_before_time = time_vector(peak_WOE_before_idx);
peak_WOE_before_sem = sem_psth_diff_before(peak_WOE_before_idx);

% Plot peak markers (only if within plot window)
if peak_AD_before_time >= plot_window(1) && peak_AD_before_time <= plot_window(2)
    plot(peak_AD_before_time, peak_AD_before_val, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    text(peak_AD_before_time + 5, peak_AD_before_val, ...
        sprintf('%.2f ± %.2f Hz', peak_AD_before_val, peak_AD_before_sem), ...
        'Color', stimcolor, 'FontSize', 10);
end

if peak_WOE_before_time >= plot_window(1) && peak_WOE_before_time <= plot_window(2)
    plot(peak_WOE_before_time, peak_WOE_before_val, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    text(peak_WOE_before_time + 5, peak_WOE_before_val, ...
        sprintf('%.2f ± %.2f Hz', peak_WOE_before_val, peak_WOE_before_sem), ...
        'Color', washcolor, 'FontSize', 10);
end

% Formatting
yline(0, '--', 'Color', [0.4 0.4 0.4]);
xlabel('Time from CSM (ms)');
ylabel('Δ Firing Rate (Hz)');
title('Avg ΔFR Before-CSM');
xlim(plot_window);
ylim([min([curve_AD_before - sem_AD_before, curve_WOE_before - sem_WOE_before]) - 5, ...
      max([curve_AD_before + sem_AD_before, curve_WOE_before + sem_WOE_before]) + 5]);
legend({'AD SEM', 'WOE SEM', 'AD Mean', 'WOE Mean'}, 'Location', 'best');

hold off;

% Print peaks info BEFORE
fprintf('\nBEFORE-CSM peaks (within %d to %d ms):\n', time_window(1), time_window(2));
fprintf('AD Before: Peak max = %.2f ± %.2f Hz at %d ms\n', peak_AD_before_val, peak_AD_before_sem, round(peak_AD_before_time));
fprintf('WOE Before: Peak min = %.2f ± %.2f Hz at %d ms\n', peak_WOE_before_val, peak_WOE_before_sem, round(peak_WOE_before_time));



% ====== AFTER Figure =======
figure('Position', [100, 100, 600, 300]); 
hold on;

% Curves and SEM (subset for plotting)
curve_AD_after = avg_psth_diff_after_AD(plot_idx);
sem_AD_after = sem_psth_diff_after_AD(plot_idx);

curve_WOE_after = avg_psth_diff_after(plot_idx);
sem_WOE_after = sem_psth_diff_after(plot_idx);

% Plot shaded error bars
fill([times_plot, fliplr(times_plot)], ...
     [curve_AD_after + sem_AD_after, fliplr(curve_AD_after - sem_AD_after)], ...
     stimcolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([times_plot, fliplr(times_plot)], ...
     [curve_WOE_after + sem_WOE_after, fliplr(curve_WOE_after - sem_WOE_after)], ...
     washcolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Plot mean curves
plot(times_plot, curve_AD_after, 'Color', stimcolor, 'LineWidth', 2);
plot(times_plot, curve_WOE_after, 'Color', washcolor, 'LineWidth', 2);

% Find peaks in peak window (-50 to 50 ms)
% AD After: max in window
[peak_AD_after_val, peak_AD_after_idx_rel] = max(avg_psth_diff_after_AD(window_idx));
peak_AD_after_idx = find(window_idx, 1, 'first') - 1 + peak_AD_after_idx_rel;
peak_AD_after_time = time_vector(peak_AD_after_idx);
peak_AD_after_sem = sem_psth_diff_after_AD(peak_AD_after_idx);

% WOE After: min in window
[peak_WOE_after_val, peak_WOE_after_idx_rel] = min(avg_psth_diff_after(window_idx));
peak_WOE_after_idx = find(window_idx, 1, 'first') - 1 + peak_WOE_after_idx_rel;
peak_WOE_after_time = time_vector(peak_WOE_after_idx);
peak_WOE_after_sem = sem_psth_diff_after(peak_WOE_after_idx);

% Plot peak markers (only if within plot window)
if peak_AD_after_time >= plot_window(1) && peak_AD_after_time <= plot_window(2)
    plot(peak_AD_after_time, peak_AD_after_val, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    text(peak_AD_after_time + 5, peak_AD_after_val, ...
        sprintf('%.2f ± %.2f Hz', peak_AD_after_val, peak_AD_after_sem), ...
        'Color', stimcolor, 'FontSize', 10);
end

if peak_WOE_after_time >= plot_window(1) && peak_WOE_after_time <= plot_window(2)
    plot(peak_WOE_after_time, peak_WOE_after_val, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    text(peak_WOE_after_time + 5, peak_WOE_after_val, ...
        sprintf('%.2f ± %.2f Hz', peak_WOE_after_val, peak_WOE_after_sem), ...
        'Color', washcolor, 'FontSize', 10);
end

% Formatting
yline(0, '--', 'Color', [0.4 0.4 0.4]);
xlabel('Time from CSM (ms)');
ylabel('Δ Firing Rate (Hz)');
title('Avg ΔFR After-CSM');
xlim(plot_window);
ylim([min([curve_AD_after - sem_AD_after, curve_WOE_after - sem_WOE_after]) - 5, ...
      max([curve_AD_after + sem_AD_after, curve_WOE_after + sem_WOE_after]) + 5]);
legend({'AD SEM', 'WOE SEM', 'AD Mean', 'WOE Mean'}, 'Location', 'best');

hold off;

% Print peaks info AFTER
fprintf('\nAFTER-CSM peaks (within %d to %d ms):\n', time_window(1), time_window(2));
fprintf('AD After: Peak max = %.2f ± %.2f Hz at %d ms\n', peak_AD_after_val, peak_AD_after_sem, round(peak_AD_after_time));
fprintf('WOE After: Peak min = %.2f ± %.2f Hz at %d ms\n', peak_WOE_after_val, peak_WOE_after_sem, round(peak_WOE_after_time));


%%
% Define plot window
plot_window = [-250, 250];
plot_idx = time_vector >= plot_window(1) & time_vector <= plot_window(2);

% Define time window to search for inflection points
time_window = [-150, 150];
window_idx = time_vector >= time_window(1) & time_vector <= time_window(2);
time_win_idx = find(window_idx);
time_win_vec = time_vector(time_win_idx);

% Init structure
inflections = struct();

% ====== BEFORE Figure =======
figure('Position', [100, 100, 600, 300]); 
hold on;

% Extract data
curve_AD_before = avg_psth_diff_before_AD(plot_idx);
sem_AD_before = sem_psth_diff_before_AD(plot_idx);
curve_WOE_before = avg_psth_diff_before(plot_idx);
sem_WOE_before = sem_psth_diff_before(plot_idx);
times_plot = time_vector(plot_idx);

% Plot shaded error
fill([times_plot, fliplr(times_plot)], ...
     [curve_AD_before + sem_AD_before, fliplr(curve_AD_before - sem_AD_before)], ...
     stimcolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([times_plot, fliplr(times_plot)], ...
     [curve_WOE_before + sem_WOE_before, fliplr(curve_WOE_before - sem_WOE_before)], ...
     washcolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Plot curves
plot(times_plot, curve_AD_before, 'Color', stimcolor, 'LineWidth', 2);
plot(times_plot, curve_WOE_before, 'Color', washcolor, 'LineWidth', 2);

% -------- AD BEFORE INFLECTIONS --------
curve_AD_win = avg_psth_diff_before_AD(time_win_idx);
sem_AD_win = sem_psth_diff_before_AD(time_win_idx);
[pks, locs_pks] = findpeaks(curve_AD_win);
[vals_neg, locs_vals] = findpeaks(-curve_AD_win); vals = -vals_neg;
inflections.AD_before.peaks = [];
inflections.AD_before.valleys = [];

for i = 1:length(pks)
    t = time_win_vec(locs_pks(i)); y = pks(i); s = sem_AD_win(locs_pks(i));
    plot(t, y, 'k^', 'MarkerFaceColor', stimcolor, 'MarkerSize', 6);
    text(t+5, y, sprintf('%.2f ± %.2f', y, s), 'Color', stimcolor, 'FontSize', 9);
    inflections.AD_before.peaks(end+1,:) = [t, y, s];
end
for i = 1:length(vals)
    t = time_win_vec(locs_vals(i)); y = vals(i); s = sem_AD_win(locs_vals(i));
    plot(t, y, 'kv', 'MarkerFaceColor', stimcolor, 'MarkerSize', 6);
    text(t+5, y, sprintf('%.2f ± %.2f', y, s), 'Color', stimcolor, 'FontSize', 9);
    inflections.AD_before.valleys(end+1,:) = [t, y, s];
end

% -------- WOE BEFORE INFLECTIONS --------
curve_WOE_win = avg_psth_diff_before(time_win_idx);
sem_WOE_win = sem_psth_diff_before(time_win_idx);
[pks, locs_pks] = findpeaks(curve_WOE_win);
[vals_neg, locs_vals] = findpeaks(-curve_WOE_win); vals = -vals_neg;
inflections.WOE_before.peaks = [];
inflections.WOE_before.valleys = [];

for i = 1:length(pks)
    t = time_win_vec(locs_pks(i)); y = pks(i); s = sem_WOE_win(locs_pks(i));
    plot(t, y, 'k^', 'MarkerFaceColor', washcolor, 'MarkerSize', 6);
    text(t+5, y, sprintf('%.2f ± %.2f', y, s), 'Color', washcolor, 'FontSize', 9);
    inflections.WOE_before.peaks(end+1,:) = [t, y, s];
end
for i = 1:length(vals)
    t = time_win_vec(locs_vals(i)); y = vals(i); s = sem_WOE_win(locs_vals(i));
    plot(t, y, 'kv', 'MarkerFaceColor', washcolor, 'MarkerSize', 6);
    text(t+5, y, sprintf('%.2f ± %.2f', y, s), 'Color', washcolor, 'FontSize', 9);
    inflections.WOE_before.valleys(end+1,:) = [t, y, s];
end

% Format plot
yline(0, '--', 'Color', [0.4 0.4 0.4]);
xline(0,'--')
%xlabel('Time from CSM (ms)');
xlabel('Time from Stim (ms)');
ylabel('Δ Firing Rate (Hz)');
title('Avg ΔFR Before-CSM');
xlim(plot_window);
ylim([min([curve_AD_before - sem_AD_before, curve_WOE_before - sem_WOE_before]) - 5, ...
      max([curve_AD_before + sem_AD_before, curve_WOE_before + sem_WOE_before]) + 5]);
legend({'AD SEM', 'WOE SEM', 'AD Mean', 'WOE Mean'}, 'Location', 'best');
hold off;

% ====== AFTER Figure =======
figure('Position', [100, 100, 600, 300]); 
hold on;

curve_AD_after = avg_psth_diff_after_AD(plot_idx);
sem_AD_after = sem_psth_diff_after_AD(plot_idx);
curve_WOE_after = avg_psth_diff_after(plot_idx);
sem_WOE_after = sem_psth_diff_after(plot_idx);

% Plot shaded error
fill([times_plot, fliplr(times_plot)], ...
     [curve_AD_after + sem_AD_after, fliplr(curve_AD_after - sem_AD_after)], ...
     stimcolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([times_plot, fliplr(times_plot)], ...
     [curve_WOE_after + sem_WOE_after, fliplr(curve_WOE_after - sem_WOE_after)], ...
     washcolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Plot curves
plot(times_plot, curve_AD_after, 'Color', stimcolor, 'LineWidth', 2);
plot(times_plot, curve_WOE_after, 'Color', washcolor, 'LineWidth', 2);

% -------- AD AFTER INFLECTIONS --------
curve_AD_win = avg_psth_diff_after_AD(time_win_idx);
sem_AD_win = sem_psth_diff_after_AD(time_win_idx);
[pks, locs_pks] = findpeaks(curve_AD_win);
[vals_neg, locs_vals] = findpeaks(-curve_AD_win); vals = -vals_neg;
inflections.AD_after.peaks = [];
inflections.AD_after.valleys = [];

for i = 1:length(pks)
    t = time_win_vec(locs_pks(i)); y = pks(i); s = sem_AD_win(locs_pks(i));
    plot(t, y, 'k^', 'MarkerFaceColor', stimcolor, 'MarkerSize', 6);
    text(t+5, y, sprintf('%.2f ± %.2f', y, s), 'Color', stimcolor, 'FontSize', 9);
    inflections.AD_after.peaks(end+1,:) = [t, y, s];
end
for i = 1:length(vals)
    t = time_win_vec(locs_vals(i)); y = vals(i); s = sem_AD_win(locs_vals(i));
    plot(t, y, 'kv', 'MarkerFaceColor', stimcolor, 'MarkerSize', 6);
    text(t+5, y, sprintf('%.2f ± %.2f', y, s), 'Color', stimcolor, 'FontSize', 9);
    inflections.AD_after.valleys(end+1,:) = [t, y, s];
end

% -------- WOE AFTER INFLECTIONS --------
curve_WOE_win = avg_psth_diff_after(time_win_idx);
sem_WOE_win = sem_psth_diff_after(time_win_idx);
[pks, locs_pks] = findpeaks(curve_WOE_win);
[vals_neg, locs_vals] = findpeaks(-curve_WOE_win); vals = -vals_neg;
inflections.WOE_after.peaks = [];
inflections.WOE_after.valleys = [];

for i = 1:length(pks)
    t = time_win_vec(locs_pks(i)); y = pks(i); s = sem_WOE_win(locs_pks(i));
    plot(t, y, 'k^', 'MarkerFaceColor', washcolor, 'MarkerSize', 6);
    text(t+5, y, sprintf('%.2f ± %.2f', y, s), 'Color', washcolor, 'FontSize', 9);
    inflections.WOE_after.peaks(end+1,:) = [t, y, s];
end
for i = 1:length(vals)
    t = time_win_vec(locs_vals(i)); y = vals(i); s = sem_WOE_win(locs_vals(i));
    plot(t, y, 'kv', 'MarkerFaceColor', washcolor, 'MarkerSize', 6);
    text(t+5, y, sprintf('%.2f ± %.2f', y, s), 'Color', washcolor, 'FontSize', 9);
    inflections.WOE_after.valleys(end+1,:) = [t, y, s];
end

% Format plot
yline(0, '--', 'Color', [0.4 0.4 0.4]);
xline(0,'--')
%xlabel('Time from CSM (ms)');
xlabel('Time from Stim (ms)');
ylabel('Δ Firing Rate (Hz)');
title('Avg ΔFR After-CSM');
xlim(plot_window);
ylim([min([curve_AD_after - sem_AD_after, curve_WOE_after - sem_WOE_after]) - 5, ...
      max([curve_AD_after + sem_AD_after, curve_WOE_after + sem_WOE_after]) + 5]);
legend({'AD SEM', 'WOE SEM', 'AD Mean', 'WOE Mean'}, 'Location', 'best');
hold off;



%%
% ====== AFTER Figure =======
figure('Position', [100, 100, 600, 300]); 
hold on;

% Curves and SEM (subset for plotting)
curve_AD_after = avg_psth_diff_after_AD(plot_idx);
sem_AD_after = sem_psth_diff_after_AD(plot_idx);

curve_WOE_after = avg_psth_diff_after(plot_idx);
sem_WOE_after = sem_psth_diff_after(plot_idx);

% Plot shaded error bars
fill([times_plot, fliplr(times_plot)], ...
     [curve_AD_after + sem_AD_after, fliplr(curve_AD_after - sem_AD_after)], ...
     stimcolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([times_plot, fliplr(times_plot)], ...
     [curve_WOE_after + sem_WOE_after, fliplr(curve_WOE_after - sem_WOE_after)], ...
     washcolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Plot mean curves
plot(times_plot, curve_AD_after, 'Color', stimcolor, 'LineWidth', 2);
plot(times_plot, curve_WOE_after, 'Color', washcolor, 'LineWidth', 2);

% Find peaks in window (-50 to 50 ms) relative to full time_vector
[peak_AD_after_val, peak_AD_after_idx_rel] = max(avg_psth_diff_after_AD(window_idx));
peak_AD_after_idx = find(window_idx, 1, 'first') - 1 + peak_AD_after_idx_rel;
peak_AD_after_time = time_vector(peak_AD_after_idx);
peak_AD_after_sem = sem_psth_diff_after_AD(peak_AD_after_idx);

[peak_WOE_after_val, peak_WOE_after_idx_rel] = min(avg_psth_diff_after(window_idx));
peak_WOE_after_idx = find(window_idx, 1, 'first') - 1 + peak_WOE_after_idx_rel;
peak_WOE_after_time = time_vector(peak_WOE_after_idx);
peak_WOE_after_sem = sem_psth_diff_after(peak_WOE_after_idx);

% Plot peak markers (only if within plot window)
if peak_AD_after_time >= plot_window(1) && peak_AD_after_time <= plot_window(2)
    plot(peak_AD_after_time, peak_AD_after_val, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    text(peak_AD_after_time + 5, peak_AD_after_val, ...
        sprintf('%.2f ± %.2f Hz', peak_AD_after_val, peak_AD_after_sem), ...
        'Color', stimcolor, 'FontSize', 10);
end

if peak_WOE_after_time >= plot_window(1) && peak_WOE_after_time <= plot_window(2)
    plot(peak_WOE_after_time, peak_WOE_after_val, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    text(peak_WOE_after_time + 5, peak_WOE_after_val, ...
        sprintf('%.2f ± %.2f Hz', peak_WOE_after_val, peak_WOE_after_sem), ...
        'Color', washcolor, 'FontSize', 10);
end

% Formatting
yline(0, '--', 'Color', [0.4 0.4 0.4]);
xlabel('Time from CSM (ms)');
ylabel('Δ Firing Rate (Hz)');
title('Avg ΔFR After-CSM');
xlim(plot_window);
ylim([min([curve_AD_after - sem_AD_after, curve_WOE_after - sem_WOE_after]) - 5, ...
      max([curve_AD_after + sem_AD_after, curve_WOE_after + sem_WOE_after]) + 5]);
legend({'AD SEM', 'WOE SEM', 'AD Mean', 'WOE Mean'}, 'Location', 'best');

hold off;

% Print peaks info AFTER
fprintf('\nAFTER-CSM peaks (within %d to %d ms):\n', time_window(1), time_window(2));
fprintf('AD After: Peak max = %.2f ± %.2f Hz at %d ms\n', peak_AD_after_val, peak_AD_after_sem, round(peak_AD_after_time));
fprintf('WOE After: Peak min = %.2f ± %.2f Hz at %d ms\n', peak_WOE_after_val, peak_WOE_after_sem, round(peak_WOE_after_time));



%%
% Define time vector and target window
time_window = [-50, 50];
%time_window = [0, 50];
window_idx = time_vector >= time_window(1) & time_vector <= time_window(2);
time_win = time_vector(window_idx);

% ====== BEFORE CONDITION ======
% 1. Find peak in AD-BL within -50 to 50
[peak_val_before_AD, local_peak_idx_before_AD] = max(avg_psth_diff_before_AD(window_idx));

% 2. Get absolute time index for that peak
global_peak_idx_before_AD = find(window_idx);  % indices of interest in full vector
peak_time_idx_before = global_peak_idx_before_AD(local_peak_idx_before_AD);
peak_time_before = time_vector(peak_time_idx_before);

% 3. Extract WOE-BL ΔFR at that time
val_WOE_before = avg_psth_diff_before(peak_time_idx_before);

% 4. Store values
before_vals = [peak_val_before_AD, val_WOE_before];

% ====== AFTER CONDITION ======
[peak_val_after_AD, local_peak_idx_after_AD] = max(avg_psth_diff_after_AD(window_idx));
global_peak_idx_after_AD = find(window_idx);
peak_time_idx_after = global_peak_idx_after_AD(local_peak_idx_after_AD);
peak_time_after = time_vector(peak_time_idx_after);

val_WOE_after = avg_psth_diff_after(peak_time_idx_after);
after_vals = [peak_val_after_AD, val_WOE_after];

% ====== PLOT COMPARISON ======
figure('Position', [100, 100, 500, 300]);
hold on;

% Bar plot or dot plot (dot plot here for clarity)
x = [1, 2];  % 1 = AD-BL, 2 = WOE-BL
bar_width = 0.35;

% Before
b1 = bar(x - 0.2, before_vals, bar_width, 'FaceColor', [0.7 0.7 0.9]);
% After
b2 = bar(x + 0.2, after_vals, bar_width, 'FaceColor', [0.9 0.7 0.9]);

% Labels and aesthetics
xticks(x);
xticklabels({'Peak AD-BL', 'WOE-BL at AD peak'});
ylabel('Δ Firing Rate (Hz)');
title('ΔFR at AD-BL peak time');
legend({'Before', 'After'}, 'Location', 'northeast');

% Annotate time of peak
text(x(1)-0.2, before_vals(1)+1, sprintf('Peak at %d ms', round(peak_time_before)), 'HorizontalAlignment', 'center', 'FontSize', 9);
text(x(1)+0.2, after_vals(1)+1, sprintf('Peak at %d ms', round(peak_time_after)), 'HorizontalAlignment', 'center', 'FontSize', 9);

ylim auto;
box off;


%%
% Define time vector and peak window
time_window = [-50, 50];
%time_window = [0, 50];
window_idx = time_vector >= time_window(1) & time_vector <= time_window(2);
time_bins_in_window = find(window_idx);  % absolute indices in time_vector

% ====== Compute differences ======
diff_AD_before = Psth_AD_CSM_before - Psth_BL_CSM_before;
diff_WOE_before = Psth_WOE_CSM_before - Psth_BL_CSM_before;

diff_AD_after = Psth_AD_CSM_after - Psth_BL_CSM_after;
diff_WOE_after = Psth_WOE_CSM_after - Psth_BL_CSM_after;

% Use minimum unit count across all datasets
n_units = min([ ...
    size(diff_AD_before,1), ...
    size(diff_WOE_before,1), ...
    size(diff_AD_after,1), ...
    size(diff_WOE_after,1)]);

% ====== BEFORE Analysis ======
AD_peaks_before = zeros(n_units, 1);
WOE_at_AD_peak_before = zeros(n_units, 1);
peak_time_idx_before = zeros(n_units, 1);

for i = 1:n_units
    % Find AD peak index (within window)
    [~, local_peak_idx] = max(diff_AD_before(i, window_idx));
    global_idx = time_bins_in_window(local_peak_idx);

    % Record values
    AD_peaks_before(i) = diff_AD_before(i, global_idx);
    WOE_at_AD_peak_before(i) = diff_WOE_before(i, global_idx);
    peak_time_idx_before(i) = global_idx;
end

% ====== AFTER Analysis ======
AD_peaks_after = zeros(n_units, 1);
WOE_at_AD_peak_after = zeros(n_units, 1);
peak_time_idx_after = zeros(n_units, 1);

for i = 1:n_units
    [~, local_peak_idx] = max(diff_AD_after(i, window_idx));
    global_idx = time_bins_in_window(local_peak_idx);

    AD_peaks_after(i) = diff_AD_after(i, global_idx);
    WOE_at_AD_peak_after(i) = diff_WOE_after(i, global_idx);
    peak_time_idx_after(i) = global_idx;
end

% ====== Plot Paired Comparisons ======
figure('Position', [100, 100, 700, 300]);

subplot(1,2,1); hold on;
for i = 1:n_units
    plot([1 2], [AD_peaks_before(i), WOE_at_AD_peak_before(i)], '-o', ...
         'Color', [0.5 0.5 0.5], 'MarkerFaceColor', 'k');
end
xlim([0.5 2.5]);
xticks([1 2]);
xticklabels({'AD-BL', 'WOE-BL'});
ylabel('ΔFR at AD peak time (Hz)');
title('Before-CSM');
box off;

subplot(1,2,2); hold on;
for i = 1:n_units
    plot([1 2], [AD_peaks_after(i), WOE_at_AD_peak_after(i)], '-o', ...
         'Color', [0.5 0.5 0.5], 'MarkerFaceColor', 'k');
end
xlim([0.5 2.5]);
xticks([1 2]);
xticklabels({'AD-BL', 'WOE-BL'});
ylabel('ΔFR at AD peak time (Hz)');
title('After-CSM');
box off;

% ====== Optional Statistical Tests ======
[p_before, ~] = signrank(AD_peaks_before, WOE_at_AD_peak_before);
[p_after, ~] = signrank(AD_peaks_after, WOE_at_AD_peak_after);

fprintf('\n=== Comparison of ΔFR at AD peak time ===\n');
fprintf('Before-CSM:  p = %.4f (WOE vs AD)\n', p_before);
fprintf('After-CSM:   p = %.4f (WOE vs AD)\n', p_after);

%%
% Time window for peak search
time_window = [-50, 50];
%time_window = [0, 50];
window_idx = time_vector >= time_window(1) & time_vector <= time_window(2);
time_bins_in_window = find(window_idx);

% % ====== BEFORE Analysis ======
% diff_AD_before = Psth_AD_CSM_before - Psth_BL_CSM_before;
% diff_WOE_before = Psth_WOE_CSM_before - Psth_BL_CSM_before;
% 
% n_units_before = size(diff_AD_before, 1);
% 
% AD_peaks_before = zeros(n_units_before, 1);
% WOE_at_AD_peak_before = zeros(n_units_before, 1);
% peak_time_idx_before = zeros(n_units_before, 1);
% 
% for i = 1:n_units_before
%     [~, local_peak_idx] = max(diff_AD_before(i, window_idx));
%     global_idx = time_bins_in_window(local_peak_idx);
%     AD_peaks_before(i) = diff_AD_before(i, global_idx);
%     WOE_at_AD_peak_before(i) = diff_WOE_before(i, global_idx);
%     peak_time_idx_before(i) = global_idx;
% end
% 
% % ====== AFTER Analysis ======
% diff_AD_after = Psth_AD_CSM_after - Psth_BL_CSM_after;
% diff_WOE_after = Psth_WOE_CSM_after - Psth_BL_CSM_after;
% 
% n_units_after = size(diff_AD_after, 1);
% 
% AD_peaks_after = zeros(n_units_after, 1);
% WOE_at_AD_peak_after = zeros(n_units_after, 1);
% peak_time_idx_after = zeros(n_units_after, 1);
% 
% for i = 1:n_units_after
%     [~, local_peak_idx] = max(diff_AD_after(i, window_idx));
%     global_idx = time_bins_in_window(local_peak_idx);
%     AD_peaks_after(i) = diff_AD_after(i, global_idx);
%     WOE_at_AD_peak_after(i) = diff_WOE_after(i, global_idx);
%     peak_time_idx_after(i) = global_idx;
% end
% ====== BEFORE Analysis ======
diff_AD_before = Psth_AD_CSM_before - Psth_BL_CSM_before;
diff_WOE_before = Psth_WOE_CSM_before - Psth_BL_CSM_before;

n_units_before = size(diff_AD_before, 1);

AD_peaks_before = zeros(n_units_before, 1);
WOE_at_AD_peak_before = zeros(n_units_before, 1);
peak_time_idx_before = zeros(n_units_before, 1);
peak_time_AD_before_ms = zeros(n_units_before, 1);  % NEW

for i = 1:n_units_before
    [~, local_peak_idx] = max(diff_AD_before(i, window_idx));
    global_idx = time_bins_in_window(local_peak_idx);
    AD_peaks_before(i) = diff_AD_before(i, global_idx);
    WOE_at_AD_peak_before(i) = diff_WOE_before(i, global_idx);
    peak_time_idx_before(i) = global_idx;
    peak_time_AD_before_ms(i) = time_vector(global_idx);  % NEW: Save actual peak time (ms)
end

% ====== AFTER Analysis ======
diff_AD_after = Psth_AD_CSM_after - Psth_BL_CSM_after;
diff_WOE_after = Psth_WOE_CSM_after - Psth_BL_CSM_after;

n_units_after = size(diff_AD_after, 1);

AD_peaks_after = zeros(n_units_after, 1);
WOE_at_AD_peak_after = zeros(n_units_after, 1);
peak_time_idx_after = zeros(n_units_after, 1);
peak_time_AD_after_ms = zeros(n_units_after, 1);  % NEW

for i = 1:n_units_after
    [~, local_peak_idx] = max(diff_AD_after(i, window_idx));
    global_idx = time_bins_in_window(local_peak_idx);
    AD_peaks_after(i) = diff_AD_after(i, global_idx);
    WOE_at_AD_peak_after(i) = diff_WOE_after(i, global_idx);
    peak_time_idx_after(i) = global_idx;
    peak_time_AD_after_ms(i) = time_vector(global_idx);  % NEW: Save actual peak time (ms)
end


% ====== Stats ======
[p_before, ~] = signrank(AD_peaks_before, WOE_at_AD_peak_before);
[p_after, ~] = signrank(AD_peaks_after, WOE_at_AD_peak_after);

mean_vals_before = [mean(AD_peaks_before), mean(WOE_at_AD_peak_before)];
sem_vals_before  = [std(AD_peaks_before) / sqrt(n_units_before), std(WOE_at_AD_peak_before) / sqrt(n_units_before)];

mean_vals_after = [mean(AD_peaks_after), mean(WOE_at_AD_peak_after)];
sem_vals_after  = [std(AD_peaks_after) / sqrt(n_units_after), std(WOE_at_AD_peak_after) / sqrt(n_units_after)];

cellcolor=rgb('Gray');
%figure('Position', [100, 100, 750, 300]);
figure('Position', [100, 100, 500, 300]);

% ------------ BEFORE --------------
subplot(1,2,1); hold on;

% Paired lines
for i = 1:n_units_before
    plot([1 2], [AD_peaks_before(i), WOE_at_AD_peak_before(i)], '-o', ...
         'MarkerSize',3,'Color', [0.5 0.5 0.5], 'MarkerFaceColor', cellcolor);
end

% Error bars
errorbar([1, 2], mean_vals_before, sem_vals_before, 'k-', ...
    'LineWidth', 2, 'CapSize', 6);

% Large mean markers
plot([1 2], mean_vals_before, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');

xlim([0.5 2.5]);
xticks([1 2]);
yline(0,'--');
xticklabels({'AD-BL', 'WOE-BL'});
ylabel('ΔFR at AD peak time (Hz)');
title('Early Cells -CSM Aligned');
box off;

text(1.5, max([AD_peaks_before; WOE_at_AD_peak_before]) + 1, ...
    sprintf('p = %.4f\nn = %d', p_before, n_units_before), ...
    'HorizontalAlignment', 'center', 'FontSize', 10);

% ------------ AFTER --------------
subplot(1,2,2); hold on;

% Paired lines
for i = 1:n_units_after
    plot([1 2], [AD_peaks_after(i), WOE_at_AD_peak_after(i)], '-o', ...
         'MarkerSize',3,'Color', [0.5 0.5 0.5], 'MarkerFaceColor', cellcolor);
end

% Error bars
errorbar([1, 2], mean_vals_after, sem_vals_after, 'k-', ...
    'LineWidth', 2, 'CapSize', 6);

% Large mean markers
plot([1 2], mean_vals_after, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');

xlim([0.5 2.5]);
xticks([1 2]);
yline(0,'--')
xticklabels({'AD-BL', 'WOE-BL'});
ylabel('ΔFR at AD peak time (Hz)');
title('Late Cells -CSM Aligned');
box off;

text(1.5, max([AD_peaks_after; WOE_at_AD_peak_after]) + 1, ...
    sprintf('p = %.4f\nn = %d', p_after, n_units_after), ...
    'HorizontalAlignment', 'center', 'FontSize', 10);

sgtitle('ΔFR at AD Peak Time (Per-Unit Analysis)');

% ====== Summary of Results ======
fprintf('\n===== SUMMARY: ΔFR at AD Peak Time =====\n');

% BEFORE
fprintf('BEFORE (Early Cells)\n');
fprintf('AD-BL:  Mean = %.2f Hz, SEM = %.2f Hz\n', mean_vals_before(1), sem_vals_before(1));
fprintf('WOE-BL: Mean = %.2f Hz, SEM = %.2f Hz\n', mean_vals_before(2), sem_vals_before(2));
fprintf('p = %.4f (Wilcoxon signed-rank test)\n\n', p_before);

% AFTER
fprintf('AFTER (Late Cells)\n');
fprintf('AD-BL:  Mean = %.2f Hz, SEM = %.2f Hz\n', mean_vals_after(1), sem_vals_after(1));
fprintf('WOE-BL: Mean = %.2f Hz, SEM = %.2f Hz\n', mean_vals_after(2), sem_vals_after(2));
fprintf('p = %.4f (Wilcoxon signed-rank test)\n', p_after);

%% 
[f_ctrl, xi_ctrl] = ksdensity(peak_time_AD_before_ms,'Bandwidth',10);

[f_ctrlfirst, xi_ctrlfirst] = ksdensity(peak_time_AD_after_ms,'Bandwidth',10);


% -----------------------------
% Plot the KDEs with Full Area Under the Curve Filled
% -----------------------------
figure; hold on;
plot(xi_ctrl, mean(f_ctrl, 1), 'k', 'LineWidth', 3); % Ctrlend (black)
plot(xi_ctrlfirst, mean(f_ctrlfirst, 1), 'b', 'LineWidth', 3); % Ctrlfirst (blue)
xlim([-60 60])
%% 
% Select data
x1 = peak_time_AD_before_ms;
x2 = peak_time_AD_after_ms;

% Use Sturges' Rule for bin number
nbins1 = ceil(log2(length(x1)) + 1);
nbins2 = ceil(log2(length(x2)) + 1);
nbins = max([nbins1, nbins2]);

% Determine common bin edges across both datasets
all_data = [x1; x2];
bin_edges = linspace(min(all_data), max(all_data), nbins+1);

% Plot overlaid histograms
figure('Position', [400, 300, 600, 350]); hold on;

histogram(x1, bin_edges, 'Normalization', 'probability', ...
    'FaceColor', rgb('DodgerBlue'), 'FaceAlpha', 0.5, 'EdgeColor', 'k');

histogram(x2, bin_edges, 'Normalization', 'probability', ...
    'FaceColor', rgb('DarkViolet'), 'FaceAlpha', 0.5, 'EdgeColor', 'k');

% Labels and formatting
xlabel('Time of Peak ΔFR (ms)');
ylabel('Probability');
legend({'Before (Early Cells)', 'After (Late Cells)'}, 'Location', 'NorthEast');
title('Overlaid Histograms of Peak ΔFR Timing (AD-BL Units)');
xlim([-50 50]);
yline(0, '--k');
box off;
%%
% Select data
x1 = peak_time_AD_before_ms;
x2 = peak_time_AD_after_ms;

% ----- Use Scott’s Rule for bin width -----
all_data = [x1; x2];
bin_width = 3.5 * std(all_data) / (length(all_data)^(1/3));
nbins = round(range(all_data) / bin_width);
bin_edges = linspace(min(all_data), max(all_data), nbins + 1);

% ----- Plot overlaid histograms -----
figure('Position', [400, 300, 600, 350]); hold on;

h1 = histogram(x1, bin_edges, 'Normalization', 'probability', ...
    'FaceColor', rgb('DodgerBlue'), 'FaceAlpha', 0.5, 'EdgeColor', 'k');

h2 = histogram(x2, bin_edges, 'Normalization', 'probability', ...
    'FaceColor', rgb('DarkViolet'), 'FaceAlpha', 0.5, 'EdgeColor', 'k');

% ----- Plot mean lines -----
mean_x1 = mean(x1);
mean_x2 = mean(x2);

plot([mean_x1 mean_x1], [0, max([h1.Values h2.Values])], '--', ...
    'Color', rgb('DodgerBlue'), 'LineWidth', 1.5);

plot([mean_x2 mean_x2], [0, max([h1.Values h2.Values])], '--', ...
    'Color', rgb('DarkViolet'), 'LineWidth', 1.5);

% ----- Labels and formatting -----
xlabel('Time of Peak ΔFR (ms)');
ylabel('Probability');
legend({'Before (Early Cells)', 'After (Late Cells)'}, 'Location', 'NorthEast');
title('Overlaid Histograms of Peak ΔFR Timing (AD-BL Units)');
xlim([-50 50]);
yline(0, '--k');
box off;

%% Compute AUC (area under curve) of ΔFR in 0–50 ms window for each unit
adcolor=rgb('DodgerBlue');
woecolor=rgb('Violet');
diff_before_WOE = Psth_WOE_CSM_before - Psth_BL_CSM_before;
diff_after_WOE  = Psth_WOE_CSM_after  - Psth_BL_CSM_after;

diff_before = Psth_AD_CSM_before - Psth_BL_CSM_before;
diff_after  = Psth_AD_CSM_after  - Psth_BL_CSM_after;

% Use trapz over 0-50 ms window (time_vector is in ms)
auc_diff_before_WOE = trapz(time_vector(time_window_idx), diff_before_WOE(:, time_window_idx), 2);
auc_diff_after_WOE  = trapz(time_vector(time_window_idx), diff_after_WOE(:, time_window_idx), 2);

auc_diff_before = trapz(time_vector(time_window_idx), diff_before(:, time_window_idx), 2);
auc_diff_after  = trapz(time_vector(time_window_idx), diff_after(:, time_window_idx), 2);

% Convert from Hz·ms to Hz·s by dividing by 1000
auc_diff_before_WOE = auc_diff_before_WOE / 1000;
auc_diff_after_WOE  = auc_diff_after_WOE / 1000;
auc_diff_before     = auc_diff_before / 1000;
auc_diff_after      = auc_diff_after / 1000;

% Combine data for plotting
auc_before_all = [auc_diff_before; auc_diff_before_WOE];
auc_after_all  = [auc_diff_after;  auc_diff_after_WOE];

% Create group labels: 1=AD, 2=WOE
group_before = [ones(size(auc_diff_before)); 2 * ones(size(auc_diff_before_WOE))];
group_after  = [ones(size(auc_diff_after));  2 * ones(size(auc_diff_after_WOE))];

% Mean and SEM per group for before
mean_auc_before = [nanmean(auc_diff_before), nanmean(auc_diff_before_WOE)];
sem_auc_before  = [nanstd(auc_diff_before)/sqrt(sum(~isnan(auc_diff_before))), ...
                   nanstd(auc_diff_before_WOE)/sqrt(sum(~isnan(auc_diff_before_WOE)))];

% Mean and SEM per group for after
mean_auc_after = [nanmean(auc_diff_after), nanmean(auc_diff_after_WOE)];
sem_auc_after  = [nanstd(auc_diff_after)/sqrt(sum(~isnan(auc_diff_after))), ...
                  nanstd(auc_diff_after_WOE)/sqrt(sum(~isnan(auc_diff_after_WOE)))];

% Plot AUC for BEFORE comparison
figure('Position', [100, 200, 250, 300]); hold on;

x_vals_before = [1, 2];
colors = {[0.85 0.3 0.3], [0.3 0.3 0.85]};

for i = 1:2
    scatter(repmat(x_vals_before(i), sum(group_before==i), 1), ...
        auc_before_all(group_before==i), 25, 'filled', ...
        'MarkerFaceColor', colors{i}, 'MarkerFaceAlpha', 0.5);
end

errorbar(x_vals_before, mean_auc_before, sem_auc_before, 'k', 'LineWidth', 1.5, 'CapSize', 10);
scatter(x_vals_before, mean_auc_before, 100, 'k', 'filled');
xlim([0.5 2.5]);
xticks(x_vals_before); xticklabels({'AD', 'WOE'});
ylabel('AUC ΔFR (Hz·s) in 0–50 ms');
title('AUC Comparison BEFORE');
yline(0, '--');
box off;

% Paired t-test BEFORE
[h_before, p_before, ci_before, stats_before] = ttest(auc_diff_before, auc_diff_before_WOE);
fprintf('Paired t-test BEFORE: p = %.4f, tstat = %.2f, df = %d\n', p_before, stats_before.tstat, stats_before.df);

% Add p-value text and asterisk to BEFORE plot
ymax_before = max(auc_before_all) + 0.1 * range(auc_before_all);
if p_before < 0.001
    sigtext_before = '***';
elseif p_before < 0.01
    sigtext_before = '**';
elseif p_before < 0.05
    sigtext_before = '*';
else
    sigtext_before = 'n.s.';
end
text(1.5, ymax_before, sprintf('p = %.3f\n%s', p_before, sigtext_before), ...
    'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');

% Plot AUC for AFTER comparison

figure('Position', [100, 200, 250, 300]); hold on;
x_vals_after = [1, 2];

for i = 1:2
    scatter(repmat(x_vals_after(i), sum(group_after==i), 1), ...
        auc_after_all(group_after==i), 25, 'filled', ...
        'MarkerFaceColor', colors{i}, 'MarkerFaceAlpha', 0.5);
end

errorbar(x_vals_after, mean_auc_after, sem_auc_after, 'k', 'LineWidth', 1.5, 'CapSize', 10);
scatter(x_vals_after, mean_auc_after, 100, 'k', 'filled');
xlim([0.5 2.5]);
xticks(x_vals_after); xticklabels({'AD', 'WOE'});
ylabel('AUC ΔFR (Hz·s) in 0–50 ms');
title('AUC Comparison AFTER');
yline(0, '--');
box off;

% Paired t-test AFTER
[h_after, p_after, ci_after, stats_after] = ttest(auc_diff_after, auc_diff_after_WOE);
fprintf('Paired t-test AFTER: p = %.4f, tstat = %.2f, df = %d\n', p_after, stats_after.tstat, stats_after.df);

% Add p-value text and asterisk to AFTER plot
ymax_after = max(auc_after_all) + 0.1 * range(auc_after_all);
if p_after < 0.001
    sigtext_after = '***';
elseif p_after < 0.01
    sigtext_after = '**';
elseif p_after < 0.05
    sigtext_after = '*';
else
    sigtext_after = 'n.s.';
end
text(1.5, ymax_after, sprintf('p = %.3f\n%s', p_after, sigtext_after), ...
    'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');

%%
% Define colors
adcolor = rgb('DodgerBlue');
woecolor = rgb('Violet');
graycolor = [0.5 0.5 0.5];
time_window_idx = find(time_vector >= 0 & time_vector < 50);

% Compute AUC (area under curve) of ΔFR in 0–50 ms window for each unit
diff_before_WOE = Psth_WOE_CSM_before - Psth_BL_CSM_before;
diff_after_WOE  = Psth_WOE_CSM_after  - Psth_BL_CSM_after;

diff_before = Psth_AD_CSM_before - Psth_BL_CSM_before;
diff_after  = Psth_AD_CSM_after  - Psth_BL_CSM_after;

% Use trapz over 0-50 ms window (time_vector is in ms)
auc_diff_before_WOE = trapz(time_vector(time_window_idx), diff_before_WOE(:, time_window_idx), 2);
auc_diff_after_WOE  = trapz(time_vector(time_window_idx), diff_after_WOE(:, time_window_idx), 2);

auc_diff_before = trapz(time_vector(time_window_idx), diff_before(:, time_window_idx), 2);
auc_diff_after  = trapz(time_vector(time_window_idx), diff_after(:, time_window_idx), 2);

% Convert from Hz·ms to Hz·s by dividing by 1000
auc_diff_before_WOE = auc_diff_before_WOE / 1000;
auc_diff_after_WOE  = auc_diff_after_WOE / 1000;
auc_diff_before     = auc_diff_before / 1000;
auc_diff_after      = auc_diff_after / 1000;

% Combine data for plotting
auc_before_all = [auc_diff_before; auc_diff_before_WOE];
auc_after_all  = [auc_diff_after;  auc_diff_after_WOE];

% Create group labels: 1=AD, 2=WOE
group_before = [ones(size(auc_diff_before)); 2 * ones(size(auc_diff_before_WOE))];
group_after  = [ones(size(auc_diff_after));  2 * ones(size(auc_diff_after_WOE))];

% Mean and SEM per group for before
mean_auc_before = [nanmean(auc_diff_before), nanmean(auc_diff_before_WOE)];
sem_auc_before  = [nanstd(auc_diff_before)/sqrt(sum(~isnan(auc_diff_before))), ...
                   nanstd(auc_diff_before_WOE)/sqrt(sum(~isnan(auc_diff_before_WOE)))];

% Mean and SEM per group for after
mean_auc_after = [nanmean(auc_diff_after), nanmean(auc_diff_after_WOE)];
sem_auc_after  = [nanstd(auc_diff_after)/sqrt(sum(~isnan(auc_diff_after))), ...
                  nanstd(auc_diff_after_WOE)/sqrt(sum(~isnan(auc_diff_after_WOE)))];

% Plot AUC for BEFORE comparison
figure('Position', [100, 200, 250, 300]); hold on;
x_vals_before = [1, 2];

% Plot gray scatter points
scatter(repmat(x_vals_before(1), sum(group_before==1), 1), ...
    auc_diff_before, 25, 'filled', ...
    'MarkerFaceColor', graycolor, 'MarkerFaceAlpha', 0.5);

scatter(repmat(x_vals_before(2), sum(group_before==2), 1), ...
    auc_diff_before_WOE, 25, 'filled', ...
    'MarkerFaceColor', graycolor, 'MarkerFaceAlpha', 0.5);

% Plot mean ± SEM and large colored dots
errorbar(x_vals_before, mean_auc_before, sem_auc_before, ...
    'k', 'LineWidth', 1.5, 'CapSize', 10);
scatter(x_vals_before(1), mean_auc_before(1), 100, 'filled', 'MarkerFaceColor', adcolor);
scatter(x_vals_before(2), mean_auc_before(2), 100, 'filled', 'MarkerFaceColor', woecolor);

xlim([0.5 2.5]);
xticks(x_vals_before); xticklabels({'AD', 'WOE'});
ylabel('AUC ΔFR (Hz·s) in 0–50 ms');
title('AUC Comparison BEFORE');
yline(0, '--');
box off;

% Paired t-test BEFORE
[h_before, p_before, ci_before, stats_before] = ttest(auc_diff_before, auc_diff_before_WOE);
fprintf('Paired t-test BEFORE: p = %.4f, tstat = %.2f, df = %d\n', ...
    p_before, stats_before.tstat, stats_before.df);

% Add p-value text and asterisk to BEFORE plot
ymax_before = max(auc_before_all) + 0.1 * range(auc_before_all);
if p_before < 0.001
    sigtext_before = '***';
elseif p_before < 0.01
    sigtext_before = '**';
elseif p_before < 0.05
    sigtext_before = '*';
else
    sigtext_before = 'n.s.';
end
text(1.5, ymax_before, sprintf('p = %.3f\n%s', p_before, sigtext_before), ...
    'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');

% Plot AUC for AFTER comparison
figure('Position', [400, 200, 250, 300]); hold on;
x_vals_after = [1, 2];

% Plot gray scatter points
scatter(repmat(x_vals_after(1), sum(group_after==1), 1), ...
    auc_diff_after, 25, 'filled', ...
    'MarkerFaceColor', graycolor, 'MarkerFaceAlpha', 0.5);

scatter(repmat(x_vals_after(2), sum(group_after==2), 1), ...
    auc_diff_after_WOE, 25, 'filled', ...
    'MarkerFaceColor', graycolor, 'MarkerFaceAlpha', 0.5);

% Plot mean ± SEM and large colored dots
errorbar(x_vals_after, mean_auc_after, sem_auc_after, ...
    'k', 'LineWidth', 1.5, 'CapSize', 10);
scatter(x_vals_after(1), mean_auc_after(1), 100, 'filled', 'MarkerFaceColor', adcolor);
scatter(x_vals_after(2), mean_auc_after(2), 100, 'filled', 'MarkerFaceColor', woecolor);

xlim([0.5 2.5]);
xticks(x_vals_after); xticklabels({'AD', 'WOE'});
ylabel('AUC ΔFR (Hz·s) in 0–50 ms');
title('AUC Comparison AFTER');
yline(0, '--');
box off;

% Paired t-test AFTER
[h_after, p_after, ci_after, stats_after] = ttest(auc_diff_after, auc_diff_after_WOE);
fprintf('Paired t-test AFTER: p = %.4f, tstat = %.2f, df = %d\n', ...
    p_after, stats_after.tstat, stats_after.df);

% Add p-value text and asterisk to AFTER plot
ymax_after = max(auc_after_all) + 0.1 * range(auc_after_all);
if p_after < 0.001
    sigtext_after = '***';
elseif p_after < 0.01
    sigtext_after = '**';
elseif p_after < 0.05
    sigtext_after = '*';
else
    sigtext_after = 'n.s.';
end
text(1.5, ymax_after, sprintf('p = %.3f\n%s', p_after, sigtext_after), ...
    'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');


%% NP data- remove entires with no input - 5 and 7

NPdata_edit = NPdata;
NPdata_edit([5,6,7]) = [];
%% Now look at kinematics
% Collect all NoCSM reaches and CSM
BLAll= struct; 
ADAll=struct;
WOEAll=struct;

CSMAll=struct;
NoCSMAll=struct;

% Combine all reachse across animals 
for i = 1:length(NPdata_edit)
 
    Reach=NPdata_edit(i).ExcitRR_BL.BL;
    BLAll(i).BL=Reach; 

    ReachAD=NPdata_edit(i).ExcitRR_AD.AD;
    ADAll(i).AD=ReachAD; 

    ReachWOE=NPdata_edit(i).ExcitRR_WOE.WO;
    WOEAll(i).WOE=ReachWOE; 
   
    clearvars Reach ReachAD ReachWOE
end

%% collect reaches and get average 
BL_stimvel=struct;
for i=1:length(BLAll)
   
    BLR=BLAll(i).BL;

    if isempty(BLR)
        C=[];
    else
        dim=7; %max outward velocity
        dir=1;
        thresh= 0.1;
       % C= AlignPeakV(NoCSMR,dim,1000);
        %C= AlignPeakV(BLR,dim,250);
        C= AlignThreshCSM_T2HF(BLR,thresh,dir,250);
        %C= AlignThreshCSMHF(BLR,thresh,dir,250);
        %C= AlignThreshCSM_T0HF(NoCSMR,thresh,dir,250);
        C=stimvel_accel(C);
        C=stimvel_accel_jerk(C);
    end
    BL_stimvel(i).stimvel=C;
end
%%
AD_stimvel=struct;
for i=1:length(ADAll)
   
    ADR=ADAll(i).AD;

    if isempty(ADR)
        C=[];
    else
        dim=7; %max outward velocity
        dir=1;
        thresh= 0.1;
       % C= AlignPeakV(NoCSMR,dim,1000);
        %C= AlignPeakV(ADR,dim,250);
        C= AlignThreshCSM_T2HF(ADR,thresh,dir,250);
        %C= AlignThreshCSM_T0HF(NoCSMR,thresh,dir,250);
        %C=  AlignThreshCSMHF(ADR,thresh,dir,250);
        C=stimvel_accel(C);
        C=stimvel_accel_jerk(C);
    end
    AD_stimvel(i).stimvel=C;
end
%%
WOE_stimvel=struct;
for i=1:length(WOEAll)
   
    WOER=WOEAll(i).WOE;

    if isempty(WOER)
        C=[];
    else
        dim=7; %max outward velocity
        dir=1;
        thresh= 0.1;
       % C= AlignPeakV(NoCSMR,dim,1000);
        %C= AlignPeakV(WOER,dim,250);
        C= AlignThreshCSM_T2HF(WOER,thresh,dir,250);
        %C= AlignThreshCSM_T0HF(NoCSMR,thresh,dir,250);
       % C= AlignThreshCSMHF(WOER,thresh,dir,250);
        C=stimvel_accel(C);
        C=stimvel_accel_jerk(C);
    end
    WOE_stimvel(i).stimvel=C;
end
%%
BL_big=struct;

for a = 1:length(BL_stimvel) 
    BL_big.stimvel(:,:,a) =BL_stimvel(a).stimvel.Average;
end
%NoReReach
BL_big.Average = nanmean(BL_big.stimvel,3);
BL_big.Std = nanstd(BL_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(BL_big.Average),1); 
for nnum = 1:length(BL_big.Average)
    samp = squeeze(BL_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
BL_big.N = N;
%%
AD_big=struct;

for a = 1:length(AD_stimvel) 
   AD_big.stimvel(:,:,a) = AD_stimvel(a).stimvel.Average;
end
%NoReReach
AD_big.Average = nanmean(AD_big.stimvel,3);
AD_big.Std = nanstd(AD_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(AD_big.Average),1); 
for nnum = 1:length(AD_big.Average)
    samp = squeeze(AD_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
AD_big.N = N;
%%
WOE_big=struct;

for a = 1:length(WOE_stimvel) 
   WOE_big.stimvel(:,:,a) = WOE_stimvel(a).stimvel.Average;
end
%NoReReach
WOE_big.Average = nanmean(WOE_big.stimvel,3);
WOE_big.Std = nanstd(WOE_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(WOE_big.Average),1); 
for nnum = 1:length(WOE_big.Average)
    samp = squeeze(WOE_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
WOE_big.N = N;
%% Stimvel Plot 
BLcolor=rgb('Black');
ADcolor=rgb('DodgerBlue');
WOcolor=rgb('Magenta');

dim=5;
StimvelPlot(BL_big,dim,BLcolor,1);
StimvelPlot(AD_big,dim,ADcolor,0);
StimvelPlot(WOE_big,dim,WOcolor,0);
yline(0,'--')
ylim('auto')
xlim([-0.25 0.25])
ylabel('Outward Velocity (cm/s)')
xlabel('Time from Peak Velocity (s)')
% f = gcf; % Get current figure handle
% f.Units = 'inches';
% pos = f.Position;
% f.Position = [pos(1:2), pos(3), pos(4)/2]; % Halve the height



%% 
diff_avg = WOE_big.Average(:,5) - BL_big.Average(:,5);
diff_avg=diff_avg*10;
time=BL_big.Average(:,1);

% Standard error of the difference (element-wise)
SE_diff = sqrt((WOE_big.Std.^2 ./ WOE_big.N) + (BL_big.Std.^2 ./ BL_big.N));
%%
figure; hold on
plot(time,diff_avg)
xlim([-0.25 0.25])
yline(0,'--')
xlabel('Time from postional threshold (ms)')
ylabel('Change in Outward Acceleration')


%%
% For example, just plot the y-dimension (2nd column)
diff_avg_y = WOE_big.Average(:,8) - BL_big.Average(:,8);
diff_avg_y=diff_avg_y*10;
SE_diff_y = sqrt((WOE_big.Std(:,8).^2 ./ WOE_big.N) + (BL_big.Std(:,8).^2 ./ BL_big.N));

time = (BL_big.Average(:,1))'; % ensure it's a column vector

% Now plot
figure('Position', [100, 100, 600, 300]);  
shadedErrorBar(time, diff_avg_y, SE_diff_y, 'lineProps', {'k', 'LineWidth', 2});
xlabel('Time (ms)');
xlim([-0.25 0.25])
ylabel('Y-Velocity Difference (CSM - NoCSM)');
title('CSM vs NoCSM Kinematic Difference (Y-Direction)');
yline(0, '--');
%%
% For example, just plot the y-dimension (2nd column)
diff_avg_y = AD_big.Average(:,5) - BL_big.Average(:,5);
diff_avg_y=diff_avg_y*10;
SE_diff_y = sqrt((AD_big.Std(:,5).^2 ./ AD_big.N) + (BL_big.Std(:,5).^2 ./ BL_big.N));

time = (BL_big.Average(:,1))'; % ensure it's a column vector

% Now plot
figure('Position', [100, 100, 600, 300]);  
shadedErrorBar(time, diff_avg_y, SE_diff_y, 'lineProps', {'k', 'LineWidth', 2});
xlabel('Time (ms)');
xlim([-0.25 0.25])
ylabel('Y-Velocity Difference (CSM - NoCSM)');
title('CSM vs NoCSM Kinematic Difference (Y-Direction)');
yline(0, '--');
%% Generate stimvel zscore
% Initialize arrays for p-values and Z-scores
parray0 = NaN(size(BL_big.stimvel, 1), 1);  % Store p-values
zarray0 =   NaN(size(BL_big.stimvel, 1), 1)   ;    % Store Z-scores (same shape as stimvel)
Zscores0= NaN(size(BL_big.stimvel, 1), 5);

for pnum = 1:length(parray0)
    % Extract trial data for C_big (control) and T0_big (experimental)
    %a4b = squeeze(WOE_big.stimvel(pnum, 3, :));  % Experimental data (T0_big)
    a4b = squeeze(AD_big.stimvel(pnum, 3, :));  % Experimental data (T0_big)
    b4b = squeeze(BL_big.stimvel(pnum,3, :));   % Control data (C_big), matching trial

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

dim = 3;
StimvelPlot(BL_big, dim, BLcolor, 1);
StimvelPlot(WOE_big, dim,WOcolor, 0);
StimvelPlot(AD_big, dim,ADcolor, 0);

% Plot reference line
yline(0, '--k');

% Axis limits and labels
ylim('auto')
%ylim([-20 20])
xlabel('Time from MaxV (s)');

ylabel('Outward Acceleration');
title('WOEvsBL')

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

dim = 2;

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


