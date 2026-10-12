%% Get data - COde for CSM FR analysis - bursters and decreases 



cd("C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\ShamAnimals\ReachandCSMMod_allCSMStogther\CSMAll_wtihMaxV")
MyFolderInfo = dir('*CSMAll.mat'); 

tic
clearvars -except MyFolderInfo
%% === INITIALIZE MASTER MATRICES ===================================== %%
NPdata = struct;

PsthMaster_CSM     = [];
PsthMaster_NoCSM   = [];

Psth_CSM_before = []; 
Psth_CSM_after = [];
Psth_NoCSM_before = []; 
Psth_NoCSM_after = [];

% === MAIN LOOP ======================================================= %%
for day = 1:length(MyFolderInfo)
    load(MyFolderInfo(day).name, ...
         'Reaches', 'data', 'filtered_modulated_indices', ...
         'neg_corr_cells','pos_corr_cells','FRBins', ...
         'FRCorrStats_edit', ...
         'ExcitRR');

    if isempty(filtered_modulated_indices)
        fprintf('Day %d: No modulated indices. Skipping...\n', day);
        continue;
    end

    data_opto_all = data.units(filtered_modulated_indices);
    timevec = data.units(filtered_modulated_indices(1)).AllCSM_pos_thresh.time_bins;

    baseline_win = find(timevec >= -1000 & timevec <= -500);
    %baseline_win = find(timevec >= -1000 & timevec <= -200);
    [~, zero_idx] = min(abs(timevec));
    %valid_peak_range = find(timevec >= -250 & timevec <= 250);
    valid_peak_range = find(timevec >= -100 & timevec <= 100);

    keep_unit = false(numel(data_opto_all), 1);
    peak_times = NaN(numel(data_opto_all), 1);

    for t = 1:numel(data_opto_all)
        psth = data_opto_all(t).max_vx.PSTH_smoothbin;
        [~, peak_idx] = max(psth);
        peak_times(t) = peak_idx;
        keep_unit(t) = ismember(peak_idx, valid_peak_range);
    end

    data_opto = data_opto_all(keep_unit);
    peak_times = peak_times(keep_unit);

    is_before = peak_times < zero_idx;
    is_after  = peak_times >= zero_idx;

    if isempty(data_opto)
        fprintf('Day %d: No valid units with peak between -250 and 250 ms. Skipping...\n', day);
        continue
    end

    % Preallocate matrices
    psth_len = numel(data_opto(1).AllCSM_pos_thresh.PSTH_smoothbin);
    csm_mat   = NaN(numel(data_opto), psth_len);
    nocsm_mat = NaN(numel(data_opto), psth_len);

    % Loop over units
    for t = 1:numel(data_opto)
        % Extract PSTHs
        % csm    = data_opto(t).AllCSM_pos_thresh.PSTH_smoothbin;
        % nocsm  = data_opto(t).NoCSM_pos_thresh.PSTH_smoothbin;

        csm    = data_opto(t).AllCSM_csmT2.PSTH_smoothbin;
        nocsm  = data_opto(t).NoCSM_csmT2.PSTH_smoothbin;
        % 
        %  csm    = data_opto(t).AllCSM_csmT0.PSTH_smoothbin;
        % nocsm  = data_opto(t).NoCSM_csmT0.PSTH_smoothbin;

        bl_ref = csm;  % z-score using the CSM PSTH


        % Z-score
        baseline_mean = mean(bl_ref(baseline_win));
        baseline_std = std(bl_ref(baseline_win));
        if baseline_std == 0
            baseline_std = 1;
        end

        csm_z    = (csm - baseline_mean) / baseline_std;
        nocsm_z  = (nocsm - baseline_mean) / baseline_std;

        % Store
        csm_mat(t,:)   = csm_z;
        nocsm_mat(t,:) = nocsm_z;
    end

    % Aggregate
    PsthMaster_CSM   = vertcat(PsthMaster_CSM,   csm_mat);
    PsthMaster_NoCSM = vertcat(PsthMaster_NoCSM, nocsm_mat);

    Psth_CSM_before     = vertcat(Psth_CSM_before,     csm_mat(is_before,:));
    Psth_CSM_after      = vertcat(Psth_CSM_after,      csm_mat(is_after,:));
    Psth_NoCSM_before   = vertcat(Psth_NoCSM_before,   nocsm_mat(is_before,:));
    Psth_NoCSM_after    = vertcat(Psth_NoCSM_after,    nocsm_mat(is_after,:));

    % Store metadata
    NPdata(day).peak_times     = peak_times;
    NPdata(day).kept_indices   = filtered_modulated_indices(keep_unit);
    NPdata(day).is_before      = is_before;
    NPdata(day).is_after       = is_after;
    NPdata(day).ExcitRR        = ExcitRR;
    NPdata(day).cells          = data;

    fprintf('Day %d: Kept %d/%d modulated units with peak between -100 and 100 ms\n', ...
        day, sum(keep_unit), numel(keep_unit));
end
%% Looking for units that are outliers
% Define time window of interest
target_start = 150;
target_end = 250;

% Get indices corresponding to that time window
idx_range = timevec >= target_start & timevec <= target_end;

% Compute per-unit variability (standard deviation) in this window
std_per_unit = std(Psth_NoCSM_before(:, idx_range), 0, 2);  % use whichever condition looks weird

% OR: max absolute difference between adjacent bins
max_jump_per_unit = max(abs(diff(Psth_NoCSM_before(:, idx_range), 1, 2)), [], 2);

% Plot to see which one is the outlier
figure;
subplot(1,2,1)
plot(std_per_unit, 'o');
title('STD per unit in 150–250ms');
xlabel('Unit'); ylabel('STD');

subplot(1,2,2)
plot(max_jump_per_unit, 'o');
title('Max jump per unit in 150–250ms');
xlabel('Unit'); ylabel('Max abs diff');

outlier_std_idx = find(std_per_unit > 2);  % pick threshold based on plot
outlier_jump_idx = find(max_jump_per_unit > 3);  % adjust if needed

% Combine or inspect
bad_units = unique([outlier_std_idx; outlier_jump_idx]);

disp('Problematic units:')
disp(bad_units)

%%
figure;
hold on;
for i = 1:length(bad_units)
    plot(timevec, Psth_NoCSM_before(bad_units(i), :));
end
xlim([-250 250])
title('Problematic Units');
xlabel('Time (ms)');
ylabel('Z-scored Firing Rate');
%% 
Color_BL_CSM      = rgb('DarkMagenta');
Color_BL_NoCSM    = rgb('Black');
% Remove from both conditions to keep them matched
Psth_CSM_before_cleaned = Psth_CSM_before;
Psth_NoCSM_before_cleaned = Psth_NoCSM_before;
%bad_units=[18; 25]; %Badunis identified from pos thresh crossing
%bad_units=[7;18;25];
% bad_units=[3];
bad_units=[2;7;9;10;12;13;14;16;17;18;21;22;24;25];
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
num_bins = size(Psth_CSM_before, 2);
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

% === PLOT: AFTER PEAK UNITS ===
%figure; hold on;

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

%xline(csmmean,'--')
xline(0,'--')
title('PSTH – After Peak Firing');
xlabel('Time from POs Thresh(Stim) (ms)');

ylabel('Firing rate (zscore)');
xlim([-250 250])
ylim([0 18])
box off;



%% 
Color_BL_CSM      = rgb('DarkMagenta');
Color_BL_NoCSM    = rgb('Black');

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


%% Define colors
csmmean=81.5;
Color_BL_CSM      = rgb('SlateBlue');
Color_BL_NoCSM    = rgb('Black');


% Compute average PSTHs and SEMs — BEFORE PEAK
avg_before_BL_CSM    = nanmean(Psth_CSM_before, 1);
avg_before_BL_NoCSM  = nanmean(Psth_NoCSM_before, 1);


sem_before_BL_CSM    = nanstd(Psth_CSM_before, 0, 1) / sqrt(size(Psth_CSM_before, 1));
sem_before_BL_NoCSM  = nanstd(Psth_NoCSM_before, 0, 1) / sqrt(size(Psth_NoCSM_before, 1));


% Compute average PSTHs and SEMs — AFTER PEAK
avg_after_BL_CSM    = nanmean(Psth_CSM_after, 1);
avg_after_BL_NoCSM  = nanmean(Psth_NoCSM_after, 1);


sem_after_BL_CSM    = nanstd(Psth_CSM_after, 0, 1) / sqrt(size(Psth_CSM_after, 1));
sem_after_BL_NoCSM  = nanstd(Psth_NoCSM_after, 0, 1) / sqrt(size(Psth_NoCSM_after, 1));


% Time vector
num_bins = size(Psth_CSM_before, 2);
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


xline(csmmean,'--')
xline(0,'--')
title('PSTH – Before Peak Firing');

xlabel('Time from Pos Thresh (stim)(ms)');
ylabel('Firing rate (zscore)');
xlim([-250 250])
ylim([0 12])
box off;

% === PLOT: AFTER PEAK UNITS ===
%figure; hold on;

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

xline(csmmean,'--')
xline(0,'--')
title('PSTH – After Peak Firing');
xlabel('Time from POs Thresh(Stim) (ms)');

ylabel('Firing rate (zscore)');
xlim([-250 250])
ylim([0 16])
box off;


%% Diff for before and after csm vs nocsm
beforecolor=rgb('OrangeRed');
aftercolor=rgb('DarkBlue');


% Compute average PSTHs (mean across units)
% diff_before=PsthMaster_WOE_CSM-PsthMaster_BL_CSM;
% diff_after=PsthMaster_WOE_NoCSM-PsthMaster_BL_NoCSM;
% 
%diff_before=Psth_CSM_before-Psth_NoCSM_before;
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
%target_x = 106; % from peak and valley data T2 aligned
target_x = 116; % from peak and valley data T2 aligned
%target_x = 15; % from peak and valley data T2 aligned
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

xlim([-100 250])
%xlim([-200 150])
%ylim([-3 5])
xline(106)
xline(116)
%xline(15)
%ylim([-25 25])
text(target_x + 5, y_after, sprintf('%.2f', y_after), 'Color', aftercolor, 'FontSize', 10);
xlabel('Time from CSM(ms)');
ylabel('Firing rate (Hz)');
title('Av. Diff in FR WOE-BL CSM for before and after');
%title('Av. Diff in FR AD-BL CSM for before and after');
%% calculate the p value 
z = 3.38;
p = 2 * (1 - normcdf(abs(z)));

%% NP data- remove entires with no input - 5 and 7

NPdata_edit = NPdata;
NPdata_edit([5,6, 7]) = [];
%% Now look at kinematics
% Collect all NoCSM reaches and CSM
BLAll= struct; 
ADAll=struct;
WOEAll=struct;
% Combine all reachse across animals 
for i = 1:length(NPdata_edit)
 
    Reach=NPdata_edit(i).ExcitRR_BL.ReReaches;
    BLAll(i).BL=Reach; 

    ReachAD=NPdata_edit(i).ExcitRR_AD.ReReaches;
    ADAll(i).AD=ReachAD; 

    ReachWOE=NPdata_edit(i).ExcitRR_WOE.ReReaches;
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
      
        %C= AlignPeakV(BLR,dim,250);
        %C= AlignThresh(BLR,thresh,dir,250);
        C= AlignThreshCSMHF(BLR,thresh,dir,250);
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
        %C= AlignPeakV(ADR,dim,250);
        %C= AlignStim(ADR,250);
        C= AlignThreshCSMHF(ADR,thresh,dir,250);
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
      
        %C= AlignPeakV(WOER,dim,250);
       % C= AlignThresh(WOER,thresh,dir,250);
        C= AlignThreshCSMHF(WOER,thresh,dir,250);
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
x1=0 ;x2=0.05;
y_upper=30; y_lower=-40;

dim=5;

figure('Position', [100, 100, 600, 300]);  hold on
StimvelPlot(BL_big,dim,BLcolor,0);
StimvelPlot(AD_big,dim,ADcolor,0);
StimvelPlot(WOE_big,dim,WOcolor,0);
yline(0,'--')
xline(0,'--')
ylim('auto')
%xline(csmmean_sec,':')
hold on 
% Fill the area between the vertical lines
%fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
xlim([-0.25 0.25])
ylabel('Outward Velocity (cm/s)')
%xlabel('Time from Stim (s)')
xlabel('Time from CSM  (s)')
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
x1=0 ;x2=0.05;
y_upper=200; y_lower=-200;

figure('Position', [100, 100, 600, 300]);  
shadedErrorBar(time, diff_avg_y, SE_diff_y, 'lineProps', {'k', 'LineWidth', 2});
hold on 
% Fill the area between the vertical lines
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
csmmean_sec=0.0815;
xline(csmmean_sec,':')
xlabel('Time from Stim (ms)');
xlim([-0.25 0.25])
ylabel('Delta- WOE-BL');
title('WOE-BL Kin Diff');
yline(0, '--');
%%
% For example, just plot the y-dimension (2nd column)
% diff_avg_y = AD_big.Average(:,5) - BL_big.Average(:,5);
% diff_avg_y=diff_avg_y*10;
% SE_diff_y = sqrt((AD_big.Std(:,5).^2 ./ AD_big.N) + (BL_big.Std(:,5).^2 ./ BL_big.N));

diff_avg_y = AD_big.Average(:,8) - BL_big.Average(:,8);
diff_avg_y=diff_avg_y*10;
SE_diff_y = sqrt((AD_big.Std(:,8).^2 ./ AD_big.N) + (BL_big.Std(:,8).^2 ./ BL_big.N));

time = (BL_big.Average(:,1))'; % ensure it's a column vector
x1=0 ;x2=0.05;
y_upper=200; y_lower=-200;
% Now plot
figure('Position', [100, 100, 600, 300]);  
shadedErrorBar(time, diff_avg_y, SE_diff_y, 'lineProps', {'k', 'LineWidth', 2});
hold on 
% Fill the area between the vertical lines
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
csmmean_sec=0.0815;
xline(csmmean_sec,':')
xlabel('Time from Stim (ms)');
xlim([-0.25 0.25])
ylabel('Y-Velocity Difference (AD - BL)');
title('AD-BL Kin Diff');
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


