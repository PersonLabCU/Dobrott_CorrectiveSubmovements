clear; clc;

cd("C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\ShamAnimals\ReachandCSMMod_allCSMStogther\CSMAll_wtihMaxV")
MyFolderInfo = dir('*CSMAll.mat');

%%
total_NC = 0;

for day = 1:length(MyFolderInfo)

    fprintf('\n========== SESSION %d / %d ==========\n', day, length(MyFolderInfo));
    load(MyFolderInfo(day).name, 'data');

    % ---- find a neurontype-like field name robustly ----
    fn = fieldnames(data.units);


    ntField = 'neurontype';


    % ---- extract values and normalize (string, trim whitespace, uppercase) ----
    raw = {data.units.(ntField)};
    types = upper(strtrim(string(raw)));



    nNC = sum(types == "NC");
    fprintf('Session %d: %d NC neurons\n', day, nNC);

    total_NC = total_NC + nNC;

end

fprintf('\nTOTAL NC neurons across all sessions: %d\n', total_NC);

%%
% ===================== NC ONLY: MASTER PSTH (max_vx) =====================
cd("C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\ShamAnimals\ReachandCSMMod_allCSMStogther\CSMAll_wtihMaxV");
MyFolderInfo = dir('*CSMAll.mat');

tic
clearvars -except MyFolderInfo

PsthMasterNC_MaxV = [];     % <-- ONE master PSTH (NC only)
NC_total = 0;

for day = 1:length(MyFolderInfo)

    fprintf('\n========== SESSION %d / %d ==========\n', day, length(MyFolderInfo));
    load(MyFolderInfo(day).name, 'Reaches', 'data', 'filtered_modulated_indices','ExcitRR');

 

    % ---------------------- FIND NC UNITS (full session) ----------------------
    ntField = 'neurontype';  % (you confirmed this is correct)
    raw  = {data.units.(ntField)};
    types = upper(strtrim(string(raw)));

    NC_all_idx = find(types == "NC");

    NC_keep_idx = NC_all_idx;


    % Now this matches your style: "data_opto_all = data.units(filtered_modulated_indices)"
    data_opto_all = data.units(NC_keep_idx);

    % ---------------------- TIME VECTOR + PEAK WINDOW ------------------------
    timevec = data_opto_all(1).max_vx.time_bins;
    %timevec = data_opto_all(1).end_point.time_bins;
    valid_peak_range = find(timevec >= -1000 & timevec <= 1000);

    % Gate: keep units with peak in window
    keep_unit = false(length(data_opto_all), 1);
    peak_idx_all = NaN(length(data_opto_all), 1);

    for t = 1:length(data_opto_all)
        psth = data_opto_all(t).max_vx.PSTH_smoothbin;
        [~, pk] = max(psth);
        peak_idx_all(t) = pk;
        if ismember(pk, valid_peak_range)
            keep_unit(t) = true;
        end
    end

    data_opto = data_opto_all(keep_unit);

    if isempty(data_opto)
        fprintf('Day %d: No NC units with peaks in window. Skipping.\n', day);
        continue
    end

    % ---------------------- BUILD MASTER PSTH (NC only, max_vx) --------------
    nCells = length(data_opto);
    nBins  = length(data_opto(1).max_vx.PSTH_smoothbin);

    psth_master_nc = NaN(nCells, nBins);

    for t = 1:nCells
        psth = data_opto(t).max_vx.PSTH_smoothbin;
        %psth = data_opto(t).end_point.PSTH_smoothbin;
        if max(psth) > 0
            psth_master_nc(t,:) = psth / max(psth);   % normalize to max (your style)
        else
            psth_master_nc(t,:) = NaN(1,nBins);
        end
    end

    % Concatenate across days
    PsthMasterNC_MaxV = vertcat(PsthMasterNC_MaxV, psth_master_nc);

    NC_total = NC_total + nCells;
    fprintf('Day %d: Added %d NC units (running total %d)\n', day, nCells, NC_total);

    clear data Reaches filtered_modulated_indices ExcitRR
end

fprintf('\nTOTAL NC units added to Master PSTH: %d\n', size(PsthMasterNC_MaxV,1));
toc
%%
% ---------------------- HEATMAP ----------------------
if ~isempty(PsthMasterNC_MaxV)
    figure('Color','w','Position',[150 150 800 650]);
    imagesc(timevec, 1:size(PsthMasterNC_MaxV,1), PsthMasterNC_MaxV);
    axis tight;
    xlabel('Time from max V (ms)');
    ylabel('NC units');
    title(sprintf('NC Master PSTH aligned to max V (n=%d)', size(PsthMasterNC_MaxV,1)));
    colorbar;
end

%%
% ---------------------- SORT BY PEAK TIME ----------------------
% (uses the same normalized PSTH you already built)
[~, peak_idx] = max(PsthMasterNC_MaxV, [], 2);        % peak bin per unit
[~, sort_idx] = sort(peak_idx, 'ascend');             % early -> late
PsthMasterNC_MaxV_sorted = PsthMasterNC_MaxV(sort_idx,:);

% ---------------------- HEATMAP (SORTED) ----------------------
figure('Color','w','Position',[150 150 800 650]);
imagesc(timevec, 1:size(PsthMasterNC_MaxV_sorted,1), PsthMasterNC_MaxV_sorted);
axis tight;
xline(-250,'r')
xline(250,'r')
xlabel('Time from max V (ms)');
ylabel('NC units (sorted by peak time)');
title(sprintf('NC Master PSTH aligned to max V (sorted by peak time) | n=%d', size(PsthMasterNC_MaxV_sorted,1)));
colormap(slanCM('GnBu'));
colorbar;

%% looking for units that are also reach modulated : 
%%
% ===================== NC ONLY: MASTER PSTH (max_vx) =====================
cd("C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\ShamAnimals\ReachandCSMMod_allCSMStogther\CSMAll_wtihMaxV");
MyFolderInfo = dir('*CSMAll.mat');

tic
clearvars -except MyFolderInfo

PsthMasterNC_MaxV = [];          % master PSTH (NC only)
IsModulatedMaster = [];          % parallel logical vector: 1 = reach modulated
UnitDayMaster = [];              % optional: store day per unit
OrigUnitIdxMaster = [];          % optional: original unit index in session

NC_total = 0;

for day = 1:length(MyFolderInfo)

    fprintf('\n========== SESSION %d / %d ==========\n', day, length(MyFolderInfo));
    load(MyFolderInfo(day).name, 'Reaches', 'data', 'filtered_modulated_indices', 'ExcitRR');

    % ---------------------- FIND NC UNITS (full session) ----------------------
    ntField = 'neurontype';
    raw  = {data.units.(ntField)};
    types = upper(strtrim(string(raw)));

    NC_all_idx = find(types == "NC");
    NC_keep_idx = NC_all_idx;

    % full NC unit structs for this session
    data_opto_all = data.units(NC_keep_idx);

    % ---------------------- TIME VECTOR + PEAK WINDOW ----------------------
    timevec = data_opto_all(1).max_vx.time_bins;
    valid_peak_range = find(timevec >= -1000 & timevec <= 1000);

    % ---------------------- GATE: keep units with peak in window ----------------------
    keep_unit = false(length(data_opto_all), 1);
    peak_idx_all = NaN(length(data_opto_all), 1);

    for t = 1:length(data_opto_all)
        psth = data_opto_all(t).max_vx.PSTH_smoothbin;
        [~, pk] = max(psth);
        peak_idx_all(t) = pk;

        if ismember(pk, valid_peak_range)
            keep_unit(t) = true;
        end
    end

    % apply peak gate
    data_opto = data_opto_all(keep_unit);
    kept_orig_idx = NC_keep_idx(keep_unit);   % original session indices of kept NC units

    if isempty(data_opto)
        fprintf('Day %d: No NC units with peaks in window. Skipping.\n', day);
        continue
    end

    % ---------------------- WHICH KEPT NC UNITS ARE REACH MODULATED? ----------------------
    % filtered_modulated_indices are assumed to refer to indices in data.units
    is_modulated_kept = ismember(kept_orig_idx, filtered_modulated_indices);

    % ---------------------- BUILD MASTER PSTH (NC only, max_vx) ----------------------
    nCells = length(data_opto);
    nBins  = length(data_opto(1).max_vx.PSTH_smoothbin);

    psth_master_nc = NaN(nCells, nBins);

    for t = 1:nCells
        psth = data_opto(t).max_vx.PSTH_smoothbin;

        if max(psth) > 0
            psth_master_nc(t,:) = psth / max(psth);   % normalize to max
        else
            psth_master_nc(t,:) = NaN(1,nBins);
        end
    end

    % ---------------------- CONCATENATE ACROSS DAYS ----------------------
    PsthMasterNC_MaxV = vertcat(PsthMasterNC_MaxV, psth_master_nc);
    IsModulatedMaster = vertcat(IsModulatedMaster, is_modulated_kept(:));
    UnitDayMaster     = vertcat(UnitDayMaster, repmat(day, nCells, 1));
    OrigUnitIdxMaster = vertcat(OrigUnitIdxMaster, kept_orig_idx(:));

    NC_total = NC_total + nCells;
    fprintf('Day %d: Added %d NC units (running total %d)\n', day, nCells, NC_total);
    fprintf('        %d / %d kept NC units are reach modulated\n', sum(is_modulated_kept), nCells);

    clear data Reaches filtered_modulated_indices ExcitRR
end

fprintf('\nTOTAL NC units added to Master PSTH: %d\n', size(PsthMasterNC_MaxV,1));
fprintf('TOTAL reach-modulated NC units in master PSTH: %d\n', sum(IsModulatedMaster));
toc

%%
% ---------------------- HEATMAP (UNSORTED) ----------------------
if ~isempty(PsthMasterNC_MaxV)
    figure('Color','w','Position',[150 150 850 650]);
    imagesc(timevec, 1:size(PsthMasterNC_MaxV,1), PsthMasterNC_MaxV);
    axis tight;
    set(gca,'YDir','normal');
    xlabel('Time from max V (ms)');
    ylabel('NC units');
    title(sprintf('NC Master PSTH aligned to max V (n=%d)', size(PsthMasterNC_MaxV,1)));
    colormap(slanCM('GnBu'));
    colorbar;
    hold on;

    % mark reach-modulated rows with black dots just left of heatmap
    mod_rows = find(IsModulatedMaster);
    x_mark = min(timevec) - 0.03*range(timevec);

    plot(repmat(x_mark, size(mod_rows)), mod_rows, 'k.', 'MarkerSize', 12);

    % expand xlim so the dots are visible
    xlim([x_mark - 0.02*range(timevec), max(timevec)]);
end

%%
% ---------------------- SORT BY PEAK TIME ----------------------
[~, peak_idx] = max(PsthMasterNC_MaxV, [], 2);        % peak bin per unit
[~, sort_idx] = sort(peak_idx, 'ascend');             % early -> late

PsthMasterNC_MaxV_sorted = PsthMasterNC_MaxV(sort_idx,:);
IsModulatedMaster_sorted = IsModulatedMaster(sort_idx);
peak_idx_sorted          = peak_idx(sort_idx);
peak_time_sorted         = timevec(peak_idx_sorted);

%%
% ---------------------- HEATMAP (SORTED) ----------------------
figure('Color','w','Position',[150 150 900 650]);
imagesc(timevec, 1:size(PsthMasterNC_MaxV_sorted,1), PsthMasterNC_MaxV_sorted);
axis tight;
set(gca,'YDir','normal');
hold on;

xline(-250,'r','LineWidth',1.5);
xline(250,'r','LineWidth',1.5);

xlabel('Time from max V (ms)');
ylabel('NC units (sorted by peak time)');
title(sprintf('NC Master PSTH aligned to max V (sorted by peak time) | n=%d', ...
    size(PsthMasterNC_MaxV_sorted,1)));
colormap(slanCM('GnBu'));
colorbar;

% ---------------------- OVERLAY MODULATED UNITS ----------------------
% Option 1: black dots to the left of rows that are reach modulated
mod_rows_sorted = find(IsModulatedMaster_sorted);
x_mark = min(timevec) - 0.03*range(timevec);

plot(repmat(x_mark, size(mod_rows_sorted)), mod_rows_sorted, 'k.', 'MarkerSize', 12);

% expand xlim so marker column is visible
xlim([x_mark - 0.02*range(timevec), max(timevec)]);

% ---------------------- OPTIONAL: plot peak times ----------------------
% plot(peak_time_sorted, 1:length(peak_time_sorted), 'k.', 'MarkerSize', 6);


% ---------------------- HEATMAP (SORTED) ----------------------
figure('Color','w','Position',[150 150 900 650]);
imagesc(timevec, 1:size(PsthMasterNC_MaxV_sorted,1), PsthMasterNC_MaxV_sorted);
axis tight;
set(gca,'YDir','reverse');   % <-- earliest units now at the TOP
hold on;

xline(-250,'r','LineWidth',1.5);
xline(250,'r','LineWidth',1.5);

xlabel('Time from max V (ms)');
ylabel('NC units (sorted by peak time)');
title(sprintf('NC Master PSTH aligned to max V (sorted by peak time) | n=%d', ...
    size(PsthMasterNC_MaxV_sorted,1)));

colormap(slanCM('GnBu'));
colorbar;

% ---------------------- OVERLAY MODULATED UNITS ----------------------
mod_rows_sorted = find(IsModulatedMaster_sorted);
x_mark = min(timevec) - 0.03*range(timevec);

plot(repmat(x_mark, size(mod_rows_sorted)), mod_rows_sorted, 'k.', 'MarkerSize', 12);

xlim([x_mark - 0.02*range(timevec), max(timevec)]);


%%
% ===================== NC ONLY: MASTER PSTH (max_vx) =====================
% Keeps NC units only
% Additional gating:
%   (1) peak must fall within [-1000, 1000] ms
%   (2) unit must be active between [-250, 250] ms around max V
%       defined here as PSTH reaching at least 25% of that unit's max
% Also tracks which kept units are in filtered_modulated_indices
% Final heatmap is sorted by peak time with earliest peak units at the TOP

cd("C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\ShamAnimals\ReachandCSMMod_allCSMStogther\CSMAll_wtihMaxV");
MyFolderInfo = dir('*CSMAll.mat');

tic
clearvars -except MyFolderInfo

% ---------------------- MASTER ARRAYS ----------------------
PsthMasterNC_MaxV   = [];
IsModulatedMaster   = [];
UnitDayMaster       = [];
OrigUnitIdxMaster   = [];
PeakTimeMaster      = [];

NC_total = 0;

% ---------------------- USER SETTINGS ----------------------
peak_window_ms      = [-1000 1000];   % peak must fall here
active_window_ms    = [-250 250];     % must be active here
activity_frac_thr   = 0.25;           % active if PSTH >= 25% of unit max in active window

for day = 1:length(MyFolderInfo)

    fprintf('\n========== SESSION %d / %d ==========\n', day, length(MyFolderInfo));

    load(MyFolderInfo(day).name, 'Reaches', 'data', 'filtered_modulated_indices', 'ExcitRR');

    % ---------------------- FIND NC UNITS ----------------------
    ntField = 'neurontype';
    raw   = {data.units.(ntField)};
    types = upper(strtrim(string(raw)));

    NC_all_idx  = find(types == "NC");
    NC_keep_idx = NC_all_idx;

    if isempty(NC_keep_idx)
        fprintf('Day %d: No NC units found. Skipping.\n', day);
        clear data Reaches filtered_modulated_indices ExcitRR
        continue
    end

    data_opto_all = data.units(NC_keep_idx);

    % ---------------------- TIME VECTOR ----------------------
    timevec = data_opto_all(1).max_vx.time_bins;

    peak_window_idx   = find(timevec >= peak_window_ms(1)   & timevec <= peak_window_ms(2));
    active_window_idx = find(timevec >= active_window_ms(1) & timevec <= active_window_ms(2));

    if isempty(peak_window_idx) || isempty(active_window_idx)
        error('Peak or active window is empty. Check timevec and requested windows.');
    end

    % ---------------------- GATING ----------------------
    nUnits_all        = length(data_opto_all);
    keep_unit         = false(nUnits_all, 1);
    peak_idx_all      = NaN(nUnits_all, 1);
    peak_time_all     = NaN(nUnits_all, 1);
    active_flag_all   = false(nUnits_all, 1);

    for t = 1:nUnits_all

        psth = data_opto_all(t).max_vx.PSTH_smoothbin;

        if isempty(psth) || all(isnan(psth))
            continue
        end

        [pk_val, pk] = max(psth);
        peak_idx_all(t)  = pk;
        peak_time_all(t) = timevec(pk);

        % test 1: peak falls within desired peak window
        peak_in_range = ismember(pk, peak_window_idx);

        % test 2: unit is active between -250 and 250 ms
        psth_win = psth(active_window_idx);

        if pk_val > 0 && ~all(isnan(psth_win))
            active_near_maxV = any(psth_win >= activity_frac_thr * pk_val);
        else
            active_near_maxV = false;
        end

        active_flag_all(t) = active_near_maxV;

        % final keep rule
        keep_unit(t) = peak_in_range && active_near_maxV;
    end

    % apply gating
    data_opto       = data_opto_all(keep_unit);
    kept_orig_idx   = NC_keep_idx(keep_unit);
    kept_peak_times = peak_time_all(keep_unit);

    if isempty(data_opto)
        fprintf('Day %d: No NC units passed gates. Skipping.\n', day);
        clear data Reaches filtered_modulated_indices ExcitRR
        continue
    end

    % ---------------------- WHICH KEPT UNITS ARE REACH MODULATED? ----------------------
    is_modulated_kept = ismember(kept_orig_idx, filtered_modulated_indices);

    % ---------------------- BUILD SESSION PSTH MASTER ----------------------
    nCells = length(data_opto);
    nBins  = length(data_opto(1).max_vx.PSTH_smoothbin);

    psth_master_nc = NaN(nCells, nBins);

    for t = 1:nCells
        psth = data_opto(t).max_vx.PSTH_smoothbin;

        if max(psth) > 0
            psth_master_nc(t,:) = psth ./ max(psth);   % normalize each unit to its own max
        else
            psth_master_nc(t,:) = NaN(1, nBins);
        end
    end

    % ---------------------- CONCATENATE ACROSS DAYS ----------------------
    PsthMasterNC_MaxV = vertcat(PsthMasterNC_MaxV, psth_master_nc);
    IsModulatedMaster = vertcat(IsModulatedMaster, is_modulated_kept(:));
    UnitDayMaster     = vertcat(UnitDayMaster, repmat(day, nCells, 1));
    OrigUnitIdxMaster = vertcat(OrigUnitIdxMaster, kept_orig_idx(:));
    PeakTimeMaster    = vertcat(PeakTimeMaster, kept_peak_times(:));

    NC_total = NC_total + nCells;

    fprintf('Day %d: Added %d NC units (running total %d)\n', day, nCells, NC_total);
    fprintf('        Reach-modulated among kept: %d / %d\n', sum(is_modulated_kept), nCells);

    clear data Reaches filtered_modulated_indices ExcitRR
end

fprintf('\nTOTAL NC units added to Master PSTH: %d\n', size(PsthMasterNC_MaxV,1));
fprintf('TOTAL reach-modulated NC units in master PSTH: %d\n', sum(IsModulatedMaster));
toc

%%
% ===================== UNSORTED HEATMAP =====================
if ~isempty(PsthMasterNC_MaxV)

    figure('Color','w','Position',[150 150 900 650]);
    imagesc(timevec, 1:size(PsthMasterNC_MaxV,1), PsthMasterNC_MaxV);
    axis tight;
    set(gca,'YDir','reverse');   % first rows at top
    hold on;

    xline(active_window_ms(1), 'r', 'LineWidth', 1.5);
    xline(active_window_ms(2), 'r', 'LineWidth', 1.5);

    xlabel('Time from max V (ms)');
    ylabel('NC units');
    title(sprintf('NC Master PSTH aligned to max V | n=%d', size(PsthMasterNC_MaxV,1)));
    colormap(slanCM('GnBu'));
    colorbar;

    % mark reach-modulated units on left
    mod_rows = find(IsModulatedMaster);
    x_mark = min(timevec) - 0.03 * range(timevec);

    plot(repmat(x_mark, size(mod_rows)), mod_rows, 'k.', 'MarkerSize', 12);
    xlim([x_mark - 0.02 * range(timevec), max(timevec)]);
end

%%
% ===================== SORTED HEATMAP =====================
% Red dots = units with peak between -250 and 250 ms
% Black dots = all other units

% define which sorted units are "included"
peak_in_center_sorted = PeakTimeMaster_sorted >= -250 & PeakTimeMaster_sorted <= 250;

figure('Color','w','Position',[150 150 900 650]);
imagesc(timevec, 1:size(PsthMasterNC_MaxV_sorted,1), PsthMasterNC_MaxV_sorted);
axis tight;
set(gca,'YDir','reverse');   % earliest peak units at TOP
hold on;

% vertical reference lines
xline(active_window_ms(1), 'r', 'LineWidth', 1.5);
xline(active_window_ms(2), 'r', 'LineWidth', 1.5);

% labels
xlabel('Time from max V (ms)');
ylabel('NC units (sorted by peak time)');
title(sprintf('NC Master PSTH aligned to max V (sorted by peak time) | n=%d', ...
    size(PsthMasterNC_MaxV_sorted,1)));

% colormap
colormap(slanCM('GnBu'));
colorbar;

% ---------------------- LEFT-SIDE DOT ANNOTATION ----------------------
x_mark = min(timevec) - 0.03 * range(timevec);

rows_red   = find(peak_in_center_sorted);    % included units
rows_black = find(~peak_in_center_sorted);   % all others

% black dots = not included
plot(repmat(x_mark, numel(rows_black), 1), rows_black, ...
    '.', 'Color', 'k', 'MarkerSize', 12);

% red dots = included
plot(repmat(x_mark, numel(rows_red), 1), rows_red, ...
    '.', 'Color', 'r', 'MarkerSize', 14);

% expand x-axis so dots are visible
xlim([x_mark - 0.02 * range(timevec), max(timevec)]);

% ---------------------- OPTIONAL: PEAK TIME DOTS ON HEATMAP ----------------------
% plot(PeakTimeMaster_sorted, 1:length(PeakTimeMaster_sorted), ...
%     '.', 'Color', [0 0 0], 'MarkerSize', 6);

% ---------------------- SUMMARY ----------------------
fprintf('\nSorted heatmap summary:\n');
fprintf('  Total kept NC units: %d\n', size(PsthMasterNC_MaxV_sorted,1));
fprintf('  Units with peak between -250 and 250 ms: %d\n', sum(peak_in_center_sorted));
fprintf('  Units outside that window: %d\n', sum(~peak_in_center_sorted));

%%
% ===================== SORT BY PEAK TIME =====================
[~, peak_idx] = max(PsthMasterNC_MaxV, [], 2);
[~, sort_idx] = sort(peak_idx, 'ascend');   % early peak -> late peak

PsthMasterNC_MaxV_sorted = PsthMasterNC_MaxV(sort_idx, :);
IsModulatedMaster_sorted = IsModulatedMaster(sort_idx);
PeakTimeMaster_sorted    = PeakTimeMaster(sort_idx);
UnitDayMaster_sorted     = UnitDayMaster(sort_idx);
OrigUnitIdxMaster_sorted = OrigUnitIdxMaster(sort_idx);

%%
% ===================== SORTED HEATMAP =====================
figure('Color','w','Position',[150 150 900 650]);
imagesc(timevec, 1:size(PsthMasterNC_MaxV_sorted,1), PsthMasterNC_MaxV_sorted);
axis tight;
set(gca,'YDir','reverse');   % earliest peak units at TOP
hold on;

xline(active_window_ms(1), 'r', 'LineWidth', 1.5);
xline(active_window_ms(2), 'r', 'LineWidth', 1.5);

xlabel('Time from max V (ms)');
ylabel('NC units (sorted by peak time)');
title(sprintf('NC Master PSTH aligned to max V (sorted by peak time) | n=%d', ...
    size(PsthMasterNC_MaxV_sorted,1)));
colormap(slanCM('GnBu'));
colorbar;

% ---------------------- MARK REACH-MODULATED UNITS ----------------------
mod_rows_sorted = find(IsModulatedMaster_sorted);
x_mark = min(timevec) - 0.03 * range(timevec);

plot(repmat(x_mark, size(mod_rows_sorted)), mod_rows_sorted, 'k.', 'MarkerSize', 12);
xlim([x_mark - 0.02 * range(timevec), max(timevec)]);

% ---------------------- OPTIONAL: PLOT PEAK TIMES ----------------------
% plot(PeakTimeMaster_sorted, 1:length(PeakTimeMaster_sorted), 'k.', 'MarkerSize', 6);

%
% ===================== OPTIONAL: BOX UNITS WITH PEAKS IN -250 TO 250 =====================
peak_in_center_sorted = PeakTimeMaster_sorted >= -250 & PeakTimeMaster_sorted <= 250;

d = diff([0; peak_in_center_sorted; 0]);
run_starts = find(d == 1);
run_ends   = find(d == -1) - 1;

for i = 1:length(run_starts)
    y1 = run_starts(i) - 0.5;
    y2 = run_ends(i)   + 0.5;

    rectangle('Position', [timevec(1), y1, timevec(end)-timevec(1), y2-y1], ...
              'EdgeColor', [0 0 0], ...
              'LineWidth', 1.5, ...
              'LineStyle', '-');
end

%
% ===================== OPTIONAL SUMMARY =====================
fprintf('\nSorted heatmap summary:\n');
fprintf('  Total kept NC units: %d\n', size(PsthMasterNC_MaxV_sorted,1));
fprintf('  Reach-modulated kept NC units: %d\n', sum(IsModulatedMaster_sorted));
fprintf('  Units with peak between -250 and 250 ms: %d\n', sum(peak_in_center_sorted));
%% 
% BLOCK STIM

% ===================== NC ONLY: MASTER PSTH (max_vx) =====================
% ---------------------- LOAD SESSIONS -----------------------
cd("C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\StimAnimals\BlockSeparated_CSMvsNoCSM_inducedalignment");
MyFolderInfo = dir('*Block_CSMvsNoCSM.mat');


tic
clearvars -except MyFolderInfo

PsthMasterNC_MaxV = [];     % <-- ONE master PSTH (NC only)
NC_total = 0;

for day = 1:length(MyFolderInfo)

    fprintf('\n========== SESSION %d / %d ==========\n', day, length(MyFolderInfo));
    load(MyFolderInfo(day).name, 'Reaches', 'data', 'filtered_modulated_indices');

 

    % ---------------------- FIND NC UNITS (full session) ----------------------
    ntField = 'neurontype';  % (you confirmed this is correct)
    raw  = {data.units.(ntField)};
    types = upper(strtrim(string(raw)));

    NC_all_idx = find(types == "NC");

    NC_keep_idx = NC_all_idx;


    % Now this matches your style: "data_opto_all = data.units(filtered_modulated_indices)"
    data_opto_all = data.units(NC_keep_idx);

    % ---------------------- TIME VECTOR + PEAK WINDOW ------------------------
    timevec = data_opto_all(1).Ctrl_max_vx.time_bins;
    %timevec = data_opto_all(1).end_point.time_bins;
    valid_peak_range = find(timevec >= -1000 & timevec <= 1000);

    % Gate: keep units with peak in window
    keep_unit = false(length(data_opto_all), 1);
    peak_idx_all = NaN(length(data_opto_all), 1);

    for t = 1:length(data_opto_all)
        psth = data_opto_all(t).Ctrl_max_vx.PSTH_smoothbin;
        [~, pk] = max(psth);
        peak_idx_all(t) = pk;
        if ismember(pk, valid_peak_range)
            keep_unit(t) = true;
        end
    end

    data_opto = data_opto_all(keep_unit);

    if isempty(data_opto)
        fprintf('Day %d: No NC units with peaks in window. Skipping.\n', day);
        continue
    end

    % ---------------------- BUILD MASTER PSTH (NC only, max_vx) --------------
    nCells = length(data_opto);
    nBins  = length(data_opto(1).Ctrl_max_vx.PSTH_smoothbin);

    psth_master_nc = NaN(nCells, nBins);

    for t = 1:nCells
        %psth = data_opto(t).max_vx.PSTH_smoothbin;
        psth = data_opto_all(t).Ctrl_max_vx.PSTH_smoothbin;
        psth = data_opto_all(t).Ctrl_endpt.PSTH_smoothbin;
        %psth = data_opto(t).end_point.PSTH_smoothbin;
        if max(psth) > 0
            psth_master_nc(t,:) = psth / max(psth);   % normalize to max (your style)
        else
            psth_master_nc(t,:) = NaN(1,nBins);
        end
    end

    % Concatenate across days
    PsthMasterNC_MaxV = vertcat(PsthMasterNC_MaxV, psth_master_nc);

    NC_total = NC_total + nCells;
    fprintf('Day %d: Added %d NC units (running total %d)\n', day, nCells, NC_total);

    clear data Reaches filtered_modulated_indices ExcitRR
end

fprintf('\nTOTAL NC units added to Master PSTH: %d\n', size(PsthMasterNC_MaxV,1));
toc
%%
% ---------------------- HEATMAP ----------------------
if ~isempty(PsthMasterNC_MaxV)
    figure('Color','w','Position',[150 150 800 650]);
    imagesc(timevec, 1:size(PsthMasterNC_MaxV,1), PsthMasterNC_MaxV);
    axis tight;
    xlabel('Time from max V (ms)');
    ylabel('NC units');
    title(sprintf('NC Master PSTH aligned to max V (n=%d)', size(PsthMasterNC_MaxV,1)));
    colorbar;
end

%%
% ---------------------- SORT BY PEAK TIME ----------------------
% (uses the same normalized PSTH you already built)
[~, peak_idx] = max(PsthMasterNC_MaxV, [], 2);        % peak bin per unit
[~, sort_idx] = sort(peak_idx, 'ascend');             % early -> late
PsthMasterNC_MaxV_sorted = PsthMasterNC_MaxV(sort_idx,:);

% ---------------------- HEATMAP (SORTED) ----------------------
figure('Color','w','Position',[150 150 800 650]);
imagesc(timevec, 1:size(PsthMasterNC_MaxV_sorted,1), PsthMasterNC_MaxV_sorted);
axis tight;
xlabel('Time from max V (ms)');
ylabel('NC units (sorted by peak time)');
title(sprintf('NC Master PSTH aligned to max V (sorted by peak time) | n=%d', size(PsthMasterNC_MaxV_sorted,1)));
colormap(slanCM('GnBu'));
colorbar;