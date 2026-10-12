

%% ============================================================
% integrated_xcorr_and_psth_classification_NEG_LAG_ONLY_NO_BA.m
% SHAM: *CSMAll.mat
%
% DOES:
%   (1) Compute IFR vs velocity cross-correlation across lags
%   (2) NEGATIVE-LAG ONLY selection: choose lag<0 with maximum |medianCorr|
%       -> category = negLag_posCorr or negLag_negCorr (by sign of chosen corr)
%   (3) KEEP ONLY units whose max_vx PSTH peak is within valid_peak_range_ms
%   (4) Compute Z-scored PSTHs for CSM and NoCSM (baseline from CSM baseline_win_def)
%   (5) Aggregate PSTHs into 2 groups (NO before/after designation)
%   (6) Plot: (A) peak lag vs corr scatter with labels, (B) group PSTHs
%
% EDIT THESE IF NEEDED:
%   - CSM_PSTH_FIELD / NoCSM_PSTH_FIELD
%   - baseline_win_def, valid_peak_range_ms, plot_win
% ============================================================

clear; clc;

% ---------------------- LOAD SESSIONS -----------------------
cd("C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\ShamAnimals\ReachandCSMMod_allCSMStogther\CSMAll_wtihMaxV")
MyFolderInfo = dir('*CSMAll.mat');

%% ---------------------- ANALYSIS PARAMS ---------------------
rasterField   = 'max_vx';
kinVar        = 7;     % aligning to max outward velocity
gauss_sigma   = 10;
resample_fs   = 200;
win_pre_ms    = 100;
win_post_ms   = 150;
maxLag_ms     = 100;
lag_step      = 5;

baseline_win_def     = [-1000 -500]; % ms (for PSTH zscore baseline)
%baseline_win_def     = [-500 -250]; % ms (for PSTH zscore baseline)
%valid_peak_range_ms  = [-100 100];   % ms (keep units whose max_vx PSTH peak falls here)
valid_peak_range_ms  = [-250 250];   % ms (keep units whose max_vx PSTH peak falls here)

% % Choose which PSTHs to aggregate (EDIT)
CSM_PSTH_FIELD   = 'AllCSM_pos_thresh'; % e.g. 'AllCSM_pos_thresh' or 'AllCSM_csmT2'
NoCSM_PSTH_FIELD = 'NoCSM_pos_thresh';  % e.g. 'NoCSM_pos_thresh' or 'NoCSM_csmT2'
% 
% CSM_PSTH_FIELD   = 'AllCSM_csmT2'; % e.g. 'AllCSM_pos_thresh' or 'AllCSM_csmT2'
% NoCSM_PSTH_FIELD = 'NoCSM_csmT2';  % e.g. 'NoCSM_pos_thresh' or 'NoCSM_csmT2'

% CSM_PSTH_FIELD   = 'AllCSM_csmT0'; % e.g. 'AllCSM_pos_thresh' or 'AllCSM_csmT2'
% NoCSM_PSTH_FIELD = 'NoCSM_csmT0';  % e.g. 'NoCSM_pos_thresh' or 'NoCSM_csmT2'

% Plot window for PSTHs
plot_win = [-500 500];

% ---------------------- MASTER STRUCTS ----------------------
XCorrData = struct;
NPdata    = struct;

group_names = {'negLag_posCorr','negLag_negCorr'};
AllGroups  = struct();
for g = 1:numel(group_names)
    gn = group_names{g};
    AllGroups.(gn).CSM   = [];
    AllGroups.(gn).NoCSM = [];
end

%% ====================== MAIN SESSION LOOP ======================
for day = 1:length(MyFolderInfo)

    fprintf('\n========== SESSION %d / %d ==========\n', day, length(MyFolderInfo));
    load(MyFolderInfo(day).name, 'Reaches','data','filtered_modulated_indices','ExcitRR');

    if isempty(filtered_modulated_indices)
        fprintf(' Day %d: no modulated units. Skipping.\n', day);
        continue;
    end

    % optional: save kinematics structs
    NPdata(day).NoCSM     = ExcitRR.NoCSM;
    NPdata(day).ReReaches = ExcitRR.ReReaches;

    session_units = filtered_modulated_indices;
    nUnits = length(session_units);

    % timevec for PSTHs
    timevec = data.units(session_units(1)).max_vx.time_bins; % ms
    baseline_win     = find(timevec >= baseline_win_def(1)    & timevec <= baseline_win_def(2));
    valid_peak_range = find(timevec >= valid_peak_range_ms(1) & timevec <= valid_peak_range_ms(2));

    kept_units_session = [];
    negLag_posCorr_units = [];
    negLag_negCorr_units = [];

    % ------------------ LOOP ACROSS UNITS -------------------
    for uu = 1:nUnits

        unitIdx = session_units(uu);
        fprintf('  Analyzing UNIT %d...\n', unitIdx);

        % ----------------- LOAD SPIKE RASTER -----------------
        spikes_full = data.units(unitIdx).(rasterField).Raster;
        [nTrials, nTime] = size(spikes_full);
        full_time_1ms = (-1000:1000)';

        % ----------------- COMPUTE IFR -----------------------
        IFR_full = zeros(nTime, nTrials);
        kern_len = round(6*gauss_sigma);
        gk = fspecial('gaussian',[1 kern_len], gauss_sigma);
        gk = gk / sum(gk);

        for tr = 1:nTrials
            sp_idx = find(spikes_full(tr,:) > 0);
            if numel(sp_idx) < 2, continue; end

            spike_times = full_time_1ms(sp_idx);
            isis = diff(spike_times)/1000;
            inst_rate = 1 ./ isis;

            rate_vec = zeros(nTime,1);
            for s = 1:length(isis)
                interval_idx = sp_idx(s):sp_idx(s+1)-1;
                rate_vec(interval_idx) = inst_rate(s);
            end

            IFR_full(:,tr) = conv(rate_vec, gk, 'same');
        end

        % Extract analysis window
        win_idx = full_time_1ms >= -win_pre_ms & full_time_1ms <= win_post_ms;
        IFR_1ms = IFR_full(win_idx, :);
        ifr_time_1ms = full_time_1ms(win_idx);

        % ----------------- RESAMPLE IFR -> 200 Hz ----------------
        time_200Hz = (-win_pre_ms:1000/resample_fs:win_post_ms)';
        nSamples200 = length(time_200Hz);

        IFR = nan(nSamples200, nTrials);
        for tr = 1:nTrials
            IFR(:,tr) = interp1(ifr_time_1ms, IFR_1ms(:,tr), time_200Hz, 'linear', 0);
        end

        % ----------------- LOAD KINEMATICS ----------------------
        %Reaches_use = ExcitRR.ReReaches; % matches your previous workflow
        Reaches_use=Reaches;
        C = AlignPeakV_1ms_edit(Reaches_use, kinVar, win_pre_ms, win_post_ms);
        C = stimvel_accel(C);

        vel_1ms_all  = squeeze(C.stimvel(:,5,:)); % velocity
        vel_time_1ms = linspace(-win_pre_ms, win_post_ms, size(vel_1ms_all,1))';

        Vel = nan(nSamples200, nTrials);
        for tr = 1:nTrials
            Vel(:,tr) = interp1(vel_time_1ms, vel_1ms_all(:,tr), time_200Hz, 'linear', NaN);
        end

        % ---------------- SLIDING WINDOW CROSS-CORRELATION --------
        lags_ms = (-maxLag_ms:lag_step:maxLag_ms)';
        nLags = length(lags_ms);
        xcorr_trials = nan(nLags, nTrials);

        for li = 1:nLags
            lag = lags_ms(li);
            shift = round(lag / (1000/resample_fs));

            for tr = 1:nTrials
                fr = IFR(:,tr);
                v  = Vel(:,tr);

                if shift > 0
                    fr_seg = fr(1+shift:end);
                    v_seg  = v(1:end-shift);
                elseif shift < 0
                    fr_seg = fr(1:end+shift);
                    v_seg  = v(1-shift:end);
                else
                    fr_seg = fr;
                    v_seg  = v;
                end

                if length(fr_seg) < 5, continue; end

                fr_sd = std(fr_seg);
                v_sd  = std(v_seg);
                if fr_sd==0 || v_sd==0 || isnan(fr_sd) || isnan(v_sd), continue; end

                fr_z = (fr_seg - mean(fr_seg)) ./ fr_sd;
                v_z  = (v_seg  - mean(v_seg))  ./ v_sd;

                if any(isnan(fr_z)) || any(isnan(v_z)), continue; end
                xcorr_trials(li,tr) = mean(fr_z .* v_z);
            end
        end

        medianCorr = nanmedian(xcorr_trials,2);

        % =========================================================
        % NEGATIVE-LAG ONLY: choose lag<0 with maximum |medianCorr|
        % + quadratic interpolation around selected point
        % =========================================================
        negMask = (lags_ms < 0);
        negIdx  = find(negMask);

        negLag_peakLag = NaN;
        negLag_peakVal = NaN;
        negLag_peakAbs = NaN;
        cat_label      = 'negLag_unknown';

        if ~isempty(negIdx) && any(~isnan(medianCorr(negIdx)))

            [~, kRel] = max(abs(medianCorr(negIdx)));
            idxStar = negIdx(kRel);

            if idxStar > 1 && idxStar < length(medianCorr) && all(~isnan(medianCorr(idxStar-1:idxStar+1)))
                x = lags_ms(idxStar-1:idxStar+1);
                y = medianCorr(idxStar-1:idxStar+1);
                p = polyfit(x, y, 2);

                lag_hat = -p(2)/(2*p(1));
                lag_hat = max(min(lag_hat, max(x)), min(x)); % clamp
                val_hat = polyval(p, lag_hat);

                negLag_peakLag = lag_hat;
                negLag_peakVal = val_hat;
            else
                negLag_peakLag = lags_ms(idxStar);
                negLag_peakVal = medianCorr(idxStar);
            end

            negLag_peakAbs = abs(negLag_peakVal);

            if negLag_peakVal >= 0
                cat_label = 'negLag_posCorr';
                negLag_posCorr_units(end+1) = unitIdx;
            else
                cat_label = 'negLag_negCorr';
                negLag_negCorr_units(end+1) = unitIdx;
            end
        end

        % ---------------- KEEP FILTER (PSTH PEAK IN WINDOW) ----------------
        psth_ref = data.units(unitIdx).max_vx.PSTH_smoothbin;
        [~, peak_idx] = max(psth_ref);

        keep = ismember(peak_idx, valid_peak_range);

        % ---------------- SAVE UNIT-LEVEL RESULTS ----------------
        XCorrData(day).units(uu).unitIdx        = unitIdx;
        XCorrData(day).units(uu).lags_ms        = lags_ms;
        XCorrData(day).units(uu).xcorr_trials   = xcorr_trials;
        XCorrData(day).units(uu).medianCorr     = medianCorr;

        XCorrData(day).units(uu).negLag_peakLag = negLag_peakLag;
        XCorrData(day).units(uu).negLag_peakVal = negLag_peakVal;
        XCorrData(day).units(uu).negLag_peakAbs = negLag_peakAbs;
        XCorrData(day).units(uu).category       = cat_label;

        XCorrData(day).units(uu).keep           = keep;
        XCorrData(day).units(uu).psth_peak_idx  = peak_idx;
        XCorrData(day).units(uu).psth_peak_time = timevec(peak_idx);

        if ~keep
            continue; % skip PSTH aggregation for non-kept units
        end

        kept_units_session(end+1) = unitIdx;

        % ---------------- PSTHs (Z-SCORED) ----------------
        csm   = data.units(unitIdx).(CSM_PSTH_FIELD).PSTH_smoothbin;
        nocsm = data.units(unitIdx).(NoCSM_PSTH_FIELD).PSTH_smoothbin;

        baseline_mean = mean(csm(baseline_win));
        baseline_std  = std(csm(baseline_win));
        if baseline_std == 0 || isnan(baseline_std), baseline_std = 1; end

        csm_z   = (csm   - baseline_mean) / baseline_std;
        nocsm_z = (nocsm - baseline_mean) / baseline_std;

        % ---------------- AGGREGATE BY 2-GROUP CATEGORY ----------------
        if ismember(cat_label, group_names)
            AllGroups.(cat_label).CSM   = vertcat(AllGroups.(cat_label).CSM,   csm_z);
            AllGroups.(cat_label).NoCSM = vertcat(AllGroups.(cat_label).NoCSM, nocsm_z);
        end

    end % unit loop

    % session-level lists
    XCorrData(day).kept_units             = kept_units_session;
    XCorrData(day).negLag_posCorr_units   = negLag_posCorr_units;
    XCorrData(day).negLag_negCorr_units   = negLag_negCorr_units;

    fprintf(' Session %d kept units=%d | negLag_posCorr=%d | negLag_negCorr=%d\n', ...
        day, numel(kept_units_session), numel(negLag_posCorr_units), numel(negLag_negCorr_units));
end

fprintf('\n================ ALL SESSIONS COMPLETE ================\n');

%% ============================================================
% PLOT A: Peak lag vs peak corr (ONLY kept units) + session/unit labels
% ============================================================

allLag  = [];
allCorr = [];
allCat  = {};
allDay  = [];
allUnit = [];

for day = 1:numel(XCorrData)
    if ~isfield(XCorrData(day),'units') || isempty(XCorrData(day).units), continue; end

    for uu = 1:numel(XCorrData(day).units)
        u = XCorrData(day).units(uu);

        if ~isfield(u,'keep') || ~u.keep, continue; end
        if ~isfield(u,'negLag_peakLag') || ~isfield(u,'negLag_peakVal') || ~isfield(u,'category'), continue; end

        lag  = u.negLag_peakLag;
        corr = u.negLag_peakVal;

        if isnan(lag) || isnan(corr), continue; end
        if lag >= 0, continue; end

        allLag(end+1,1)  = lag;
        allCorr(end+1,1) = corr;
        allCat{end+1,1}  = u.category;
        allDay(end+1,1)  = day;
        allUnit{end+1,1} = num2str(u.unitIdx);
    end
end

figure('Position',[350 250 980 560]); hold on;
ms = 42;

idxPos = strcmp(allCat,'negLag_posCorr');
idxNeg = strcmp(allCat,'negLag_negCorr');

h1 = scatter(allLag(idxPos), allCorr(idxPos), ms, 'filled');
h2 = scatter(allLag(idxNeg), allCorr(idxNeg), ms, 'filled');

xline(0,'k--','LineWidth',1);
yline(0,'k--','LineWidth',1);
xlabel('Peak lag (ms) [negative lags only]');
ylabel('Peak correlation (median z·z)');
title('Neg-lag max-|corr| peak per unit (kept units only)');
legend([h1 h2], {'negLag\_posCorr','negLag\_negCorr'}, 'Location','best');
grid on; box on;

dx = 0.6;     % ms
dy = 0.015;   % corr
for i = 1:numel(allLag)
    text(allLag(i)+dx, allCorr(i)+dy, sprintf('S%d:U%s', allDay(i), allUnit{i}), ...
        'FontSize',8, 'Color',[0.25 0.25 0.25], 'Interpreter','none');
end



%% ============================================================
% PLOT A: Peak lag vs peak corr (ONLY kept units)
% Publication-style, SteelBlue (posCorr) + DarkGreen (negCorr)
% NO labels, NO grid
% AXES (actual plot) is ~half-height of the figure
% ============================================================

allLag  = [];
allCorr = [];
allCat  = {};

for day = 1:numel(XCorrData)
    if ~isfield(XCorrData(day),'units') || isempty(XCorrData(day).units), continue; end

    for uu = 1:numel(XCorrData(day).units)
        u = XCorrData(day).units(uu);

        if ~isfield(u,'keep') || ~u.keep, continue; end
        if ~isfield(u,'negLag_peakLag') || ~isfield(u,'negLag_peakVal') || ~isfield(u,'category'), continue; end

        lag  = u.negLag_peakLag;
        corr = u.negLag_peakVal;

        if isnan(lag) || isnan(corr), continue; end
        if lag >= 0, continue; end

        allLag(end+1,1)  = lag;
        allCorr(end+1,1) = corr;
        allCat{end+1,1}  = u.category;
    end
end

% Colors (requires rgb.m in path)
Color_pos = rgb('SteelBlue');
Color_neg = rgb('DarkGreen');

% ---- Figure (keep whatever window size you want)
fig = figure('Color','w','Units','pixels','Position',[550 350 650 650]);

% ---- Axes: make the ACTUAL plot ~half the figure height
ax = axes(fig);
%set(ax, 'Position', [0.15 0.33 0.75 0.45]);  % [left bottom width height] (normalized)
set(ax, 'Position', [0.15 0.33 0.62 0.45]); 
hold(ax,'on');

ms  = 75;

idxPos = strcmp(allCat,'negLag_posCorr');
idxNeg = strcmp(allCat,'negLag_negCorr');

h1 = scatter(ax, allLag(idxPos), allCorr(idxPos), ms, ...
    'MarkerFaceColor', Color_pos, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceAlpha', 0.85, 'MarkerEdgeAlpha', 0.95, 'LineWidth', 0.75);

h2 = scatter(ax, allLag(idxNeg), allCorr(idxNeg), ms, ...
    'MarkerFaceColor', Color_neg, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceAlpha', 0.85, 'MarkerEdgeAlpha', 0.95, 'LineWidth', 0.75);

% Reference lines
xline(ax, 0, '--', 'Color', [0 0 0], 'LineWidth', 1.0, 'Alpha', 0.7);
yline(ax, 0, '--', 'Color', [0 0 0], 'LineWidth', 1.0, 'Alpha', 0.7);

% Limits
ylim(ax, [-1 1]);
if ~isempty(allLag)
    xpad = 0.03 * range(allLag);
    if xpad == 0, xpad = 1; end
    xlim(ax, [min(allLag)-xpad, max(allLag)+xpad]);
end

% Labels
xlabel(ax, 'Peak lag (ms) [negative lags only]');
ylabel(ax, 'Peak correlation (median z·z)');

% --- REMOVE GRID COMPLETELY ---
grid(ax,'off');

% Publication styling
ax.LineWidth  = 1.0;
ax.TickDir    = 'out';
ax.TickLength = [0.012 0.012];
ax.FontName   = 'Arial';
ax.FontSize   = 12;
ax.Box        = 'on';
ax.Layer      = 'top';

% IMPORTANT: do NOT force square, otherwise height changes won't apply
axis(ax,'normal');

% Legend only
lgd = legend(ax, [h1 h2], {'negLag\_posCorr','negLag\_negCorr'}, ...
    'Location','best', 'Box','off');
lgd.FontSize = 11;

% ===================== ADD MARGINAL HISTOGRAMS (TOP + RIGHT) =====================
% --- pull the final limits from the main axis
xl = xlim(ax);
yl = ylim(ax);

% --- main axis position (normalized)
pos = ax.Position;              % [L B W H]
gap  = 0.015;                   % small gap between main plot and hist axes
histH = 0.16;                   % height of TOP hist axis
histW = 0.16;                   % width  of RIGHT hist axis

% --- create TOP histogram axis (above main axis, outside it)
axTop = axes(fig);
axTop.Position = [pos(1), pos(2)+pos(4)+gap, pos(3), histH];
hold(axTop,'on');

% --- create RIGHT histogram axis (to the right of main axis, outside it)
axRight = axes(fig);
axRight.Position = [pos(1)+pos(3)+gap, pos(2), histW, pos(4)];
hold(axRight,'on');

% --- choose bins (feel free to tweak)
nBinsLag  = 30;
nBinsCorr = 30;

% --- TOP: distribution of lag (overlapping)
edgesLag = linspace(xl(1), xl(2), nBinsLag+1);

histogram(axTop, allLag(idxPos), edgesLag, ...
    'FaceColor', Color_pos, 'EdgeColor', 'k', 'FaceAlpha', 0.45);

histogram(axTop, allLag(idxNeg), edgesLag, ...
    'FaceColor', Color_neg, 'EdgeColor', 'k', 'FaceAlpha', 0.45);

xlim(axTop, xl);

% --- RIGHT: distribution of corr (overlapping; horizontal orientation)
edgesCorr = linspace(yl(1), yl(2), nBinsCorr+1);

histogram(axRight, allCorr(idxPos), edgesCorr, ...
    'FaceColor', Color_pos, 'EdgeColor', 'k', 'FaceAlpha', 0.45, 'Orientation','horizontal');

histogram(axRight, allCorr(idxNeg), edgesCorr, ...
    'FaceColor', Color_neg, 'EdgeColor', 'k', 'FaceAlpha', 0.45, 'Orientation','horizontal');

ylim(axRight, yl);

% --- styling: keep them close, clean, and "off" the main plot
set([axTop axRight], 'Color','w', 'Box','off', 'TickDir','out', ...
    'FontName','Arial', 'FontSize', 11, 'LineWidth', 1.0);
grid(axTop,'off');  grid(axRight,'off');

% remove redundant ticks/labels
axTop.XTickLabel   = [];
axTop.XLabel.String = '';
axTop.YLabel.String = '';

axRight.YTickLabel   = [];
axRight.YLabel.String = '';
axRight.XLabel.String = '';

% subtle: only show the outer axes lines you want
axTop.XAxis.Visible   = 'off';   % hides top hist x-axis line/ticks
axRight.YAxis.Visible = 'off';   % hides right hist y-axis line/ticks

% keep limits locked & aligned
linkaxes([ax axTop],'x');
linkaxes([ax axRight],'y');

% make sure the main axis stays on top visually
uistack(ax,'top');


%% ============================================================
% PSTH plotting code UPDATED for CURRENT analysis (NO before/after)
% Current analysis produces ONLY:
%   AllGroups.negLag_posCorr.CSM
%   AllGroups.negLag_posCorr.NoCSM
%   AllGroups.negLag_negCorr.CSM
%   AllGroups.negLag_negCorr.NoCSM
%
% This block mirrors your old style (mean+SEM plots + delta plots),
% but uses the NEW fields and makes NO assumptions about unit #9 etc.
% ============================================================

% -------------------- GET MASTER PSTHs --------------------
PSTH_poscorr_neglag_CSM   = AllGroups.negLag_posCorr.CSM;
PSTH_poscorr_neglag_NoCSM = AllGroups.negLag_posCorr.NoCSM;

PSTH_negcorr_neglag_CSM   = AllGroups.negLag_negCorr.CSM;
PSTH_negcorr_neglag_NoCSM = AllGroups.negLag_negCorr.NoCSM;

% -------------------- SETUP --------------------
Color_CSM   = rgb('DarkMagenta');
Color_NoCSM = rgb('Black');

DeltaColor  = rgb('DarkGreen');
AltDeltaColor = rgb('SteelBlue');

% robust timevec grab:
% use the timevec already defined in your analysis (recommended),
% otherwise fall back to a unit that definitely exists
if exist('timevec','var') && ~isempty(timevec)
    full_timevec = timevec;
else
    % fallback: find first non-empty session/unit
    found = false;
    for d = 1:numel(XCorrData)
        if isfield(XCorrData(d),'units') && ~isempty(XCorrData(d).units)
            u0 = XCorrData(d).units(1).unitIdx;
            full_timevec = data.units(u0).AllCSM_pos_thresh.time_bins; %#ok<NASGU>
            found = true;
            break
        end
    end
    if ~found
        error('Could not determine PSTH time vector (timevec/full_timevec not found).');
    end
end

% Indices for plotting window
time_idx    = find(full_timevec >= -250 & full_timevec <= 250);
%time_idx    = find(full_timevec >= -500 & full_timevec <= 500);
time_vector = full_timevec(time_idx);

% Helper for mean + SEM (matches your previous pattern)
mean_sem = @(mat) deal( ...
    nanmean(mat(:,time_idx),1), ...
    nanstd(mat(:,time_idx),0,1) ./ sqrt(max(1,size(mat,1))) );

% Helper for delta mean + SEM
delta_sem = @(csm,no) deal( ...
    nanmean(csm(:,time_idx) - no(:,time_idx),1), ...
    nanstd(csm(:,time_idx) - no(:,time_idx),0,1) ./ sqrt(max(1,size(csm,1))) );

%% ==================== PLOTTING: MEAN +/- SEM ====================

% 1) NEG LAG, POS CORR
if ~isempty(PSTH_poscorr_neglag_CSM) && ~isempty(PSTH_poscorr_neglag_NoCSM)

    [mCSM, semCSM] = mean_sem(PSTH_poscorr_neglag_CSM);
    [mNo,  semNo ] = mean_sem(PSTH_poscorr_neglag_NoCSM);

    figure('Name','NEG LAG, POS CORR','Position',[100,100,600,300]); hold on;
    plot(time_vector, mCSM,'Color',Color_CSM,'LineWidth',2);
    fill([time_vector fliplr(time_vector)], ...
         [mCSM+semCSM fliplr(mCSM-semCSM)], ...
         Color_CSM,'FaceAlpha',0.2,'EdgeColor','none');

    plot(time_vector, mNo,'Color',Color_NoCSM,'LineWidth',2);
    fill([time_vector fliplr(time_vector)], ...
         [mNo+semNo fliplr(mNo-semNo)], ...
         Color_NoCSM,'FaceAlpha',0.2,'EdgeColor','none');

    xline(0,'--k');
    xlabel('Time from pos. thresh (ms)'); ylabel('Z-scored firing rate');
    %xlim([-250 250]);
    xlim('auto')
    title('Negative Lag, Positive Correlation');
    box off;

else
    fprintf('Skipping NEG LAG POS CORR mean plot (empty group).\n');
end

% 2) NEG LAG, NEG CORR
if ~isempty(PSTH_negcorr_neglag_CSM) && ~isempty(PSTH_negcorr_neglag_NoCSM)

    [mCSM, semCSM] = mean_sem(PSTH_negcorr_neglag_CSM);
    [mNo,  semNo ] = mean_sem(PSTH_negcorr_neglag_NoCSM);

    figure('Name','NEG LAG, NEG CORR','Position',[100,450,600,300]); hold on;
    plot(time_vector, mCSM,'Color',Color_CSM,'LineWidth',2);
    fill([time_vector fliplr(time_vector)], ...
         [mCSM+semCSM fliplr(mCSM-semCSM)], ...
         Color_CSM,'FaceAlpha',0.2,'EdgeColor','none');

    plot(time_vector, mNo,'Color',Color_NoCSM,'LineWidth',2);
    fill([time_vector fliplr(time_vector)], ...
         [mNo+semNo fliplr(mNo-semNo)], ...
         Color_NoCSM,'FaceAlpha',0.2,'EdgeColor','none');

    xline(0,'--k');
    xlabel('Time from pos. thresh (ms)'); ylabel('Z-scored firing rate');
    %xlim([-250 250]);
    xlim('auto')
    title('Negative Lag, Negative Correlation');
    box off;

else
    fprintf('Skipping NEG LAG NEG CORR mean plot (empty group).\n');
end

%% ==================== DELTA PLOTS: (CSM - NoCSM) ====================

% 3) Δ NEG LAG, POS CORR
if ~isempty(PSTH_poscorr_neglag_CSM) && ~isempty(PSTH_poscorr_neglag_NoCSM)

    [mdelta, semdelta] = delta_sem(PSTH_poscorr_neglag_CSM, PSTH_poscorr_neglag_NoCSM);

    figure('Name','Δ NEG LAG POS CORR','Position',[800,100,600,300]); hold on;
    plot(time_vector, mdelta, 'Color', AltDeltaColor, 'LineWidth', 2);
    fill([time_vector fliplr(time_vector)], ...
         [mdelta+semdelta fliplr(mdelta-semdelta)], ...
         AltDeltaColor,'FaceAlpha',0.2,'EdgeColor','none');

    xline(0,'--k'); yline(0,'--k');
    yline(-2,':'); yline(2,':');
    %xlabel('Time from pos. thresh (ms)'); 
    %xlabel('Time from csmT2 (ms)'); 
    xlabel('Time from csmT0 (ms)'); 
    ylabel('Δ firing rate (CSM - NoCSM)');
    %xlim([-250 250]);
    title('Δ Negative Lag, Positive Correlation');
    box off;

else
    fprintf('Skipping Δ NEG LAG POS CORR plot (empty group).\n');
end

% 4) Δ NEG LAG, NEG CORR
if ~isempty(PSTH_negcorr_neglag_CSM) && ~isempty(PSTH_negcorr_neglag_NoCSM)

    [mdelta, semdelta] = delta_sem(PSTH_negcorr_neglag_CSM, PSTH_negcorr_neglag_NoCSM);

    figure('Name','Δ NEG LAG NEG CORR','Position',[800,450,600,300]); hold on;
    plot(time_vector, mdelta, 'Color', DeltaColor, 'LineWidth', 2);
    fill([time_vector fliplr(time_vector)], ...
         [mdelta+semdelta fliplr(mdelta-semdelta)], ...
         DeltaColor,'FaceAlpha',0.2,'EdgeColor','none');

    xline(0,'--k'); yline(0,'--k');
    yline(-2,':'); yline(2,':');
    %xlabel('Time from pos. thresh (ms)'); 
    %xlabel('Time from csmT2 (ms)'); 
    xlabel('Time from csmT0 (ms)'); 
    ylabel('Δ firing rate (CSM - NoCSM)');
    %xlim([-250 250]);
    title('Δ Negative Lag, Negative Correlation');
    box off;

else
    fprintf('Skipping Δ NEG LAG NEG CORR plot (empty group).\n');
end


%% ==================== DELTA PLOTS: (CSM - NoCSM) ON SAME FIGURE ====================

figure('Name','Δ NEG LAG COMPARISON','Position',[800,200,700,350]); 
hold on;

plotted_any = false;

% ----- Δ NEG LAG, POS CORR -----
if ~isempty(PSTH_poscorr_neglag_CSM) && ~isempty(PSTH_poscorr_neglag_NoCSM)

    [mdelta_pos, semdelta_pos] = delta_sem(PSTH_poscorr_neglag_CSM, PSTH_poscorr_neglag_NoCSM);

    plot(time_vector, mdelta_pos, 'Color', AltDeltaColor, 'LineWidth', 2);
    fill([time_vector fliplr(time_vector)], ...
         [mdelta_pos+semdelta_pos fliplr(mdelta_pos-semdelta_pos)], ...
         AltDeltaColor,'FaceAlpha',0.2,'EdgeColor','none');

    plotted_any = true;

else
    fprintf('Skipping Δ NEG LAG POS CORR plot (empty group).\n');
end

% ----- Δ NEG LAG, NEG CORR -----
if ~isempty(PSTH_negcorr_neglag_CSM) && ~isempty(PSTH_negcorr_neglag_NoCSM)

    [mdelta_neg, semdelta_neg] = delta_sem(PSTH_negcorr_neglag_CSM, PSTH_negcorr_neglag_NoCSM);

    plot(time_vector, mdelta_neg, 'Color', DeltaColor, 'LineWidth', 2);
    fill([time_vector fliplr(time_vector)], ...
         [mdelta_neg+semdelta_neg fliplr(mdelta_neg-semdelta_neg)], ...
         DeltaColor,'FaceAlpha',0.2,'EdgeColor','none');

    plotted_any = true;

else
    fprintf('Skipping Δ NEG LAG NEG CORR plot (empty group).\n');
end

% ----- Finalize only if something was plotted -----
if plotted_any
    xline(0,'--k'); 
    yline(0,'--k');
    yline(-2,':'); 
    yline(2,':');

    xline(-55)
    xline(87)

    %xlabel('Time from csmT0 (ms)');
    %xlabel('Time from csmT2 (ms)');
    xlabel('Time from threshold crossing(ms)');
    ylabel('Δ firing rate (CSM - NoCSM)');
    %xlim([-250 250]);
    ylim([-3 4])


    title('Δ Negative Lag: Pos Corr vs Neg Corr');
    legend({'NegLag PosCorr','','NegLag NegCorr',''}, 'Location','best');

    box off;
else
    close(gcf);
    fprintf('No delta plots generated (both groups empty).\n');
end

%% Now group them togeher 
%% ============================================================
% PSTH plotting code UPDATED for CURRENT analysis (NO before/after)
% + NEW: ALL-KEPT grouped PSTH (posCorr + negCorr together)
% + NEW: ALL-KEPT DELTA (CSM - NoCSM)
% ============================================================

% -------------------- GET MASTER PSTHs --------------------
PSTH_poscorr_neglag_CSM   = AllGroups.negLag_posCorr.CSM;
PSTH_poscorr_neglag_NoCSM = AllGroups.negLag_posCorr.NoCSM;

PSTH_negcorr_neglag_CSM   = AllGroups.negLag_negCorr.CSM;
PSTH_negcorr_neglag_NoCSM = AllGroups.negLag_negCorr.NoCSM;

% -------------------- NEW: ALL KEPT (BOTH CATS TOGETHER) --------------------
PSTH_allkept_CSM   = vertcat(PSTH_poscorr_neglag_CSM,   PSTH_negcorr_neglag_CSM);
PSTH_allkept_NoCSM = vertcat(PSTH_poscorr_neglag_NoCSM, PSTH_negcorr_neglag_NoCSM);

% -------------------- SETUP --------------------
Color_CSM   = rgb('DarkMagenta');
Color_NoCSM = rgb('Black');

DeltaColor     = rgb('DarkGreen');   % negCorr delta
AltDeltaColor  = rgb('SteelBlue');   % posCorr delta
DeltaColor_All = rgb('Gray');        % all-kept delta (neutral)

% robust timevec grab:
if exist('timevec','var') && ~isempty(timevec)
    full_timevec = timevec;
else
    found = false;
    for d = 1:numel(XCorrData)
        if isfield(XCorrData(d),'units') && ~isempty(XCorrData(d).units)
            u0 = XCorrData(d).units(1).unitIdx;
            full_timevec = data.units(u0).AllCSM_pos_thresh.time_bins; %#ok<NASGU>
            found = true;
            break
        end
    end
    if ~found
        error('Could not determine PSTH time vector (timevec/full_timevec not found).');
    end
end

% Indices for plotting window
time_idx    = find(full_timevec >= -500 & full_timevec <= 500);
%time_idx    = find(full_timevec >= -250 & full_timevec <= 250);
time_vector = full_timevec(time_idx);

% Helper for mean + SEM
mean_sem = @(mat) deal( ...
    nanmean(mat(:,time_idx),1), ...
    nanstd(mat(:,time_idx),0,1) ./ sqrt(max(1,size(mat,1))) );

% Helper for delta mean + SEM
delta_sem = @(csm,no) deal( ...
    nanmean(csm(:,time_idx) - no(:,time_idx),1), ...
    nanstd(csm(:,time_idx) - no(:,time_idx),0,1) ./ sqrt(max(1,size(csm,1))) );

%% ==================== NEW: ALL-KEPT MEAN +/- SEM ====================
if ~isempty(PSTH_allkept_CSM) && ~isempty(PSTH_allkept_NoCSM)

    [mCSM_all, semCSM_all] = mean_sem(PSTH_allkept_CSM);
    [mNo_all,  semNo_all ] = mean_sem(PSTH_allkept_NoCSM);

    figure('Name','ALL KEPT (POS+NEG CORR)','Position',[100,760,650,300]); hold on;

    plot(time_vector, mCSM_all,'Color',Color_CSM,'LineWidth',2);
    fill([time_vector fliplr(time_vector)], ...
         [mCSM_all+semCSM_all fliplr(mCSM_all-semCSM_all)], ...
         Color_CSM,'FaceAlpha',0.2,'EdgeColor','none');

    plot(time_vector, mNo_all,'Color',Color_NoCSM,'LineWidth',2);
    fill([time_vector fliplr(time_vector)], ...
         [mNo_all+semNo_all fliplr(mNo_all-semNo_all)], ...
         Color_NoCSM,'FaceAlpha',0.2,'EdgeColor','none');

    xline(0,'--k');
    xlabel('Time (ms)'); ylabel('Z-scored firing rate');
    title(sprintf('ALL KEPT units (N=%d): CSM vs NoCSM', size(PSTH_allkept_CSM,1)));
    box off;

else
    fprintf('Skipping ALL-KEPT mean plot (empty).\n');
end

%% ==================== NEW: ALL-KEPT DELTA (CSM - NoCSM) ====================
if ~isempty(PSTH_allkept_CSM) && ~isempty(PSTH_allkept_NoCSM)

    [mdelta_all, semdelta_all] = delta_sem(PSTH_allkept_CSM, PSTH_allkept_NoCSM);

    figure('Name','Δ ALL KEPT (CSM - NoCSM)','Position',[600,400,650,300]); hold on;

    plot(time_vector, mdelta_all, 'Color', DeltaColor_All, 'LineWidth', 2);
    fill([time_vector fliplr(time_vector)], ...
         [mdelta_all+semdelta_all fliplr(mdelta_all-semdelta_all)], ...
         DeltaColor_All,'FaceAlpha',0.25,'EdgeColor','none');

    xline(0,'--k'); yline(0,'--k');
    yline(-2,':'); yline(2,':');
    xlim([-250 250])
    ylim([-3 3])
    xlabel('Time (ms)'); ylabel('Δ firing rate (CSM - NoCSM)');
    title(sprintf('Δ ALL KEPT units (N=%d)', size(PSTH_allkept_CSM,1)));
    box off;

else
    fprintf('Skipping Δ ALL-KEPT plot (empty).\n');
end

%% ==================== PLOTTING: MEAN +/- SEM (PER CATEGORY) ====================

% 1) NEG LAG, POS CORR
if ~isempty(PSTH_poscorr_neglag_CSM) && ~isempty(PSTH_poscorr_neglag_NoCSM)

    [mCSM, semCSM] = mean_sem(PSTH_poscorr_neglag_CSM);
    [mNo,  semNo ] = mean_sem(PSTH_poscorr_neglag_NoCSM);

    figure('Name','NEG LAG, POS CORR','Position',[100,100,600,300]); hold on;
    plot(time_vector, mCSM,'Color',Color_CSM,'LineWidth',2);
    fill([time_vector fliplr(time_vector)], ...
         [mCSM+semCSM fliplr(mCSM-semCSM)], ...
         Color_CSM,'FaceAlpha',0.2,'EdgeColor','none');

    plot(time_vector, mNo,'Color',Color_NoCSM,'LineWidth',2);
    fill([time_vector fliplr(time_vector)], ...
         [mNo+semNo fliplr(mNo-semNo)], ...
         Color_NoCSM,'FaceAlpha',0.2,'EdgeColor','none');

    xline(0,'--k');
    xlabel('Time from pos. thresh (ms)'); ylabel('Z-scored firing rate');
    title('Negative Lag, Positive Correlation');
    box off;

else
    fprintf('Skipping NEG LAG POS CORR mean plot (empty group).\n');
end

% 2) NEG LAG, NEG CORR
if ~isempty(PSTH_negcorr_neglag_CSM) && ~isempty(PSTH_negcorr_neglag_NoCSM)

    [mCSM, semCSM] = mean_sem(PSTH_negcorr_neglag_CSM);
    [mNo,  semNo ] = mean_sem(PSTH_negcorr_neglag_NoCSM);

    figure('Name','NEG LAG, NEG CORR','Position',[100,450,600,300]); hold on;
    plot(time_vector, mCSM,'Color',Color_CSM,'LineWidth',2);
    fill([time_vector fliplr(time_vector)], ...
         [mCSM+semCSM fliplr(mCSM-semCSM)], ...
         Color_CSM,'FaceAlpha',0.2,'EdgeColor','none');

    plot(time_vector, mNo,'Color',Color_NoCSM,'LineWidth',2);
    fill([time_vector fliplr(time_vector)], ...
         [mNo+semNo fliplr(mNo-semNo)], ...
         Color_NoCSM,'FaceAlpha',0.2,'EdgeColor','none');

    xline(0,'--k');
    xlabel('Time from pos. thresh (ms)'); ylabel('Z-scored firing rate');
    title('Negative Lag, Negative Correlation');
    box off;

else
    fprintf('Skipping NEG LAG NEG CORR mean plot (empty group).\n');
end

%% ==================== DELTA PLOTS: (CSM - NoCSM) ====================

% 3) Δ NEG LAG, POS CORR
if ~isempty(PSTH_poscorr_neglag_CSM) && ~isempty(PSTH_poscorr_neglag_NoCSM)

    [mdelta, semdelta] = delta_sem(PSTH_poscorr_neglag_CSM, PSTH_poscorr_neglag_NoCSM);

    figure('Name','Δ NEG LAG POS CORR','Position',[800,100,600,300]); hold on;
    plot(time_vector, mdelta, 'Color', AltDeltaColor, 'LineWidth', 2);
    fill([time_vector fliplr(time_vector)], ...
         [mdelta+semdelta fliplr(mdelta-semdelta)], ...
         AltDeltaColor,'FaceAlpha',0.2,'EdgeColor','none');

    xline(0,'--k'); yline(0,'--k');
    yline(-2,':'); yline(2,':');
    xlabel('Time (ms)');
    ylabel('Δ firing rate (CSM - NoCSM)');
    title('Δ Negative Lag, Positive Correlation');
    box off;

else
    fprintf('Skipping Δ NEG LAG POS CORR plot (empty group).\n');
end

% 4) Δ NEG LAG, NEG CORR
if ~isempty(PSTH_negcorr_neglag_CSM) && ~isempty(PSTH_negcorr_neglag_NoCSM)

    [mdelta, semdelta] = delta_sem(PSTH_negcorr_neglag_CSM, PSTH_negcorr_neglag_NoCSM);

    figure('Name','Δ NEG LAG NEG CORR','Position',[800,450,600,300]); hold on;
    plot(time_vector, mdelta, 'Color', DeltaColor, 'LineWidth', 2);
    fill([time_vector fliplr(time_vector)], ...
         [mdelta+semdelta fliplr(mdelta-semdelta)], ...
         DeltaColor,'FaceAlpha',0.2,'EdgeColor','none');

    xline(0,'--k'); yline(0,'--k');
    yline(-2,':'); yline(2,':');
    xlabel('Time (ms)');
    ylabel('Δ firing rate (CSM - NoCSM)');
    title('Δ Negative Lag, Negative Correlation');
    box off;

else
    fprintf('Skipping Δ NEG LAG NEG CORR plot (empty group).\n');
end
