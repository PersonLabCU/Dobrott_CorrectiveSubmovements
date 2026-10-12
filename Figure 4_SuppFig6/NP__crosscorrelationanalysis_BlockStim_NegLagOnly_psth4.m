

% ============================================================
% integrated_BlockStim_NEG_LAG_xcorr__CSMmod_DELTAZshuffle__PSTHavg.m
%
% BLOCK STIM SESSIONS (*Block_CSMvsNoCSM.mat)
%
% GOALS (your spec):
%   (1) XCORR categorization (NEGATIVE-LAG ONLY):
%         - choose lag < 0 with maximum |medianCorr|
%         - category = negLag_posCorr or negLag_negCorr (sign of chosen corr)
%
%   (2) CSM modulation classification:
%         - sustained deltaZ (z-score) + UNIT-LEVEL TRIAL-LABEL SHUFFLE p-value
%         - store per-unit: sig_any / sig_inc / sig_dec / sig_bi / p_unit etc.
%
%   (3) Average PSTH outputs:
%         - for ALL KEPT units (peak within [-250 250] ms around max_vx)
%         - for KEPT units that are significantly modulated during CSM (sig_any)
%         - (optionally also split by the two negLag categories)
%
% NOTES:
%   - “baseline” vs “stimulated” here means: baseline_psth_field vs stim_psth_field
%   - z-scored PSTHs use baseline trials’ baseline window (baseline_win_def)
%   - peak-gate uses peakgate_field PSTH (typically BLCSM_max_vx)
% ============================================================

clear; clc;

% ---------------------- LOAD SESSIONS -----------------------
cd("C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\StimAnimals\BlockSeparated_CSMvsNoCSM_inducedalignment");
MyFolderInfo = dir('*Block_CSMvsNoCSM.mat');

%% ---------------------- USER SETTINGS -----------------------
% ----- XCORR settings -----
rasterField   = 'max_vx';
kinVar        = 7;
gauss_sigma   = 10;
resample_fs   = 200;
win_pre_ms    = 100;
win_post_ms   = 150;
maxLag_ms     = 100;
lag_step      = 5;

% ----- PSTH peak-gate + zscore baseline window -----
valid_peak_range_ms  = [-250 250];     % KEEP if peak time (peakgate_field) is within this


baseline_win_def     = [-1000 -500];    % baseline window for PSTH z-scoring (ms)
plot_win             = [-500 500];     % for plotting later if you want

% ----- Choose fields (BASELINE vs STIMULATED) -----
% Peak gate usually done on baseline (BL) max_vx
peakgate_field       = 'BLCSM_max_vx';          % for keep gate + xcorr raster source

% % PSTHs you want to average (baseline vs stimulated)
% baseline_psth_field  = 'BLCSM_csmT2';      % BASELINE PSTH field
% stim_psth_field      = 'ADCSM_csmT2';      % STIMULATED PSTH field
% wash_psth_field      = 'WOECSM_csmT2';      % STIMULATED PSTH field


% PSTHs you want to average (baseline vs stimulated)
baseline_psth_field  = 'BLCSM_pos_thresh';      % BASELINE PSTH field
stim_psth_field      = 'ADCSM_pos_thresh';      % STIMULATED PSTH field
wash_psth_field      = 'WOECSM_pos_thresh';      % STIMULATED PSTH field
washall_psth_field      = 'WOCSM_pos_thresh';      % STIMULATED PSTH field
norm_psth_field      = 'Ctrl_pos_thresh';
% e.g. swap to WOE:
% stim_psth_field    = 'WOECSM_pos_thresh';

% ----- CSM modulation (deltaZ sustained + shuffle) -----
csm_epoch_win          = [-50 100];   % ms
csm_baseline_win       = [-500 -250]; % ms (baseline for pooled baseline noise)
min_trials_per_cond    = 5;
min_spikes_total_epoch = 5;
min_spikes_per_cond    = 3;

z_thresh        = 2;
min_consec_bins = 2;

min_base_std    = 1;    % Hz SD floor
bin_ms          = 10;
step_ms         = 5;

do_unit_shuffle  = true;
nPerm_unit       = 1000;
alpha_unit       = 0.05;

print_unit_qc    = false;
print_session_qc = true;

% ---------------------- MASTER STRUCTS ----------------------
XCorrData = struct;
NPdata    = struct;

group_names = {'negLag_posCorr','negLag_negCorr'};

% PSTH aggregation containers (ALL kept; and MOD kept)
Agg = struct;
Agg.AllKept.baseline = [];   % [nUnits x nTime]
Agg.AllKept.stim     = [];
Agg.AllKept.wash     = [];
Agg.ModKept.baseline = [];
Agg.ModKept.stim     = [];
Agg.ModKept.wash     = [];

% Optional: also store by NEG-LAG category
for g = 1:numel(group_names)
    gn = group_names{g};
    Agg.ByGroup.(gn).AllKept.baseline = [];
    Agg.ByGroup.(gn).AllKept.stim     = [];
    Agg.ByGroup.(gn).AllKept.wash     = [];
    Agg.ByGroup.(gn).ModKept.baseline = [];
    Agg.ByGroup.(gn).ModKept.stim     = [];
    Agg.ByGroup.(gn).ModKept.wash     = [];
end

%% ====================== MAIN SESSION LOOP ======================
for day = 1:numel(MyFolderInfo)

    fprintf('\n========== SESSION %d / %d ==========\n', day, numel(MyFolderInfo));
    load(MyFolderInfo(day).name, ...
        'Reaches','data','ExcitRR_BL','ExcitRR_AD','ExcitRR_WOE','ExcitRR_WO','Salt_index');

    if ~exist('Salt_index','var') || isempty(Salt_index)
        fprintf(' Day %d: Salt_index empty or missing. Skipping.\n', day);
        continue;
    end

    % Optional: keep kinematics structs (if you need later)
    NPdata(day).BLCSM   = ExcitRR_BL.ReReaches;
    NPdata(day).BLNoCSM = ExcitRR_BL.NoReReach;
    NPdata(day).ADCSM   = ExcitRR_AD.ReReaches;
    NPdata(day).WOCSM   = ExcitRR_WO.ReReaches;
    NPdata(day).WOECSM  = ExcitRR_WOE.ReReaches;

    session_units = Salt_index(:);
    nUnits = numel(session_units);

    % time vector (from peakgate field)
    timevec = data.units(session_units(1)).(peakgate_field).time_bins; % ms
    baseline_win     = find(timevec >= baseline_win_def(1)    & timevec <= baseline_win_def(2));
    valid_peak_range = find(timevec >= valid_peak_range_ms(1) & timevec <= valid_peak_range_ms(2));

    % session bookkeeping
    kept_units_session = [];
    mod_units_session  = [];
    negLag_posCorr_units = [];
    negLag_negCorr_units = [];

    % session skip counters (mod test)
    skip = struct;
    skip.peak_gate_fail   = 0;
    skip.missing_rasters  = 0;
    skip.empty_windows    = 0;
    skip.min_trials       = 0;
    skip.min_spikes       = 0;
    skip.ok               = 0;

    for uu = 1:nUnits

        unitIdx = session_units(uu);
        fprintf('  Analyzing UNIT %d...\n', unitIdx);

        % ---------------- KEEP GATE (peak within valid window) ----------------
        if ~isfield(data.units(unitIdx), peakgate_field) || ~isfield(data.units(unitIdx).(peakgate_field),'PSTH_smoothbin')
            XCorrData(day).units(uu).unitIdx = unitIdx;
            XCorrData(day).units(uu).keep   = false;
            skip.peak_gate_fail = skip.peak_gate_fail + 1;
            continue;
        end

        psth_peakgate = data.units(unitIdx).(peakgate_field).PSTH_smoothbin;
        [~, peak_idx] = max(psth_peakgate);

        keep = ismember(peak_idx, valid_peak_range);

        XCorrData(day).units(uu).unitIdx        = unitIdx;
        XCorrData(day).units(uu).keep           = keep;
        XCorrData(day).units(uu).psth_peak_idx  = peak_idx;
        XCorrData(day).units(uu).psth_peak_time = timevec(peak_idx);

        if ~keep
            skip.peak_gate_fail = skip.peak_gate_fail + 1;
            continue;
        end

        kept_units_session(end+1,1) = unitIdx;

        % ---------------- XCORR (NEGATIVE-LAG ONLY selection) ----------------
        % Load raster for IFR (use peakgate_field Raster; typically BLCSM_max_vx)
        if ~isfield(data.units(unitIdx).(peakgate_field),'Raster')
            XCorrData(day).units(uu).negLag_peakLag = NaN;
            XCorrData(day).units(uu).negLag_peakVal = NaN;
            XCorrData(day).units(uu).category       = 'negLag_unknown';
        else
            spikes_full = data.units(unitIdx).(peakgate_field).Raster;  % trials x time (1ms)
            [nTrials, nTime] = size(spikes_full);

            full_time_1ms = (-1000:1000)';
            if nTime ~= numel(full_time_1ms)
                if nTime < numel(full_time_1ms)
                    full_time_1ms = full_time_1ms(1:nTime);
                else
                    full_time_1ms = linspace(-1000,1000,nTime)';
                end
            end

            % IFR
            IFR_full = zeros(nTime, nTrials);
            kern_len = max(1, round(6*gauss_sigma));
            gk = fspecial('gaussian',[1 kern_len], gauss_sigma);
            gk = gk / sum(gk);

            for tr = 1:nTrials
                sp_idx = find(spikes_full(tr,:) > 0);
                if numel(sp_idx) < 2, continue; end

                spike_times = full_time_1ms(sp_idx);
                isis = diff(spike_times)/1000;
                inst_rate = 1 ./ isis;

                rate_vec = zeros(nTime,1);
                for s = 1:numel(isis)
                    interval_idx = sp_idx(s):sp_idx(s+1)-1;
                    rate_vec(interval_idx) = inst_rate(s);
                end

                IFR_full(:,tr) = conv(rate_vec, gk, 'same');
            end

            win_idx = full_time_1ms >= -win_pre_ms & full_time_1ms <= win_post_ms;
            IFR_1ms = IFR_full(win_idx,:);
            ifr_time_1ms = full_time_1ms(win_idx);

            time_200Hz   = (-win_pre_ms:1000/resample_fs:win_post_ms)'; % ms
            nSamples200  = numel(time_200Hz);
            IFR = nan(nSamples200, nTrials);
            for tr = 1:nTrials
                IFR(:,tr) = interp1(ifr_time_1ms, IFR_1ms(:,tr), time_200Hz, 'linear', 0);
            end

            % Kinematics (baseline reaches for xcorr: BL)
            Reaches_use = ExcitRR_BL.ReReaches;
            C = AlignPeakV_1ms_edit(Reaches_use, kinVar, win_pre_ms, win_post_ms);
            C = stimvel_accel(C);

            vel_1ms_all  = squeeze(C.stimvel(:,5,:)); % velocity
            vel_time_1ms = linspace(-win_pre_ms, win_post_ms, size(vel_1ms_all,1))';

            nTrials_vel = size(vel_1ms_all,2);
            nTrials_use = min(size(IFR,2), nTrials_vel);
            IFR = IFR(:,1:nTrials_use);

            Vel = nan(nSamples200, nTrials_use);
            for tr = 1:nTrials_use
                Vel(:,tr) = interp1(vel_time_1ms, vel_1ms_all(:,tr), time_200Hz, 'linear', NaN);
            end

            lags_ms = (-maxLag_ms:lag_step:maxLag_ms)';
            nLags   = numel(lags_ms);
            xcorr_trials = nan(nLags, nTrials_use);

            for li = 1:nLags
                lag   = lags_ms(li);
                shift = round(lag / (1000/resample_fs));

                for tr = 1:nTrials_use
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

                    if numel(fr_seg) < 5, continue; end
                    if all(isnan(fr_seg)) || all(isnan(v_seg)), continue; end

                    fr_sd = std(fr_seg,'omitnan');
                    v_sd  = std(v_seg, 'omitnan');
                    if fr_sd < 1e-9 || v_sd < 1e-9, continue; end

                    fr_z = (fr_seg - mean(fr_seg,'omitnan')) ./ fr_sd;
                    v_z  = (v_seg  - mean(v_seg, 'omitnan')) ./ v_sd;

                    if any(isnan(fr_z)) || any(isnan(v_z)), continue; end
                    xcorr_trials(li,tr) = mean(fr_z .* v_z, 'omitnan');
                end
            end

            medianCorr = nanmedian(xcorr_trials,2);

            % NEGATIVE-LAG ONLY choice (max |medianCorr| among lags<0) + quadratic interpolation
            negMask = (lags_ms < 0);
            negIdx  = find(negMask);

            negLag_peakLag = NaN;
            negLag_peakVal = NaN;
            cat_label      = 'negLag_unknown';

            if ~isempty(negIdx) && any(~isnan(medianCorr(negIdx)))
                [~, kRel] = max(abs(medianCorr(negIdx)));
                idxStar = negIdx(kRel);

                if idxStar > 1 && idxStar < numel(medianCorr) && all(~isnan(medianCorr(idxStar-1:idxStar+1)))
                    x = lags_ms(idxStar-1:idxStar+1);
                    y = medianCorr(idxStar-1:idxStar+1);
                    p = polyfit(x, y, 2);
                    lag_hat = -p(2)/(2*p(1));
                    lag_hat = max(min(lag_hat, max(x)), min(x));
                    val_hat = polyval(p, lag_hat);
                    negLag_peakLag = lag_hat;
                    negLag_peakVal = val_hat;
                else
                    negLag_peakLag = lags_ms(idxStar);
                    negLag_peakVal = medianCorr(idxStar);
                end

                if negLag_peakVal >= 0
                    cat_label = 'negLag_posCorr';
                    negLag_posCorr_units(end+1,1) = unitIdx;
                else
                    cat_label = 'negLag_negCorr';
                    negLag_negCorr_units(end+1,1) = unitIdx;
                end
            end

            XCorrData(day).units(uu).lags_ms        = lags_ms;
            XCorrData(day).units(uu).medianCorr     = medianCorr;
            XCorrData(day).units(uu).xcorr_trials   = xcorr_trials;
            XCorrData(day).units(uu).negLag_peakLag = negLag_peakLag;
            XCorrData(day).units(uu).negLag_peakVal = negLag_peakVal;
            XCorrData(day).units(uu).category       = cat_label;
        end

        % % ---------------- PSTH Z-SCORES (baseline vs stimulated) ----------------
        % % Pull PSTHs (must exist)
        % if ~isfield(data.units(unitIdx), baseline_psth_field) || ~isfield(data.units(unitIdx).(baseline_psth_field),'PSTH_smoothbin') || ...
        %    ~isfield(data.units(unitIdx), stim_psth_field)     || ~isfield(data.units(unitIdx).(stim_psth_field),    'PSTH_smoothbin') || ...
        %    ~isfield(data.units(unitIdx), wash_psth_field)    || ~isfield(data.units(unitIdx).(wash_psth_field),    'PSTH_smoothbin')
        %     XCorrData(day).units(uu).PSTH = struct('skip',true,'reason','missing_psth_fields');
        % else
        %     psth_base = data.units(unitIdx).(baseline_psth_field).PSTH_smoothbin;
        %     psth_stim = data.units(unitIdx).(stim_psth_field).PSTH_smoothbin;
        %     psth_wash = data.units(unitIdx).(wash_psth_field).PSTH_smoothbin;
        % 
        % 
        %     mu0 = mean(psth_base(baseline_win),'omitnan');
        %     sd0 = std(psth_base(baseline_win),0,'omitnan');
        %     if ~isfinite(sd0) || sd0 == 0, sd0 = 1; end
        % 
        %     base_z = (psth_base - mu0) ./ sd0;
        %     stim_z = (psth_stim - mu0) ./ sd0;
        %     wash_z = (psth_wash - mu0) ./ sd0;
        % 
        %     XCorrData(day).units(uu).PSTH.skip   = false;
        %     XCorrData(day).units(uu).PSTH.base_z = base_z(:)';
        %     XCorrData(day).units(uu).PSTH.stim_z = stim_z(:)';
        %     XCorrData(day).units(uu).PSTH.wash_z = wash_z(:)';
        % 
        %     % global aggregation: ALL kept
        %     Agg.AllKept.baseline = vertcat(Agg.AllKept.baseline, base_z(:)');
        %     Agg.AllKept.stim     = vertcat(Agg.AllKept.stim,     stim_z(:)');
        %     Agg.AllKept.wash     = vertcat(Agg.AllKept.wash,     wash_z(:)');
        % 
        %     % optional: by negLag category (if assigned)
        %     if isfield(XCorrData(day).units(uu),'category')
        %         gn = XCorrData(day).units(uu).category;
        %         if ismember(gn, group_names)
        %             Agg.ByGroup.(gn).AllKept.baseline = vertcat(Agg.ByGroup.(gn).AllKept.baseline, base_z(:)');
        %             Agg.ByGroup.(gn).AllKept.stim     = vertcat(Agg.ByGroup.(gn).AllKept.stim,     stim_z(:)');
        %             Agg.ByGroup.(gn).AllKept.wash     = vertcat(Agg.ByGroup.(gn).AllKept.wash,     wash_z(:)');
        %         end
        %     end
        % end

        % ---------------- PSTH Z-SCORES (baseline vs stimulated) ----------------
        % Only change here: use PSTH_smooth instead of PSTH_smoothbin
        % Assumes PSTH_smooth is 1x1971 spanning -1000 to 1000 ms

        % if ~isfield(data.units(unitIdx), baseline_psth_field) || ~isfield(data.units(unitIdx).(baseline_psth_field),'PSTH_smooth') || ...
        %         ~isfield(data.units(unitIdx), stim_psth_field)     || ~isfield(data.units(unitIdx).(stim_psth_field),    'PSTH_smooth') || ...
        %         ~isfield(data.units(unitIdx), wash_psth_field)     || ~isfield(data.units(unitIdx).(wash_psth_field),    'PSTH_smooth')
        % 
        %     XCorrData(day).units(uu).PSTH = struct('skip',true,'reason','missing_psth_fields');
        % 
        % else
        %     psth_norm=  data.units(unitIdx).(norm_psth_field).PSTH_smooth;
        %     psth_base = data.units(unitIdx).(baseline_psth_field).PSTH_smooth;
        %     psth_stim = data.units(unitIdx).(stim_psth_field).PSTH_smooth;
        %     psth_wash = data.units(unitIdx).(wash_psth_field).PSTH_smooth;
        %     psth_washall= data.units(unitIdx).(washall_psth_field).PSTH_smooth;

        %Recomputing PSTH
        raster_norm    = data.units(unitIdx).(norm_psth_field).Raster;
        raster_base    = data.units(unitIdx).(baseline_psth_field).Raster;
        raster_stim    = data.units(unitIdx).(stim_psth_field).Raster;
        raster_wash    = data.units(unitIdx).(wash_psth_field).Raster;
        raster_washall = data.units(unitIdx).(washall_psth_field).Raster;

        bin_size = 10;
        smoothing_bins = 3; % originally 3 

        % --- compute PSTHs ---
        [psth_raw_norm,    psth_smooth_norm,    time_bins] = computePSTH(raster_norm,    bin_size, smoothing_bins);
        [psth_raw_base,    psth_smooth_base,    ~]         = computePSTH(raster_base,    bin_size, smoothing_bins);
        [psth_raw_stim,    psth_smooth_stim,    ~]         = computePSTH(raster_stim,    bin_size, smoothing_bins);
        [psth_raw_wash,    psth_smooth_wash,    ~]         = computePSTH(raster_wash,    bin_size, smoothing_bins);
        [psth_raw_washall, psth_smooth_washall, ~]         = computePSTH(raster_washall, bin_size, smoothing_bins);

        % --- row vectors ---
        psth_norm    = psth_smooth_norm(:)';
        psth_base    = psth_smooth_base(:)';
        psth_stim    = psth_smooth_stim(:)';
        psth_wash    = psth_smooth_wash(:)';
        psth_washall = psth_smooth_washall(:)';

        % --- use returned time bins directly ---
        psth_timevec = time_bins(:)';

        % baseline window for z-scoring
        baseline_win_psth = find(psth_timevec >= baseline_win_def(1) & psth_timevec <= baseline_win_def(2));

        if isempty(baseline_win_psth)
            error('Baseline window does not overlap PSTH time vector.');
        end

        % --- z-score relative to BASELINE PSTH baseline window ---
        mu0 = mean(psth_base(baseline_win_psth), 'omitnan');
        sd0 = std(psth_base(baseline_win_psth), 0, 'omitnan');
        if ~isfinite(sd0) || sd0 == 0
            sd0 = 1;
        end

        % Optional alternative references (kept here if you want them later)
        mu0wash = mean(psth_washall(baseline_win_psth), 'omitnan');
        sd0wash = std(psth_washall(baseline_win_psth), 0, 'omitnan');
        if ~isfinite(sd0wash) || sd0wash == 0
            sd0wash = 1;
        end

        mu0stim = mean(psth_stim(baseline_win_psth), 'omitnan');
        sd0stim = std(psth_stim(baseline_win_psth), 0, 'omitnan');
        if ~isfinite(sd0stim) || sd0stim == 0
            sd0stim = 1;
        end

        base_z = (psth_base - mu0) ./ sd0;
        stim_z = (psth_stim - mu0) ./ sd0;
        wash_z = (psth_wash - mu0) ./ sd0;

        % base_z =psth_base ;
        % stim_z = psth_stim;
        % wash_z = psth_wash ;


        % If instead you want condition-specific z-scoring, use these:
        % stim_z = (psth_stim - mu0stim) ./ sd0stim;
        % wash_z = (psth_wash - mu0wash) ./ sd0wash;

        XCorrData(day).units(uu).PSTH.skip           = false;
        XCorrData(day).units(uu).PSTH.timevec        = psth_timevec;
        XCorrData(day).units(uu).PSTH.bin_size       = bin_size;
        XCorrData(day).units(uu).PSTH.smoothing_bins = smoothing_bins;

        XCorrData(day).units(uu).PSTH.base_raw       = psth_raw_base(:)';
        XCorrData(day).units(uu).PSTH.stim_raw       = psth_raw_stim(:)';
        XCorrData(day).units(uu).PSTH.wash_raw       = psth_raw_wash(:)';

        XCorrData(day).units(uu).PSTH.base_smooth    = psth_base;
        XCorrData(day).units(uu).PSTH.stim_smooth    = psth_stim;
        XCorrData(day).units(uu).PSTH.wash_smooth    = psth_wash;

        XCorrData(day).units(uu).PSTH.base_z         = base_z(:)';
        XCorrData(day).units(uu).PSTH.stim_z         = stim_z(:)';
        XCorrData(day).units(uu).PSTH.wash_z         = wash_z(:)';

        % global aggregation: ALL kept
        Agg.AllKept.baseline = vertcat(Agg.AllKept.baseline, base_z(:)');
        Agg.AllKept.stim     = vertcat(Agg.AllKept.stim,     stim_z(:)');
        Agg.AllKept.wash     = vertcat(Agg.AllKept.wash,     wash_z(:)');

        % optional: by negLag category (if assigned)
        if isfield(XCorrData(day).units(uu), 'category')
            gn = XCorrData(day).units(uu).category;
            if ismember(gn, group_names)
                Agg.ByGroup.(gn).AllKept.baseline = vertcat(Agg.ByGroup.(gn).AllKept.baseline, base_z(:)');
                Agg.ByGroup.(gn).AllKept.stim     = vertcat(Agg.ByGroup.(gn).AllKept.stim,     stim_z(:)');
                Agg.ByGroup.(gn).AllKept.wash     = vertcat(Agg.ByGroup.(gn).AllKept.wash,     wash_z(:)');
            end
        end


        % ---------------- CSM MODULATION (deltaZ sustained + shuffle) ----------------
        % Need Rasters for baseline + stim for deltaZ (these should match your *_csmT2 style fields)
        % If your modulation should be computed on csmT2 regardless of pos_thresh PSTH averaging,
        % set these explicitly here:
        bl_mod_field   = 'BLCSM_csmT2';
        stim_mod_field = 'ADCSM_csmT2';   % or 'WOECSM_csmT2'

        need_bl = isfield(data.units(unitIdx), bl_mod_field)   && isfield(data.units(unitIdx).(bl_mod_field),  'Raster');
        need_st = isfield(data.units(unitIdx), stim_mod_field) && isfield(data.units(unitIdx).(stim_mod_field),'Raster');

        if ~need_bl || ~need_st
            XCorrData(day).units(uu).CSM_zmod = struct('skip',true,'reason','missing_required_rasters');
            skip.missing_rasters = skip.missing_rasters + 1;
        else
            BL_raster   = data.units(unitIdx).(bl_mod_field).Raster;      % trials x time (1ms)
            Stim_raster = data.units(unitIdx).(stim_mod_field).Raster;

            t = (-1000:1000)';  % assumed
            nT = min([size(BL_raster,2), size(Stim_raster,2), numel(t)]);
            t = t(1:nT);
            BL_raster   = BL_raster(:,1:nT);
            Stim_raster = Stim_raster(:,1:nT);

            base_idx = (t >= csm_baseline_win(1)) & (t <= csm_baseline_win(2));
            csm_idx  = (t >= csm_epoch_win(1))    & (t <= csm_epoch_win(2));

            if ~any(base_idx) || ~any(csm_idx)
                XCorrData(day).units(uu).CSM_zmod = struct('skip',true,'reason','empty_base_or_csm_window');
                skip.empty_windows = skip.empty_windows + 1;
            else
                nBL   = size(BL_raster,1);
                nStim = size(Stim_raster,1);

                if (nBL < min_trials_per_cond) || (nStim < min_trials_per_cond)
                    XCorrData(day).units(uu).CSM_zmod = struct('skip',true,'reason','min_trials_per_cond');
                    skip.min_trials = skip.min_trials + 1;
                else
                    BL_spk_epoch_total   = sum(BL_raster(:,csm_idx),   'all');
                    Stim_spk_epoch_total = sum(Stim_raster(:,csm_idx), 'all');

                    if (BL_spk_epoch_total + Stim_spk_epoch_total) < min_spikes_total_epoch || ...
                       (BL_spk_epoch_total < min_spikes_per_cond) || (Stim_spk_epoch_total < min_spikes_per_cond)
                        XCorrData(day).units(uu).CSM_zmod = struct('skip',true,'reason','min_spikes_epoch');
                        skip.min_spikes = skip.min_spikes + 1;
                    else
                        % --- binning setup ---
                        bin_s = bin_ms/1000;
                        bin_starts = csm_epoch_win(1):step_ms:(csm_epoch_win(2)-bin_ms);
                        nBins = numel(bin_starts);

                        base_s = sum(base_idx)/1000;

                        % pooled baseline noise
                        BL_base_FR_trials   = sum(BL_raster(:,base_idx),2)   ./ base_s;
                        Stim_base_FR_trials = sum(Stim_raster(:,base_idx),2) ./ base_s;

                        base_all = [BL_base_FR_trials; Stim_base_FR_trials];
                        mu_base  = mean(base_all, 'omitnan');
                        sd_base  = std(base_all,  0, 'omitnan');
                        if ~isfinite(sd_base) || sd_base < min_base_std
                            sd_base = min_base_std;
                        end

                        FR_BL   = nan(nBins,1);
                        FR_Stim = nan(nBins,1);
                        FR_all_trials = nan(nBL+nStim, nBins);

                        for b = 1:nBins
                            b0 = bin_starts(b);
                            b_idx = (t >= b0) & (t < (b0 + bin_ms));
                            if ~any(b_idx), continue; end

                            BL_bin_FR_trials   = sum(BL_raster(:,b_idx),2)   ./ bin_s;
                            Stim_bin_FR_trials = sum(Stim_raster(:,b_idx),2) ./ bin_s;

                            FR_BL(b)   = mean(BL_bin_FR_trials,   'omitnan');
                            FR_Stim(b) = mean(Stim_bin_FR_trials, 'omitnan');

                            FR_all_trials(:,b) = [Stim_bin_FR_trials; BL_bin_FR_trials]; % [Stim; BL]
                        end

                        zBL   = (FR_BL   - mu_base) ./ sd_base;
                        zStim = (FR_Stim - mu_base) ./ sd_base;
                        deltaZ = zStim - zBL;

                        above = isfinite(deltaZ) & (deltaZ >=  z_thresh);
                        below = isfinite(deltaZ) & (deltaZ <= -z_thresh);

                        % max consecutive WITHOUT defining a function
                        dA = diff([0; above(:); 0]);
                        runS = find(dA==1); runE = find(dA==-1)-1;
                        if isempty(runS), max_inc_run = 0; else, max_inc_run = max(runE-runS+1); end

                        dB = diff([0; below(:); 0]);
                        runS = find(dB==1); runE = find(dB==-1)-1;
                        if isempty(runS), max_dec_run = 0; else, max_dec_run = max(runE-runS+1); end

                        sig_inc_obs = (max_inc_run >= min_consec_bins);
                        sig_dec_obs = (max_dec_run >= min_consec_bins);
                        sig_bi_obs  = sig_inc_obs & sig_dec_obs;

                        run_stat_obs = max(max_inc_run, max_dec_run);

                        dz_peak_pos = max(deltaZ, [], 'omitnan');
                        dz_peak_neg = min(deltaZ, [], 'omitnan');
                        dz_mean     = mean(deltaZ, 'omitnan');

                        if (sig_inc_obs || sig_dec_obs)
                            if abs(dz_peak_pos) >= abs(dz_peak_neg)
                                dom_dir = +1; dom_val = dz_peak_pos;
                            else
                                dom_dir = -1; dom_val = dz_peak_neg;
                            end
                        else
                            dom_dir = 0;
                            dom_val = dz_peak_pos;
                            if ~isfinite(dom_val), dom_val = dz_peak_neg; end
                        end

                        frac_bins_over = mean(isfinite(deltaZ) & (abs(deltaZ) >= z_thresh), 'omitnan');

                        % unit-level shuffle p-value on run_stat
                        labels0 = [true(nStim,1); false(nBL,1)]; % order matches FR_all_trials: [Stim; BL]
                        p_unit = NaN;

                        if do_unit_shuffle
                            null_stats = zeros(nPerm_unit,1);
                            for pp = 1:nPerm_unit
                                permLab = labels0(randperm(numel(labels0)));

                                mStim = mean(FR_all_trials(permLab,:),  1, 'omitnan');
                                mBL   = mean(FR_all_trials(~permLab,:), 1, 'omitnan');

                                zStim_p = (mStim(:) - mu_base) ./ sd_base;
                                zBL_p   = (mBL(:)   - mu_base) ./ sd_base;
                                dZp     = zStim_p - zBL_p;

                                above_p = isfinite(dZp) & (dZp >=  z_thresh);
                                below_p = isfinite(dZp) & (dZp <= -z_thresh);

                                dAp = diff([0; above_p(:); 0]);
                                s = find(dAp==1); e = find(dAp==-1)-1;
                                if isempty(s), r_inc = 0; else, r_inc = max(e-s+1); end

                                dBp = diff([0; below_p(:); 0]);
                                s = find(dBp==1); e = find(dBp==-1)-1;
                                if isempty(s), r_dec = 0; else, r_dec = max(e-s+1); end

                                null_stats(pp) = max(r_inc, r_dec);
                            end

                            p_unit = (1 + sum(null_stats >= run_stat_obs)) / (nPerm_unit + 1);
                            sig_any = (run_stat_obs >= min_consec_bins) && (p_unit < alpha_unit);
                        else
                            sig_any = (run_stat_obs >= min_consec_bins);
                        end

                        sig_inc = sig_inc_obs && sig_any;
                        sig_dec = sig_dec_obs && sig_any;
                        sig_bi  = sig_inc && sig_dec;

                        if print_unit_qc
                            fprintf('    [zQC] sd_base=%.3fHz | frac|dz|>=%.1f=%.2f | runs inc=%d dec=%d | runStat=%d | p=%.4f | sig=%d\n', ...
                                sd_base, z_thresh, frac_bins_over, max_inc_run, max_dec_run, run_stat_obs, p_unit, sig_any);
                        end

                        XCorrData(day).units(uu).CSM_zmod.skip            = false;
                        XCorrData(day).units(uu).CSM_zmod.method          = 'pooled_baseline_noise + sustained_run + unit_shuffle';
                        XCorrData(day).units(uu).CSM_zmod.bl_field        = bl_mod_field;
                        XCorrData(day).units(uu).CSM_zmod.stim_field      = stim_mod_field;

                        XCorrData(day).units(uu).CSM_zmod.z_thresh        = z_thresh;
                        XCorrData(day).units(uu).CSM_zmod.min_consec_bins = min_consec_bins;
                        XCorrData(day).units(uu).CSM_zmod.min_base_std_Hz = min_base_std;

                        XCorrData(day).units(uu).CSM_zmod.bin_ms          = bin_ms;
                        XCorrData(day).units(uu).CSM_zmod.step_ms         = step_ms;
                        XCorrData(day).units(uu).CSM_zmod.bin_starts_ms   = bin_starts(:);

                        XCorrData(day).units(uu).CSM_zmod.mu_base_Hz      = mu_base;
                        XCorrData(day).units(uu).CSM_zmod.sd_base_Hz      = sd_base;
                        XCorrData(day).units(uu).CSM_zmod.frac_bins_over  = frac_bins_over;

                        XCorrData(day).units(uu).CSM_zmod.deltaZ          = deltaZ(:);
                        XCorrData(day).units(uu).CSM_zmod.deltaZ_mean     = dz_mean;
                        XCorrData(day).units(uu).CSM_zmod.deltaZ_max      = dz_peak_pos;
                        XCorrData(day).units(uu).CSM_zmod.deltaZ_min      = dz_peak_neg;

                        XCorrData(day).units(uu).CSM_zmod.max_inc_run     = max_inc_run;
                        XCorrData(day).units(uu).CSM_zmod.max_dec_run     = max_dec_run;
                        XCorrData(day).units(uu).CSM_zmod.run_stat_obs    = run_stat_obs;

                        XCorrData(day).units(uu).CSM_zmod.do_unit_shuffle = do_unit_shuffle;
                        XCorrData(day).units(uu).CSM_zmod.nPerm_unit      = nPerm_unit;
                        XCorrData(day).units(uu).CSM_zmod.alpha_unit      = alpha_unit;
                        XCorrData(day).units(uu).CSM_zmod.p_unit          = p_unit;

                        XCorrData(day).units(uu).CSM_zmod.sig_inc         = sig_inc;
                        XCorrData(day).units(uu).CSM_zmod.sig_dec         = sig_dec;
                        XCorrData(day).units(uu).CSM_zmod.sig_bi          = sig_bi;
                        XCorrData(day).units(uu).CSM_zmod.sig_any         = sig_any;

                        XCorrData(day).units(uu).CSM_zmod.dom_dir         = dom_dir;
                        XCorrData(day).units(uu).CSM_zmod.dom_val         = dom_val;

                        if sig_any
                            mod_units_session(end+1,1) = unitIdx;

                            % If PSTH z exists (from earlier), add to MOD averages
                            if isfield(XCorrData(day).units(uu),'PSTH') && isfield(XCorrData(day).units(uu).PSTH,'skip') && ~XCorrData(day).units(uu).PSTH.skip
                                Agg.ModKept.baseline = vertcat(Agg.ModKept.baseline, XCorrData(day).units(uu).PSTH.base_z);
                                Agg.ModKept.stim     = vertcat(Agg.ModKept.stim,     XCorrData(day).units(uu).PSTH.stim_z);
                                Agg.ModKept.wash     = vertcat(Agg.ModKept.wash,     XCorrData(day).units(uu).PSTH.wash_z);

                                if isfield(XCorrData(day).units(uu),'category')
                                    gn = XCorrData(day).units(uu).category;
                                    if ismember(gn, group_names)
                                        Agg.ByGroup.(gn).ModKept.baseline = vertcat(Agg.ByGroup.(gn).ModKept.baseline, XCorrData(day).units(uu).PSTH.base_z);
                                        Agg.ByGroup.(gn).ModKept.stim     = vertcat(Agg.ByGroup.(gn).ModKept.stim,     XCorrData(day).units(uu).PSTH.stim_z);
                                        Agg.ByGroup.(gn).ModKept.wash     = vertcat(Agg.ByGroup.(gn).ModKept.wash,     XCorrData(day).units(uu).PSTH.wash_z);
                                    end
                                end
                            end
                        end

                        skip.ok = skip.ok + 1;
                    end
                end
            end
        end

    end % unit loop

    % session-level lists
    XCorrData(day).kept_units           = kept_units_session;
    XCorrData(day).mod_units            = mod_units_session;
    XCorrData(day).negLag_posCorr_units = negLag_posCorr_units;
    XCorrData(day).negLag_negCorr_units = negLag_negCorr_units;

    fprintf(' Session %d: kept=%d | mod(sig_any)=%d | negLag_posCorr=%d | negLag_negCorr=%d\n', ...
        day, numel(kept_units_session), numel(mod_units_session), numel(negLag_posCorr_units), numel(negLag_negCorr_units));

    if print_session_qc
        fprintf('   [SKIP] peakGate=%d | missingRasters=%d | emptyWin=%d | minTrials=%d | minSpikes=%d | OK=%d\n', ...
            skip.peak_gate_fail, skip.missing_rasters, skip.empty_windows, skip.min_trials, skip.min_spikes, skip.ok);
    end
end

fprintf('\n================ ALL SESSIONS COMPLETE ================\n');

%% ===================== WATERFALL / ROW PLOT: PEAK TIMING RELATIVE TO MAX OUTWARD V =====================
% What it does:
%   - Uses your stored keep-gated peak time from peakgate_field (already relative to max_vx)
%   - Makes one row per kept unit
%   - Draws a short vertical tick at that unit's peak time
%   - Colors by neg-lag xcorr category:
%         negLag_posCorr = rgb('SteelBlue')
%         negLag_negCorr = rgb('DarkGreen')
%
% Put this AFTER your session loop finishes (after "ALL SESSIONS COMPLETE")

% ---- collect peaks + categories across all sessions ----
peak_ms = [];
catStr  = {};
unitID  = [];
dayID   = [];

for day = 1:numel(XCorrData)
    if ~isfield(XCorrData(day),'units') || isempty(XCorrData(day).units), continue; end
    for uu = 1:numel(XCorrData(day).units)
        if ~isfield(XCorrData(day).units(uu),'keep') || ~XCorrData(day).units(uu).keep
            continue
        end

        % peak time relative to max outward velocity (because peakgate_field is *_max_vx)
        if ~isfield(XCorrData(day).units(uu),'psth_peak_time') || ~isfinite(XCorrData(day).units(uu).psth_peak_time)
            continue
        end

        peak_ms(end+1,1) = XCorrData(day).units(uu).psth_peak_time;

        if isfield(XCorrData(day).units(uu),'category')
            catStr{end+1,1} = XCorrData(day).units(uu).category;
        else
            catStr{end+1,1} = 'negLag_unknown';
        end

        if isfield(XCorrData(day).units(uu),'unitIdx')
            unitID(end+1,1) = XCorrData(day).units(uu).unitIdx;
        else
            unitID(end+1,1) = NaN;
        end
        dayID(end+1,1) = day;
    end
end

if isempty(peak_ms)
    warning('No kept units with peak times found (peak_ms empty).');
else
    % ---- split groups ----
    isPos = strcmp(catStr,'negLag_posCorr');
    isNeg = strcmp(catStr,'negLag_negCorr');

    % ---- sorting (optional): sort within each group by peak time ----
    % If you prefer: sort all together by peak_ms, just replace the block below.
    idxPos = find(isPos); [~,sp] = sort(peak_ms(idxPos)); idxPos = idxPos(sp);
    idxNeg = find(isNeg); [~,sn] = sort(peak_ms(idxNeg)); idxNeg = idxNeg(sn);
    idxOther = find(~isPos & ~isNeg); [~,so] = sort(peak_ms(idxOther)); idxOther = idxOther(so);

    plotOrder = [idxPos; idxNeg; idxOther];

    peak_plot = peak_ms(plotOrder);
    cat_plot  = catStr(plotOrder);

    nU = numel(plotOrder);
    y  = (1:nU)';

    % ---- colors ----
    cPos = rgb('SteelBlue');   % user request
    cNeg = rgb('DarkGreen');   % user request
    cOth = [0.6 0.6 0.6];      % fallback for unknowns

    % ---- make plot ----
    %figure('Color','w','Position',[200 100 420 650]); hold on
    figure('Color','w','Position',[200 100 650 400]); hold on

    % draw faint horizontal guides (optional)
    for i = 1:nU
        plot([-250 250],[i i],'-','Color',[0.92 0.92 0.92],'LineWidth',0.5);
    end

    % draw peak ticks
    tickHalfHeight = 0.35;  % thickness of each row tick
    for i = 1:nU
        if strcmp(cat_plot{i},'negLag_posCorr')
            col = cPos;
        elseif strcmp(cat_plot{i},'negLag_negCorr')
            col = cNeg;
        else
            col = cOth;
        end

        x = peak_plot(i);
        plot([x x],[y(i)-tickHalfHeight y(i)+tickHalfHeight],'-','Color',col,'LineWidth',2.0);
    end

    % reference line at 0 ms (aligned max_vx)
    xline(0,'k-','LineWidth',1);

    % axes
    xlim(valid_peak_range_ms);   % uses your existing [-250 250] range
    ylim([0 nU+1]);
    set(gca,'YDir','reverse');   % like a waterfall raster (top = first)
    xlabel('Peak time (ms) relative to max outward velocity');
    ylabel('Unit (sorted)');

    % minimalist y ticks (optional)
    yticks([]);

    % legend proxy lines
    p1 = plot(nan,nan,'-','Color',cPos,'LineWidth',2);
    p2 = plot(nan,nan,'-','Color',cNeg,'LineWidth',2);
    p3 = plot(nan,nan,'-','Color',cOth,'LineWidth',2);
    legend([p1 p2 p3],{'negLag\_posCorr','negLag\_negCorr','other/unknown'},'Location','southoutside','Box','off');

    title(sprintf('Peak timing by unit (N=%d kept units)', nU));

    % ---- optional: print counts ----
    fprintf('[PEAK ROW PLOT] kept=%d | posCorr=%d | negCorr=%d | other=%d\n', ...
        nU, sum(strcmp(cat_plot,'negLag_posCorr')), sum(strcmp(cat_plot,'negLag_negCorr')), sum(~strcmp(cat_plot,'negLag_posCorr') & ~strcmp(cat_plot,'negLag_negCorr')));

    % ---- optional: save editable vector ----
    % exportgraphics(gcf,'PeakTiming_byUnit.pdf','ContentType','vector');
end
%% ===================== WATERFALL / ROW PLOT: PEAK TIMING SORTED BY PEAK (COLOR BY CORR SIGN) =====================
% Copy/paste this whole block AFTER your session loop finishes.
% Requirements:
%   - Uses XCorrData(day).units(uu).keep == true
%   - Uses XCorrData(day).units(uu).psth_peak_time (already relative to max_vx since peakgate_field is *_max_vx)
%   - Colors by XCorrData(day).units(uu).category:
%         'negLag_posCorr' -> rgb('SteelBlue')
%         'negLag_negCorr' -> rgb('DarkGreen')
%   - Sorts ALL kept units by peak time so pos/neg units are intermingled.

% ---- collect kept-unit peak times + categories across all sessions ----
peak_ms = [];
catStr  = {};
unitID  = [];
dayID   = [];

for day = 1:numel(XCorrData)

    if ~isfield(XCorrData(day),'units') || isempty(XCorrData(day).units)
        continue
    end

    for uu = 1:numel(XCorrData(day).units)

        % keep gate
        if ~isfield(XCorrData(day).units(uu),'keep') || ~XCorrData(day).units(uu).keep
            continue
        end

        % peak time (ms) relative to max outward velocity
        if ~isfield(XCorrData(day).units(uu),'psth_peak_time') || ~isfinite(XCorrData(day).units(uu).psth_peak_time)
            continue
        end

        peak_ms(end+1,1) = XCorrData(day).units(uu).psth_peak_time;

        % category for color
        if isfield(XCorrData(day).units(uu),'category') && ~isempty(XCorrData(day).units(uu).category)
            catStr{end+1,1} = XCorrData(day).units(uu).category;
        else
            catStr{end+1,1} = 'negLag_unknown';
        end

        % IDs (optional bookkeeping)
        if isfield(XCorrData(day).units(uu),'unitIdx')
            unitID(end+1,1) = XCorrData(day).units(uu).unitIdx;
        else
            unitID(end+1,1) = NaN;
        end
        dayID(end+1,1) = day;

    end
end

if isempty(peak_ms)
    warning('No kept units with peak times found. Nothing to plot.');
else
    % ---- sort ALL units by PEAK TIME (intermingle pos/neg) ----
    [~, order] = sort(peak_ms, 'ascend');
    peak_plot  = peak_ms(order);
    cat_plot   = catStr(order);
    unitIDplot = unitID(order);
    dayIDplot  = dayID(order);

    nU = numel(order);
    y  = (1:nU)';

    % ---- colors ----
    cPos = rgb('SteelBlue');    % positively correlated
    cNeg = rgb('DarkGreen');    % negatively correlated
    cOth = [0.60 0.60 0.60];    % unknown/other

    % ---- make figure ----
    figure('Color','w','Position',[200 100 360 650]); hold on

    % faint row guides
    for i = 1:nU
        plot([valid_peak_range_ms(1) valid_peak_range_ms(2)], [i i], '-', ...
            'Color',[0.92 0.92 0.92], 'LineWidth', 2);
    end

    % peak tick per unit
    tickHalfHeight = 0.35;
    for i = 1:nU

        if strcmp(cat_plot{i}, 'negLag_posCorr')
            col = cPos;
        elseif strcmp(cat_plot{i}, 'negLag_negCorr')
            col = cNeg;
        else
            col = cOth;
        end

        x = peak_plot(i);
        plot([x x], [y(i)-tickHalfHeight y(i)+tickHalfHeight], '-', ...
            'Color', col, 'LineWidth', 2.0);
    end

    % reference at 0 ms (max_vx alignment)
    xline(0,'k-','LineWidth',1);

    % axes formatting
    xlim(valid_peak_range_ms);   % uses your existing [-250 250]
    ylim([0 nU+1]);
    set(gca,'YDir','reverse');   % top = earliest peak (since order is ascending)
    xlabel('Peak time (ms) relative to max outward velocity');
    ylabel('Unit (sorted by peak time)');
    yticks([]);

    % legend
    p1 = plot(nan,nan,'-','Color',cPos,'LineWidth',2);
    p2 = plot(nan,nan,'-','Color',cNeg,'LineWidth',2);
    p3 = plot(nan,nan,'-','Color',cOth,'LineWidth',2);
    legend([p1 p2 p3], {'negLag\_posCorr','negLag\_negCorr','other/unknown'}, ...
        'Location','southoutside','Box','off');

    title(sprintf('Peak timing by unit (sorted by peak time; N=%d)', nU));

    % ---- optional: print counts ----
    nPos = sum(strcmp(cat_plot,'negLag_posCorr'));
    nNeg = sum(strcmp(cat_plot,'negLag_negCorr'));
    nOth = nU - nPos - nNeg;
    fprintf('[PEAK SORTED PLOT] kept=%d | posCorr=%d | negCorr=%d | other=%d\n', nU, nPos, nNeg, nOth);

    % ---- optional: export ----
    % exportgraphics(gcf,'PeakTiming_sortedByPeakTime.pdf','ContentType','vector');

    % ---- optional: table in workspace for debugging ----
    PeakPlotTable = table(dayIDplot, unitIDplot, peak_plot, string(cat_plot), ...
        'VariableNames', {'day','unitIdx','peakTime_ms','category'});
end
%% ===================== PEAK TIME SUMMARY STATS =====================
% Uses PeakPlotTable created above

if exist('PeakPlotTable','var') && ~isempty(PeakPlotTable)

    % all kept units
    all_peak = PeakPlotTable.peakTime_ms;
    all_peak = all_peak(isfinite(all_peak));

    mean_all = mean(all_peak,'omitnan');
    sd_all   = std(all_peak,'omitnan');
    sem_all  = sd_all / sqrt(sum(~isnan(all_peak)));
    n_all    = sum(~isnan(all_peak));

    % pos corr
    idxPos = PeakPlotTable.category == "negLag_posCorr";
    pos_peak = PeakPlotTable.peakTime_ms(idxPos);
    pos_peak = pos_peak(isfinite(pos_peak));

    mean_pos = mean(pos_peak,'omitnan');
    sd_pos   = std(pos_peak,'omitnan');
    sem_pos  = sd_pos / sqrt(sum(~isnan(pos_peak)));
    n_pos    = sum(~isnan(pos_peak));

    % neg corr
    idxNeg = PeakPlotTable.category == "negLag_negCorr";
    neg_peak = PeakPlotTable.peakTime_ms(idxNeg);
    neg_peak = neg_peak(isfinite(neg_peak));

    mean_neg = mean(neg_peak,'omitnan');
    sd_neg   = std(neg_peak,'omitnan');
    sem_neg  = sd_neg / sqrt(sum(~isnan(neg_peak)));
    n_neg    = sum(~isnan(neg_peak));

    % print to command window
    fprintf('\n================ PEAK TIME SUMMARY ================\n');
    fprintf('All kept units     : mean = %.2f ms | SD = %.2f | SEM = %.2f | N = %d\n', ...
        mean_all, sd_all, sem_all, n_all);
    fprintf('negLag_posCorr     : mean = %.2f ms | SD = %.2f | SEM = %.2f | N = %d\n', ...
        mean_pos, sd_pos, sem_pos, n_pos);
    fprintf('negLag_negCorr     : mean = %.2f ms | SD = %.2f | SEM = %.2f | N = %d\n', ...
        mean_neg, sd_neg, sem_neg, n_neg);

    % optional statistical comparison between groups
    if ~isempty(pos_peak) && ~isempty(neg_peak)
        [~,p_peak,~,stats_peak] = ttest2(pos_peak, neg_peak);
        fprintf('posCorr vs negCorr : t(%d) = %.3f, p = %.5f\n', ...
            stats_peak.df, stats_peak.tstat, p_peak);
    end
    fprintf('===================================================\n\n');

    % optional summary table in workspace
    PeakSummary = table( ...
        ["AllKept"; "negLag_posCorr"; "negLag_negCorr"], ...
        [n_all; n_pos; n_neg], ...
        [mean_all; mean_pos; mean_neg], ...
        [sd_all; sd_pos; sd_neg], ...
        [sem_all; sem_pos; sem_neg], ...
        'VariableNames', {'Group','N','Mean_peak_ms','SD_peak_ms','SEM_peak_ms'});

end
%% ---------------------- FINAL AVERAGES (OPTIONAL) ----------------------
% At this point you have:
%   Agg.AllKept.baseline / Agg.AllKept.stim
%   Agg.ModKept.baseline / Agg.ModKept.stim
% plus Agg.ByGroup.(negLag_posCorr / negLag_negCorr).(AllKept/ModKept).(baseline/stim)
%
% Example quick means (z-score PSTHs):
if ~isempty(Agg.AllKept.baseline)
    mean_all_base = mean(Agg.AllKept.baseline,1,'omitnan');
    mean_all_stim = mean(Agg.AllKept.stim,    1,'omitnan');
    fprintf('ALL KEPT: n=%d\n', size(Agg.AllKept.baseline,1));
end
if ~isempty(Agg.ModKept.baseline)
    mean_mod_base = mean(Agg.ModKept.baseline,1,'omitnan');
    mean_mod_stim = mean(Agg.ModKept.stim,    1,'omitnan');
    fprintf('MOD KEPT (sig_any): n=%d\n', size(Agg.ModKept.baseline,1));
end


%% ============================================================
% PLOT DELTAS PER NEG-LAG CATEGORY
%   - Categories: negLag_posCorr, negLag_negCorr
%   - Two subsets overlaid on SAME axes:
%       (A) ALL kept units (peak-gated)
%       (B) CSM-modulated kept units (sig_any)
%   - Two delta traces per subset:
%       Stim - Baseline   : DodgerBlue
%       Washout - Baseline: Magenta
%
% ASSUMES you have these matrices already collected (rows=units, cols=time):
%   Agg.ByGroup.(gn).AllKept.baseline
%   Agg.ByGroup.(gn).AllKept.stim
%   Agg.ByGroup.(gn).AllKept.wo
%   Agg.ByGroup.(gn).ModKept.baseline
%   Agg.ByGroup.(gn).ModKept.stim
%   Agg.ByGroup.(gn).ModKept.wo
%
% And a time vector (ms): timevec
% ============================================================

% ---- colors ----
DodgerBlue = rgb('DodgerBlue');
Magenta    = rgb('Magenta');

% ---- categories ----
cats = {'negLag_posCorr','negLag_negCorr'};

% ---- optional plot window ----
plot_win = [-500 500];
tmask = (timevec >= plot_win(1)) & (timevec <= plot_win(2));
tplot = timevec(tmask);

figure('Position',[200 200 1100 450]); 
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

for ci = 1:numel(cats)
    gn = cats{ci};
    nexttile; hold on;

    % -------- ALL KEPT --------
    B_all  = Agg.ByGroup.(gn).AllKept.baseline;
    S_all  = Agg.ByGroup.(gn).AllKept.stim;
    W_all  = Agg.ByGroup.(gn).AllKept.wash;

    % -------- MOD KEPT (sig_any) --------
    B_mod  = Agg.ByGroup.(gn).ModKept.baseline;
    S_mod  = Agg.ByGroup.(gn).ModKept.stim;
    W_mod  = Agg.ByGroup.(gn).ModKept.wash;

    % Defensive: handle empty groups gracefully
    nAll = size(B_all,1);
    nMod = size(B_mod,1);

    % Compute deltas (mean across units)
    if nAll > 0
        dStim_all = mean(S_all - B_all, 1, 'omitnan');
        dWO_all   = mean(W_all - B_all, 1, 'omitnan');
        plot(tplot, dStim_all(tmask), 'LineWidth', 2, 'Color', DodgerBlue);
        plot(tplot, dWO_all(tmask),   'LineWidth', 2, 'Color', Magenta);
    end

    if nMod > 0
        dStim_mod = mean(S_mod - B_mod, 1, 'omitnan');
        dWO_mod   = mean(W_mod - B_mod, 1, 'omitnan');
        plot(tplot, dStim_mod(tmask), '--', 'LineWidth', 2, 'Color', DodgerBlue);
        plot(tplot, dWO_mod(tmask),   '--', 'LineWidth', 2, 'Color', Magenta);
    end

    % Axes cosmetics
    yline(0,'k-');
    xline(0,'k-');
    xlim(plot_win);
    xlabel('Time (ms)');
    ylabel('\Delta (Condition - Baseline) [z]');
    title(sprintf('%s | ALL kept n=%d | CSM-mod n=%d', gn, nAll, nMod), 'Interpreter','none');

    % Legend (only once per subplot; adapt to what exists)
    leg = {};
    if nAll > 0
        leg = [leg, {'Stim - BL (all kept)','WO - BL (all kept)'}];
    end
    if nMod > 0
        leg = [leg, {'Stim - BL (CSM-mod)','WO - BL (CSM-mod)'}];
    end
    if ~isempty(leg)
        legend(leg, 'Location','best','Box','off');
    end

    set(gca,'TickDir','out');
end


%% ============================================================
% DELTA PLOTS (Stim-BL and WO-BL) BY NEG-LAG CATEGORY
%   - Separate figures for:
%       (1) ALL KEPT units
%       (2) CSM-MODULATED KEPT units (sig_any)
%   - Two traces per figure:
%       Stim - Baseline   : DodgerBlue
%       WO   - Baseline   : Magenta
%   - Plot range: -250 to 250 ms
%
% ASSUMES you already have:
%   timevec_global (or timevec)
%   Agg.ByGroup.(cat).AllKept.baseline / stim / wo     (rows=units, cols=time)
%   Agg.ByGroup.(cat).ModKept.baseline / stim / wo
% ============================================================

% ---- colors ----
Color_StimDelta = rgb('DodgerBlue');
Color_WODelta   = rgb('Magenta');

% ---- time window ----
full_timevec = timevec;   % change to timevec if that's what you use
%full_timevec = linspace(-1000,1000,size(Agg.AllKept.baseline,2));
%time_idx     = find(full_timevec >= -250 & full_timevec <= 250);
time_idx     = find(full_timevec >= -150 & full_timevec <= 250);
%time_idx     = find(full_timevec >= -500 & full_timevec <= 500);
time_vector  = full_timevec(time_idx);

% ---- helper: delta mean+sem ----
delta_sem = @(A,B) deal( ...
    nanmean(A(:,time_idx) - B(:,time_idx), 1), ...
    nanstd( A(:,time_idx) - B(:,time_idx), 0, 1) ./ sqrt(max(1,size(A,1))) );

group_names = {'negLag_posCorr','negLag_negCorr'};

%% ===================== (1) ALL KEPT =====================
for ci = 1:numel(group_names)
    cat = group_names{ci};

     BL = Agg.ByGroup.(cat).AllKept.baseline;
     ST = Agg.ByGroup.(cat).AllKept.stim;
     WO = Agg.ByGroup.(cat).AllKept.wash;

    %BL = Agg.ByGroup.(cat).ModKept.baseline;
   % ST = Agg.ByGroup.(cat).ModKept.stim;
    %WO = Agg.ByGroup.(cat).ModKept.wash;

    if isempty(BL) || isempty(ST) || isempty(WO)
        fprintf('AllKept %s empty -> skipping.\n', cat);
        continue;
    end

    % ---- deltas ----
    [mStim, semStim] = delta_sem(ST, BL);
    [mWO,   semWO  ] = delta_sem(WO, BL);
   

    figure('Name',['ALL-KEPT DELTAS ' cat],'Position',[100,100,650,300]); hold on;

    plot(time_vector, mStim, 'Color', Color_StimDelta, 'LineWidth', 2);
    fill([time_vector fliplr(time_vector)], ...
         [mStim+semStim fliplr(mStim-semStim)], ...
         Color_StimDelta, 'FaceAlpha',0.20, 'EdgeColor','none');

    plot(time_vector, mWO, 'Color', Color_WODelta, 'LineWidth', 2);
    fill([time_vector fliplr(time_vector)], ...
         [mWO+semWO fliplr(mWO-semWO)], ...
         Color_WODelta, 'FaceAlpha',0.20, 'EdgeColor','none');

    xline(0,'--k'); yline(0,'--k');
    yline(-2,':'); yline(2,':');
    xline(-10)
    % ylim([-5 5])
    % yticks(-66:2:6)

    %xlabel('Time from max velocity (ms)');   % edit label if you want (csmT2 etc.)
    xlabel('Time from threshold crossing (ms)');   % edit label if you want (csmT2 etc.)
    ylabel('\Delta Z (Condition - Baseline)');
    xlim([-250 250]);
    title(sprintf('ALL-KEPT Δ: %s (N=%d)', cat, size(BL,1)), 'Interpreter','none');
    legend({'Stim-BL','', 'WO-BL',''}, 'Location','best','Box','off'); % blank entries hide fills
    box off;
end


%% Correct time axis 

%% ============================================================
% DELTA PLOTS (Stim-BL and WO-BL) BY NEG-LAG CATEGORY
%   - Corrected for conv(...,'valid') shrinkage
% ============================================================

% ---- colors ----
Color_StimDelta = rgb('DodgerBlue');
Color_WODelta   = rgb('Magenta');

% ---- kernel settings used in Psth1 ----
smoother   = 5;
kernel_len = smoother * 6;   % 30

% ---- actual PSTH length after valid convolution ----
psth_len = size(Agg.AllKept.baseline, 2);

% ---- infer ORIGINAL raster length ----
orig_len = psth_len + kernel_len - 1;

% ---- rebuild ORIGINAL time axis using the inferred length ----
orig_timevec = linspace(-1000, 1000, orig_len);

% ---- compute correct valid-convolution time axis ----
full_timevec = zeros(1, psth_len);
for i = 1:psth_len
    full_timevec(i) = mean(orig_timevec(i:i+kernel_len-1));
end

% ---- sanity check ----
if length(full_timevec) ~= psth_len
    error('Time vector length (%d) does not match PSTH length (%d).', ...
        length(full_timevec), psth_len);
end

% ---- time window ----
time_idx    = find(full_timevec >= -250 & full_timevec <= 250);
time_vector = full_timevec(time_idx);

group_names = {'negLag_posCorr','negLag_negCorr'};

% ===================== ALL KEPT =====================
for ci = 1:numel(group_names)
    cat = group_names{ci};

    BL = Agg.ByGroup.(cat).AllKept.baseline(:, time_idx);
    ST = Agg.ByGroup.(cat).AllKept.stim(:,     time_idx);
    WO = Agg.ByGroup.(cat).AllKept.wash(:,     time_idx);

    if isempty(BL) || isempty(ST) || isempty(WO)
        fprintf('AllKept %s empty -> skipping.\n', cat);
        continue;
    end

    % ---- paired deltas ----
    dStim = ST - BL;
    dWO   = WO - BL;

    % ---- mean ----
    mStim = mean(dStim, 1, 'omitnan');
    mWO   = mean(dWO,   1, 'omitnan');

    % ---- SEM ----
    nStim = sum(~isnan(dStim), 1);
    nWO   = sum(~isnan(dWO),   1);

    semStim = std(dStim, 0, 1, 'omitnan') ./ sqrt(nStim);
    semWO   = std(dWO,   0, 1, 'omitnan') ./ sqrt(nWO);

    mStim(nStim==0)   = NaN;
    semStim(nStim==0) = NaN;
    mWO(nWO==0)       = NaN;
    semWO(nWO==0)     = NaN;
    

    % ---- y-limits ----
    yvals = [mStim+semStim, mStim-semStim, mWO+semWO, mWO-semWO];
    ymax = max(abs(yvals), [], 'omitnan');
    if isempty(ymax) || ~isfinite(ymax) || ymax == 0
        ymax = 1;
    end
    ymax = ymax * 1.1;

    % ---- plot ----
    figure('Name',['ALL-KEPT DELTAS ' cat],'Position',[100,100,650,300]);
    hold on;

    fill([time_vector fliplr(time_vector)], ...
         [mStim+semStim fliplr(mStim-semStim)], ...
         Color_StimDelta, 'FaceAlpha',0.20, 'EdgeColor','none');

    fill([time_vector fliplr(time_vector)], ...
         [mWO+semWO fliplr(mWO-semWO)], ...
         Color_WODelta, 'FaceAlpha',0.20, 'EdgeColor','none');

    plot(time_vector, mStim, 'Color', Color_StimDelta, 'LineWidth', 2);
    plot(time_vector, mWO,   'Color', Color_WODelta,   'LineWidth', 2);

    xline(0,'--k');
    yline(0,'--k');
    xline(-10,'-k');

    yline(-2,':'); yline(2,':');

    xlim([-250 250]);
    ylim([-ymax ymax]);

    xlabel('Time from threshold crossing (ms)');
    ylabel('\Delta firing rate (Hz)');
    title(sprintf('ALL-KEPT ΔFR: %s (N=%d)', cat, size(BL,1)), 'Interpreter','none');

    legend({'Stim SEM','','WO SEM','','Stim-BL','WO-BL'}, ...
           'Location','best','Box','off');

    box off;
    set(gca,'TickDir','out');
end
%%
% ===================== ALL KEPT =====================
PeakMetrics = struct;

% ---- multiple peak settings ----
nPeaksShow  = 3;    % number of top peaks to display
minPeakDist = 5;   % minimum distance between peaks (in samples)

for ci = 1:numel(group_names)
    cat = group_names{ci};

    BL = Agg.ByGroup.(cat).AllKept.baseline(:, time_idx);
    ST = Agg.ByGroup.(cat).AllKept.stim(:,     time_idx);
    WO = Agg.ByGroup.(cat).AllKept.wash(:,     time_idx);

    if isempty(BL) || isempty(ST) || isempty(WO)
        fprintf('AllKept %s empty -> skipping.\n', cat);
        continue;
    end

    % ---- paired deltas ----
    dStim = ST - BL;
    dWO   = WO - BL;

    % ---- mean ----
    mStim = mean(dStim, 1, 'omitnan');
    mWO   = mean(dWO,   1, 'omitnan');

    % ---- SEM ----
    nStim = sum(~isnan(dStim), 1);
    nWO   = sum(~isnan(dWO),   1);

    semStim = std(dStim, 0, 1, 'omitnan') ./ sqrt(nStim);
    semWO   = std(dWO,   0, 1, 'omitnan') ./ sqrt(nWO);

    mStim(nStim==0)   = NaN;
    semStim(nStim==0) = NaN;
    mWO(nWO==0)       = NaN;
    semWO(nWO==0)     = NaN;

    % ---- find multiple peaks on |mean trace| ----
    stim_for_peaks = abs(mStim);
    wo_for_peaks   = abs(mWO);

    stim_for_peaks(~isfinite(stim_for_peaks)) = NaN;
    wo_for_peaks(~isfinite(wo_for_peaks))     = NaN;

    % replace NaNs temporarily so findpeaks does not choke
    stim_tmp = stim_for_peaks;
    wo_tmp   = wo_for_peaks;
    stim_tmp(~isfinite(stim_tmp)) = -inf;
    wo_tmp(~isfinite(wo_tmp))     = -inf;

    [pksStim, locsStim] = findpeaks(stim_tmp, 'MinPeakDistance', minPeakDist);
    [pksWO,   locsWO]   = findpeaks(wo_tmp,   'MinPeakDistance', minPeakDist);

    % fallback in case findpeaks finds nothing
    if isempty(locsStim)
        [~, locsStim] = max(stim_tmp);
        pksStim = stim_tmp(locsStim);
    end
    if isempty(locsWO)
        [~, locsWO] = max(wo_tmp);
        pksWO = wo_tmp(locsWO);
    end

    % sort peaks by amplitude descending
    [~, ordStim] = sort(pksStim, 'descend');
    [~, ordWO]   = sort(pksWO,   'descend');

    locsStim = locsStim(ordStim);
    pksStim  = pksStim(ordStim);

    locsWO = locsWO(ordWO);
    pksWO  = pksWO(ordWO);

    % keep only top N
    nStimShow = min(nPeaksShow, numel(locsStim));
    nWOShow   = min(nPeaksShow, numel(locsWO));

    locsStim = locsStim(1:nStimShow);
    pksStim  = pksStim(1:nStimShow);

    locsWO = locsWO(1:nWOShow);
    pksWO  = pksWO(1:nWOShow);

    % store peak metrics
    PeakMetrics.AllKept.(cat).Stim.locs     = locsStim;
    PeakMetrics.AllKept.(cat).Stim.time_ms  = time_vector(locsStim);
    PeakMetrics.AllKept.(cat).Stim.mean     = mStim(locsStim);
    PeakMetrics.AllKept.(cat).Stim.sem      = semStim(locsStim);
    PeakMetrics.AllKept.(cat).Stim.absmean  = abs(mStim(locsStim));

    PeakMetrics.AllKept.(cat).WO.locs       = locsWO;
    PeakMetrics.AllKept.(cat).WO.time_ms    = time_vector(locsWO);
    PeakMetrics.AllKept.(cat).WO.mean       = mWO(locsWO);
    PeakMetrics.AllKept.(cat).WO.sem        = semWO(locsWO);
    PeakMetrics.AllKept.(cat).WO.absmean    = abs(mWO(locsWO));

    % print to command window
    fprintf('\n[%s | AllKept]\n', cat);

    fprintf('Stim-BL top peaks:\n');
    for k = 1:nStimShow
        fprintf('  Peak %d: time = %.1f ms | mean = %.3f Hz | SEM = %.3f | abs = %.3f\n', ...
            k, time_vector(locsStim(k)), mStim(locsStim(k)), semStim(locsStim(k)), abs(mStim(locsStim(k))));
    end

    fprintf('WO-BL top peaks:\n');
    for k = 1:nWOShow
        fprintf('  Peak %d: time = %.1f ms | mean = %.3f Hz | SEM = %.3f | abs = %.3f\n', ...
            k, time_vector(locsWO(k)), mWO(locsWO(k)), semWO(locsWO(k)), abs(mWO(locsWO(k))));
    end

    % ---- y-limits ----
    yvals = [mStim+semStim, mStim-semStim, mWO+semWO, mWO-semWO];
    ymax = max(abs(yvals), [], 'omitnan');
    if isempty(ymax) || ~isfinite(ymax) || ymax == 0
        ymax = 1;
    end
    ymax = ymax * 1.1;

    % ---- plot ----
    figure('Name',['ALL-KEPT DELTAS ' cat],'Position',[100,100,650,300]);
    hold on;

    fill([time_vector fliplr(time_vector)], ...
         [mStim+semStim fliplr(mStim-semStim)], ...
         Color_StimDelta, 'FaceAlpha',0.20, 'EdgeColor','none');

    fill([time_vector fliplr(time_vector)], ...
         [mWO+semWO fliplr(mWO-semWO)], ...
         Color_WODelta, 'FaceAlpha',0.20, 'EdgeColor','none');

    plot(time_vector, mStim, 'Color', Color_StimDelta, 'LineWidth', 2);
    plot(time_vector, mWO,   'Color', Color_WODelta,   'LineWidth', 2);

    % ---- mark Stim peaks ----
    for k = 1:nStimShow
        xpk = time_vector(locsStim(k));
        ypk = mStim(locsStim(k));

        plot(xpk, ypk, 'o', ...
            'MarkerFaceColor', Color_StimDelta, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize', 6);

        text(xpk, ypk, sprintf('  %d) %.1f ms, %.2f', k, xpk, ypk), ...
            'Color', Color_StimDelta, ...
            'FontSize', 8, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle');
    end

    % ---- mark WO peaks ----
    for k = 1:nWOShow
        xpk = time_vector(locsWO(k));
        ypk = mWO(locsWO(k));

        plot(xpk, ypk, 'o', ...
            'MarkerFaceColor', Color_WODelta, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize', 6);

        text(xpk, ypk, sprintf('  %d) %.1f ms, %.2f', k, xpk, ypk), ...
            'Color', Color_WODelta, ...
            'FontSize', 8, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle');
    end

    xline(0,'--k');
    yline(0,'--k');
    xline(-10,'-k');

    yline(-2,':');
    yline(2,':');

    xlim([-250 250]);
    ylim([-ymax ymax]);

    xlabel('Time from threshold crossing (ms)');
    ylabel('\Delta firing rate (Hz)');
    title(sprintf('ALL-KEPT ΔFR: %s (N=%d)', cat, size(BL,1)), 'Interpreter','none');

    legend({'Stim SEM','','WO SEM','','Stim-BL','WO-BL'}, ...
           'Location','best','Box','off');

    box off;
    set(gca,'TickDir','out');
end

%% ============================================================
% DELTA PLOTS (Stim-BL and WO-BL) BY NEG-LAG CATEGORY
%   - Corrected for conv(...,'valid') shrinkage
%   - Finds multiple top peaks on |delta trace|
%   - Computes z = mean / SEM at each timepoint
%   - Computes two-tailed p value from z score
%   - Displays peak time, mean, z, and p on plot
% ============================================================

clear PeakMetrics

% ---- colors ----
Color_StimDelta = rgb('DodgerBlue');
Color_WODelta   = rgb('Magenta');

% ---- kernel settings used in Psth1 ----
smoother   = 5;
kernel_len = smoother * 6;   % 30

% ---- actual PSTH length after valid convolution ----
psth_len = size(Agg.AllKept.baseline, 2);

% ---- infer ORIGINAL raster length ----
orig_len = psth_len + kernel_len - 1;

% ---- rebuild ORIGINAL time axis using inferred length ----
orig_timevec = linspace(-1000, 1000, orig_len);

% ---- compute correct valid-convolution time axis ----
full_timevec = zeros(1, psth_len);
for ii = 1:psth_len
    full_timevec(ii) = mean(orig_timevec(ii:ii+kernel_len-1));
end

% ---- sanity check ----
if length(full_timevec) ~= psth_len
    error('Time vector length (%d) does not match PSTH length (%d).', ...
        length(full_timevec), psth_len);
end

% ---- time window ----
time_idx    = find(full_timevec >= -250 & full_timevec <= 250);
time_vector = full_timevec(time_idx);

% ---- groups ----
group_names = {'negLag_posCorr','negLag_negCorr'};

% ---- peak settings ----
nPeaksShow  = 3;    % number of top peaks to display
minPeakDist = 5;   % minimum spacing between peaks in samples

% ============================================================
% ALL KEPT
% ============================================================
for ci = 1:numel(group_names)

    cat = group_names{ci};

    BL = Agg.ByGroup.(cat).AllKept.baseline(:, time_idx);
    ST = Agg.ByGroup.(cat).AllKept.stim(:,     time_idx);
    WO = Agg.ByGroup.(cat).AllKept.wash(:,     time_idx);

    % BL = Agg.ByGroup.(cat).ModKept.baseline(:, time_idx);
    % ST = Agg.ByGroup.(cat).ModKept.stim(:,     time_idx);
    % WO = Agg.ByGroup.(cat).ModKept.wash(:,     time_idx);

    if isempty(BL) || isempty(ST) || isempty(WO)
        fprintf('AllKept %s empty -> skipping.\n', cat);
        continue
    end

    % ---------------- paired deltas ----------------
    dStim = ST - BL;
    dWO   = WO - BL;

    % ---------------- mean traces ----------------
    mStim = mean(dStim, 1, 'omitnan');
    mWO   = mean(dWO,   1, 'omitnan');

    % ---------------- SEM traces ----------------
    nStim = sum(~isnan(dStim), 1);
    nWO   = sum(~isnan(dWO),   1);

    semStim = std(dStim, 0, 1, 'omitnan') ./ sqrt(nStim);
    semWO   = std(dWO,   0, 1, 'omitnan') ./ sqrt(nWO);

    mStim(nStim==0)   = NaN;
    semStim(nStim==0) = NaN;
    mWO(nWO==0)       = NaN;
    semWO(nWO==0)     = NaN;

    % ---------------- z scores + two-tailed p values ----------------
    % keeps this consistent with your earlier z-based approach
    zStim = mStim ./ semStim;
    zWO   = mWO   ./ semWO;

    pStim = 2 * (1 - normcdf(abs(zStim), 0, 1));
    pWO   = 2 * (1 - normcdf(abs(zWO),   0, 1));

    zStim(~isfinite(zStim)) = NaN;
    zWO(~isfinite(zWO))     = NaN;
    pStim(~isfinite(pStim)) = NaN;
    pWO(~isfinite(pWO))     = NaN;

    % ---------------- find multiple top peaks ----------------
    stimAbs = abs(mStim);
    woAbs   = abs(mWO);

    stimTmp = stimAbs;
    woTmp   = woAbs;

    stimTmp(~isfinite(stimTmp)) = -inf;
    woTmp(~isfinite(woTmp))     = -inf;

    [pksStim, locsStim] = findpeaks(stimTmp, 'MinPeakDistance', minPeakDist);
    [pksWO,   locsWO]   = findpeaks(woTmp,   'MinPeakDistance', minPeakDist);

    % fallback if no local peaks are found
    if isempty(locsStim)
        [~, locsStim] = max(stimTmp);
        pksStim = stimTmp(locsStim);
    end
    if isempty(locsWO)
        [~, locsWO] = max(woTmp);
        pksWO = woTmp(locsWO);
    end

    % sort by descending peak size
    [~, ordStim] = sort(pksStim, 'descend');
    [~, ordWO]   = sort(pksWO,   'descend');

    locsStim = locsStim(ordStim);
    pksStim  = pksStim(ordStim);

    locsWO = locsWO(ordWO);
    pksWO  = pksWO(ordWO);

    % keep top N
    nStimShow = min(nPeaksShow, numel(locsStim));
    nWOShow   = min(nPeaksShow, numel(locsWO));

    locsStim = locsStim(1:nStimShow);
    locsWO   = locsWO(1:nWOShow);

    % ---------------- store metrics ----------------
    PeakMetrics.AllKept.(cat).Stim.time_ms = time_vector(locsStim);
    PeakMetrics.AllKept.(cat).Stim.mean    = mStim(locsStim);
    PeakMetrics.AllKept.(cat).Stim.sem     = semStim(locsStim);
    PeakMetrics.AllKept.(cat).Stim.z       = zStim(locsStim);
    PeakMetrics.AllKept.(cat).Stim.p       = pStim(locsStim);

    PeakMetrics.AllKept.(cat).WO.time_ms   = time_vector(locsWO);
    PeakMetrics.AllKept.(cat).WO.mean      = mWO(locsWO);
    PeakMetrics.AllKept.(cat).WO.sem       = semWO(locsWO);
    PeakMetrics.AllKept.(cat).WO.z         = zWO(locsWO);
    PeakMetrics.AllKept.(cat).WO.p         = pWO(locsWO);

    % ---------------- print to command window ----------------
    fprintf('\n[%s | AllKept]\n', cat);

    fprintf('Stim-BL top peaks:\n');
    for k = 1:nStimShow
        fprintf('  Peak %d: time = %.1f ms | mean = %.3f Hz | SEM = %.3f | z = %.3f | p = %.5g\n', ...
            k, time_vector(locsStim(k)), mStim(locsStim(k)), semStim(locsStim(k)), ...
            zStim(locsStim(k)), pStim(locsStim(k)));
    end

    fprintf('WO-BL top peaks:\n');
    for k = 1:nWOShow
        fprintf('  Peak %d: time = %.1f ms | mean = %.3f Hz | SEM = %.3f | z = %.3f | p = %.5g\n', ...
            k, time_vector(locsWO(k)), mWO(locsWO(k)), semWO(locsWO(k)), ...
            zWO(locsWO(k)), pWO(locsWO(k)));
    end

    % ---------------- y limits ----------------
    yvals = [mStim+semStim, mStim-semStim, mWO+semWO, mWO-semWO];
    ymax = max(abs(yvals), [], 'omitnan');
    if isempty(ymax) || ~isfinite(ymax) || ymax == 0
        ymax = 1;
    end
    ymax = ymax * 1.12;

    % ---------------- plot ----------------
    figure('Name',['ALL-KEPT DELTAS ' cat],'Position',[100,100,700,320]); 
    hold on

    fill([time_vector fliplr(time_vector)], ...
         [mStim+semStim fliplr(mStim-semStim)], ...
         Color_StimDelta, 'FaceAlpha',0.20, 'EdgeColor','none');

    fill([time_vector fliplr(time_vector)], ...
         [mWO+semWO fliplr(mWO-semWO)], ...
         Color_WODelta, 'FaceAlpha',0.20, 'EdgeColor','none');

    plot(time_vector, mStim, 'Color', Color_StimDelta, 'LineWidth', 2);
    plot(time_vector, mWO,   'Color', Color_WODelta,   'LineWidth', 2);

    % ---------------- mark Stim peaks ----------------
    for k = 1:nStimShow
        xpk = time_vector(locsStim(k));
        ypk = mStim(locsStim(k));

        plot(xpk, ypk, 'o', ...
            'MarkerFaceColor', Color_StimDelta, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize', 6);

        text(xpk, ypk, ...
            sprintf('  %d) %.0f ms, %.2f\\nz=%.2f, p=%.3g', ...
            k, xpk, ypk, zStim(locsStim(k)), pStim(locsStim(k))), ...
            'Color', Color_StimDelta, ...
            'FontSize', 8, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle');
    end

    % ---------------- mark WO peaks ----------------
    for k = 1:nWOShow
        xpk = time_vector(locsWO(k));
        ypk = mWO(locsWO(k));

        plot(xpk, ypk, 'o', ...
            'MarkerFaceColor', Color_WODelta, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize', 6);

        text(xpk, ypk, ...
            sprintf('  %d) %.0f ms, %.2f\\nz=%.2f, p=%.3g', ...
            k, xpk, ypk, zWO(locsWO(k)), pWO(locsWO(k))), ...
            'Color', Color_WODelta, ...
            'FontSize', 8, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle');
    end

    % ---------------- reference lines ----------------
    xline(0,'--k');
    yline(0,'--k');
    xline(-10,'-k');
    yline(-2,':k');
    yline( 2,':k');

    % ---------------- axes ----------------
    xlim([-250 250]);
    ylim([-ymax ymax]);

    xlabel('Time from threshold crossing (ms)');
    ylabel('\Delta firing rate (Hz)');
    title(sprintf('ALL-KEPT \\DeltaFR: %s (N=%d)', cat, size(BL,1)), 'Interpreter','none');

    legend({'Stim SEM','','WO SEM','','Stim-BL','WO-BL'}, ...
        'Location','best','Box','off');

    box off
    set(gca,'TickDir','out')
end

%% ============================================================
% DELTA PLOTS (Stim-BL and WO-BL) BY NEG-LAG CATEGORY
%   - Corrected for conv(...,'valid') shrinkage
%   - Uses delta traces (Stim-BL and WO-BL)
%   - Finds multiple top peaks on:
%         mean(delta)
%         mean(delta) + SEM
%         mean(delta) - SEM
%   - Reports z and two-tailed p at the MEAN-delta peaks
% ============================================================

clear PeakMetrics

% ---- colors ----
Color_StimDelta = rgb('DodgerBlue');
Color_WODelta   = rgb('Magenta');

% ---- kernel settings used in Psth1 ----
smoother   = 5;
kernel_len = smoother * 6;   % 30

% ---- actual PSTH length after valid convolution ----
psth_len = size(Agg.AllKept.baseline, 2);

% ---- infer ORIGINAL raster length ----
orig_len = psth_len + kernel_len - 1;

% ---- rebuild ORIGINAL time axis using inferred length ----
orig_timevec = linspace(-1000, 1000, orig_len);

% ---- compute correct valid-convolution time axis ----
full_timevec = zeros(1, psth_len);
for ii = 1:psth_len
    full_timevec(ii) = mean(orig_timevec(ii:ii+kernel_len-1));
end

% ---- sanity check ----
if length(full_timevec) ~= psth_len
    error('Time vector length (%d) does not match PSTH length (%d).', ...
        length(full_timevec), psth_len);
end

% ---- time window ----
time_idx    = find(full_timevec >= -250 & full_timevec <= 250);
time_vector = full_timevec(time_idx);

% ---- groups ----
group_names = {'negLag_posCorr','negLag_negCorr'};

% ---- peak settings ----
nPeaksShow  = 3;
minPeakDist = 5;

for ci = 1:numel(group_names)

    cat = group_names{ci};

    % BL = Agg.ByGroup.(cat).AllKept.baseline(:, time_idx);
    % ST = Agg.ByGroup.(cat).AllKept.stim(:,     time_idx);
    % WO = Agg.ByGroup.(cat).AllKept.wash(:,     time_idx);

    BL = Agg.ByGroup.(cat).ModKept.baseline(:, time_idx);
    ST = Agg.ByGroup.(cat).ModKept.stim(:,     time_idx);
    WO = Agg.ByGroup.(cat).ModKept.wash(:,     time_idx);

    if isempty(BL) || isempty(ST) || isempty(WO)
        fprintf('AllKept %s empty -> skipping.\n', cat);
        continue
    end

    % ============================================================
    % paired deltas
    % ============================================================
    dStim = ST - BL;
    dWO   = WO - BL;

    % ============================================================
    % mean delta + SEM
    % ============================================================
    mStim = mean(dStim, 1, 'omitnan');
    mWO   = mean(dWO,   1, 'omitnan');

    nStim = sum(~isnan(dStim), 1);
    nWO   = sum(~isnan(dWO),   1);

    semStim = std(dStim, 0, 1, 'omitnan') ./ sqrt(nStim);
    semWO   = std(dWO,   0, 1, 'omitnan') ./ sqrt(nWO);

    mStim(nStim==0)   = NaN;
    semStim(nStim==0) = NaN;
    mWO(nWO==0)       = NaN;
    semWO(nWO==0)     = NaN;

    % ---- delta +/- SEM ----
    upStim = mStim + semStim;
    loStim = mStim - semStim;

    upWO = mWO + semWO;
    loWO = mWO - semWO;

    % ============================================================
    % z and two-tailed p from mean/SEM (same style as before)
    % ============================================================
    zStim = mStim ./ semStim;
    zWO   = mWO   ./ semWO;

    pStim = 2 * (1 - normcdf(abs(zStim), 0, 1));
    pWO   = 2 * (1 - normcdf(abs(zWO),   0, 1));

    zStim(~isfinite(zStim)) = NaN;
    zWO(~isfinite(zWO))     = NaN;
    pStim(~isfinite(pStim)) = NaN;
    pWO(~isfinite(pWO))     = NaN;

    % ============================================================
    % peak finder helper
    % ============================================================
    getTopPeaksAbs = @(trace) local_getTopPeaksAbs(trace, nPeaksShow, minPeakDist);

    % ---- Stim peaks from mean / plus / minus ----
    [locMeanStim, ~] = getTopPeaksAbs(mStim);
    [locUpStim,   ~] = getTopPeaksAbs(upStim);
    [locLoStim,   ~] = getTopPeaksAbs(loStim);

    % ---- WO peaks from mean / plus / minus ----
    [locMeanWO, ~] = getTopPeaksAbs(mWO);
    [locUpWO,   ~] = getTopPeaksAbs(upWO);
    [locLoWO,   ~] = getTopPeaksAbs(loWO);

    % ============================================================
    % store metrics
    % ============================================================
    PeakMetrics.AllKept.(cat).Stim.mean.time_ms = time_vector(locMeanStim);
    PeakMetrics.AllKept.(cat).Stim.mean.value   = mStim(locMeanStim);
    PeakMetrics.AllKept.(cat).Stim.mean.sem     = semStim(locMeanStim);
    PeakMetrics.AllKept.(cat).Stim.mean.z       = zStim(locMeanStim);
    PeakMetrics.AllKept.(cat).Stim.mean.p       = pStim(locMeanStim);

    PeakMetrics.AllKept.(cat).Stim.plus.time_ms = time_vector(locUpStim);
    PeakMetrics.AllKept.(cat).Stim.plus.value   = upStim(locUpStim);

    PeakMetrics.AllKept.(cat).Stim.minus.time_ms = time_vector(locLoStim);
    PeakMetrics.AllKept.(cat).Stim.minus.value   = loStim(locLoStim);

    PeakMetrics.AllKept.(cat).WO.mean.time_ms = time_vector(locMeanWO);
    PeakMetrics.AllKept.(cat).WO.mean.value   = mWO(locMeanWO);
    PeakMetrics.AllKept.(cat).WO.mean.sem     = semWO(locMeanWO);
    PeakMetrics.AllKept.(cat).WO.mean.z       = zWO(locMeanWO);
    PeakMetrics.AllKept.(cat).WO.mean.p       = pWO(locMeanWO);

    PeakMetrics.AllKept.(cat).WO.plus.time_ms = time_vector(locUpWO);
    PeakMetrics.AllKept.(cat).WO.plus.value   = upWO(locUpWO);

    PeakMetrics.AllKept.(cat).WO.minus.time_ms = time_vector(locLoWO);
    PeakMetrics.AllKept.(cat).WO.minus.value   = loWO(locLoWO);

    % ============================================================
    % print
    % ============================================================
    fprintf('\n[%s | AllKept delta]\n', cat);

    fprintf('Stim-BL mean peaks:\n');
    for k = 1:numel(locMeanStim)
        fprintf('  Peak %d: time = %.1f ms | mean = %.3f | SEM = %.3f | z = %.3f | p = %.5g\n', ...
            k, time_vector(locMeanStim(k)), mStim(locMeanStim(k)), semStim(locMeanStim(k)), ...
            zStim(locMeanStim(k)), pStim(locMeanStim(k)));
    end
    fprintf('Stim-BL +SEM peaks:\n');
    for k = 1:numel(locUpStim)
        fprintf('  Peak %d: time = %.1f ms | mean+SEM = %.3f\n', ...
            k, time_vector(locUpStim(k)), upStim(locUpStim(k)));
    end
    fprintf('Stim-BL -SEM peaks:\n');
    for k = 1:numel(locLoStim)
        fprintf('  Peak %d: time = %.1f ms | mean-SEM = %.3f\n', ...
            k, time_vector(locLoStim(k)), loStim(locLoStim(k)));
    end

    fprintf('WO-BL mean peaks:\n');
    for k = 1:numel(locMeanWO)
        fprintf('  Peak %d: time = %.1f ms | mean = %.3f | SEM = %.3f | z = %.3f | p = %.5g\n', ...
            k, time_vector(locMeanWO(k)), mWO(locMeanWO(k)), semWO(locMeanWO(k)), ...
            zWO(locMeanWO(k)), pWO(locMeanWO(k)));
    end
    fprintf('WO-BL +SEM peaks:\n');
    for k = 1:numel(locUpWO)
        fprintf('  Peak %d: time = %.1f ms | mean+SEM = %.3f\n', ...
            k, time_vector(locUpWO(k)), upWO(locUpWO(k)));
    end
    fprintf('WO-BL -SEM peaks:\n');
    for k = 1:numel(locLoWO)
        fprintf('  Peak %d: time = %.1f ms | mean-SEM = %.3f\n', ...
            k, time_vector(locLoWO(k)), loWO(locLoWO(k)));
    end

    % ============================================================
    % y limits
    % ============================================================
    yvals = [upStim, loStim, upWO, loWO];
    ymax = max(abs(yvals), [], 'omitnan');
    if isempty(ymax) || ~isfinite(ymax) || ymax == 0
        ymax = 1;
    end
    ymax = ymax * 1.12;

    % ============================================================
    % plot
    % ============================================================
    figure('Name',['ALL-KEPT DELTAS ' cat],'Position',[100,100,760,340]);
    hold on

    % SEM fills
    fill([time_vector fliplr(time_vector)], ...
         [upStim fliplr(loStim)], ...
         Color_StimDelta, 'FaceAlpha',0.20, 'EdgeColor','none');

    fill([time_vector fliplr(time_vector)], ...
         [upWO fliplr(loWO)], ...
         Color_WODelta, 'FaceAlpha',0.20, 'EdgeColor','none');

    % mean traces
    plot(time_vector, mStim, 'Color', Color_StimDelta, 'LineWidth', 2);
    plot(time_vector, mWO,   'Color', Color_WODelta,   'LineWidth', 2);

    % ---------------- mark Stim mean peaks ----------------
    for k = 1:numel(locMeanStim)
        xpk = time_vector(locMeanStim(k));
        ypk = mStim(locMeanStim(k));

        plot(xpk, ypk, 'o', ...
            'MarkerFaceColor', Color_StimDelta, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize', 6);

        text(xpk, ypk, ...
            sprintf('  S%d) %.0f, %.2f\\nz=%.2f, p=%.3g', ...
            k, xpk, ypk, zStim(locMeanStim(k)), pStim(locMeanStim(k))), ...
            'Color', Color_StimDelta, ...
            'FontSize', 8, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle');
    end

    % ---------------- mark Stim +SEM peaks ----------------
    for k = 1:numel(locUpStim)
        xpk = time_vector(locUpStim(k));
        ypk = upStim(locUpStim(k));
        plot(xpk, ypk, '^', ...
            'MarkerFaceColor', Color_StimDelta, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize', 5);
    end

    % ---------------- mark Stim -SEM peaks ----------------
    for k = 1:numel(locLoStim)
        xpk = time_vector(locLoStim(k));
        ypk = loStim(locLoStim(k));
        plot(xpk, ypk, 'v', ...
            'MarkerFaceColor', Color_StimDelta, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize', 5);
    end

    % ---------------- mark WO mean peaks ----------------
    for k = 1:numel(locMeanWO)
        xpk = time_vector(locMeanWO(k));
        ypk = mWO(locMeanWO(k));

        plot(xpk, ypk, 'o', ...
            'MarkerFaceColor', Color_WODelta, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize', 6);

        text(xpk, ypk, ...
            sprintf('  W%d) %.0f, %.2f\\nz=%.2f, p=%.3g', ...
            k, xpk, ypk, zWO(locMeanWO(k)), pWO(locMeanWO(k))), ...
            'Color', Color_WODelta, ...
            'FontSize', 8, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle');
    end

    % ---------------- mark WO +SEM peaks ----------------
    for k = 1:numel(locUpWO)
        xpk = time_vector(locUpWO(k));
        ypk = upWO(locUpWO(k));
        plot(xpk, ypk, '^', ...
            'MarkerFaceColor', Color_WODelta, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize', 5);
    end

    % ---------------- mark WO -SEM peaks ----------------
    for k = 1:numel(locLoWO)
        xpk = time_vector(locLoWO(k));
        ypk = loWO(locLoWO(k));
        plot(xpk, ypk, 'v', ...
            'MarkerFaceColor', Color_WODelta, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize', 5);
    end

    % reference lines
    xline(0,'--k');
    yline(0,'--k');
    xline(-10,'-k');
    yline(-2,':k');
    yline( 2,':k');

    xlim([-250 250]);
    ylim([-ymax ymax]);

    xlabel('Time from threshold crossing (ms)');
    ylabel('\Delta firing rate (Hz)');
    title(sprintf('ALL-KEPT \\DeltaFR: %s (N=%d)', cat, size(BL,1)), 'Interpreter','none');

    legend({'Stim SEM','','WO SEM','','Stim-BL','WO-BL'}, ...
        'Location','best','Box','off');

    box off
    set(gca,'TickDir','out')
end


%% ============================================================
% DELTA PLOTS (Stim-BL and WO-BL) BY NEG-LAG CATEGORY
%   - Corrected for conv(...,'valid') shrinkage
%   - Uses delta traces (Stim-BL and WO-BL)
%   - Finds multiple top peaks on mean(delta)
%   - Reports SEM, z, and two-tailed p at those peak times
%   - Marks mean, mean+SEM, and mean-SEM at the same peak time
% ============================================================

clear PeakMetrics

% ---- colors ----
Color_StimDelta = rgb('DodgerBlue');
Color_WODelta   = rgb('Magenta');

% ---- kernel settings used in Psth1 ----
smoother   = 5;
kernel_len = smoother * 6;   % 30

% ---- actual PSTH length after valid convolution ----
psth_len = size(Agg.AllKept.baseline, 2);

% ---- infer ORIGINAL raster length ----
orig_len = psth_len + kernel_len - 1;

% ---- rebuild ORIGINAL time axis using inferred length ----
orig_timevec = linspace(-1000, 1000, orig_len);

% ---- compute correct valid-convolution time axis ----
full_timevec = zeros(1, psth_len);
for ii = 1:psth_len
    full_timevec(ii) = mean(orig_timevec(ii:ii+kernel_len-1));
end

% ---- sanity check ----
if length(full_timevec) ~= psth_len
    error('Time vector length (%d) does not match PSTH length (%d).', ...
        length(full_timevec), psth_len);
end

% ---- time window ----
time_idx    = find(full_timevec >= -250 & full_timevec <= 250);
time_vector = full_timevec(time_idx);

% ---- groups ----
group_names = {'negLag_posCorr','negLag_negCorr'};

% ---- peak settings ----
nPeaksShow  = 3;
minPeakDist = 5;

for ci = 1:numel(group_names)

    cat = group_names{ci};

    BL = Agg.ByGroup.(cat).AllKept.baseline(:, time_idx);
    ST = Agg.ByGroup.(cat).AllKept.stim(:,     time_idx);
    WO = Agg.ByGroup.(cat).AllKept.wash(:,     time_idx);

    if isempty(BL) || isempty(ST) || isempty(WO)
        fprintf('AllKept %s empty -> skipping.\n', cat);
        continue
    end

    % ============================================================
    % paired deltas
    % ============================================================
    dStim = ST - BL;
    dWO   = WO - BL;

    % ============================================================
    % mean delta + SEM
    % ============================================================
    mStim = mean(dStim, 1, 'omitnan');
    mWO   = mean(dWO,   1, 'omitnan');

    nStim = sum(~isnan(dStim), 1);
    nWO   = sum(~isnan(dWO),   1);

    semStim = std(dStim, 0, 1, 'omitnan') ./ sqrt(nStim);
    semWO   = std(dWO,   0, 1, 'omitnan') ./ sqrt(nWO);

    mStim(nStim==0)   = NaN;
    semStim(nStim==0) = NaN;
    mWO(nWO==0)       = NaN;
    semWO(nWO==0)     = NaN;

    % ---- delta +/- SEM ----
    upStim = mStim + semStim;
    loStim = mStim - semStim;

    upWO = mWO + semWO;
    loWO = mWO - semWO;

    % ============================================================
    % z and two-tailed p from mean/SEM
    % ============================================================
    zStim = mStim ./ semStim;
    zWO   = mWO   ./ semWO;

    pStim = 2 * (1 - normcdf(abs(zStim), 0, 1));
    pWO   = 2 * (1 - normcdf(abs(zWO),   0, 1));

    zStim(~isfinite(zStim)) = NaN;
    zWO(~isfinite(zWO))     = NaN;
    pStim(~isfinite(pStim)) = NaN;
    pWO(~isfinite(pWO))     = NaN;

    % ============================================================
    % find multiple top peaks on MEAN delta only
    % ============================================================
    getTopPeaksAbs = @(trace) local_getTopPeaksAbs(trace, nPeaksShow, minPeakDist);

    [locMeanStim, ~] = getTopPeaksAbs(mStim);
    [locMeanWO,   ~] = getTopPeaksAbs(mWO);

    % ============================================================
    % store metrics
    % ============================================================
    PeakMetrics.AllKept.(cat).Stim.time_ms  = time_vector(locMeanStim);
    PeakMetrics.AllKept.(cat).Stim.mean     = mStim(locMeanStim);
    PeakMetrics.AllKept.(cat).Stim.sem      = semStim(locMeanStim);
    PeakMetrics.AllKept.(cat).Stim.plusSEM  = upStim(locMeanStim);
    PeakMetrics.AllKept.(cat).Stim.minusSEM = loStim(locMeanStim);
    PeakMetrics.AllKept.(cat).Stim.z        = zStim(locMeanStim);
    PeakMetrics.AllKept.(cat).Stim.p        = pStim(locMeanStim);

    PeakMetrics.AllKept.(cat).WO.time_ms    = time_vector(locMeanWO);
    PeakMetrics.AllKept.(cat).WO.mean       = mWO(locMeanWO);
    PeakMetrics.AllKept.(cat).WO.sem        = semWO(locMeanWO);
    PeakMetrics.AllKept.(cat).WO.plusSEM    = upWO(locMeanWO);
    PeakMetrics.AllKept.(cat).WO.minusSEM   = loWO(locMeanWO);
    PeakMetrics.AllKept.(cat).WO.z          = zWO(locMeanWO);
    PeakMetrics.AllKept.(cat).WO.p          = pWO(locMeanWO);

    % ============================================================
    % print
    % ============================================================
    fprintf('\n[%s | AllKept delta]\n', cat);

    fprintf('Stim-BL mean peaks:\n');
    for k = 1:numel(locMeanStim)
        fprintf('  Peak %d: time = %.1f ms | mean = %.3f | SEM = %.3f | mean+SEM = %.3f | mean-SEM = %.3f | z = %.3f | p = %.5g\n', ...
            k, time_vector(locMeanStim(k)), mStim(locMeanStim(k)), semStim(locMeanStim(k)), ...
            upStim(locMeanStim(k)), loStim(locMeanStim(k)), zStim(locMeanStim(k)), pStim(locMeanStim(k)));
    end

    fprintf('WO-BL mean peaks:\n');
    for k = 1:numel(locMeanWO)
        fprintf('  Peak %d: time = %.1f ms | mean = %.3f | SEM = %.3f | mean+SEM = %.3f | mean-SEM = %.3f | z = %.3f | p = %.5g\n', ...
            k, time_vector(locMeanWO(k)), mWO(locMeanWO(k)), semWO(locMeanWO(k)), ...
            upWO(locMeanWO(k)), loWO(locMeanWO(k)), zWO(locMeanWO(k)), pWO(locMeanWO(k)));
    end

    % ============================================================
    % y limits
    % ============================================================
    yvals = [upStim, loStim, upWO, loWO];
    ymax = max(abs(yvals), [], 'omitnan');
    if isempty(ymax) || ~isfinite(ymax) || ymax == 0
        ymax = 1;
    end
    ymax = ymax * 1.12;

    % ============================================================
    % plot
    % ============================================================
    figure('Name',['ALL-KEPT DELTAS ' cat],'Position',[100,100,760,340]);
    hold on

    % SEM fills
    fill([time_vector fliplr(time_vector)], ...
         [upStim fliplr(loStim)], ...
         Color_StimDelta, 'FaceAlpha',0.20, 'EdgeColor','none');

    fill([time_vector fliplr(time_vector)], ...
         [upWO fliplr(loWO)], ...
         Color_WODelta, 'FaceAlpha',0.20, 'EdgeColor','none');

    % mean traces
    plot(time_vector, mStim, 'Color', Color_StimDelta, 'LineWidth', 2);
    plot(time_vector, mWO,   'Color', Color_WODelta,   'LineWidth', 2);

    % ---------------- mark Stim peaks ----------------
    for k = 1:numel(locMeanStim)
        xpk = time_vector(locMeanStim(k));

        y_mean  = mStim(locMeanStim(k));
        y_plus  = upStim(locMeanStim(k));
        y_minus = loStim(locMeanStim(k));
        thisSEM = semStim(locMeanStim(k));

        % mean marker
        plot(xpk, y_mean, 'o', ...
            'MarkerFaceColor', Color_StimDelta, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize', 6);

        % +SEM / -SEM markers at same timepoint
        plot(xpk, y_plus, '^', ...
            'MarkerFaceColor', Color_StimDelta, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize', 5);

        plot(xpk, y_minus, 'v', ...
            'MarkerFaceColor', Color_StimDelta, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize', 5);

        text(xpk, y_mean, ...
            sprintf('  S%d) %.0f, %.2f\\nSEM=%.2f, z=%.2f, p=%.3g', ...
            k, xpk, y_mean, thisSEM, zStim(locMeanStim(k)), pStim(locMeanStim(k))), ...
            'Color', Color_StimDelta, ...
            'FontSize', 8, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle');
    end

    % ---------------- mark WO peaks ----------------
    for k = 1:numel(locMeanWO)
        xpk = time_vector(locMeanWO(k));

        y_mean  = mWO(locMeanWO(k));
        y_plus  = upWO(locMeanWO(k));
        y_minus = loWO(locMeanWO(k));
        thisSEM = semWO(locMeanWO(k));

        % mean marker
        plot(xpk, y_mean, 'o', ...
            'MarkerFaceColor', Color_WODelta, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize', 6);

        % +SEM / -SEM markers at same timepoint
        plot(xpk, y_plus, '^', ...
            'MarkerFaceColor', Color_WODelta, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize', 5);

        plot(xpk, y_minus, 'v', ...
            'MarkerFaceColor', Color_WODelta, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize', 5);

        text(xpk, y_mean, ...
            sprintf('  W%d) %.0f, %.2f\\nSEM=%.2f, z=%.2f, p=%.3g', ...
            k, xpk, y_mean, thisSEM, zWO(locMeanWO(k)), pWO(locMeanWO(k))), ...
            'Color', Color_WODelta, ...
            'FontSize', 8, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle');
    end

    % reference lines
    xline(0,'--k');
    yline(0,'--k');
    xline(-10,'-k');
    yline(-2,':k');
    yline( 2,':k');

    xlim([-250 250]);
    ylim([-ymax ymax]);

    xlabel('Time from threshold crossing (ms)');
    ylabel('\Delta firing rate (Hz)');
    title(sprintf('ALL-KEPT \\DeltaFR: %s (N=%d)', cat, size(BL,1)), 'Interpreter','none');

    legend({'Stim SEM','','WO SEM','','Stim-BL','WO-BL'}, ...
        'Location','best','Box','off');

    box off
    set(gca,'TickDir','out')
end


%% ============================================================
% local function
% ============================================================
function [locsOut, pksOut] = local_getTopPeaksAbs(traceIn, nPeaksShow, minPeakDist)

traceTmp = abs(traceIn);
traceTmp(~isfinite(traceTmp)) = -inf;

[pks, locs] = findpeaks(traceTmp, 'MinPeakDistance', minPeakDist);

if isempty(locs)
    [~, locs] = max(traceTmp);
    pks = traceTmp(locs);
end

[~, ord] = sort(pks, 'descend');
locs = locs(ord);
pks  = pks(ord);

nShow = min(nPeaksShow, numel(locs));
locsOut = locs(1:nShow);
pksOut  = pks(1:nShow);

end
%% DELTA PLOTS (Stim-BL and WO-BL) BY NEG-LAG CATEGORY - Smoothed 
% USED IN THE PAPER 
% ============================================================

% ---- colors ----
Color_StimDelta = rgb('DodgerBlue');
Color_WODelta   = rgb('Magenta');

% ---- time window ----
%timevec= XCorrData(1).units(6).PSTH.timevec;

full_timevec = timevec; 
%full_timevec = linspace(-1000,1000,size(Agg.AllKept.baseline,2));
time_idx     = find(full_timevec >= -250 & full_timevec <= 250);
time_vector  = full_timevec(time_idx);

% ---- optional light smoothing for DISPLAY ONLY ----
smooth_bins = 3;   % try 5, 7, 9, or 11

smooth_plot = @(x) smoothdata(x, 'gaussian', smooth_bins);

% ---- helper: delta mean+sem ----
delta_sem = @(A,B) deal( ...
    nanmean(A(:,time_idx) - B(:,time_idx), 1), ...
    nanstd( A(:,time_idx) - B(:,time_idx), 0, 1) ./ sqrt(max(1,size(A,1))) );

group_names = {'negLag_posCorr','negLag_negCorr'};

%% ===================== (1) ALL KEPT =====================
for ci = 1:numel(group_names)
    cat = group_names{ci};
    % 
    BL = Agg.ByGroup.(cat).AllKept.baseline;
    ST = Agg.ByGroup.(cat).AllKept.stim;
    WO = Agg.ByGroup.(cat).AllKept.wash;
    % BL = Agg.ByGroup.(cat).ModKept.baseline;
    % ST = Agg.ByGroup.(cat).ModKept.stim;
    % WO = Agg.ByGroup.(cat).ModKept.wash;

    if isempty(BL) || isempty(ST) || isempty(WO)
        fprintf('AllKept %s empty -> skipping.\n', cat);
        continue;
    end

    % ---- raw deltas ----
    [mStim, semStim] = delta_sem(ST, BL);
    [mWO,   semWO  ] = delta_sem(WO, BL);

    % ---- smooth for plotting only ----
    mStim_plot   = smooth_plot(mStim);
    semStim_plot = smooth_plot(semStim);

    mWO_plot     = smooth_plot(mWO);
    semWO_plot   = smooth_plot(semWO);

    figure('Name',['ALL-KEPT DELTAS ' cat],'Position',[100,100,650,300]); 
    hold on;

    % SEM first so lines sit on top
    fill([time_vector fliplr(time_vector)], ...
         [mStim_plot + semStim_plot, fliplr(mStim_plot - semStim_plot)], ...
         Color_StimDelta, 'FaceAlpha',0.20, 'EdgeColor','none');

    fill([time_vector fliplr(time_vector)], ...
         [mWO_plot + semWO_plot, fliplr(mWO_plot - semWO_plot)], ...
         Color_WODelta, 'FaceAlpha',0.20, 'EdgeColor','none');

    plot(time_vector, mStim_plot, 'Color', Color_StimDelta, 'LineWidth', 2);
    plot(time_vector, mWO_plot,   'Color', Color_WODelta,   'LineWidth', 2);

    xline(0,'--k'); 
    yline(0,'--k');
    yline(-2,':'); 
    yline(2,':');
    xline(-10)
    xline(81.5)
    
     %  ylim([-5 5])
     % yticks(-6:2:6);
   ylim([-10 8]);
    yticks(-10:2:10);
    

    xlabel('Time from threshold crossing (ms)');
    ylabel('\Delta Z (Condition - Baseline)');
    xlim([-150 250]);

    title(sprintf('ALL-KEPT Δ: %s (N=%d)', cat, size(BL,1)), 'Interpreter','none');
    legend({'Stim-BL','WO-BL'}, 'Location','best','Box','off');
    box off;
    set(gca,'TickDir','out');
end

%% Combined plot with accel
%% ============================================================
% COMBINED PLOT:
%   LEFT AXIS  = Neural Stim-BL deltas for:
%                  (1) negLag_posCorr
%                  (2) negLag_negCorr
%   RIGHT AXIS = Acceleration delta overlaid
%
% NOTES:
%   - Uses your smoothed neural display traces
%   - Replace colAcc with the correct acceleration column
%   - Replace ACC_ST and ACC_BL with the behavior structs you want
%     (e.g., AD_big and BL_big, or Stim_big and BL_big)
% ============================================================

% -------------------- COLORS --------------------
Color_pos = rgb('SteelBlue');
Color_neg = rgb('DarkGreen');
Color_acc = [0 0 0];

% -------------------- SETTINGS ------------------
xlim_ms = [-260 260];

smooth_bins = 3;   % display only

% --- choose behavior structs for acceleration delta ---
ACC_ST = AD_big;   % <-- replace if needed
ACC_BL = BL_big;   % <-- replace if needed

% --- choose acceleration column ---
colAcc = 5;        % <-- REPLACE with correct acceleration column

% -------------------- TIME WINDOW ----------------
% full_timevec = timevec;
% time_idx     = find(full_timevec >= xlim_ms(1) & full_timevec <= xlim_ms(2));
% time_vector  = full_timevec(time_idx);

% % % ---- time window ----
% full_timevec = timevec; 
% % %full_timevec = linspace(-1000,1000,size(Agg.AllKept.baseline,2));
%  time_idx     = find(full_timevec >= -250 & full_timevec <= 250);
%  time_vector  = full_timevec(time_idx);

% xlim_ms = [-250 250];
% psth_timevec = time_bins(:)' - bin_size/2;
% 


timevec = XCorrData(1).units(6).PSTH.timevec(:)';
time_idx = find(timevec >= xlim_ms(1) & timevec <= xlim_ms(2));
time_vector = timevec(time_idx);


smooth_plot = @(x) smoothdata(x, 'gaussian', smooth_bins);

% % helper: paired delta mean+sem
% delta_sem = @(A,B) deal( ...
%     nanmean(A(:,time_idx) - B(:,time_idx), 1), ...
%     nanstd(A(:,time_idx) - B(:,time_idx), 0, 1) ./ sqrt(sum(~isnan(A(:,time_idx) - B(:,time_idx)),1)) );

delta_sem = @(A,B) deal( ...
    nanmean(A(:,time_idx) - B(:,time_idx), 1), ...
    nanstd(A(:,time_idx) - B(:,time_idx), 0, 1) ./ ...
    sqrt(sum(~isnan(A(:,time_idx) - B(:,time_idx)),1)) );
%
% ============================================================
% PREP NEURAL DATA: Stim-BL only, both categories on same plot
% ============================================================
cats   = {'negLag_posCorr','negLag_negCorr'};
colors = {Color_pos, Color_neg};
labels = {'+Corr Stim-BL','-Corr Stim-BL'};

neural_mean = cell(1,numel(cats));
neural_sem  = cell(1,numel(cats));
neural_n    = nan(1,numel(cats));
keep_cat    = false(1,numel(cats));

for ci = 1:numel(cats)
    cat = cats{ci};

    BL = Agg.ByGroup.(cat).AllKept.baseline;
    ST = Agg.ByGroup.(cat).AllKept.stim;

    if isempty(BL) || isempty(ST)
        fprintf('AllKept %s empty -> skipping.\n', cat);
        continue;
    end

    [mStim, semStim] = delta_sem(ST, BL);

    neural_mean{ci} = smooth_plot(mStim);
    neural_sem{ci}  = smooth_plot(semStim);
    neural_n(ci)    = size(BL,1);
    keep_cat(ci)    = true;
end

% ============================================================
% PREP ACCELERATION DELTA
% ============================================================
diff_acc = ACC_ST.Average(:,colAcc) - ACC_BL.Average(:,colAcc);
diff_acc = diff_acc * 10;

SE_acc = sqrt( ...
    (ACC_ST.Std(:,colAcc).^2 ./ ACC_ST.N) + ...
    (ACC_BL.Std(:,colAcc).^2 ./ ACC_BL.N) );

time_s      = ACC_BL.Average(:,1);
time_ms_acc = time_s(:)' * 1000;

diff_acc = diff_acc(:)';
SE_acc   = SE_acc(:)';

acc_keep = (time_ms_acc >= xlim_ms(1)) & (time_ms_acc <= xlim_ms(2));
time_acc = time_ms_acc(acc_keep);
mAcc     = diff_acc(acc_keep);
semAcc   = SE_acc(acc_keep);

% ============================================================
% FIGURE
% ============================================================
fig = figure('Color','w','Position',[100,100,650,300]);
set(fig,'Renderer','painters');

% ---------------- RIGHT AXIS: ACCEL ----------------
ax2 = axes('Parent',fig);
hold(ax2,'on');
box(ax2,'off');

fill(ax2, ...
    [time_acc fliplr(time_acc)], ...
    [mAcc + semAcc, fliplr(mAcc - semAcc)], ...
    Color_acc, ...
    'FaceAlpha',0.15, ...
    'EdgeColor','none');

plot(ax2, time_acc, mAcc, 'Color', Color_acc, 'LineWidth', 2);

yline(ax2, 0, '--k', 'LineWidth', 1);
xline(ax2, 0, '--k', 'LineWidth', 1);

ax2.YAxisLocation = 'right';
ax2.XLim          = [-250 250];

ax2.XColor        = 'k';
ax2.YColor        = 'k';
ax2.LineWidth     = 1;
ax2.FontSize      = 11;
ax2.TickDir       = 'out';
ax2.XTickLabel    = [];
ax2.XLabel.String = '';
ylabel(ax2, '\Delta acceleration');

% ---------------- LEFT AXIS: NEURAL ----------------
ax1 = axes('Parent',fig);
hold(ax1,'on');
box(ax1,'off');

legend_handles = [];
legend_labels  = {};

for ci = 1:numel(cats)
    if ~keep_cat(ci)
        continue;
    end

    mPlot   = neural_mean{ci};
    semPlot = neural_sem{ci};

    fill(ax1, ...
        [time_vector fliplr(time_vector)], ...
        [mPlot + semPlot, fliplr(mPlot - semPlot)], ...
        colors{ci}, ...
        'FaceAlpha',0.20, ...
        'EdgeColor','none');

    h = plot(ax1, time_vector, mPlot, ...
        'Color', colors{ci}, 'LineWidth', 2);

    legend_handles(end+1) = h; %#ok<SAGROW>
    legend_labels{end+1}  = labels{ci}; %#ok<SAGROW>
end

xline(ax1, 0, '--k', 'LineWidth', 1);
yline(ax1, 0, '--k', 'LineWidth', 1);
yline(ax1, -2, ':k', 'LineWidth', 1);
yline(ax1,  2, ':k', 'LineWidth', 1);
xline(-10)
xline(81.5)

ax1.Color     = 'none';
ax1.XLim      = [-250 250];
ax1.YLim= ([-5 5]);
ax1.YTick= -6:2:6;
ax1.XColor    = 'k';
ax1.YColor    = 'k';
ax1.LineWidth = 1;
ax1.FontSize  = 11;
ax1.TickDir   = 'out';
xlabel(ax1, 'Time from threshold crossing (ms)');
ylabel(ax1, '\Delta Z (Stim - BL)');

% ---------------- MATCH AXES ----------------
ax1.Position = ax2.Position;
linkaxes([ax1 ax2],'x');

drawnow;

yl1 = ylim(ax1);
M1 = max(abs(yl1));
if ~isfinite(M1) || M1 == 0
    M1 = 1;
end
ylim(ax1, [-M1 M1]);

yl2 = ylim(ax2);
M2 = max(abs(yl2));
if ~isfinite(M2) || M2 == 0
    M2 = 1;
end
ylim(ax2, [-M2 M2]);

ax1.Position = ax2.Position;

if ~isempty(legend_handles)
    legend(ax1, legend_handles, legend_labels, ...
        'Location','best', 'Box','off');
end

title(ax1, 'Stim-BL neural deltas (+Corr / -Corr) with acceleration overlay');

% ---------------- OPTIONAL EXPORT ----------------
% print(fig, fullfile(output_folder,'StimBL_posNegCorr_withAccel.pdf'), '-dpdf', '-painters');

% ============================================================
% GET PEAKS FOR +Corr AND -Corr UNITS
% Uses the smoothed Stim-BL neural traces already computed above
% ============================================================

fprintf('\n==================== PEAK SUMMARY ====================\n');

for ci = 1:numel(cats)
    if ~keep_cat(ci)
        fprintf('%s: no data\n', cats{ci});
        continue;
    end

    trace = neural_mean{ci};   % smoothed Stim-BL trace
    t     = time_vector;       % matching time vector (ms)

    % --- positive peak ---
    [peakMaxVal, idxMax] = max(trace);
    peakMaxTime = t(idxMax);

    % --- negative peak ---
    [peakMinVal, idxMin] = min(trace);
    peakMinTime = t(idxMin);

    % --- absolute peak (largest magnitude, regardless of sign) ---
    [~, idxAbs] = max(abs(trace));
    peakAbsVal  = trace(idxAbs);
    peakAbsTime = t(idxAbs);

    fprintf('\nCategory: %s\n', cats{ci});
    fprintf('  Positive peak:  %.3f at %.1f ms\n', peakMaxVal, peakMaxTime);
    fprintf('  Negative peak:  %.3f at %.1f ms\n', peakMinVal, peakMinTime);
    fprintf('  Absolute peak:  %.3f at %.1f ms\n', peakAbsVal, peakAbsTime);
end
%% just peaks 
%% ============================================================
% POPULATION PEAKS FROM AVERAGE TRACE (50 to 250 ms)
% ============================================================

win = [50 250];   % ms
win_idx = time_vector >= win(1) & time_vector <= win(2);
t_win   = time_vector(win_idx);

fprintf('\n==================== POPULATION PEAKS (50–250 ms) ====================\n');

for ci = 1:numel(cats)
    if ~keep_cat(ci)
        fprintf('%s: no data\n', cats{ci});
        continue;
    end

    trace_win = neural_mean{ci}(win_idx);

    % positive peak
    [peakVal, idxPeak] = max(trace_win);
    peakTime = t_win(idxPeak);

    % negative peak
    [troughVal, idxTrough] = min(trace_win);
    troughTime = t_win(idxTrough);

    fprintf('\nCategory: %s\n', cats{ci});
    fprintf('  Positive peak: %.3f at %.1f ms\n', peakVal, peakTime);
    fprintf('  Negative peak: %.3f at %.1f ms\n', troughVal, troughTime);
end
%% ============================================================
% ALL KEPT: Single plot with AD-BL deltas for BOTH categories
%   - negLag_posCorr  -> SteelBlue
%   - negLag_negCorr  -> DarkGreen
% ============================================================

Color_pos = rgb('SteelBlue');
Color_neg = rgb('DarkGreen');

figure('Name','ALL-KEPT AD-BL (both categories)','Position',[100,100,650,300]); 
hold on;

cats = {'negLag_posCorr','negLag_negCorr'};
colors = {Color_pos, Color_neg};

for ci = 1:numel(cats)

    cat = cats{ci};

    BL = Agg.ByGroup.(cat).AllKept.baseline;
    ST = Agg.ByGroup.(cat).AllKept.stim;

    % BL = Agg.ByGroup.(cat).ModKept.baseline;
    % ST = Agg.ByGroup.(cat).ModKept.stim;

    if isempty(BL) || isempty(ST)
        fprintf('AllKept %s empty -> skipping.\n', cat);
        continue;
    end

    % ---- deltas (AD - BL) ----
    [mStim, semStim] = delta_sem(ST, BL);

    % plot mean
    plot(time_vector, mStim, 'Color', colors{ci}, 'LineWidth', 2);

    % SEM shading
    fill([time_vector fliplr(time_vector)], ...
         [mStim+semStim fliplr(mStim-semStim)], ...
         colors{ci}, 'FaceAlpha',0.20, 'EdgeColor','none');

end

% common cosmetics
xline(0,'--k'); 
yline(0,'--k');
yline(-2,':'); 
yline(2,':');
xline(81.5)

xlabel('Time from threshold crossing (ms)');
ylabel('\Delta Z (AD - BL)');

xlim([-250 250]);
ylim([-8 8])

title('ALL-KEPT Δ (AD - BL) : posCorr vs negCorr');

legend({'posCorr AD-BL','','negCorr AD-BL',''}, ...
       'Location','best','Box','off');

box off;

%%
% For example, just plot the y-dimension (2nd column)
diff_avg_y = AD_big.Average(:,5) - BL_big.Average(:,5);
diff_avg_y=diff_avg_y*10;
SE_diff_y = sqrt((AD_big.Std(:,5).^2 ./ AD_big.N) + (BL_big.Std(:,5).^2 ./ BL_big.N));

% diff_avg_y = RR_big.Average(:,5) - NRR_big.Average(:,5);
% diff_avg_y=diff_avg_y*10;
% SE_diff_y = sqrt((RR_big.Std(:,5).^2 ./ RR_big.N) + (NRR_big.Std(:,5).^2 ./ NRR_big.N));

time = (BL_big.Average(:,1))'; % ensure it's a column vector

% Now plot
figure('Position', [100, 100, 600, 200]);  
shadedErrorBar(time, diff_avg_y, SE_diff_y, 'lineProps', {'k', 'LineWidth', 2});
xlabel('Time from threshold crossing (s)');
xlim([-0.25 0.25])
%xlim('auto')
ylim([-8 8])
ylabel('Y-Velocity Difference (CSM - NoCSM)');
title('CSM vs NoCSM Kinematic Difference (Y-Direction)');
yline(0, '--');
xline(-0.055)
xline(0.087)
xline(0,'--')

%% COMBINED

%% ============================================================
% COMBINED DELTA PLOT: NEURAL (left y) + VELOCITY (right y)
% AD-BL neural deltas for BOTH categories + AD-BL behavior velocity
% Zero aligned across both y-axes
% Illustrator-editable (painters renderer)
% ============================================================

% -------------------------- COLORS --------------------------
Color_pos = rgb('SteelBlue');
Color_neg = rgb('DarkGreen');
Color_beh = [0 0 0];

% -------------------------- SETTINGS ------------------------
xlim_ms = [-250 250];

% behavior event lines in seconds -> ms
xlines_behavior_s  = [-0.055, 0.0815, 0];
xlines_behavior_ms = xlines_behavior_s * 1000;

do_neural_ref_lines = true;
colY = 5;   % velocity column

% ------------------ PREP BEHAVIOR DELTA ---------------------
% AD - BL behavior delta
diff_avg_y = AD_big.Average(:,colY) - BL_big.Average(:,colY);
diff_avg_y = diff_avg_y * 10;

SE_diff_y = sqrt( ...
    (AD_big.Std(:,colY).^2 ./ AD_big.N) + ...
    (BL_big.Std(:,colY).^2 ./ BL_big.N) );

time_s  = BL_big.Average(:,1);
time_ms = time_s(:)' * 1000;

diff_avg_y = diff_avg_y(:)';
SE_diff_y  = SE_diff_y(:)';

% window behavior to x-lims
beh_keep  = (time_ms >= xlim_ms(1)) & (time_ms <= xlim_ms(2));
time_ms_b = time_ms(beh_keep);
diff_b    = diff_avg_y(beh_keep);
SE_b      = SE_diff_y(beh_keep);

% -------------------------- FIGURE --------------------------
fig = figure('Name','ALL-KEPT AD-BL Neural + Velocity', ...
             'Color','w', ...
             'Position',[120,120,650,320]);
set(fig,'Renderer','painters');

% ============================================================
% AX2: BEHAVIOR / VELOCITY (RIGHT Y-AXIS) FIRST
% ============================================================
ax2 = axes('Parent',fig);
hold(ax2,'on');
box(ax2,'off');

% SEM patch
fill(ax2, ...
    [time_ms_b fliplr(time_ms_b)], ...
    [diff_b + SE_b, fliplr(diff_b - SE_b)], ...
    Color_beh, ...
    'FaceAlpha',0.18, ...
    'EdgeColor','none');

% mean trace
plot(ax2, time_ms_b, diff_b, 'Color', Color_beh, 'LineWidth', 2);

% behavior reference lines
yline(ax2, 0, '--k', 'LineWidth', 1);
for xl = xlines_behavior_ms
    xline(ax2, xl, '--k', 'LineWidth', 1);
end

% behavior axis formatting
ax2.YAxisLocation = 'right';
ax2.XLim          = xlim_ms;
ax2.XColor        = 'k';
ax2.YColor        = 'k';
ax2.LineWidth     = 1;
ax2.FontSize      = 11;
ax2.TickDir       = 'out';
ax2.XTickLabel    = [];
ax2.XLabel.String = '';
ylabel(ax2, 'Y-velocity \Delta (AD - BL)');

% ============================================================
% AX1: NEURAL (LEFT Y-AXIS) SECOND
% ============================================================
ax1 = axes('Parent',fig);
hold(ax1,'on');
box(ax1,'off');

cats   = {'negLag_posCorr','negLag_negCorr'};
colors = {Color_pos, Color_neg};

legend_handles = [];
legend_labels  = {};

for ci = 1:numel(cats)

    cat = cats{ci};

    BL = Agg.ByGroup.(cat).AllKept.baseline;
    ST = Agg.ByGroup.(cat).AllKept.stim;

    if isempty(BL) || isempty(ST)
        fprintf('AllKept %s empty -> skipping.\n', cat);
        continue;
    end

    % AD - BL neural delta
    [mStim, semStim] = delta_sem(ST, BL);

    % SEM patch
    fill(ax1, ...
        [time_vector fliplr(time_vector)], ...
        [mStim + semStim, fliplr(mStim - semStim)], ...
        colors{ci}, ...
        'FaceAlpha',0.20, ...
        'EdgeColor','none');

    % mean trace
    h = plot(ax1, time_vector, mStim, ...
        'Color', colors{ci}, 'LineWidth', 2);

    legend_handles(end+1) = h; %#ok<SAGROW>

    if strcmp(cat,'negLag_posCorr')
        legend_labels{end+1} = 'posCorr AD-BL'; %#ok<SAGROW>
    else
        legend_labels{end+1} = 'negCorr AD-BL'; %#ok<SAGROW>
    end
end

% neural reference lines
if do_neural_ref_lines
    xline(ax1, 0, '--k', 'LineWidth', 1);
    xline(ax1, 81.5, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1);
    yline(ax1, 0, '--k', 'LineWidth', 1);
    yline(ax1, -2, ':k', 'LineWidth', 1);
    yline(ax1,  2, ':k', 'LineWidth', 1);
end

% neural axis formatting
ax1.Color     = 'none';
ax1.XLim      = xlim_ms;
ax1.XColor    = 'k';
ax1.YColor    = 'k';
ax1.LineWidth = 1;
ax1.FontSize  = 11;
ax1.TickDir   = 'out';
xlabel(ax1, 'Time from threshold crossing (ms)');
ylabel(ax1, '\Delta Z (AD - BL)');

% match positions exactly and link x
ax1.Position = ax2.Position;
linkaxes([ax1 ax2], 'x');

% ============================================================
% FORCE ZERO TO BE AT THE SAME HEIGHT ON BOTH AXES
% ============================================================
drawnow;

yl1 = ylim(ax1);
M1  = max(abs(yl1));
if M1 == 0
    M1 = 1;
end
ylim(ax1, [-M1 M1]);

yl2 = ylim(ax2);
M2  = max(abs(yl2));
if M2 == 0
    M2 = 1;
end
ylim(ax2, [-M2 M2]);

% optional fixed velocity ylim if you want exactly what you had before:
% ylim(ax2, [-8 8]);

% re-apply exact matching position
ax1.Position = ax2.Position;

% legend
if ~isempty(legend_handles)
    legend(ax1, legend_handles, legend_labels, ...
        'Location','best', ...
        'Box','off');
end

title(ax1, 'ALL-KEPT \Delta Neural (AD-BL, left) + \Delta Velocity (AD-BL, right)');

% ---------------- OPTIONAL EXPORT ----------------
% if ~exist(output_folder,'dir'); mkdir(output_folder); end
% outFile = fullfile(output_folder, 'AllKept_ADminusBL_NeuralPlusVelocity.pdf');
% print(fig, outFile, '-dpdf', '-painters');
%% ============================================================
% NET ACTIVITY: ALL KEPT (AD - BL) combining posCorr and negCorr as "vectors"
%
% Goal:
%   - Compute the AD-BL delta trace for posCorr units and negCorr units
%   - Then compute a NET trace that reflects "vector" combination:
%         NET = (posCorr delta) - (negCorr delta)
%     (i.e., opposite-signed populations cancel)
%
% Also plots:
%   - posCorr delta (SteelBlue)
%   - negCorr delta (DarkGreen)
%   - NET (black)
%
% Assumes:
%   - Agg.ByGroup.negLag_posCorr.AllKept.baseline / .stim exist
%   - Agg.ByGroup.negLag_negCorr.AllKept.baseline / .stim exist
%   - time_vector exists (same length as PSTHs)
%   - delta_sem(ST,BL) returns [meanDelta, semDelta] for (ST - BL)
% ============================================================

Color_pos = rgb('SteelBlue');
Color_neg = rgb('DarkGreen');
Color_net = [0 0 0];

figure('Name','ALL-KEPT NET (AD-BL): pos/neg + net','Position',[100,100,750,330]);
hold on;

% --------------------- pull arrays ---------------------
BL_pos = Agg.ByGroup.negLag_posCorr.AllKept.baseline;
ST_pos = Agg.ByGroup.negLag_posCorr.AllKept.stim;

BL_neg = Agg.ByGroup.negLag_negCorr.AllKept.baseline;
ST_neg = Agg.ByGroup.negLag_negCorr.AllKept.stim;

% --------------------- sanity checks ---------------------
if isempty(BL_pos) || isempty(ST_pos)
    error('AllKept posCorr arrays are empty.');
end
if isempty(BL_neg) || isempty(ST_neg)
    error('AllKept negCorr arrays are empty.');
end

% --------------------- deltas per group ---------------------
[mPos, semPos] = delta_sem(ST_pos, BL_pos);  % AD - BL for posCorr
[mNeg, semNeg] = delta_sem(ST_neg, BL_neg);  % AD - BL for negCorr

% --------------------- NET vector combination ---------------------
% Interpret negCorr population as contributing opposite direction to posCorr
mNet = mPos - mNeg;

% SEM propagation assuming independence:
% var(net) = var(pos) + var(neg)
semNet = sqrt(semPos.^2 + semNeg.^2);

% --------------------- plot posCorr ---------------------
plot(time_vector, mPos, 'Color', Color_pos, 'LineWidth', 2);
fill([time_vector fliplr(time_vector)], ...
     [mPos+semPos fliplr(mPos-semPos)], ...
     Color_pos, 'FaceAlpha',0.20, 'EdgeColor','none');

% --------------------- plot negCorr ---------------------
plot(time_vector, mNeg, 'Color', Color_neg, 'LineWidth', 2);
fill([time_vector fliplr(time_vector)], ...
     [mNeg+semNeg fliplr(mNeg-semNeg)], ...
     Color_neg, 'FaceAlpha',0.20, 'EdgeColor','none');

% --------------------- plot NET ---------------------
plot(time_vector, mNet, 'Color', Color_net, 'LineWidth', 2.5);
fill([time_vector fliplr(time_vector)], ...
     [mNet+semNet fliplr(mNet-semNet)], ...
     Color_net, 'FaceAlpha',0.10, 'EdgeColor','none');

% --------------------- cosmetics ---------------------
xline(0,'--k');
yline(0,'--k');
yline(-2,':'); yline(2,':');

% keep your extra reference if you want
xline(81.5);

xlabel('Time from threshold crossing (ms)');
ylabel('\Delta Z (AD - BL)');
xlim([-250 250]);

title('ALL-KEPT \Delta (AD - BL): posCorr, negCorr, and NET = pos - neg');

legend({'posCorr mean','posCorr SEM','negCorr mean','negCorr SEM','NET mean','NET SEM'}, ...
       'Location','best','Box','off');

box off;
%% ============================================================
% BASELINE NET (difference of means) within time_idx window
% time_vector is already full_timevec(time_idx)
% ============================================================

Color_pos = rgb('SteelBlue');
Color_neg = rgb('DarkGreen');
Color_net = [0 0 0];

figure('Name','ALL-KEPT BASELINE NET (windowed): mean(+Corr)-mean(-Corr)','Position',[100,100,750,330]);
hold on;

% --------------------- pull arrays ---------------------
BL_pos = Agg.ByGroup.negLag_posCorr.AllKept.baseline;   % [nPos x Tfull]
BL_neg = Agg.ByGroup.negLag_negCorr.AllKept.baseline;   % [nNeg x Tfull]

% --------------------- sanity checks ---------------------
if isempty(BL_pos) || isempty(BL_neg)
    error('Baseline arrays are empty (BL_pos or BL_neg).');
end
if ~exist('time_idx','var') || isempty(time_idx)
    error('time_idx missing/empty.');
end
if ~exist('time_vector','var') || isempty(time_vector)
    error('time_vector missing/empty.');
end

% --------------------- compute within window only ---------------------
A = BL_pos(:, time_idx);
B = BL_neg(:, time_idx);

% means
mPos = mean(A, 1, 'omitnan');
mNeg = mean(B, 1, 'omitnan');

% SEMs (across units)
nPos   = sum(~isnan(A), 1);
nNeg   = sum(~isnan(B), 1);

semPos = std(A, 0, 1, 'omitnan') ./ sqrt(max(nPos,1));
semNeg = std(B, 0, 1, 'omitnan') ./ sqrt(max(nNeg,1));

% NET = difference of means
mNet   = mPos - mNeg;
semNet = sqrt(semPos.^2 + semNeg.^2);   % independence propagation

% --------------------- plot posCorr ---------------------
plot(time_vector, mPos, 'Color', Color_pos, 'LineWidth', 2);
fill([time_vector fliplr(time_vector)], ...
     [mPos+semPos fliplr(mPos-semPos)], ...
     Color_pos, 'FaceAlpha',0.20, 'EdgeColor','none');

% --------------------- plot negCorr ---------------------
plot(time_vector, mNeg, 'Color', Color_neg, 'LineWidth', 2);
fill([time_vector fliplr(time_vector)], ...
     [mNeg+semNeg fliplr(mNeg-semNeg)], ...
     Color_neg, 'FaceAlpha',0.20, 'EdgeColor','none');

% --------------------- plot NET ---------------------
plot(time_vector, mNet, 'Color', Color_net, 'LineWidth', 2.5);
fill([time_vector fliplr(time_vector)], ...
     [mNet+semNet fliplr(mNet-semNet)], ...
     Color_net, 'FaceAlpha',0.10, 'EdgeColor','none');

% --------------------- cosmetics ---------------------
xline(0,'--k');
yline(0,'--k');
yline(-2,':'); yline(2,':');
xline(81.5);

xlabel('Time from threshold crossing (ms)');
ylabel('\Delta Z (Baseline: mean +Corr - mean -Corr)');
xlim([-250 250]);

title('ALL-KEPT BASELINE (windowed): +Corr, -Corr, NET = mean(+Corr) - mean(-Corr)');

legend({'posCorr mean','posCorr SEM','negCorr mean','negCorr SEM','NET mean','NET SEM'}, ...
       'Location','best','Box','off');

box off;
%% ============================================================
% STIM NET (difference of means) within time_idx window
% time_vector is already full_timevec(time_idx)
% ============================================================

Color_pos = rgb('SteelBlue');
Color_neg = rgb('DarkGreen');
Color_net = [0 0 0];

figure('Name','ALL-KEPT STIM NET (windowed): mean(+Corr)-mean(-Corr)','Position',[100,100,750,330]);
hold on;

% --------------------- pull arrays ---------------------
ST_pos = Agg.ByGroup.negLag_posCorr.AllKept.stim;   % [nPos x Tfull]
ST_neg = Agg.ByGroup.negLag_negCorr.AllKept.stim;   % [nNeg x Tfull]

% --------------------- sanity checks ---------------------
if isempty(ST_pos) || isempty(ST_neg)
    error('Stim arrays are empty (ST_pos or ST_neg).');
end
if ~exist('time_idx','var') || isempty(time_idx)
    error('time_idx missing/empty.');
end
if ~exist('time_vector','var') || isempty(time_vector)
    error('time_vector missing/empty.');
end

% --------------------- compute within window only ---------------------
A = ST_pos(:, time_idx);
B = ST_neg(:, time_idx);

% means
mPos = mean(A, 1, 'omitnan');
mNeg = mean(B, 1, 'omitnan');

% SEMs (across units)
nPos   = sum(~isnan(A), 1);
nNeg   = sum(~isnan(B), 1);

semPos = std(A, 0, 1, 'omitnan') ./ sqrt(max(nPos,1));
semNeg = std(B, 0, 1, 'omitnan') ./ sqrt(max(nNeg,1));

% NET = difference of means
mNet   = mPos - mNeg;
semNet = sqrt(semPos.^2 + semNeg.^2);   % independence propagation

% --------------------- plot posCorr ---------------------
plot(time_vector, mPos, 'Color', Color_pos, 'LineWidth', 2);
fill([time_vector fliplr(time_vector)], ...
     [mPos+semPos fliplr(mPos-semPos)], ...
     Color_pos, 'FaceAlpha',0.20, 'EdgeColor','none');

% --------------------- plot negCorr ---------------------
plot(time_vector, mNeg, 'Color', Color_neg, 'LineWidth', 2);
fill([time_vector fliplr(time_vector)], ...
     [mNeg+semNeg fliplr(mNeg-semNeg)], ...
     Color_neg, 'FaceAlpha',0.20, 'EdgeColor','none');

% --------------------- plot NET ---------------------
plot(time_vector, mNet, 'Color', Color_net, 'LineWidth', 2.5);
fill([time_vector fliplr(time_vector)], ...
     [mNet+semNet fliplr(mNet-semNet)], ...
     Color_net, 'FaceAlpha',0.10, 'EdgeColor','none');

% --------------------- cosmetics ---------------------
xline(0,'--k');
yline(0,'--k');
yline(-2,':'); yline(2,':');
xline(81.5);

xlabel('Time from threshold crossing (ms)');
ylabel('\Delta Z (Stim: mean +Corr - mean -Corr)');
xlim([-250 250]);

title('ALL-KEPT STIM (windowed): +Corr, -Corr, NET = mean(+Corr) - mean(-Corr)');

legend({'posCorr mean','posCorr SEM','negCorr mean','negCorr SEM','NET mean','NET SEM'}, ...
       'Location','best','Box','off');

box off;

%% ============================================================
% WO NET (difference of means) within time_idx window
% time_vector is already full_timevec(time_idx)
% ============================================================

Color_pos = rgb('SteelBlue');
Color_neg = rgb('DarkGreen');
Color_net = [0 0 0];

figure('Name','ALL-KEPT WO NET (windowed): mean(+Corr)-mean(-Corr)','Position',[100,100,750,330]);
hold on;

% --------------------- pull arrays ---------------------
WO_pos = Agg.ByGroup.negLag_posCorr.AllKept.wash;   % [nPos x Tfull]
WO_neg = Agg.ByGroup.negLag_negCorr.AllKept.wash;   % [nNeg x Tfull]

% --------------------- sanity checks ---------------------
if isempty(WO_pos) || isempty(WO_neg)
    error('WO arrays are empty (WO_pos or WO_neg).');
end
if ~exist('time_idx','var') || isempty(time_idx)
    error('time_idx missing/empty.');
end
if ~exist('time_vector','var') || isempty(time_vector)
    error('time_vector missing/empty.');
end

% --------------------- compute within window only ---------------------
A = WO_pos(:, time_idx);
B = WO_neg(:, time_idx);

% means
mPos = mean(A, 1, 'omitnan');
mNeg = mean(B, 1, 'omitnan');

% SEMs (across units)
nPos   = sum(~isnan(A), 1);
nNeg   = sum(~isnan(B), 1);

semPos = std(A, 0, 1, 'omitnan') ./ sqrt(max(nPos,1));
semNeg = std(B, 0, 1, 'omitnan') ./ sqrt(max(nNeg,1));

% NET = difference of means
mNet   = mPos - mNeg;
semNet = sqrt(semPos.^2 + semNeg.^2);   % independence propagation

% --------------------- plot posCorr ---------------------
plot(time_vector, mPos, 'Color', Color_pos, 'LineWidth', 2);
fill([time_vector fliplr(time_vector)], ...
     [mPos+semPos fliplr(mPos-semPos)], ...
     Color_pos, 'FaceAlpha',0.20, 'EdgeColor','none');

% --------------------- plot negCorr ---------------------
plot(time_vector, mNeg, 'Color', Color_neg, 'LineWidth', 2);
fill([time_vector fliplr(time_vector)], ...
     [mNeg+semNeg fliplr(mNeg-semNeg)], ...
     Color_neg, 'FaceAlpha',0.20, 'EdgeColor','none');

% --------------------- plot NET ---------------------
plot(time_vector, mNet, 'Color', Color_net, 'LineWidth', 2.5);
fill([time_vector fliplr(time_vector)], ...
     [mNet+semNet fliplr(mNet-semNet)], ...
     Color_net, 'FaceAlpha',0.10, 'EdgeColor','none');

% --------------------- cosmetics ---------------------
xline(0,'--k');
yline(0,'--k');
yline(-2,':'); yline(2,':');
xline(81.5);

xlabel('Time from threshold crossing (ms)');
ylabel('\Delta Z (WO: mean +Corr - mean -Corr)');
xlim([-250 250]);

title('ALL-KEPT WO (windowed): +Corr, -Corr, NET = mean(+Corr) - mean(-Corr)');

legend({'posCorr mean','posCorr SEM','negCorr mean','negCorr SEM','NET mean','NET SEM'}, ...
       'Location','best','Box','off');

box off;
%% ============================================================
% NET OVERLAY: Baseline (black), Stim (DodgerBlue), Washout (magenta)
%   NET = mean(+Corr) - mean(-Corr)   (within time_idx window)
%
% Assumes:
%   - Agg.ByGroup.negLag_posCorr.AllKept.(baseline/stim/wash) exist
%   - Agg.ByGroup.negLag_negCorr.AllKept.(baseline/stim/wash) exist
%   - time_idx + time_vector already defined (windowed)
% ============================================================

Color_BL = [0 0 0];                 % black
Color_ST = rgb('DodgerBlue');       % dodger blue
Color_WO = rgb('Magenta');      % magenta

figure('Name','NET Overlay: BL vs ST vs WO','Position',[100,100,780,330]);
hold on;

% --------------------- sanity checks ---------------------
if ~exist('time_idx','var') || isempty(time_idx),     error('time_idx missing/empty.'); end
if ~exist('time_vector','var') || isempty(time_vector), error('time_vector missing/empty.'); end

% --------------------- helper: NET from pos/neg arrays ---------------------
net_from_groups = @(POS,NEG) deal( ...
    mean(POS(:,time_idx),1,'omitnan') - mean(NEG(:,time_idx),1,'omitnan'), ...
    sqrt( ...
        ( std(POS(:,time_idx),0,1,'omitnan') ./ sqrt(max(1,sum(~isnan(POS(:,time_idx)),1))) ).^2 + ...
        ( std(NEG(:,time_idx),0,1,'omitnan') ./ sqrt(max(1,sum(~isnan(NEG(:,time_idx)),1))) ).^2 ) ...
    ) ...

% --------------------- pull arrays ---------------------
BL_pos = Agg.ByGroup.negLag_posCorr.AllKept.baseline;
BL_neg = Agg.ByGroup.negLag_negCorr.AllKept.baseline;

ST_pos = Agg.ByGroup.negLag_posCorr.AllKept.stim;
ST_neg = Agg.ByGroup.negLag_negCorr.AllKept.stim;

WO_pos = Agg.ByGroup.negLag_posCorr.AllKept.wash;   % <-- change to .wo or .washout if needed
WO_neg = Agg.ByGroup.negLag_negCorr.AllKept.wash;   % <-- change to .wo or .washout if needed

if isempty(BL_pos) || isempty(BL_neg), error('Baseline arrays empty.'); end
if isempty(ST_pos) || isempty(ST_neg), error('Stim arrays empty.'); end
if isempty(WO_pos) || isempty(WO_neg), error('Washout arrays empty.'); end

% --------------------- compute NET traces ---------------------
[mBL, semBL] = net_from_groups(BL_pos, BL_neg);
[mST, semST] = net_from_groups(ST_pos, ST_neg);
[mWO, semWO] = net_from_groups(WO_pos, WO_neg);

% --------------------- plot NET only ---------------------
plot(time_vector, mBL, 'Color', Color_BL, 'LineWidth', 2.5);
fill([time_vector fliplr(time_vector)], ...
     [mBL+semBL fliplr(mBL-semBL)], ...
     Color_BL, 'FaceAlpha',0.10, 'EdgeColor','none');

plot(time_vector, mST, 'Color', Color_ST, 'LineWidth', 2.5);
fill([time_vector fliplr(time_vector)], ...
     [mST+semST fliplr(mST-semST)], ...
     Color_ST, 'FaceAlpha',0.12, 'EdgeColor','none');

plot(time_vector, mWO, 'Color', Color_WO, 'LineWidth', 2.5);
fill([time_vector fliplr(time_vector)], ...
     [mWO+semWO fliplr(mWO-semWO)], ...
     Color_WO, 'FaceAlpha',0.12, 'EdgeColor','none');

% --------------------- cosmetics ---------------------
xline(0,'--k');
yline(0,'--k');
xline(81.5);

xlabel('Time from threshold crossing (ms)');
ylabel('\Delta Z NET = mean(+Corr) - mean(-Corr)');
xlim([-250 250]);

title('NET activity overlay: Baseline (black), Stim (DodgerBlue), Washout (magenta)');

legend({'BL NET','BL SEM','ST NET','ST SEM','WO NET','WO SEM'}, ...
       'Location','best','Box','off');

box off;
%% ============================================================
% POPULATION VECTOR NET ACTIVITY (AD - BL)  [dimension-safe]
% ============================================================

do_bootstrap = true;
nBoot        = 2000;
alpha        = 0.05;

Color_pos = rgb('SteelBlue');
Color_neg = rgb('DarkGreen');
Color_net = [0 0 0];

% ------------------- pull data -------------------
BL_pos = Agg.ByGroup.negLag_posCorr.AllKept.baseline;
AD_pos = Agg.ByGroup.negLag_posCorr.AllKept.stim;

BL_neg = Agg.ByGroup.negLag_negCorr.AllKept.baseline;
AD_neg = Agg.ByGroup.negLag_negCorr.AllKept.stim;

if isempty(BL_pos) || isempty(AD_pos), error('posCorr AllKept baseline/stim empty'); end
if isempty(BL_neg) || isempty(AD_neg), error('negCorr AllKept baseline/stim empty'); end

% ------------------- per-unit deltas -------------------
dPos = AD_pos - BL_pos;   % [nPos x nT]
dNeg = AD_neg - BL_neg;   % [nNeg x nT]

nPos = size(dPos,1);
nNeg = size(dNeg,1);
nT   = size(dPos,2);

% ------------------- choose / build a time axis that MATCHES nT -------------------
% Priority:
%   1) time_vector if it matches
%   2) timevec (from your earlier PSTH structs) if it matches
%   3) fallback to linear axis over [-250 250] (or your xlim) with nT samples

if exist('time_vector','var') && ~isempty(time_vector) && numel(time_vector)==nT
    tv = time_vector(:)';   % force row
elseif exist('timevec','var') && ~isempty(timevec) && numel(timevec)==nT
    tv = timevec(:)';       % force row
else
    % fallback: assume your intended plot window is [-250 250]
    tv = linspace(-250, 250, nT);
end

% If still mismatched (e.g., dNeg has different nT), hard-trim to min length
nT2 = min([numel(tv), size(dPos,2), size(dNeg,2)]);
tv  = tv(1:nT2);
dPos = dPos(:,1:nT2);
dNeg = dNeg(:,1:nT2);

% ------------------- population vector construction -------------------
Vpos = +1 .* dPos;     % signed contributions
Vneg = -1 .* dNeg;

Vall = [Vpos; Vneg];   % [N x nT2]
N    = size(Vall,1);

net_sum  = sum(Vall, 1, 'omitnan');
net_mean = mean(Vall, 1, 'omitnan');
net_sem  = std(Vall, 0, 1, 'omitnan') ./ sqrt(N);

% ------------------- bootstrap CI for net_mean (optional) -------------------
net_ci_lo = nan(1, size(Vall,2));
net_ci_hi = nan(1, size(Vall,2));

if do_bootstrap
    boot_means = nan(nBoot, size(Vall,2));
    for b = 1:nBoot
        idx = randi(N, [N 1]);
        boot_means(b,:) = mean(Vall(idx,:), 1, 'omitnan');
    end
    loq = 100*(alpha/2);
    hiq = 100*(1 - alpha/2);
    net_ci_lo = prctile(boot_means, loq, 1);
    net_ci_hi = prctile(boot_means, hiq, 1);
end

% ------------------- group means for context (unsigned AD-BL) -------------------
mPos = mean(dPos,1,'omitnan');     % 1 x nT2
mNeg = mean(dNeg,1,'omitnan');     % 1 x nT2

% force row vectors (plot-safe)
tv       = tv(:)';
mPos     = mPos(:)';
mNeg     = mNeg(:)';
net_mean = net_mean(:)';
net_sem  = net_sem(:)';
net_ci_lo = net_ci_lo(:)';
net_ci_hi = net_ci_hi(:)';

% ------------------- plot -------------------
figure('Name','Population Vector NET (AD-BL)','Position',[100,100,820,360]); hold on;

plot(tv, mPos, 'Color', Color_pos, 'LineWidth', 1.75);
plot(tv, mNeg, 'Color', Color_neg, 'LineWidth', 1.75);

plot(tv, net_mean, 'Color', Color_net, 'LineWidth', 2.75);

% shading: bootstrap CI if usable; else SEM
useCI = do_bootstrap && any(isfinite(net_ci_lo)) && any(isfinite(net_ci_hi));
if useCI
    ok = isfinite(tv) & isfinite(net_ci_lo) & isfinite(net_ci_hi);
    fill([tv(ok) fliplr(tv(ok))], ...
         [net_ci_hi(ok) fliplr(net_ci_lo(ok))], ...
         Color_net, 'FaceAlpha',0.12, 'EdgeColor','none');
else
    fill([tv fliplr(tv)], ...
         [net_mean+net_sem fliplr(net_mean-net_sem)], ...
         Color_net, 'FaceAlpha',0.12, 'EdgeColor','none');
end

xline(0,'--k');
yline(0,'--k');

% keep your preferred xlim if it matches your alignment
xlim([-250 250]);

xlabel('Time (ms) relative to alignment event');
ylabel('\Delta Z (AD - BL)');
title(sprintf('Population Vector NET: mean(dir_i*(AD-BL)) | N=%d (pos=%d, neg=%d)', N, nPos, nNeg));

legend({'posCorr mean (AD-BL)','negCorr mean (AD-BL)','NET mean','NET CI/SEM'}, ...
       'Location','best','Box','off');
box off;

% ------------------- outputs -------------------
PopVec = struct();
PopVec.dPos      = dPos;
PopVec.dNeg      = dNeg;
PopVec.Vall      = Vall;
PopVec.net_sum   = net_sum;
PopVec.net_mean  = net_mean;
PopVec.net_sem   = net_sem;
PopVec.net_ci_lo = net_ci_lo;
PopVec.net_ci_hi = net_ci_hi;
PopVec.time_ms   = tv;

fprintf('[PopVec] nPos=%d | nNeg=%d | N=%d | nT=%d\n', nPos, nNeg, N, numel(tv));

%% ============================================================
% NET ACTIVITY: ALL KEPT (WO - BL) combining posCorr and negCorr as "vectors"
%
% Goal:
%   - Compute the WO-BL delta trace for posCorr units and negCorr units
%   - Then compute a NET trace that reflects "vector" combination:
%         NET = (posCorr delta) - (negCorr delta)
%
% Uses:
%   - Agg.ByGroup.(cat).AllKept.baseline  (BL)
%   - Agg.ByGroup.(cat).AllKept.wash      (WO)   <-- NOTE: uses .wash, not .stim
%
% Assumes:
%   - time_vector exists (same length as PSTHs); if not, use tv handling elsewhere
%   - delta_sem(ST,BL) returns [meanDelta, semDelta] for (ST - BL)
% ============================================================

Color_pos = rgb('SteelBlue');
Color_neg = rgb('DarkGreen');
Color_net = [0 0 0];

figure('Name','ALL-KEPT NET (WO-BL): pos/neg + net','Position',[100,100,750,330]);
hold on;

% --------------------- pull arrays ---------------------
BL_pos = Agg.ByGroup.negLag_posCorr.AllKept.baseline;
WO_pos = Agg.ByGroup.negLag_posCorr.AllKept.wash;     % <- WO

BL_neg = Agg.ByGroup.negLag_negCorr.AllKept.baseline;
WO_neg = Agg.ByGroup.negLag_negCorr.AllKept.wash;     % <- WO

% --------------------- sanity checks ---------------------
if isempty(BL_pos) || isempty(WO_pos)
    error('AllKept posCorr arrays are empty (baseline or wash).');
end
if isempty(BL_neg) || isempty(WO_neg)
    error('AllKept negCorr arrays are empty (baseline or wash).');
end

% --------------------- deltas per group ---------------------
[mPos, semPos] = delta_sem(WO_pos, BL_pos);  % WO - BL for posCorr
[mNeg, semNeg] = delta_sem(WO_neg, BL_neg);  % WO - BL for negCorr

% --------------------- NET vector combination ---------------------
mNet  = mPos - mNeg;
semNet = sqrt(semPos.^2 + semNeg.^2);

% --------------------- plot posCorr ---------------------
plot(time_vector, mPos, 'Color', Color_pos, 'LineWidth', 2);
fill([time_vector fliplr(time_vector)], ...
     [mPos+semPos fliplr(mPos-semPos)], ...
     Color_pos, 'FaceAlpha',0.20, 'EdgeColor','none');

% --------------------- plot negCorr ---------------------
plot(time_vector, mNeg, 'Color', Color_neg, 'LineWidth', 2);
fill([time_vector fliplr(time_vector)], ...
     [mNeg+semNeg fliplr(mNeg-semNeg)], ...
     Color_neg, 'FaceAlpha',0.20, 'EdgeColor','none');

% --------------------- plot NET ---------------------
plot(time_vector, mNet, 'Color', Color_net, 'LineWidth', 2.5);
fill([time_vector fliplr(time_vector)], ...
     [mNet+semNet fliplr(mNet-semNet)], ...
     Color_net, 'FaceAlpha',0.10, 'EdgeColor','none');

% --------------------- cosmetics ---------------------
xline(0,'--k');
yline(0,'--k');
yline(-2,':'); yline(2,':');

xline(81.5);

xlabel('Time from threshold crossing (ms)');
ylabel('\Delta Z (WO - BL)');
xlim([-250 250]);

title('ALL-KEPT \Delta (WO - BL): posCorr, negCorr, and NET = pos - neg');

legend({'posCorr mean','posCorr SEM','negCorr mean','negCorr SEM','NET mean','NET SEM'}, ...
       'Location','best','Box','off');

box off;
%% ===================== (2) MODULATED KEPT (sig_any) =====================
for ci = 1:numel(group_names)
    cat = group_names{ci};

    BL = Agg.ByGroup.(cat).ModKept.baseline;
    ST = Agg.ByGroup.(cat).ModKept.stim;
    WO = Agg.ByGroup.(cat).ModKept.wash;

    if isempty(BL) || isempty(ST) || isempty(WO)
        fprintf('ModKept %s empty -> skipping.\n', cat);
        continue;
    end

    % ---- deltas ----
    [mStim, semStim] = delta_sem(ST, BL);
    [mWO,   semWO  ] = delta_sem(WO, BL);

    figure('Name',['CSM-MOD-KEPT DELTAS ' cat],'Position',[100,100,650,300]); hold on;

    plot(time_vector, mStim, 'Color', Color_StimDelta, 'LineWidth', 2);
    fill([time_vector fliplr(time_vector)], ...
         [mStim+semStim fliplr(mStim-semStim)], ...
         Color_StimDelta, 'FaceAlpha',0.20, 'EdgeColor','none');

    plot(time_vector, mWO, 'Color', Color_WODelta, 'LineWidth', 2);
    fill([time_vector fliplr(time_vector)], ...
         [mWO+semWO fliplr(mWO-semWO)], ...
         Color_WODelta, 'FaceAlpha',0.20, 'EdgeColor','none');

    xline(0,'--k'); yline(0,'--k');
    yline(-2,':'); yline(2,':');

    %xlabel('Time from max velocity (ms)');   % edit label if you want (csmT2 etc.)
    xlabel('Time from threshold crossing (ms)');   % edit label if you want (csmT2 etc.)
    ylabel('\Delta Z (Condition - Baseline)');
    xlim([-250 250]);
    title(sprintf('CSM-MOD-KEPT Δ: %s (N=%d)', cat, size(BL,1)), 'Interpreter','none');
    legend({'Stim-BL','', 'WO-BL',''}, 'Location','best','Box','off');
    box off;
end
%% ============================================================
% CSM-MODULATED KEPT: Single plot with AD-BL deltas for BOTH categories
%   - negLag_posCorr  -> SteelBlue
%   - negLag_negCorr  -> DarkGreen
% ============================================================

Color_pos = rgb('SteelBlue');
Color_neg = rgb('DarkGreen');

figure('Name','CSM-MOD-KEPT AD-BL (both categories)','Position',[100,100,650,300]); 
hold on;

cats = {'negLag_posCorr','negLag_negCorr'};
colors = {Color_pos, Color_neg};

for ci = 1:numel(cats)

    cat = cats{ci};

    BL = Agg.ByGroup.(cat).ModKept.baseline;
    ST = Agg.ByGroup.(cat).ModKept.stim;

    if isempty(BL) || isempty(ST)
        fprintf('ModKept %s empty -> skipping.\n', cat);
        continue;
    end

    % ---- deltas (AD - BL) ----
    [mStim, semStim] = delta_sem(ST, BL);

    % plot mean
    plot(time_vector, mStim, 'Color', colors{ci}, 'LineWidth', 2);

    % SEM shading
    fill([time_vector fliplr(time_vector)], ...
         [mStim+semStim fliplr(mStim-semStim)], ...
         colors{ci}, 'FaceAlpha',0.20, 'EdgeColor','none');

end

% common cosmetics
xline(0,'--k'); 
yline(0,'--k');
yline(-2,':'); 
yline(2,':');
xline(81.5)

xlabel('Time from threshold crossing (ms)');
ylabel('\Delta Z (AD - BL)');

xlim([-250 250]);

title('CSM-MOD-KEPT Δ (AD - BL) : posCorr vs negCorr');

legend({'posCorr AD-BL','','negCorr AD-BL',''}, ...
       'Location','best','Box','off');

box off;

%% ============================================================
% CSM-MODULATED KEPT: NET ACTIVITY (AD - BL)
%   - posCorr delta (SteelBlue)
%   - negCorr delta (DarkGreen)
%   - NET = pos - neg (Black)
% ============================================================

Color_pos = rgb('SteelBlue');
Color_neg = rgb('DarkGreen');
Color_net = [0 0 0];

figure('Name','CSM-MOD-KEPT NET (AD-BL)','Position',[100,100,700,320]);
hold on;

% --------------------- pull arrays ---------------------
BL_pos = Agg.ByGroup.negLag_posCorr.ModKept.baseline;
AD_pos = Agg.ByGroup.negLag_posCorr.ModKept.stim;

BL_neg = Agg.ByGroup.negLag_negCorr.ModKept.baseline;
AD_neg = Agg.ByGroup.negLag_negCorr.ModKept.stim;

% --------------------- sanity checks ---------------------
if isempty(BL_pos) || isempty(AD_pos)
    error('ModKept posCorr arrays are empty.');
end
if isempty(BL_neg) || isempty(AD_neg)
    error('ModKept negCorr arrays are empty.');
end

% --------------------- deltas ---------------------
[mPos, semPos] = delta_sem(AD_pos, BL_pos);   % AD - BL (pos)
[mNeg, semNeg] = delta_sem(AD_neg, BL_neg);   % AD - BL (neg)

% --------------------- NET (vector difference) ---------------------
mNet   = mPos - mNeg;
semNet = sqrt(semPos.^2 + semNeg.^2);

% --------------------- plot posCorr ---------------------
plot(time_vector, mPos, 'Color', Color_pos, 'LineWidth', 2);
fill([time_vector fliplr(time_vector)], ...
     [mPos+semPos fliplr(mPos-semPos)], ...
     Color_pos, 'FaceAlpha',0.20, 'EdgeColor','none');

% --------------------- plot negCorr ---------------------
plot(time_vector, mNeg, 'Color', Color_neg, 'LineWidth', 2);
fill([time_vector fliplr(time_vector)], ...
     [mNeg+semNeg fliplr(mNeg-semNeg)], ...
     Color_neg, 'FaceAlpha',0.20, 'EdgeColor','none');

% --------------------- plot NET ---------------------
plot(time_vector, mNet, 'Color', Color_net, 'LineWidth', 2.5);
fill([time_vector fliplr(time_vector)], ...
     [mNet+semNet fliplr(mNet-semNet)], ...
     Color_net, 'FaceAlpha',0.10, 'EdgeColor','none');

% --------------------- cosmetics ---------------------
xline(0,'--k');
yline(0,'--k');
yline(-2,':'); 
yline(2,':');
xline(81.5);

xlabel('Time from threshold crossing (ms)');
ylabel('\Delta Z (AD - BL)');
xlim([-250 250]);

title('CSM-MOD-KEPT \Delta (AD - BL): posCorr, negCorr, and NET = pos - neg');

legend({'posCorr mean','posCorr SEM','negCorr mean','negCorr SEM','NET mean','NET SEM'}, ...
       'Location','best','Box','off');

box off;

%% ============================================================
% CSM-MODULATED KEPT: NET ACTIVITY (WO - BL)
%   - posCorr delta (SteelBlue)
%   - negCorr delta (DarkGreen)
%   - NET = pos - neg (Black)
% ============================================================

Color_pos = rgb('SteelBlue');
Color_neg = rgb('DarkGreen');
Color_net = [0 0 0];

figure('Name','CSM-MOD-KEPT NET (WO-BL)','Position',[100,100,700,320]);
hold on;

% --------------------- pull arrays ---------------------
BL_pos = Agg.ByGroup.negLag_posCorr.ModKept.baseline;
WO_pos = Agg.ByGroup.negLag_posCorr.ModKept.wash;   % <- WO

BL_neg = Agg.ByGroup.negLag_negCorr.ModKept.baseline;
WO_neg = Agg.ByGroup.negLag_negCorr.ModKept.wash;   % <- WO

% --------------------- sanity checks ---------------------
if isempty(BL_pos) || isempty(WO_pos)
    error('ModKept posCorr arrays are empty (baseline or wash).');
end
if isempty(BL_neg) || isempty(WO_neg)
    error('ModKept negCorr arrays are empty (baseline or wash).');
end

% --------------------- deltas ---------------------
[mPos, semPos] = delta_sem(WO_pos, BL_pos);   % WO - BL (pos)
[mNeg, semNeg] = delta_sem(WO_neg, BL_neg);   % WO - BL (neg)

% --------------------- NET (vector difference) ---------------------
mNet   = mPos - mNeg;
semNet = sqrt(semPos.^2 + semNeg.^2);

% --------------------- plot posCorr ---------------------
plot(time_vector, mPos, 'Color', Color_pos, 'LineWidth', 2);
fill([time_vector fliplr(time_vector)], ...
     [mPos+semPos fliplr(mPos-semPos)], ...
     Color_pos, 'FaceAlpha',0.20, 'EdgeColor','none');

% --------------------- plot negCorr ---------------------
plot(time_vector, mNeg, 'Color', Color_neg, 'LineWidth', 2);
fill([time_vector fliplr(time_vector)], ...
     [mNeg+semNeg fliplr(mNeg-semNeg)], ...
     Color_neg, 'FaceAlpha',0.20, 'EdgeColor','none');

% --------------------- plot NET ---------------------
plot(time_vector, mNet, 'Color', Color_net, 'LineWidth', 2.5);
fill([time_vector fliplr(time_vector)], ...
     [mNet+semNet fliplr(mNet-semNet)], ...
     Color_net, 'FaceAlpha',0.10, 'EdgeColor','none');

% --------------------- cosmetics ---------------------
xline(0,'--k');
yline(0,'--k');
yline(-2,':');
yline(2,':');
xline(81.5);

xlabel('Time from threshold crossing (ms)');
ylabel('\Delta Z (WO - BL)');
xlim([-250 250]);

title('CSM-MOD-KEPT \Delta (WO - BL): posCorr, negCorr, and NET = pos - neg');

legend({'posCorr mean','posCorr SEM','negCorr mean','negCorr SEM','NET mean','NET SEM'}, ...
       'Location','best','Box','off');

box off;
%% ============================================================
% FIG: NEG-LAG ONLY — Epoch DELTAZ scatter by correlation sign (shuffle-based)
%   - NOT significant: gray
%   - sig +corr (negLag_posCorr): SteelBlue
%   - sig -corr (negLag_negCorr): DarkGreen
%
% Uses (per unit):
%   U.keep == true
%   U.category in {'negLag_posCorr','negLag_negCorr'}
%   U.CSM_zmod.skip == false
%   U.CSM_zmod.sig_any
%   U.CSM_zmod.deltaZ_mean   (or swap to U.CSM_zmod.dom_val if you prefer)
%
% Drop AFTER main script (XCorrData exists).


% -------- settings --------
groups = {'negLag_posCorr','negLag_negCorr'};

if exist('rgb','file')==2 || exist('rgb','builtin')==5
    c_gray = rgb('Gray');
    c_pos  = rgb('SteelBlue');
    c_neg  = rgb('DarkGreen');
else
    c_gray = [0.6 0.6 0.6];
    c_pos  = [70 130 180]/255;  % steelblue-ish
    c_neg  = [0 100 0]/255;     % darkgreen-ish
end

jitter_amt = 0.25;

% -------- collect per-unit deltaZ-shuffle results --------
Unit = struct('cat',{},'day',{},'unitIdx',{});
dEp  = [];   % epoch effect size (deltaZ_mean OR dom_val)
sigE = [];   % sig_any

for day = 1:numel(XCorrData)
    if ~isfield(XCorrData(day),'units') || isempty(XCorrData(day).units), continue; end

    for uu = 1:numel(XCorrData(day).units)
        U = XCorrData(day).units(uu);

        if ~isfield(U,'keep') || ~U.keep, continue; end
        if ~isfield(U,'category') || ~ismember(U.category, groups), continue; end

        if ~isfield(U,'CSM_zmod') || ~isfield(U.CSM_zmod,'skip') || U.CSM_zmod.skip, continue; end
        if ~isfield(U.CSM_zmod,'sig_any'), continue; end

        % choose what to plot on Y:
        %   (A) mean deltaZ across bins in epoch
        if ~isfield(U.CSM_zmod,'deltaZ_mean') || ~isfinite(U.CSM_zmod.deltaZ_mean), continue; end
        yval = U.CSM_zmod.deltaZ_mean;

        %   (B) OR: dominant signed peak (uncomment to use)
        % if ~isfield(U.CSM_zmod,'dom_val') || ~isfinite(U.CSM_zmod.dom_val), continue; end
        % yval = U.CSM_zmod.dom_val;

        Unit(end+1).cat   = U.category;
        Unit(end).day     = day;
        Unit(end).unitIdx = U.unitIdx;

        dEp(end+1,1)  = yval;
        sigE(end+1,1) = logical(U.CSM_zmod.sig_any);
    end
end

if isempty(dEp)
    warning('No units found matching: kept + negLag_(pos/neg)Corr + valid CSM_zmod.');
else
    fprintf('Plotting %d kept units with CSM_zmod results (neg-lag only).\n', numel(dEp));
end

% -------- map categories to x positions --------
x = zeros(numel(Unit),1);
for i = 1:numel(Unit)
    x(i) = find(strcmp(Unit(i).cat, groups));
end
jitter = (rand(numel(x),1)-0.5) * jitter_amt;
xj = x + jitter;

idxSig = logical(sigE);
idxNo  = ~idxSig;

idxSigPos = idxSig & (x==1);  % negLag_posCorr
idxSigNeg = idxSig & (x==2);  % negLag_negCorr

%% -------- plot --------
figure('Color','w','Position',[120 120 330 420]); hold on;

sz_sig = 55;
sz_nos = 45;

% non-sig (gray)
scatter(xj(idxNo), dEp(idxNo), sz_nos, 'filled', ...
    'MarkerFaceColor', c_gray, 'MarkerFaceAlpha', 0.55, 'MarkerEdgeColor','none');

% sig posCorr (steel blue)
scatter(xj(idxSigPos), dEp(idxSigPos), sz_sig, 'filled', ...
    'MarkerFaceColor', c_pos, 'MarkerFaceAlpha', 0.95, 'MarkerEdgeColor','k');

% sig negCorr (dark green)
scatter(xj(idxSigNeg), dEp(idxSigNeg), sz_sig, 'filled', ...
    'MarkerFaceColor', c_neg, 'MarkerFaceAlpha', 0.95, 'MarkerEdgeColor','k');

xlim([0.5 numel(groups)+0.5]);
xticks(1:numel(groups));
xticklabels({'negLag posCorr','negLag negCorr'});
xtickangle(20);

yline(0,'k-','LineWidth',1);
ylabel('\DeltaZ (Stim - BL) during CSM epoch');
title('Epoch effect size (deltaZ shuffle); significant neurons colored by corr sign');

% optional y-limits (uncomment/edit)
% ylim([-3 3])

% group means (SIGNIFICANT ONLY) + counts
for g = 1:numel(groups)
    m_all = (x==g);
    m_sig = m_all & idxSig;

    if any(m_sig)
        mu = mean(dEp(m_sig),'omitnan');
        plot(g, mu, 'kd','MarkerFaceColor','k','MarkerSize',8);

        text(g, mu, sprintf('  %d sig / %d total', sum(m_sig), sum(m_all)), ...
            'FontWeight','bold','VerticalAlignment','bottom');
    else
        if any(m_all)
            text(g, max(dEp(m_all),[],'omitnan'), sprintf('  0 sig / %d total', sum(m_all)), ...
                'FontWeight','bold','VerticalAlignment','bottom');
        end
    end
end

legend({'not sig','sig +corr','sig -corr','mean (sig only)'}, 'Location','best', 'Box','off');
set(gca,'LineWidth',1,'FontSize',12,'TickDir','out');
box off;


%% ============================================================
% FIG: p-value QC by NEG-LAG category (strip plot, log y)
%   - Categories: negLag_posCorr, negLag_negCorr
%   - Plots per-unit shuffle p-values: U.CSM_zmod.p_unit
%   - Only KEPT units (peak-gated) and with valid CSM_zmod (skip==false)
%
% Drop AFTER main script (XCorrData exists).
% ============================================================

cats = {'negLag_posCorr','negLag_negCorr'};
nCats = numel(cats);

% ---- collect p_unit by category ----
p_by_cat = cell(1,nCats);
n_by_cat = zeros(1,nCats);

for day = 1:numel(XCorrData)
    if ~isfield(XCorrData(day),'units') || isempty(XCorrData(day).units), continue; end

    for uu = 1:numel(XCorrData(day).units)
        U = XCorrData(day).units(uu);

        % kept only
        if ~isfield(U,'keep') || ~U.keep, continue; end

        % category filter
        if ~isfield(U,'category') || ~ismember(U.category, cats), continue; end

        % must have CSM_zmod p-values
        if ~isfield(U,'CSM_zmod') || ~isfield(U.CSM_zmod,'skip') || U.CSM_zmod.skip, continue; end
        if ~isfield(U.CSM_zmod,'p_unit') || ~isfinite(U.CSM_zmod.p_unit), continue; end

        ci = find(strcmp(U.category, cats));
        p_by_cat{ci}(end+1,1) = U.CSM_zmod.p_unit;
        n_by_cat(ci) = n_by_cat(ci) + 1;
    end
end

%% ---- plot ----
figure('Position',[120 120 330 420],'Color','w'); hold on;

jitter_amt = 0.18;

for i = 1:nCats
    p = p_by_cat{i};
    p = p(isfinite(p));
    if isempty(p), continue; end

    x = i + jitter_amt*(rand(size(p)) - 0.5);
    plot(x, p, '.', 'MarkerSize',20);
end

% alpha line (only meaningful if do_unit_shuffle=true)
if exist('alpha_unit','var') && ~isempty(alpha_unit)
    yline(alpha_unit, '--', sprintf('alpha=%.3g', alpha_unit), 'LineWidth',1.5);
else
    yline(0.05, '--', 'alpha=0.05', 'LineWidth',1.5);
end

xlim([0.5 nCats+0.5]);
xticks(1:nCats);
xticklabels({'negLag posCorr','negLag negCorr'});
xtickangle(20);

ylabel('p_{unit} (shuffle)');
title(sprintf('Unit-level shuffle p-values by NEG-LAG category (KEPT only)  |  N=[%d %d]', ...
    n_by_cat(1), n_by_cat(2)));

set(gca,'YScale','log');
ylim([1e-4 1]);   % tweak if you want
grid on; box on;
set(gca,'TickDir','out','LineWidth',1,'FontSize',12);

%% ============================================================
% PLOT A (BLOCK STIM): Peak lag vs peak corr (NEG-LAG ONLY, KEPT units)
% Publication-style, SteelBlue (posCorr) + DarkGreen (negCorr)
% NO labels, NO grid
% AXES (actual plot) is ~half-height of the figure
%
% Requires XCorrData from your BlockStim script:
%   U.keep
%   U.negLag_peakLag
%   U.negLag_peakVal
%   U.category in {'negLag_posCorr','negLag_negCorr'}
%% ============================================================

allLag  = [];
allCorr = [];
allCat  = {};

cats_use = {'negLag_posCorr','negLag_negCorr'};

for day = 1:numel(XCorrData)
    if ~isfield(XCorrData(day),'units') || isempty(XCorrData(day).units), continue; end

    for uu = 1:numel(XCorrData(day).units)
        u = XCorrData(day).units(uu);

        if ~isfield(u,'keep') || ~u.keep, continue; end
        if ~isfield(u,'negLag_peakLag') || ~isfield(u,'negLag_peakVal') || ~isfield(u,'category'), continue; end
        if ~ismember(u.category, cats_use), continue; end

        lag  = u.negLag_peakLag;
        corr = u.negLag_peakVal;

        if ~isfinite(lag) || ~isfinite(corr), continue; end
        if lag >= 0, continue; end  % negative lags only

        allLag(end+1,1)  = lag;
        allCorr(end+1,1) = corr;
        allCat{end+1,1}  = u.category;
    end
end

% Colors (requires rgb.m in path; else fall back)
if exist('rgb','file')==2 || exist('rgb','builtin')==5
    Color_pos = rgb('SteelBlue');
    Color_neg = rgb('DarkGreen');
else
    Color_pos = [70 130 180]/255;
    Color_neg = [0 100 0]/255;
end

% ---- Figure
fig = figure('Color','w','Units','pixels','Position',[350 150 650 650]);

% ---- Axes: make the ACTUAL plot ~half the figure height
ax = axes(fig);
set(ax, 'Position', [0.15 0.33 0.75 0.45]);  % [left bottom width height] (normalized)
hold(ax,'on');

ms  = 75;

idxPos = strcmp(allCat,'negLag_posCorr');
idxNeg = strcmp(allCat,'negLag_negCorr');

if any(idxPos)
    h1 = scatter(ax, allLag(idxPos), allCorr(idxPos), ms, ...
        'MarkerFaceColor', Color_pos, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceAlpha', 0.85, 'MarkerEdgeAlpha', 0.95, 'LineWidth', 0.75);
else
    h1 = scatter(ax, NaN, NaN, ms, 'filled', 'MarkerFaceColor', Color_pos, 'MarkerEdgeColor','k');
end

if any(idxNeg)
    h2 = scatter(ax, allLag(idxNeg), allCorr(idxNeg), ms, ...
        'MarkerFaceColor', Color_neg, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceAlpha', 0.85, 'MarkerEdgeAlpha', 0.95, 'LineWidth', 0.75);
else
    h2 = scatter(ax, NaN, NaN, ms, 'filled', 'MarkerFaceColor', Color_neg, 'MarkerEdgeColor','k');
end

% Reference lines
xline(ax, 0, '--', 'Color', [0 0 0], 'LineWidth', 1.0, 'Alpha', 0.7);
yline(ax, 0, '--', 'Color', [0 0 0], 'LineWidth', 1.0, 'Alpha', 0.7);

% Limits
ylim(ax, [-1 1]);
xlim(ax, [-maxLag_ms 0]);  % uses your script param; safe default for block stim

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

axis(ax,'normal'); % do NOT force square

% Legend only
lgd = legend(ax, [h1 h2], {'negLag\_posCorr','negLag\_negCorr'}, ...
    'Location','best', 'Box','off');
lgd.FontSize = 11;

% Optional: add N text below plot
text(ax, 0.02, -0.25, sprintf('N pos=%d | N neg=%d', sum(idxPos), sum(idxNeg)), ...
    'Units','normalized','FontSize',11,'Color',[0 0 0]);
