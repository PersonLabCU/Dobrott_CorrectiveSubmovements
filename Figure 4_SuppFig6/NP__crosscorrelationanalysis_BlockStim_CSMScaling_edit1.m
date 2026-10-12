
% ============================================================
% integrated_BlockStim_NEG_LAG_xcorr__CSMmod_DELTAZshuffle__PSTHavg__SCALING.m
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
%   (4) NEW: SCALING ANALYSIS (ONLY ADDITIVE; DOES NOT CHANGE YOUR EXISTING LOGIC)
%         - For CSM-modulated units, test whether trial-by-trial neural response
%           scales with trial-by-trial CSM magnitude.
%         - Neural window can be shifted by the unit’s neg-lag peakLag.
%         - Outputs: long-format table + per-unit correlations + mixed-effects models
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

% Threshold for AlignThreshCSMHF_1ms_edit
posThresh = 0.1;     % <-- set to whatever your pos_thresh alignment uses
alignDir  = 1;        % 1 = outreach


% ----- PSTH peak-gate + zscore baseline window -----
valid_peak_range_ms  = [-250 250];     % KEEP if peak time (peakgate_field) is within this
baseline_win_def     = [-1000 -500];   % baseline window for PSTH z-scoring (ms)
plot_win             = [-500 500];     % for plotting later if you want

% ----- Choose fields (BASELINE vs STIMULATED) -----
peakgate_field       = 'BLCSM_max_vx';          % for keep gate + xcorr raster source

% PSTHs you want to average (baseline vs stimulated)
baseline_psth_field  = 'BLCSM_pos_thresh';      % BASELINE PSTH field
stim_psth_field      = 'ADCSM_pos_thresh';      % STIMULATED PSTH field
wash_psth_field      = 'WOECSM_pos_thresh';     % WASH PSTH field

% ----- CSM modulation (deltaZ sustained + shuffle) -----
%csm_epoch_win          = [-70 200];   % ms
csm_epoch_win          = [-50 100];   % ms used in the other analyses 
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

 %% ---------------------- NEW: SCALING SETTINGS -----------------------
% % (Additive; does not change your existing xcorr + modulation logic)
% do_scaling = true;
% 
% % Use only CSM-modulated units (sig_any) for scaling
% scaling_use_only_sig_any = true;
% 
% % Behavior epoch used to compute CSM magnitude (relative to alignment event)
% beh_epoch_win = csm_epoch_win;  % e.g. [-50 100] ms
% 
% % Neural epoch for per-trial response:
% % If true, neural window = beh_epoch_win + negLag_peakLag (peakLag is <0)
% scale_shift_neural_by_peakLag = true;
% 
% % Optional padding around neural window
% pad_ms = 0;  % e.g., 10 if you want slack
% 
% % Choose CSM magnitude metric from aligned kinematics in beh_epoch_win:
% % 'peakVel' | 'peakAbsVel' | 'aucVel' | 'peakAccel'
% csmMag_metric = 'peakVel';
% 
% % Dimension index in stimvel/stimaccel (you used 5 above for outward)
% kin_dim = 5;
% 
% % Minimum trials to compute within-unit correlation
% min_trials_for_unit_corr = 8;
% 
% % Kinematics alignment window for scaling extraction (ms)
% % Keep reasonably wide so epoch windows always exist.
% scale_align_pre_ms  = 500;
% scale_align_post_ms = 500;

%% ---------------------- NEW: SCALING SETTINGS -----------------------
do_scaling = true;
scaling_use_only_sig_any = true;

% --- FIXED kinematic window for CSM magnitude ---
%beh_epoch_win = [0 150];     % <-- fixed window for kinematics (ms) % good neg corr units 
%beh_epoch_win = [0 100]; 
%beh_epoch_win = [0 75];     % <-- fixed window for kinematics (ms) % good neg corr units 
%beh_epoch_win = [25 100];
beh_epoch_win = [-25 100];

% --- NEURAL window uses peakLag shift relative to beh_epoch_win ---
scale_shift_neural_by_peakLag = true;
pad_ms = 0;                    % optional padding around neural window

csmMag_metric = 'peakVel';
kin_dim = 5;

min_trials_for_unit_corr = 8;

scale_align_pre_ms  = 500;
scale_align_post_ms = 500;

%
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

% NEW: Scaling long-format table (trial-level)
Scale = struct();
Scale.rows = table();  % appended as we go

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

        % ---------------- PSTH Z-SCORES (baseline vs stimulated) ----------------
        % Pull PSTHs (must exist)
        if ~isfield(data.units(unitIdx), baseline_psth_field) || ~isfield(data.units(unitIdx).(baseline_psth_field),'PSTH_smoothbin') || ...
           ~isfield(data.units(unitIdx), stim_psth_field)     || ~isfield(data.units(unitIdx).(stim_psth_field),    'PSTH_smoothbin') || ...
           ~isfield(data.units(unitIdx), wash_psth_field)     || ~isfield(data.units(unitIdx).(wash_psth_field),    'PSTH_smoothbin')
            XCorrData(day).units(uu).PSTH = struct('skip',true,'reason','missing_psth_fields');
        else
            psth_base = data.units(unitIdx).(baseline_psth_field).PSTH_smoothbin;
            psth_stim = data.units(unitIdx).(stim_psth_field).PSTH_smoothbin;
            psth_wash = data.units(unitIdx).(wash_psth_field).PSTH_smoothbin;

            mu0 = mean(psth_base(baseline_win),'omitnan');
            sd0 = std(psth_base(baseline_win),0,'omitnan');
            if ~isfinite(sd0) || sd0 == 0, sd0 = 1; end

            base_z = (psth_base - mu0) ./ sd0;
            stim_z = (psth_stim - mu0) ./ sd0;
            wash_z = (psth_wash - mu0) ./ sd0;

            XCorrData(day).units(uu).PSTH.skip   = false;
            XCorrData(day).units(uu).PSTH.base_z = base_z(:)';
            XCorrData(day).units(uu).PSTH.stim_z = stim_z(:)';
            XCorrData(day).units(uu).PSTH.wash_z = wash_z(:)';

            % global aggregation: ALL kept
            Agg.AllKept.baseline = vertcat(Agg.AllKept.baseline, base_z(:)');
            Agg.AllKept.stim     = vertcat(Agg.AllKept.stim,     stim_z(:)');
            Agg.AllKept.wash     = vertcat(Agg.AllKept.wash,     wash_z(:)');

            % optional: by negLag category (if assigned)
            if isfield(XCorrData(day).units(uu),'category')
                gn = XCorrData(day).units(uu).category;
                if ismember(gn, group_names)
                    Agg.ByGroup.(gn).AllKept.baseline = vertcat(Agg.ByGroup.(gn).AllKept.baseline, base_z(:)');
                    Agg.ByGroup.(gn).AllKept.stim     = vertcat(Agg.ByGroup.(gn).AllKept.stim,     stim_z(:)');
                    Agg.ByGroup.(gn).AllKept.wash     = vertcat(Agg.ByGroup.(gn).AllKept.wash,     wash_z(:)');
                end
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

                        % ===================== NEW: SCALING EXTRACTION =====================
                        % Only additive; uses unit’s peakLag (optional) and trial-by-trial rasters/kinematics.
                        % NEW BEHAVIOR:
                        %   - KINEMATICS uses FIXED beh_epoch_win (NOT shifted)
                        %   - NEURAL window is shifted by pkLag (optionally) but uses SAME beh_epoch_win edges
                        %   - NEURAL response is PEAK FR (Hz) within the neural window
                        % USING PEAK LAG AS THE WINDOW FOR THE NEURAL DATA
                        if do_scaling
                            if ~scaling_use_only_sig_any || sig_any

                                % Pull unit peakLag and category
                                pkLag = NaN;
                                if isfield(XCorrData(day).units(uu),'negLag_peakLag')
                                    pkLag = XCorrData(day).units(uu).negLag_peakLag;
                                end

                                if scale_shift_neural_by_peakLag && ~isfinite(pkLag)
                                    % can't shift if no peakLag
                                else
                                    cat = "negLag_unknown";
                                    if isfield(XCorrData(day).units(uu),'category')
                                        cat = string(XCorrData(day).units(uu).category);
                                    end

                                    % Use Stim condition (MATCHES stim_mod_field)
                                    % If you switch stim_mod_field to WOE, also switch this:
                                    Reaches_scale = ExcitRR_AD.ReReaches;
                                    Raster_scale  = Stim_raster;  % from stim_mod_field above

                                    % Kinematics aligned to SAME EVENT as raster_scale
                                    [vel_1ms, tK_ms, accel_1ms] = getCSMAlignedKinematics_1ms( ...
                                        Reaches_scale, scale_align_pre_ms, scale_align_post_ms, kinVar, kin_dim, posThresh, alignDir);

                                    % Trial matching
                                    nTrN = size(Raster_scale,1);
                                    nTrK = size(vel_1ms,2);
                                    nTr  = min(nTrN, nTrK);

                                    if nTr >= min_trials_per_cond
                                        Raster_scale = Raster_scale(1:nTr,:);
                                        vel_1ms      = vel_1ms(:,1:nTr);
                                        if ~isempty(accel_1ms), accel_1ms = accel_1ms(:,1:nTr); end

                                        % ============================================================
                                        % OLD CODE (kept for reference)
                                        % Behavior window idx (kinematics)
                                        % beh_idx = (tK_ms >= beh_epoch_win(1)) & (tK_ms <= beh_epoch_win(2));
                                        %
                                        % % Neural time base (same assumption used in modulation)
                                        % tN_ms = (-1000:1000)';
                                        % nTn   = min(size(Raster_scale,2), numel(tN_ms));
                                        % tN_ms = tN_ms(1:nTn);
                                        % Raster_scale = Raster_scale(:,1:nTn);
                                        %
                                        % % Neural window: shift by peakLag if requested
                                        % if scale_shift_neural_by_peakLag
                                        %     neu_win = [beh_epoch_win(1)+pkLag-pad_ms, beh_epoch_win(2)+pkLag+pad_ms];
                                        % else
                                        %     neu_win = [beh_epoch_win(1)-pad_ms, beh_epoch_win(2)+pad_ms];
                                        % end
                                        % neu_idx = (tN_ms >= neu_win(1)) & (tN_ms <= neu_win(2));
                                        % ============================================================

                                        % ============================================================
                                        % NEW CODE INSERTED HERE:
                                        %   FIXED kinematic window (beh_epoch_win) + peakLag-shifted neural window
                                        % ============================================================

                                        % Behavior window idx (FIXED; never shifted)
                                        beh_idx = (tK_ms >= beh_epoch_win(1)) & (tK_ms <= beh_epoch_win(2));

                                        if any(beh_idx)

                                            % Neural time base (same assumption used in modulation)
                                            tN_ms = (-1000:1000)';
                                            nTn   = min(size(Raster_scale,2), numel(tN_ms));
                                            tN_ms = tN_ms(1:nTn);
                                            Raster_scale = Raster_scale(:,1:nTn);

                                            % Neural window (SHIFTED by peakLag; kinematics stays FIXED)
                                            if scale_shift_neural_by_peakLag && isfinite(pkLag)
                                                neu_win = [beh_epoch_win(1)+pkLag-pad_ms, beh_epoch_win(2)+pkLag+pad_ms];
                                            else
                                                neu_win = [beh_epoch_win(1)-pad_ms, beh_epoch_win(2)+pad_ms];
                                            end
                                            neu_idx = (tN_ms >= neu_win(1)) & (tN_ms <= neu_win(2));

                                            if any(neu_idx)

                                                % ============================================================
                                                % OLD NEURAL RESPONSE (MEAN FR) (kept for reference)
                                                % win_s = sum(neu_idx)/1000;
                                                % neuralFR = sum(Raster_scale(:,neu_idx),2) ./ win_s;  % [nTr x 1]
                                                % ============================================================

                                                % ============================================================
                                                % NEW NEURAL RESPONSE: PEAK FIRING RATE (Hz) WITHIN neu_idx
                                                % (Gaussian smoothing on 1ms raster)
                                                % ============================================================
                                                if ~exist('gauss_sigma','var') || isempty(gauss_sigma) || ~isfinite(gauss_sigma)
                                                    gauss_sigma = 10;
                                                end
                                                kern_len = max(1, round(6*gauss_sigma));              % points (ms)
                                                gk = fspecial('gaussian',[1 kern_len], gauss_sigma);
                                                gk = gk / sum(gk);

                                                neuralFR = nan(nTr,1);
                                                for tr = 1:nTr
                                                    sp    = double(Raster_scale(tr,:));   % spikes per ms (0/1)
                                                    fr_ms = conv(sp, gk, 'same');         % spikes/ms
                                                    fr_hz = fr_ms * 1000;                 % Hz
                                                    neuralFR(tr) = max(fr_hz(neu_idx), [], 'omitnan');  % PEAK in window
                                                end

                                                % ===================== CSM MAGNITUDE PER TRIAL (FIXED beh_idx) =====================
                                                switch lower(csmMag_metric)
                                                    case 'peakvel'
                                                        csmMag = max(vel_1ms(beh_idx,:),[],1)';          % [nTr x 1]
                                                    case 'peakabsvel'
                                                        csmMag = max(abs(vel_1ms(beh_idx,:)),[],1)';     % [nTr x 1]
                                                    case 'aucvel'
                                                        dt = median(diff(tK_ms))/1000;
                                                        csmMag = trapz(vel_1ms(beh_idx,:))*dt;
                                                        csmMag = csmMag(:);
                                                    case 'peakaccel'
                                                        if isempty(accel_1ms)
                                                            csmMag = NaN(nTr,1);
                                                        else
                                                            csmMag = max(accel_1ms(beh_idx,:),[],1)';    % [nTr x 1]
                                                        end
                                                    otherwise
                                                        error('Unknown csmMag_metric: %s', csmMag_metric);
                                                end

                                                % Build chunk table
                                                trIdx = (1:nTr)';
                                                Tchunk = table();
                                                Tchunk.day      = repmat(day, nTr, 1);
                                                Tchunk.unitIdx  = repmat(unitIdx, nTr, 1);
                                                Tchunk.trial    = trIdx;
                                                Tchunk.category = repmat(cat, nTr, 1);
                                                Tchunk.dom_dir  = repmat(dom_dir, nTr, 1);
                                                Tchunk.peakLag  = repmat(pkLag, nTr, 1);

                                                % store BOTH windows (helps QC)
                                                Tchunk.behWin1  = repmat(beh_epoch_win(1), nTr, 1);
                                                Tchunk.behWin2  = repmat(beh_epoch_win(2), nTr, 1);
                                                Tchunk.neuWin1  = repmat(neu_win(1), nTr, 1);
                                                Tchunk.neuWin2  = repmat(neu_win(2), nTr, 1);

                                                Tchunk.neuralFR = neuralFR;   % peak FR (Hz)
                                                Tchunk.csmMag   = csmMag;

                                                % Remove NaNs
                                                bad = ~isfinite(Tchunk.neuralFR) | ~isfinite(Tchunk.csmMag);
                                                Tchunk(bad,:) = [];

                                                if ~isempty(Tchunk)
                                                    Scale.rows = [Scale.rows; Tchunk];
                                                end
                                            end
                                        end
                                        % ============================================================
                                        % END NEW INSERT
                                        % ============================================================
                                    end
                                end
                            end
                        end
                        % =================== END NEW: SCALING EXTRACTION ===================

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
%%
% ---------- Count SIG_ANY units by category from XCorrData ----------
nSig_pos = 0; nSig_neg = 0;
allSig_units_pos = []; allSig_units_neg = [];

for day = 1:numel(XCorrData)
    if ~isfield(XCorrData(day),'units'), continue; end
    for uu = 1:numel(XCorrData(day).units)
        u = XCorrData(day).units(uu);
        if ~isfield(u,'keep') || ~u.keep, continue; end
        if ~isfield(u,'CSM_zmod') || ~isfield(u.CSM_zmod,'skip') || u.CSM_zmod.skip, continue; end
        if ~isfield(u.CSM_zmod,'sig_any') || ~u.CSM_zmod.sig_any, continue; end
        if ~isfield(u,'category'), continue; end

        if strcmp(u.category,'negLag_posCorr')
            nSig_pos = nSig_pos + 1;
            allSig_units_pos(end+1) = u.unitIdx;
        elseif strcmp(u.category,'negLag_negCorr')
            nSig_neg = nSig_neg + 1;
            allSig_units_neg(end+1) = u.unitIdx; 
        end
    end
end

fprintf('\nSIG_ANY counts from modulation:\n  posCorr=%d | negCorr=%d\n', nSig_pos, nSig_neg);

% ---------- Count UNIQUE units that actually made it into Scale.rows ----------
if exist('Scale','var') && isfield(Scale,'rows') && ~isempty(Scale.rows)
    scaleUnits = unique(Scale.rows.unitIdx);
    fprintf('Unique units in Scale.rows = %d\n', numel(scaleUnits));

    % How many sig_any units are missing from Scale.rows?
    miss_pos = setdiff(unique(allSig_units_pos(:)), scaleUnits);
    miss_neg = setdiff(unique(allSig_units_neg(:)), scaleUnits);

    fprintf('Sig_any units missing from Scale.rows:\n  posCorr missing=%d | negCorr missing=%d\n', ...
        numel(miss_pos), numel(miss_neg));
end
%% ===================== POST: SCALING ANALYSIS =====================
if do_scaling && isfield(Scale,'rows') && ~isempty(Scale.rows)

    T = Scale.rows;

    % --------- make sure category is string (handles char/cellstr) ----------
    if iscell(T.category)
        T.category = string(T.category);
    elseif ischar(T.category)
        T.category = string(cellstr(T.category));
    else
        T.category = string(T.category);
    end

    % Within-unit z-scoring (trial-to-trial covariation; removes unit scale)
    T.neuralFR_z = NaN(height(T),1);
    T.csmMag_z   = NaN(height(T),1);

    uList = unique(T.unitIdx);

    unitStats = table('Size',[0 9], ...
        'VariableTypes',{'double','double','double','double','double','double','double','double','string'}, ...
        'VariableNames',{'unitIdx','nTrials','rho_spear','p_spear','r_pear','p_pear','slope','intercept','category'});

    for i = 1:numel(uList)
        u  = uList(i);
        ix = (T.unitIdx == u);

        x = T.neuralFR(ix);
        y = T.csmMag(ix);

        if sum(isfinite(x) & isfinite(y)) < min_trials_for_unit_corr
            continue
        end

        % within-unit z
        sx = std(x,0,'omitnan'); if ~isfinite(sx) || sx==0, sx = 1; end
        sy = std(y,0,'omitnan'); if ~isfinite(sy) || sy==0, sy = 1; end
        T.neuralFR_z(ix) = (x - mean(x,'omitnan')) ./ sx;
        T.csmMag_z(ix)   = (y - mean(y,'omitnan')) ./ sy;

        % correlations
        [rho, pR] = corr(x, y, 'type','Spearman','rows','complete');
        [rP,  pP] = corr(x, y, 'type','Pearson', 'rows','complete');

        % linear fit within unit (y = a + b*x)
        P = polyfit(x, y, 1);
        slope = P(1);
        intercept = P(2);

        % --------- string-safe "mode" for category ----------
        cats = string(T.category(ix));
        cats = cats(~ismissing(cats));
        if isempty(cats)
            catU = "negLag_unknown";
        else
            uCats = unique(cats);
            counts = zeros(numel(uCats),1);
            for k = 1:numel(uCats)
                counts(k) = sum(cats == uCats(k));
            end
            [~, imax] = max(counts);
            catU = uCats(imax);
        end

        unitStats = [unitStats; {u, sum(isfinite(x)&isfinite(y)), rho, pR, rP, pP, slope, intercept, catU}]; %#ok<AGROW>
    end

    % Mixed effects: pooled trial-level scaling, random intercept per unit
    lme_raw = [];
    lme_z   = [];
    try
        lme_raw = fitlme(T, 'csmMag ~ neuralFR + (1|unitIdx)');
    catch
        lme_raw = [];
    end
    try
        lme_z   = fitlme(T, 'csmMag_z ~ neuralFR_z + (1|unitIdx)');
    catch
        lme_z = [];
    end

    % Figures
    figure('Color','w','Position',[200 200 900 380]);

    subplot(1,2,1); hold on;
    scatter(T.neuralFR, T.csmMag, 12, 'filled', 'MarkerFaceAlpha',0.25);
    xlabel('Neural response (Hz)'); ylabel('CSM magnitude');
    title('All trials (CSM-mod units)');

    % guard against all-NaN / constant values
    if all(~isnan(T.neuralFR)) && all(~isnan(T.csmMag)) && numel(unique(T.neuralFR(~isnan(T.neuralFR)))) > 1
        xg = linspace(min(T.neuralFR,[],'omitnan'), max(T.neuralFR,[],'omitnan'), 100);
        Pg = polyfit(T.neuralFR, T.csmMag, 1);
        plot(xg, polyval(Pg,xg), 'k-', 'LineWidth',2);
    end

    subplot(1,2,2); hold on;
    scatter(T.neuralFR_z, T.csmMag_z, 12, 'filled', 'MarkerFaceAlpha',0.25);
    xlabel('Neural response (z within unit)'); ylabel('CSM magnitude (z within unit)');
    title('Within-unit scaling');

    if all(~isnan(T.neuralFR_z)) && all(~isnan(T.csmMag_z)) && numel(unique(T.neuralFR_z(~isnan(T.neuralFR_z)))) > 1
        xg2 = linspace(min(T.neuralFR_z,[],'omitnan'), max(T.neuralFR_z,[],'omitnan'), 100);
        Pg2 = polyfit(T.neuralFR_z, T.csmMag_z, 1);
        plot(xg2, polyval(Pg2,xg2), 'k-', 'LineWidth',2);
    end

    % Print summaries
    fprintf('\n===== SCALING SUMMARY =====\n');
    fprintf('Total observations (unit x trial rows): %d\n', height(T));
    fprintf('Units with >=%d trials: %d\n', min_trials_for_unit_corr, height(unitStats));
    if ~isempty(unitStats)
        fprintf('Median Spearman rho across units: %.3f\n', median(unitStats.rho_spear,'omitnan'));
        fprintf('Units with Spearman p<0.05: %d / %d\n', sum(unitStats.p_spear<0.05), height(unitStats));
    end

    if ~isempty(lme_raw)
        disp('--- LME (raw): csmMag ~ neuralFR + (1|unitIdx) ---');
        disp(lme_raw);
    else
        disp('--- LME (raw) FAILED (fitlme unavailable or error). ---');
    end
    if ~isempty(lme_z)
        disp('--- LME (within-unit z): csmMag_z ~ neuralFR_z + (1|unitIdx) ---');
        disp(lme_z);
    else
        disp('--- LME (z) FAILED (fitlme unavailable or error). ---');
    end

    % Store outputs
    Scale.table     = T;
    Scale.unitStats = unitStats;
    Scale.lme_raw   = lme_raw;
    Scale.lme_z     = lme_z;

else
    fprintf('\n[SCALING] No scaling rows collected.\n');
    fprintf('Most common cause: the kinematics alignment in getCSMAlignedKinematics_1ms() is still set to AlignPeakV.\n');
    fprintf('Replace the ONE marked line in the helper with your CSM-alignment function for stim_mod_field.\n');
end


%% ===================== POST: SCALING ANALYSIS =====================
if do_scaling && isfield(Scale,'rows') && ~isempty(Scale.rows)

    T = Scale.rows;

    % --------- ensure category is STRING (handles char/cellstr/categorical) ----------
    if iscell(T.category)
        T.category = string(T.category);
    elseif ischar(T.category)
        T.category = string(cellstr(T.category));
    elseif iscategorical(T.category)
        T.category = string(T.category);
    else
        T.category = string(T.category);
    end

    % --------- define masks for negLag_posCorr vs negLag_negCorr ----------
    isPos = (T.category == "negLag_posCorr");
    isNeg = (T.category == "negLag_negCorr");

    Tpos = T(isPos,:);
    Tneg = T(isNeg,:);

    % Within-unit z-scoring (trial-to-trial covariation; removes unit scale)
    T.neuralFR_z = NaN(height(T),1);
    T.csmMag_z   = NaN(height(T),1);

    uList = unique(T.unitIdx);

    unitStats = table('Size',[0 10], ...
        'VariableTypes',{'double','double','double','double','double','double','double','double','string','string'}, ...
        'VariableNames',{'unitIdx','nTrials','rho_spear','p_spear','r_pear','p_pear','slope','intercept','category','group'});

    % ===================== PER-UNIT STATS (tag each unit as pos/neg) =====================
    for i = 1:numel(uList)
        u  = uList(i);
        ix = (T.unitIdx == u);

        x = T.neuralFR(ix);
        y = T.csmMag(ix);

        nGood = sum(isfinite(x) & isfinite(y));
        if nGood < min_trials_for_unit_corr
            continue
        end

        % within-unit z
        sx = std(x,0,'omitnan'); if ~isfinite(sx) || sx==0, sx = 1; end
        sy = std(y,0,'omitnan'); if ~isfinite(sy) || sy==0, sy = 1; end
        T.neuralFR_z(ix) = (x - mean(x,'omitnan')) ./ sx;
        T.csmMag_z(ix)   = (y - mean(y,'omitnan')) ./ sy;

        % correlations
        [rho, pR] = corr(x, y, 'type','Spearman','rows','complete');
        [rP,  pP] = corr(x, y, 'type','Pearson', 'rows','complete');

        % linear fit within unit (y = a + b*x)
        P = polyfit(x, y, 1);
        slope = P(1);
        intercept = P(2);

        % --------- string-safe "mode" for category (NO function) ----------
        cats = string(T.category(ix));
        cats = cats(~ismissing(cats));
        if isempty(cats)
            catU = "negLag_unknown";
        else
            uCats = unique(cats);
            counts = zeros(numel(uCats),1);
            for k = 1:numel(uCats)
                counts(k) = sum(cats == uCats(k));
            end
            [~, imax] = max(counts);
            catU = uCats(imax);
        end

        % group label
        if catU == "negLag_posCorr"
            grp = "posCorr";
        elseif catU == "negLag_negCorr"
            grp = "negCorr";
        else
            grp = "unknown";
        end

        unitStats = [unitStats; {u, nGood, rho, pR, rP, pP, slope, intercept, catU, grp}]; %#ok<AGROW>
    end

    % IMPORTANT: after computing within-unit z, refresh Tpos/Tneg to include *_z columns
    Tpos = T(isPos,:);
    Tneg = T(isNeg,:);

    % ===================== MIXED EFFECTS (SEPARATE BY GROUP) =====================
    lme = struct();
    lme.raw_all = []; lme.z_all = [];
    lme.raw_pos = []; lme.z_pos = [];
    lme.raw_neg = []; lme.z_neg = [];

    % All
    try, lme.raw_all = fitlme(T, 'csmMag ~ neuralFR + (1|unitIdx)'); catch, lme.raw_all = []; end
    try, lme.z_all   = fitlme(T, 'csmMag_z ~ neuralFR_z + (1|unitIdx)'); catch, lme.z_all = []; end

    % posCorr
    if ~isempty(Tpos) && numel(unique(Tpos.unitIdx)) >= 2
        try, lme.raw_pos = fitlme(Tpos, 'csmMag ~ neuralFR + (1|unitIdx)'); catch, lme.raw_pos = []; end
        try, lme.z_pos   = fitlme(Tpos, 'csmMag_z ~ neuralFR_z + (1|unitIdx)'); catch, lme.z_pos = []; end
    end

    % negCorr
    if ~isempty(Tneg) && numel(unique(Tneg.unitIdx)) >= 2
        try, lme.raw_neg = fitlme(Tneg, 'csmMag ~ neuralFR + (1|unitIdx)'); catch, lme.raw_neg = []; end
        try, lme.z_neg   = fitlme(Tneg, 'csmMag_z ~ neuralFR_z + (1|unitIdx)'); catch, lme.z_neg = []; end
    end

    % ===================== FIGURES: separate pos vs neg =====================
    figure('Color','w','Position',[200 200 1050 720]);

    % ---- RAW: posCorr ----
    subplot(2,2,1); hold on;
    scatter(Tpos.neuralFR, Tpos.csmMag, 12, 'filled', 'MarkerFaceAlpha',0.25);
    xlabel('Neural response (Hz)'); ylabel('CSM magnitude');
    title('posCorr units (raw)');
    if ~isempty(Tpos)
        xNon = Tpos.neuralFR(~isnan(Tpos.neuralFR));
        if numel(unique(xNon)) > 1
            xg = linspace(min(Tpos.neuralFR,[],'omitnan'), max(Tpos.neuralFR,[],'omitnan'), 100);
            Pg = polyfit(Tpos.neuralFR, Tpos.csmMag, 1);
            plot(xg, polyval(Pg,xg), 'k-', 'LineWidth',2);
        end
    end

    % ---- RAW: negCorr ----
    subplot(2,2,2); hold on;
    scatter(Tneg.neuralFR, Tneg.csmMag, 12, 'filled', 'MarkerFaceAlpha',0.25);
    xlabel('Neural response (Hz)'); ylabel('CSM magnitude');
    title('negCorr units (raw)');
    if ~isempty(Tneg)
        xNon = Tneg.neuralFR(~isnan(Tneg.neuralFR));
        if numel(unique(xNon)) > 1
            xg = linspace(min(Tneg.neuralFR,[],'omitnan'), max(Tneg.neuralFR,[],'omitnan'), 100);
            Pg = polyfit(Tneg.neuralFR, Tneg.csmMag, 1);
            plot(xg, polyval(Pg,xg), 'k-', 'LineWidth',2);
        end
    end

    % ---- Z within unit: posCorr ----
    subplot(2,2,3); hold on;
    scatter(Tpos.neuralFR_z, Tpos.csmMag_z, 12, 'filled', 'MarkerFaceAlpha',0.25);
    xlabel('Neural response (z within unit)'); ylabel('CSM magnitude (z within unit)');
    title('posCorr units (within-unit z)');
    if ~isempty(Tpos)
        xNon = Tpos.neuralFR_z(~isnan(Tpos.neuralFR_z));
        if numel(unique(xNon)) > 1
            xg2 = linspace(min(Tpos.neuralFR_z,[],'omitnan'), max(Tpos.neuralFR_z,[],'omitnan'), 100);
            Pg2 = polyfit(Tpos.neuralFR_z, Tpos.csmMag_z, 1);
            plot(xg2, polyval(Pg2,xg2), 'k-', 'LineWidth',2);
        end
    end

    % ---- Z within unit: negCorr ----
    subplot(2,2,4); hold on;
    scatter(Tneg.neuralFR_z, Tneg.csmMag_z, 12, 'filled', 'MarkerFaceAlpha',0.25);
    xlabel('Neural response (z within unit)'); ylabel('CSM magnitude (z within unit)');
    title('negCorr units (within-unit z)');
    if ~isempty(Tneg)
        xNon = Tneg.neuralFR_z(~isnan(Tneg.neuralFR_z));
        if numel(unique(xNon)) > 1
            xg2 = linspace(min(Tneg.neuralFR_z,[],'omitnan'), max(Tneg.neuralFR_z,[],'omitnan'), 100);
            Pg2 = polyfit(Tneg.neuralFR_z, Tneg.csmMag_z, 1);
            plot(xg2, polyval(Pg2,xg2), 'k-', 'LineWidth',2);
        end
    end

    % ===================== PRINT SUMMARIES =====================
    fprintf('\n===== SCALING SUMMARY (ALL) =====\n');
    fprintf('Total observations: %d\n', height(T));
    fprintf('Units with >=%d trials: %d\n', min_trials_for_unit_corr, height(unitStats));
    if ~isempty(unitStats)
        fprintf('Median Spearman rho (all units): %.3f\n', median(unitStats.rho_spear,'omitnan'));
        fprintf('Units with Spearman p<0.05 (all): %d / %d\n', sum(unitStats.p_spear<0.05), height(unitStats));
    end

    unitPos = unitStats(unitStats.group=="posCorr",:);
    unitNeg = unitStats(unitStats.group=="negCorr",:);

    fprintf('\n===== SCALING SUMMARY (posCorr units) =====\n');
    fprintf('Units: %d\n', height(unitPos));
    if ~isempty(unitPos)
        fprintf('Median Spearman rho: %.3f\n', median(unitPos.rho_spear,'omitnan'));
        fprintf('Units with Spearman p<0.05: %d / %d\n', sum(unitPos.p_spear<0.05), height(unitPos));
    end

    fprintf('\n===== SCALING SUMMARY (negCorr units) =====\n');
    fprintf('Units: %d\n', height(unitNeg));
    if ~isempty(unitNeg)
        fprintf('Median Spearman rho: %.3f\n', median(unitNeg.rho_spear,'omitnan'));
        fprintf('Units with Spearman p<0.05: %d / %d\n', sum(unitNeg.p_spear<0.05), height(unitNeg));
    end

    % Print LMEs
    if ~isempty(lme.raw_pos)
        disp('--- LME posCorr (raw): csmMag ~ neuralFR + (1|unitIdx) ---'); disp(lme.raw_pos);
    else
        disp('--- LME posCorr (raw) FAILED / insufficient data ---');
    end
    if ~isempty(lme.raw_neg)
        disp('--- LME negCorr (raw): csmMag ~ neuralFR + (1|unitIdx) ---'); disp(lme.raw_neg);
    else
        disp('--- LME negCorr (raw) FAILED / insufficient data ---');
    end
    if ~isempty(lme.z_pos)
        disp('--- LME posCorr (z): csmMag_z ~ neuralFR_z + (1|unitIdx) ---'); disp(lme.z_pos);
    else
        disp('--- LME posCorr (z) FAILED / insufficient data ---');
    end
    if ~isempty(lme.z_neg)
        disp('--- LME negCorr (z): csmMag_z ~ neuralFR_z + (1|unitIdx) ---'); disp(lme.z_neg);
    else
        disp('--- LME negCorr (z) FAILED / insufficient data ---');
    end

    % Store outputs
    Scale.table     = T;
    Scale.unitStats = unitStats;
    Scale.lme       = lme;
    Scale.Tpos      = Tpos;
    Scale.Tneg      = Tneg;

else
    fprintf('\n[SCALING] No scaling rows collected.\n');
    fprintf('Most common cause: kinematics alignment helper not matched to stim_mod_field.\n');
end

%% ===================== PER-UNIT CORRELATION PLOTS + SUMMARY =====================
% Uses unitStats.rho_spear (one rho per unit). Also splits by posCorr/negCorr.

if exist('unitStats','var') && ~isempty(unitStats)

    % Keep only pos/neg corr units (drop unknown)
    unitPos = unitStats(unitStats.group=="posCorr",:);
    unitNeg = unitStats(unitStats.group=="negCorr",:);

    % --- bootstrap CI helper (no functions; inline) ---
    nBoot = 5000;
    alphaCI = 0.05;

    % Prepare figure
    figure('Color','w','Position',[200 200 900 380]);

    % ===================== PANEL A: per-unit rho scatter =====================
    subplot(1,2,1); hold on;

    % x positions
    xPos = ones(height(unitPos),1) * 1;
    xNeg = ones(height(unitNeg),1) * 2;

    % jitter
    jit = 0.10;
    xPos = xPos + (rand(size(xPos))-0.5)*2*jit;
    xNeg = xNeg + (rand(size(xNeg))-0.5)*2*jit;

    % points
    scatter(xPos, unitPos.rho_spear, 35, 'filled', 'MarkerFaceAlpha',0.6);
    scatter(xNeg, unitNeg.rho_spear, 35, 'filled', 'MarkerFaceAlpha',0.6);

    % group means + boot CI
    % ---- pos ----
    if ~isempty(unitPos)
        vals = unitPos.rho_spear;
        vals = vals(isfinite(vals));
        if ~isempty(vals)
            bootMean = zeros(nBoot,1);
            n = numel(vals);
            for b = 1:nBoot
                samp = vals(randi(n,n,1));
                bootMean(b) = mean(samp,'omitnan');
            end
            mPos = mean(vals,'omitnan');
            ciPos = prctile(bootMean, [100*alphaCI/2, 100*(1-alphaCI/2)]);
            plot([0.85 1.15],[mPos mPos],'k-','LineWidth',2);
            plot([1 1],[ciPos(1) ciPos(2)],'k-','LineWidth',2);
        end
    end

    % ---- neg ----
    if ~isempty(unitNeg)
        vals = unitNeg.rho_spear;
        vals = vals(isfinite(vals));
        if ~isempty(vals)
            bootMean = zeros(nBoot,1);
            n = numel(vals);
            for b = 1:nBoot
                samp = vals(randi(n,n,1));
                bootMean(b) = mean(samp,'omitnan');
            end
            mNeg = mean(vals,'omitnan');
            ciNeg = prctile(bootMean, [100*alphaCI/2, 100*(1-alphaCI/2)]);
            plot([1.85 2.15],[mNeg mNeg],'k-','LineWidth',2);
            plot([2 2],[ciNeg(1) ciNeg(2)],'k-','LineWidth',2);
        end
    end

    yline(0,'--');
    xlim([0.5 2.5]);
    xticks([1 2]);
    xticklabels({'posCorr','negCorr'});
    ylabel('Per-unit Spearman \rho (neuralFR vs CSM magnitude)');
    title('Per-unit scaling correlations');

    % ===================== PANEL B: histogram / density-ish =====================
    subplot(1,2,2); hold on;
    if ~isempty(unitPos)
        histogram(unitPos.rho_spear, 'BinWidth',0.1, 'Normalization','count');
    end
    if ~isempty(unitNeg)
        histogram(unitNeg.rho_spear, 'BinWidth',0.1, 'Normalization','count');
    end
    xline(0,'--');
    xlabel('Per-unit Spearman \rho');
    ylabel('Probability');
    title('Distribution of per-unit \rho');

    % ===================== PRINT SUMMARY STATS =====================
    fprintf('\n===== PER-UNIT CORRELATION SUMMARY =====\n');

    if ~isempty(unitPos)
        vals = unitPos.rho_spear; vals = vals(isfinite(vals));
        fprintf('posCorr units: n=%d | mean rho=%.3f | median rho=%.3f\n', ...
            numel(vals), mean(vals,'omitnan'), median(vals,'omitnan'));
        try
            p0 = signrank(vals, 0); % test median != 0
            fprintf('  signrank vs 0: p=%.4g\n', p0);
        catch
        end
    else
        fprintf('posCorr units: n=0\n');
    end

    if ~isempty(unitNeg)
        vals = unitNeg.rho_spear; vals = vals(isfinite(vals));
        fprintf('negCorr units: n=%d | mean rho=%.3f | median rho=%.3f\n', ...
            numel(vals), mean(vals,'omitnan'), median(vals,'omitnan'));
        try
            p0 = signrank(vals, 0);
            fprintf('  signrank vs 0: p=%.4g\n', p0);
        catch
        end
    else
        fprintf('negCorr units: n=0\n');
    end

    % Optional: compare pos vs neg distributions (rank-sum)
    try
        v1 = unitPos.rho_spear; v1 = v1(isfinite(v1));
        v2 = unitNeg.rho_spear; v2 = v2(isfinite(v2));
        if ~isempty(v1) && ~isempty(v2)
            pComp = ranksum(v1, v2);
            fprintf('pos vs neg (ranksum): p=%.4g\n', pComp);
        end
    catch
    end

    % Store
    Scale.unitPos = unitPos;
    Scale.unitNeg = unitNeg;

end

%% ===================== PER-UNIT CORRELATION PLOTS + SUMMARY =====================
% EDITED TO USE PEARSON CORRELATION (r) INSTEAD OF SPEARMAN (rho)
% Assumes unitStats has one row per unit and contains:
%   - unitStats.group  ("posCorr" or "negCorr")
%   - unitStats.r_pear  (Pearson r for that unit; one value per unit)
% If your column is named differently, change rCol below.

if exist('unitStats','var') && ~isempty(unitStats)

    % ---- NEW: set which column holds Pearson r ----
    rCol = "r_pear";                 % <-- CHANGE THIS if needed
    % OLD (Spearman):
    % rCol = "rho_spear";

    % Keep only pos/neg corr units (drop unknown)
    unitPos = unitStats(unitStats.group=="posCorr",:);
    unitNeg = unitStats(unitStats.group=="negCorr",:);

    % --- bootstrap CI helper (no functions; inline) ---
    nBoot = 5000;
    alphaCI = 0.05;

    % Prepare figure
    figure('Color','w','Position',[200 200 900 380]);

    % ===================== PANEL A: per-unit r scatter =====================
    subplot(1,2,1); hold on;

    % x positions
    xPos = ones(height(unitPos),1) * 1;
    xNeg = ones(height(unitNeg),1) * 2;

    % jitter
    jit = 0.10;
    xPos = xPos + (rand(size(xPos))-0.5)*2*jit;
    xNeg = xNeg + (rand(size(xNeg))-0.5)*2*jit;

    % points
    scatter(xPos, unitPos.(rCol), 35, 'filled', 'MarkerFaceAlpha',0.6);
    scatter(xNeg, unitNeg.(rCol), 35, 'filled', 'MarkerFaceAlpha',0.6);

    % group means + boot CI
    % ---- pos ----
    if ~isempty(unitPos)
        vals = unitPos.(rCol);
        vals = vals(isfinite(vals));
        if ~isempty(vals)
            bootMean = zeros(nBoot,1);
            n = numel(vals);
            for b = 1:nBoot
                samp = vals(randi(n,n,1));
                bootMean(b) = mean(samp,'omitnan');
            end
            mPos = mean(vals,'omitnan');
            ciPos = prctile(bootMean, [100*alphaCI/2, 100*(1-alphaCI/2)]);
            plot([0.85 1.15],[mPos mPos],'k-','LineWidth',2);
            plot([1 1],[ciPos(1) ciPos(2)],'k-','LineWidth',2);
        end
    end

    % ---- neg ----
    if ~isempty(unitNeg)
        vals = unitNeg.(rCol);
        vals = vals(isfinite(vals));
        if ~isempty(vals)
            bootMean = zeros(nBoot,1);
            n = numel(vals);
            for b = 1:nBoot
                samp = vals(randi(n,n,1));
                bootMean(b) = mean(samp,'omitnan');
            end
            mNeg = mean(vals,'omitnan');
            ciNeg = prctile(bootMean, [100*alphaCI/2, 100*(1-alphaCI/2)]);
            plot([1.85 2.15],[mNeg mNeg],'k-','LineWidth',2);
            plot([2 2],[ciNeg(1) ciNeg(2)],'k-','LineWidth',2);
        end
    end

    yline(0,'--');
    xlim([0.5 2.5]);
    xticks([1 2]);
    xticklabels({'posCorr','negCorr'});
    ylabel('Per-unit Pearson r (neuralFR vs CSM magnitude)');
    title('Per-unit scaling correlations');

    % ===================== PANEL B: histogram / density-ish =====================
    subplot(1,2,2); hold on;
    if ~isempty(unitPos)
        histogram(unitPos.(rCol), 'BinWidth',0.1, 'Normalization','count');
    end
    if ~isempty(unitNeg)
        histogram(unitNeg.(rCol), 'BinWidth',0.1, 'Normalization','count');
    end
    xline(0,'--');
    xlabel('Per-unit Pearson r');
    ylabel('Count');
    title('Distribution of per-unit r');

    % ===================== PRINT SUMMARY STATS =====================
    fprintf('\n===== PER-UNIT CORRELATION SUMMARY (PEARSON r) =====\n');

    if ~isempty(unitPos)
        vals = unitPos.(rCol); vals = vals(isfinite(vals));
        fprintf('posCorr units: n=%d | mean r=%.3f | median r=%.3f\n', ...
            numel(vals), mean(vals,'omitnan'), median(vals,'omitnan'));
        % NOTE: signrank tests median vs 0 (nonparam). Keeping it as-is.
        try
            p0 = signrank(vals, 0);
            fprintf('  signrank vs 0: p=%.4g\n', p0);
        catch
        end
    else
        fprintf('posCorr units: n=0\n');
    end

    if ~isempty(unitNeg)
        vals = unitNeg.(rCol); vals = vals(isfinite(vals));
        fprintf('negCorr units: n=%d | mean r=%.3f | median r=%.3f\n', ...
            numel(vals), mean(vals,'omitnan'), median(vals,'omitnan'));
        try
            p0 = signrank(vals, 0);
            fprintf('  signrank vs 0: p=%.4g\n', p0);
        catch
        end
    else
        fprintf('negCorr units: n=0\n');
    end

    % Optional: compare pos vs neg distributions (rank-sum)
    try
        v1 = unitPos.(rCol); v1 = v1(isfinite(v1));
        v2 = unitNeg.(rCol); v2 = v2(isfinite(v2));
        if ~isempty(v1) && ~isempty(v2)
            pComp = ranksum(v1, v2);
            fprintf('pos vs neg (ranksum): p=%.4g\n', pComp);
        end
    catch
    end

    % Store
    Scale.unitPos = unitPos;
    Scale.unitNeg = unitNeg;

end

%% FIGURE EDIT 
%% ===================== PER-UNIT CORRELATION PLOTS + SUMMARY =====================
% Mean shown as LARGE CIRCLE with bootstrap CI error bars

if exist('unitStats','var') && ~isempty(unitStats)

    rCol = "r_pear";   % Pearson column

    unitPos = unitStats(unitStats.group=="posCorr",:);
    unitNeg = unitStats(unitStats.group=="negCorr",:);

    nBoot = 5000;
    alphaCI = 0.05;

    figure('Color','w','Position',[200 200 500 380]);

    % ===================== PANEL A =====================
    subplot(1,2,1); hold on;

    % x positions
    xPos = ones(height(unitPos),1);
    xNeg = ones(height(unitNeg),1)*2;

    % jitter
    jit = 0.10;
    xPos = xPos + (rand(size(xPos))-0.5)*2*jit;
    xNeg = xNeg + (rand(size(xNeg))-0.5)*2*jit;

    % scatter points
    scatter(xPos, unitPos.(rCol), 35, 'filled', 'MarkerFaceAlpha',0.6);
    scatter(xNeg, unitNeg.(rCol), 35, 'filled', 'MarkerFaceAlpha',0.6);

    % ---------- POS MEAN + CI ----------
    if ~isempty(unitPos)
        vals = unitPos.(rCol);
        vals = vals(isfinite(vals));

        if ~isempty(vals)
            n = numel(vals);
            bootMean = zeros(nBoot,1);
            for b = 1:nBoot
                samp = vals(randi(n,n,1));
                bootMean(b) = mean(samp,'omitnan');
            end

            mPos = mean(vals,'omitnan');
            ciPos = prctile(bootMean,[100*alphaCI/2 100*(1-alphaCI/2)]);

            % LARGE MEAN DOT
            scatter(1, mPos, 140, 'k', 'filled', 'MarkerFaceAlpha',1);

            % ERROR BAR
            errorbar(1, mPos, mPos-ciPos(1), ciPos(2)-mPos, ...
                'k','LineWidth',2,'CapSize',10);
        end
    end

    % ---------- NEG MEAN + CI ----------
    if ~isempty(unitNeg)
        vals = unitNeg.(rCol);
        vals = vals(isfinite(vals));

        if ~isempty(vals)
            n = numel(vals);
            bootMean = zeros(nBoot,1);
            for b = 1:nBoot
                samp = vals(randi(n,n,1));
                bootMean(b) = mean(samp,'omitnan');
            end

            mNeg = mean(vals,'omitnan');
            ciNeg = prctile(bootMean,[100*alphaCI/2 100*(1-alphaCI/2)]);

            scatter(2, mNeg, 140, 'k', 'filled', 'MarkerFaceAlpha',1);

            errorbar(2, mNeg, mNeg-ciNeg(1), ciNeg(2)-mNeg, ...
                'k','LineWidth',2,'CapSize',10);
        end
    end

    yline(0,'--');
    xlim([0.5 2.5]);
    xticks([1 2]);
    xticklabels({'posCorr','negCorr'});
    ylabel('Per-unit Pearson r (neuralFR vs CSM magnitude)');
    title('Per-unit scaling correlations');

    % ===================== PANEL B =====================
    subplot(1,2,2); hold on;
    if ~isempty(unitPos)
        histogram(unitPos.(rCol), 'BinWidth',0.1);
    end
    if ~isempty(unitNeg)
        histogram(unitNeg.(rCol), 'BinWidth',0.1);
    end
    xline(0,'--');
    xlabel('Per-unit Pearson r');
    ylabel('Count');
    title('Distribution of per-unit r');

    % ===================== PRINT SUMMARY =====================
    fprintf('\n===== PER-UNIT CORRELATION SUMMARY (PEARSON r) =====\n');

    if ~isempty(unitPos)
        vals = unitPos.(rCol); vals = vals(isfinite(vals));
        fprintf('posCorr units: n=%d | mean r=%.3f | median r=%.3f\n', ...
            numel(vals), mean(vals,'omitnan'), median(vals,'omitnan'));
    else
        fprintf('posCorr units: n=0\n');
    end

    if ~isempty(unitNeg)
        vals = unitNeg.(rCol); vals = vals(isfinite(vals));
        fprintf('negCorr units: n=%d | mean r=%.3f | median r=%.3f\n', ...
            numel(vals), mean(vals,'omitnan'), median(vals,'omitnan'));
    else
        fprintf('negCorr units: n=0\n');
    end

    Scale.unitPos = unitPos;
    Scale.unitNeg = unitNeg;

end

%% ===================== VIS: SPAGHETTI PER-UNIT LINES + POOLED FIT + CI (pos vs neg) =====================
% Each thin line = per-unit linear fit across that unit's trials
% Thick line + shaded band = pooled linear fit with bootstrap CI
%
% Requires Scale.table columns: unitIdx, neuralFR, csmMag, category

if do_scaling && isfield(Scale,'table') && ~isempty(Scale.table)

    T = Scale.table;

    % --- ensure category is string ---
    if iscell(T.category)
        T.category = string(T.category);
    elseif ischar(T.category)
        T.category = string(cellstr(T.category));
    elseif iscategorical(T.category)
        T.category = string(T.category);
    else
        T.category = string(T.category);
    end

    % groups
    isPos = (T.category == "negLag_posCorr");
    isNeg = (T.category == "negLag_negCorr");

    % settings
    minTrialsPerUnit_line = max(5, min_trials_for_unit_corr);  % per-unit line requires this many valid trials
    nGrid = 120;               % x-grid for pooled fit
    nBoot = 1000;              % bootstrap draws for CI band
    alphaCI = 0.05;            % 95% CI
    thinAlpha = 0.18;          % per-unit line transparency
    thinLW = 1;                % per-unit line width
    thickLW = 3;               % pooled line width

    % helper inline: make panel for one group
    figure('Color','w','Position',[200 200 980 420]);

    for panel = 1:2

        if panel == 1
            Tg = T(isPos,:);
            panelTitle = "posCorr units";
        else
            Tg = T(isNeg,:);
            panelTitle = "negCorr units";
        end

        % keep finite
        Tg = Tg(isfinite(Tg.neuralFR) & isfinite(Tg.csmMag),:);

        subplot(1,2,panel); hold on;

        if isempty(Tg)
            text(0.1,0.5,'Not enough data','Units','normalized');
            axis off
            continue
        end

        % --------- define common x-grid for lines ----------
        xMin = min(Tg.csmMag,[],'omitnan');
        xMax = max(Tg.csmMag,[],'omitnan');
        if ~isfinite(xMin) || ~isfinite(xMax) || xMax==xMin
            text(0.1,0.5,'CSM magnitude has no range','Units','normalized');
            axis off
            continue
        end
        xGrid = linspace(xMin, xMax, nGrid)';

        % --------- thin per-unit lines ----------
        uList = unique(Tg.unitIdx);
        nLines = 0;
        for iU = 1:numel(uList)
            u = uList(iU);
            ix = (Tg.unitIdx == u);

            x = Tg.csmMag(ix);
            y = Tg.neuralFR(ix);

            good = isfinite(x) & isfinite(y);
            if sum(good) < minTrialsPerUnit_line
                continue
            end

            % linear fit y = a + b*x
            P = polyfit(x(good), y(good), 1);
            yHat = polyval(P, xGrid);

            % draw thin line
            h = plot(xGrid, yHat, 'k-', 'LineWidth',thinLW);
            % set transparency (MATLAB lines don't support alpha directly; use color with alpha if available)
            try
                h.Color = [0 0 0 thinAlpha];
            catch
                % older MATLAB: approximate by lighter gray
                h.Color = [0.65 0.65 0.65];
            end

            nLines = nLines + 1;
        end

        % --------- pooled fit (thick) ----------
        xAll = Tg.csmMag;
        yAll = Tg.neuralFR;

        Pg = polyfit(xAll, yAll, 1);
        yPooled = polyval(Pg, xGrid);

        % --------- bootstrap CI band around pooled fit ----------
        % resample rows (trials) WITH replacement, fit line each time, eval on xGrid
        Yboot = nan(nGrid, nBoot);
        nObs = height(Tg);
        for b = 1:nBoot
            ii = randi(nObs, nObs, 1);
            xb = xAll(ii);
            yb = yAll(ii);

            % guard: if xb has no variance, skip
            if numel(unique(xb)) < 2
                continue
            end

            Pb = polyfit(xb, yb, 1);
            Yboot(:,b) = polyval(Pb, xGrid);
        end

        lo = prctile(Yboot', 100*(alphaCI/2))';          % [nGrid x 1]
        hi = prctile(Yboot', 100*(1-alphaCI/2))';

        okBand = isfinite(lo) & isfinite(hi);
        if any(okBand)
            xx = [xGrid(okBand); flipud(xGrid(okBand))];
            yy = [lo(okBand);    flipud(hi(okBand))];
            fill(xx, yy, 'k', 'FaceAlpha',0.12, 'EdgeColor','none');
        end

        plot(xGrid, yPooled, 'k-', 'LineWidth',thickLW);

        % --------- annotate pooled rho too ----------
        rho = NaN; pR = NaN;
        if height(Tg) >= 10 && numel(unique(xAll)) > 1
            [rho, pR] = corr(xAll, yAll, 'type','Spearman', 'rows','complete');
        end

        xlabel('CSM magnitude');
        ylabel('Neural response (Hz)');
        title(sprintf('%s: nUnits(lines)=%d | nTrials=%d | pooled \\rho=%.3f (p=%.3g)', ...
            panelTitle, nLines, height(Tg), rho, pR));

        % If you want to also SEE trial points faintly, uncomment:
        % scatter(xAll, yAll, 6, 'filled', 'MarkerFaceAlpha',0.05, 'MarkerEdgeAlpha',0);

        box off
    end

else
    fprintf('\n[SPAGHETTI+POOLED] Need Scale.table first.\n');
end


%% ===================== VIS: SPAGHETTI PER-UNIT LINES + POOLED FIT + CI (pos vs neg) =====================
% Adds:
%   - pooled Pearson r and R^2
%   - per-unit lines colored by significance (Pearson p<0.05 by default)
%
% Requires Scale.table columns: unitIdx, neuralFR, csmMag, category

if do_scaling && isfield(Scale,'table') && ~isempty(Scale.table)

    T = Scale.table;

    % --- ensure category is string ---
    if iscell(T.category)
        T.category = string(T.category);
    elseif ischar(T.category)
        T.category = string(cellstr(T.category));
    elseif iscategorical(T.category)
        T.category = string(T.category);
    else
        T.category = string(T.category);
    end

    % groups
    isPos = (T.category == "negLag_posCorr");
    isNeg = (T.category == "negLag_negCorr");

    % settings
    minTrialsPerUnit_line = max(5, min_trials_for_unit_corr);  % per-unit line needs this many valid trials
    nGrid   = 120;         % x-grid for pooled fit
    nBoot   = 1000;        % bootstrap draws for CI band
    alphaCI = 0.05;        % 95% CI

    % line styling
    lw_nonsig = 1;
    lw_sig    = 1.6;
    lw_pooled = 3;

    % Colors (MATLAB default-friendly)
    col_nonsig = [0.75 0.75 0.75];  % light gray
    col_sig    = [0.10 0.10 0.10];  % near-black (use colored if you prefer)
    % e.g., to use a color for sig lines:
    % col_sig = [0.0 0.45 0.74];   % blue

    % Significance criterion: Pearson p<0.05 for per-unit lines
    alpha_unit = 0.05;

    figure('Color','w','Position',[200 200 980 420]);

    for panel = 1:2

        if panel == 1
            Tg = T(isPos,:);
            panelTitle = "posCorr units";
        else
            Tg = T(isNeg,:);
            panelTitle = "negCorr units";
        end

        % keep finite
        Tg = Tg(isfinite(Tg.neuralFR) & isfinite(Tg.csmMag),:);

        subplot(1,2,panel); hold on;

        if isempty(Tg)
            text(0.1,0.5,'Not enough data','Units','normalized');
            axis off
            continue
        end

        % common x-grid
        xMin = min(Tg.csmMag,[],'omitnan');
        xMax = max(Tg.csmMag,[],'omitnan');
        if ~isfinite(xMin) || ~isfinite(xMax) || xMax==xMin
            text(0.1,0.5,'CSM magnitude has no range','Units','normalized');
            axis off
            continue
        end
        xGrid = linspace(xMin, xMax, nGrid)';

        % --------- per-unit lines ----------
        uList = unique(Tg.unitIdx);
        nLines = 0;
        nSig   = 0;

        for iU = 1:numel(uList)
            u = uList(iU);
            ix = (Tg.unitIdx == u);

            x = Tg.csmMag(ix);
            y = Tg.neuralFR(ix);

            good = isfinite(x) & isfinite(y);
            if sum(good) < minTrialsPerUnit_line
                continue
            end

            xg = x(good);
            yg = y(good);

            % per-unit Pearson correlation p-value (for coloring)
            pUnit = NaN; rUnit = NaN;
            if numel(unique(xg)) > 1
                [rUnit, pUnit] = corr(xg, yg, 'type','Pearson', 'rows','complete');
            end
            isSig = isfinite(pUnit) && (pUnit < alpha_unit);

            % linear fit y = a + b*x
            P = polyfit(xg, yg, 1);
            yHat = polyval(P, xGrid);

            % plot
            if isSig
                plot(xGrid, yHat, '-', 'Color',col_sig, 'LineWidth',lw_sig);
                nSig = nSig + 1;
            else
                plot(xGrid, yHat, '-', 'Color',col_nonsig, 'LineWidth',lw_nonsig);
            end

            nLines = nLines + 1;
        end

        % --------- pooled fit + pooled r and R^2 ----------
        xAll = Tg.csmMag;
        yAll = Tg.neuralFR;

        % pooled Pearson r
        rPool = NaN; pPool = NaN;
        if height(Tg) >= 10 && numel(unique(xAll)) > 1
            [rPool, pPool] = corr(xAll, yAll, 'type','Pearson', 'rows','complete');
        end

        % pooled linear fit + R^2
        Pg = polyfit(xAll, yAll, 1);
        yPooled = polyval(Pg, xGrid);

        yHatAll = polyval(Pg, xAll);
        SSres = sum((yAll - yHatAll).^2, 'omitnan');
        SStot = sum((yAll - mean(yAll,'omitnan')).^2, 'omitnan');
        if ~isfinite(SSres) || ~isfinite(SStot) || SStot==0
            R2pool = NaN;
        else
            R2pool = 1 - (SSres / SStot);
        end

        % --------- bootstrap CI band around pooled fit (resample TRIALS) ----------
        Yboot = nan(nGrid, nBoot);
        nObs = height(Tg);

        for b = 1:nBoot
            ii = randi(nObs, nObs, 1);
            xb = xAll(ii);
            yb = yAll(ii);

            if numel(unique(xb)) < 2
                continue
            end

            Pb = polyfit(xb, yb, 1);
            Yboot(:,b) = polyval(Pb, xGrid);
        end

        lo = prctile(Yboot', 100*(alphaCI/2))';
        hi = prctile(Yboot', 100*(1-alphaCI/2))';

        okBand = isfinite(lo) & isfinite(hi);
        if any(okBand)
            xx = [xGrid(okBand); flipud(xGrid(okBand))];
            yy = [lo(okBand);    flipud(hi(okBand))];
            fill(xx, yy, 'k', 'FaceAlpha',0.12, 'EdgeColor','none');
        end

        % pooled line on top
        plot(xGrid, yPooled, 'k-', 'LineWidth',lw_pooled);

        % --------- labels + annotation ----------
        xlabel('CSM magnitude');
        ylabel('Neural response (Hz)');

        title(sprintf('%s: nUnits=%d (sig=%d) | nTrials=%d', panelTitle, nLines, nSig, height(Tg)));

        % Put r + R^2 in the corner
        txt = sprintf('Pooled Pearson r = %.3f (p=%.3g)\nPooled R^2 = %.3f', rPool, pPool, R2pool);
        text(0.03, 0.97, txt, 'Units','normalized', 'VerticalAlignment','top', ...
            'FontSize',10, 'BackgroundColor','w', 'Margin',4);

        box off
    end

    % Store counts if you want
    Scale.spaghetti.minTrialsPerUnit_line = minTrialsPerUnit_line;
    Scale.spaghetti.nBoot = nBoot;

else
    fprintf('\n[SPAGHETTI+POOLED] Need Scale.table first.\n');
end

%% ===================== VIS: SPAGHETTI PER-UNIT LINES + MEAN-OF-UNITS FIT + UNIT-BOOT CI (pos vs neg) =====================
% Fixes the mismatch you observed:
%   - REPLACES pooled fit (trial-weighted) with MEAN OF PER-UNIT REGRESSIONS (unit-weighted)
%   - CI band is bootstrap over UNITS (not trials)
%   - Annotates mean per-unit r (Fisher) + mean slope (+ optional bootstrap p for slope)
%
% Requires Scale.table columns: unitIdx, neuralFR, csmMag, category

if do_scaling && isfield(Scale,'table') && ~isempty(Scale.table)

    T = Scale.table;

    % --- ensure category is string ---
    if iscell(T.category)
        T.category = string(T.category);
    elseif ischar(T.category)
        T.category = string(cellstr(T.category));
    elseif iscategorical(T.category)
        T.category = string(T.category);
    else
        T.category = string(T.category);
    end

    % groups
    isPos = (T.category == "negLag_posCorr");
    isNeg = (T.category == "negLag_negCorr");

    % settings
    minTrialsPerUnit_line = max(5, min_trials_for_unit_corr);
    nGrid   = 120;
    nBoot   = 2000;     % bootstrap over UNITS
    alphaCI = 0.05;
    alpha_unit = 0.05;  % for coloring per-unit lines

    % line styling
    lw_nonsig = 1;
    lw_sig    = 1.6;
    lw_mean   = 3;

    % Colors
    col_nonsig = [0.75 0.75 0.75];
    col_sig    = [0.10 0.10 0.10];

    figure('Color','w','Position',[200 200 980 420]);

    for panel = 1:2

        if panel == 1
            Tg = T(isPos,:);
            panelTitle = "posCorrunits";
        else
            Tg = T(isNeg,:);
            panelTitle = "negCorrunits";
        end

        % keep finite
        Tg = Tg(isfinite(Tg.neuralFR) & isfinite(Tg.csmMag),:);

        subplot(1,2,panel); hold on;

        if isempty(Tg)
            text(0.1,0.5,'Not enough data','Units','normalized');
            axis off
            continue
        end

        % common x-grid (based on this group's x-range)
        xMin = min(Tg.csmMag,[],'omitnan');
        xMax = max(Tg.csmMag,[],'omitnan');
        if ~isfinite(xMin) || ~isfinite(xMax) || xMax==xMin
            text(0.1,0.5,'CSM magnitude has no range','Units','normalized');
            axis off
            continue
        end
        xGrid = linspace(xMin, xMax, nGrid)';

        % --------- per-unit lines + store each unit's fitted line ----------
        uList = unique(Tg.unitIdx);

        nLines = 0;
        nSig   = 0;

        rUnits = [];
        pUnits = [];
        slopes = [];
        ints   = [];

        YhatUnits = [];  % nGrid x nUnitsKept

        for iU = 1:numel(uList)

            u = uList(iU);
            ix = (Tg.unitIdx == u);

            x = Tg.csmMag(ix);
            y = Tg.neuralFR(ix);

            good = isfinite(x) & isfinite(y);
            if sum(good) < minTrialsPerUnit_line
                continue
            end

            xg = x(good);
            yg = y(good);

            if numel(unique(xg)) < 2
                continue
            end

            % per-unit Pearson r/p (for your sanity check + coloring)
            [rUnit, pUnit] = corr(xg, yg, 'type','Pearson', 'rows','complete');
            isSig = isfinite(pUnit) && (pUnit < alpha_unit);

            % per-unit linear fit
            P = polyfit(xg, yg, 1);      % [slope intercept]
            yHat = polyval(P, xGrid);

            % plot spaghetti
            if isSig
                plot(xGrid, yHat, '-', 'Color',col_sig, 'LineWidth',lw_sig);
                nSig = nSig + 1;
            else
                plot(xGrid, yHat, '-', 'Color',col_nonsig, 'LineWidth',lw_nonsig);
            end

            % store (ALL units that pass thresholds)
            rUnits(end+1,1) = rUnit;         %#ok<SAGROW>
            pUnits(end+1,1) = pUnit;         %#ok<SAGROW>
            slopes(end+1,1) = P(1);          %#ok<SAGROW>
            ints(end+1,1)   = P(2);          %#ok<SAGROW>
            YhatUnits(:,end+1) = yHat;       %#ok<SAGROW>

            nLines = nLines + 1;
        end

        if isempty(YhatUnits)
            text(0.1,0.5,'No units passed thresholds','Units','normalized');
            axis off
            continue
        end

        % --------- MEAN-OF-UNITS regression line (unit-weighted) ----------
        yMean = mean(YhatUnits, 2, 'omitnan');

        % --------- UNIT-bootstrap CI band (resample UNITS) ----------
        nU = size(YhatUnits,2);
        Yboot = nan(nGrid, nBoot);
        for b = 1:nBoot
            ii = randi(nU, nU, 1);
            Yboot(:,b) = mean(YhatUnits(:,ii), 2, 'omitnan');
        end

        lo = prctile(Yboot', 100*(alphaCI/2))';
        hi = prctile(Yboot', 100*(1-alphaCI/2))';

        okBand = isfinite(lo) & isfinite(hi);
        if any(okBand)
            xx = [xGrid(okBand); flipud(xGrid(okBand))];
            yy = [lo(okBand);    flipud(hi(okBand))];
            fill(xx, yy, 'k', 'FaceAlpha',0.12, 'EdgeColor','none');
        end

        % mean line on top
        plot(xGrid, yMean, 'k-', 'LineWidth',lw_mean);

        % --------- Summary stats that match what you're comparing against ----------
        % Fisher mean of per-unit r (comparable to your per-unit r plot)
        rOk = rUnits(isfinite(rUnits) & abs(rUnits) < 1);
        if ~isempty(rOk)
            rMeanFisher = tanh(mean(atanh(rOk)));
        else
            rMeanFisher = NaN;
        end

        % mean slope across units + bootstrap p(slope != 0)
        mSlope = mean(slopes,'omitnan');
        slopeBoot = nan(nBoot,1);
        for b = 1:nBoot
            ii = randi(numel(slopes), numel(slopes), 1);
            slopeBoot(b) = mean(slopes(ii),'omitnan');
        end
        pSlope = 2 * min(mean(slopeBoot >= 0), mean(slopeBoot <= 0));

        % labels
        xlabel('CSM magnitude');
        ylabel('Neural response (Hz)');
        title(sprintf('%s: nUnits=%d (sig=%d) | nTrials=%d', panelTitle, nLines, nSig, height(Tg)));

        % annotation (NO pooled r)
        txt = sprintf(['Mean per-unit r (Fisher) = %.3f\n' ...
                       'Median per-unit r = %.3f\n' ...
                       'Mean slope = %.3g (boot p=%.3g)'], ...
                       rMeanFisher, median(rUnits,'omitnan'), mSlope, pSlope);

        text(0.03, 0.97, txt, 'Units','normalized', 'VerticalAlignment','top', ...
            'FontSize',10, 'BackgroundColor','w', 'Margin',4);

        box off

        % store if desired
        Scale.spaghettiAvg.(char(panelTitle)).xGrid = xGrid; %#ok<STRNU>
        Scale.spaghettiAvg.(char(panelTitle)).yMean = yMean;
        Scale.spaghettiAvg.(char(panelTitle)).rUnits = rUnits;
        Scale.spaghettiAvg.(char(panelTitle)).pUnits = pUnits;
        Scale.spaghettiAvg.(char(panelTitle)).slopes = slopes;
        Scale.spaghettiAvg.(char(panelTitle)).ci_lo = lo;
        Scale.spaghettiAvg.(char(panelTitle)).ci_hi = hi;

    end

    Scale.spaghettiAvg.minTrialsPerUnit_line = minTrialsPerUnit_line;
    Scale.spaghettiAvg.nBoot = nBoot;

else
    fprintf('\n[SPAGHETTI+MEAN-OF-UNITS] Need Scale.table first.\n');
end

%% ===================== NORMALIZED SPAGHETTI + MEAN-OF-UNITS (pos vs neg) =====================
% Same as previous plot but WITHIN-UNIT Z-SCORE normalization
% → removes baseline & gain differences across neurons
%
% Normalization:
%   xz = (x - mean(x)) / std(x)
%   yz = (y - mean(y)) / std(y)
%
% Average line = mean of per-unit regressions in normalized space

if do_scaling && isfield(Scale,'table') && ~isempty(Scale.table)

    T = Scale.table;

    % ensure category string
    if iscell(T.category)
        T.category = string(T.category);
    elseif ischar(T.category)
        T.category = string(cellstr(T.category));
    elseif iscategorical(T.category)
        T.category = string(T.category);
    else
        T.category = string(T.category);
    end

    isPos = (T.category == "negLag_posCorr");
    isNeg = (T.category == "negLag_negCorr");

    minTrialsPerUnit_line = max(5, min_trials_for_unit_corr);
    nGrid   = 120;
    nBoot   = 2000;
    alphaCI = 0.05;
    alpha_unit = 0.05;

    lw_nonsig = 1;
    lw_sig    = 1.6;
    lw_mean   = 3;

    col_nonsig = [0.75 0.75 0.75];
    col_sig    = [0.10 0.10 0.10];

    figure('Color','w','Position',[200 200 980 420]);

    for panel = 1:2

        if panel == 1
            Tg = T(isPos,:);
            panelTitle = "posCorr units";
        else
            Tg = T(isNeg,:);
            panelTitle = "negCorr units";
        end

        Tg = Tg(isfinite(Tg.neuralFR) & isfinite(Tg.csmMag),:);

        subplot(1,2,panel); hold on;

        if isempty(Tg)
            text(0.1,0.5,'Not enough data','Units','normalized');
            axis off
            continue
        end

        % fixed normalized grid (same for all units)
        xGrid = linspace(-3,3,nGrid)';

        uList = unique(Tg.unitIdx);

        rUnits = [];
        slopes = [];
        YhatUnits = [];

        for iU = 1:numel(uList)

            u = uList(iU);
            ix = Tg.unitIdx == u;

            x = Tg.csmMag(ix);
            y = Tg.neuralFR(ix);

            good = isfinite(x) & isfinite(y);
            if sum(good) < minTrialsPerUnit_line
                continue
            end

            xg = x(good);
            yg = y(good);

            if std(xg)==0 || std(yg)==0
                continue
            end

            % ===== NORMALIZATION =====
            xz = (xg - mean(xg)) ./ std(xg);
            yz = (yg - mean(yg)) ./ std(yg);

            % per-unit r
            rU = corr(xz, yz, 'type','Pearson','rows','complete');

            % per-unit regression
            P = polyfit(xz, yz, 1);
            yHat = polyval(P, xGrid);

            % plot spaghetti
            plot(xGrid, yHat, '-', 'Color',col_nonsig, 'LineWidth',lw_nonsig);

            rUnits(end+1,1) = rU; %#ok<SAGROW>
            slopes(end+1,1) = P(1); %#ok<SAGROW>
            YhatUnits(:,end+1) = yHat; %#ok<SAGROW>

        end

        if isempty(YhatUnits)
            text(0.1,0.5,'No units passed thresholds','Units','normalized');
            axis off
            continue
        end

        % ===== MEAN LINE =====
        yMean = mean(YhatUnits,2,'omitnan');

        % ===== UNIT BOOTSTRAP CI =====
        nU = size(YhatUnits,2);
        Yboot = nan(nGrid,nBoot);
        for b = 1:nBoot
            ii = randi(nU,nU,1);
            Yboot(:,b) = mean(YhatUnits(:,ii),2,'omitnan');
        end
        lo = prctile(Yboot',100*(alphaCI/2))';
        hi = prctile(Yboot',100*(1-alphaCI/2))';

        ok = isfinite(lo) & isfinite(hi);
        if any(ok)
            fill([xGrid(ok); flipud(xGrid(ok))], ...
                 [lo(ok); flipud(hi(ok))], ...
                 'k','FaceAlpha',0.12,'EdgeColor','none');
        end

        plot(xGrid,yMean,'k-','LineWidth',lw_mean);

        % ===== SUMMARY STATS =====
        rMeanFisher = tanh(mean(atanh(rUnits)));

        xlabel('CSM magnitude (z within unit)');
        ylabel('Firing rate (z within unit)');
        title(sprintf('%s (normalized)',panelTitle));

        txt = sprintf('Mean per-unit r (Fisher) = %.3f\nMean slope = %.3f', ...
            rMeanFisher, mean(slopes));

        text(0.03,0.97,txt,'Units','normalized','VerticalAlignment','top', ...
            'BackgroundColor','w','Margin',4);

        box off

    end

else
    fprintf('\n[NORMALIZED SPAGHETTI] Need Scale.table first.\n');
end
%% ===================== CORRELATION PLOTS: CSM MAG vs FIRING RATE (pos vs neg) =====================
% Pooled trial-level scatter (with transparency) + pooled Pearson r/p.
% Requires Scale.table columns: unitIdx, neuralFR, csmMag, category

if do_scaling && isfield(Scale,'table') && ~isempty(Scale.table)

    T = Scale.table;

    % ensure category is string
    if iscell(T.category)
        T.category = string(T.category);
    elseif ischar(T.category)
        T.category = string(cellstr(T.category));
    elseif iscategorical(T.category)
        T.category = string(T.category);
    else
        T.category = string(T.category);
    end

    grpKeys   = ["negLag_posCorr","negLag_negCorr"];
    grpTitles = ["posCorr units","negCorr units"];

    figure('Color','w','Position',[200 200 980 420]);

    for g = 1:2

        Tg = T(T.category == grpKeys(g),:);
        Tg = Tg(isfinite(Tg.neuralFR) & isfinite(Tg.csmMag),:);

        subplot(1,2,g); hold on;

        if isempty(Tg)
            text(0.1,0.5,'Not enough data','Units','normalized'); axis off; continue
        end

        x = Tg.csmMag;
        y = Tg.neuralFR;

        % pooled Pearson r
        r = NaN; p = NaN;
        if numel(x) >= 10 && numel(unique(x)) > 1
            [r,p] = corr(x, y, 'type','Pearson', 'rows','complete');
        end

        % scatter (transparent)
        scatter(x, y, 12, 'k', 'filled', 'MarkerFaceAlpha',0.15, 'MarkerEdgeAlpha',0.15);

        % pooled regression line (just for visualization)
        P = polyfit(x, y, 1);
        xGrid = linspace(min(x), max(x), 120);
        yFit  = polyval(P, xGrid);
        plot(xGrid, yFit, 'k-', 'LineWidth',2.5);

        xlabel('CSM magnitude');
        ylabel('Firing rate (Hz)');
        title(sprintf('%s (n=%d)', grpTitles(g), numel(x)));

        txt = sprintf('Pooled Pearson r = %.3f\np = %.3g', r, p);
        text(0.03,0.97,txt,'Units','normalized','VerticalAlignment','top', ...
            'BackgroundColor','w','Margin',4);

        box off
    end

else
    fprintf('\n[CORR PLOTS] Need Scale.table first.\n');
end
%% ===================== WITHIN-UNIT NORMALIZED CORRELATION PLOTS (pos vs neg) =====================
% z-score x and y within each unit, then pool and plot.

if do_scaling && isfield(Scale,'table') && ~isempty(Scale.table)

    T = Scale.table;

    if iscell(T.category)
        T.category = string(T.category);
    elseif ischar(T.category)
        T.category = string(cellstr(T.category));
    elseif iscategorical(T.category)
        T.category = string(T.category);
    else
        T.category = string(T.category);
    end

    grpKeys   = ["negLag_posCorr","negLag_negCorr"];
    grpTitles = ["posCorr units","negCorr units"];
    minTrialsPerUnit = max(5, min_trials_for_unit_corr);

    figure('Color','w','Position',[200 200 980 420]);

    for g = 1:2

        Tg = T(T.category == grpKeys(g),:);
        Tg = Tg(isfinite(Tg.neuralFR) & isfinite(Tg.csmMag),:);

        subplot(1,2,g); hold on;

        if isempty(Tg)
            text(0.1,0.5,'Not enough data','Units','normalized'); axis off; continue
        end

        uList = unique(Tg.unitIdx);
        xAll = [];
        yAll = [];

        for iU = 1:numel(uList)
            u = uList(iU);
            ix = Tg.unitIdx == u;

            x = Tg.csmMag(ix);
            y = Tg.neuralFR(ix);

            good = isfinite(x) & isfinite(y);
            if sum(good) < minTrialsPerUnit
                continue
            end

            x = x(good); y = y(good);
            if std(x)==0 || std(y)==0
                continue
            end

            xz = (x - mean(x)) ./ std(x);
            yz = (y - mean(y)) ./ std(y);

            xAll = [xAll; xz];
            yAll = [yAll; yz];
        end

        r = NaN; p = NaN;
        if numel(xAll) >= 10 && numel(unique(xAll)) > 1
            [r,p] = corr(xAll, yAll, 'type','Pearson', 'rows','complete');
        end

        scatter(xAll, yAll, 10, 'k', 'filled', 'MarkerFaceAlpha',0.18, 'MarkerEdgeAlpha',0.18);

        P = polyfit(xAll, yAll, 1);
        xGrid = linspace(min(xAll), max(xAll), 120);
        yFit  = polyval(P, xGrid);
        plot(xGrid, yFit, 'k-', 'LineWidth',2.5);

        xlabel('CSM magnitude (z within unit)');
        ylabel('Firing rate (z within unit)');
        title(sprintf('%s (pooled within-unit z)', grpTitles(g)));

        txt = sprintf('Pooled Pearson r = %.3f\np = %.3g', r, p);
        text(0.03,0.97,txt,'Units','normalized','VerticalAlignment','top', ...
            'BackgroundColor','w','Margin',4);

        box off
    end

else
    fprintf('\n[WITHIN-UNIT CORR PLOTS] Need Scale.table first.\n');
end
%% ===================== VIS: PER-UNIT r DOTS + MEAN OF INDIVIDUAL REGRESSIONS (pos vs neg) =====================
% TOP: each unit's Pearson r (dot), colored by p<0.05 + Fisher-mean r line
% BOTTOM: mean of individual unit regression lines + CI via bootstrap over UNITS
%
% Requires Scale.table columns: unitIdx, neuralFR, csmMag, category

if do_scaling && isfield(Scale,'table') && ~isempty(Scale.table)

    T = Scale.table;

    % --- ensure category is string ---
    if iscell(T.category)
        T.category = string(T.category);
    elseif ischar(T.category)
        T.category = string(cellstr(T.category));
    elseif iscategorical(T.category)
        T.category = string(T.category);
    else
        T.category = string(T.category);
    end

    % groups
    grpKeys   = ["negLag_posCorr","negLag_negCorr"];   % SAFE struct field names
    grpTitles = ["posCorr units","negCorr units"];     % for plot titles only

    % settings
    minTrialsPerUnit_line = max(5, min_trials_for_unit_corr);
    nGrid   = 120;
    nBoot   = 2000;     % unit-bootstrap
    alphaCI = 0.05;     % 95% CI
    alpha_unit = 0.05;

    % styling
    lw_nonsig = 1;
    lw_mean   = 3;

    col_nonsig = [0.75 0.75 0.75];
    col_sig    = [0.10 0.10 0.10];

    figure('Color','w','Position',[200 120 980 720]);

    % init storage
    if ~isfield(Scale,'unitSummary') || ~isstruct(Scale.unitSummary)
        Scale.unitSummary = struct();
    end

    for g = 1:2

        grpKey   = grpKeys(g);
        grpTitle = grpTitles(g);

        Tg = T(T.category == grpKey,:);
        Tg = Tg(isfinite(Tg.neuralFR) & isfinite(Tg.csmMag),:);

        % ------------------ common x-grid ------------------
        if isempty(Tg)
            subplot(2,2,g); axis off; text(0.1,0.5,'Not enough data','Units','normalized');
            subplot(2,2,g+2); axis off; text(0.1,0.5,'Not enough data','Units','normalized');
            continue
        end

        xMin = min(Tg.csmMag,[],'omitnan');
        xMax = max(Tg.csmMag,[],'omitnan');
        if ~isfinite(xMin) || ~isfinite(xMax) || xMax==xMin
            subplot(2,2,g); axis off; text(0.1,0.5,'CSM magnitude has no range','Units','normalized');
            subplot(2,2,g+2); axis off; text(0.1,0.5,'CSM magnitude has no range','Units','normalized');
            continue
        end
        xGrid = linspace(xMin, xMax, nGrid)';

        % ------------------ per-unit computations ------------------
        uList = unique(Tg.unitIdx);

        rUnits    = [];
        pUnits    = [];
        slopes    = [];
        intercepts= [];
        YhatUnits = [];   % nGrid x nUnitsKept

        for iU = 1:numel(uList)

            u  = uList(iU);
            ix = (Tg.unitIdx == u);

            x = Tg.csmMag(ix);
            y = Tg.neuralFR(ix);

            good = isfinite(x) & isfinite(y);
            if sum(good) < minTrialsPerUnit_line
                continue
            end

            xg = x(good);
            yg = y(good);

            if numel(unique(xg)) < 2
                continue
            end

            [rU, pU] = corr(xg, yg, 'type','Pearson', 'rows','complete');

            P = polyfit(xg, yg, 1);
            yHat = polyval(P, xGrid);

            rUnits(end+1,1)      = rU;    %#ok<SAGROW>
            pUnits(end+1,1)      = pU;    %#ok<SAGROW>
            slopes(end+1,1)      = P(1);  %#ok<SAGROW>
            intercepts(end+1,1)  = P(2);  %#ok<SAGROW>
            YhatUnits(:,end+1)   = yHat;  %#ok<SAGROW>
        end

        nUnits = numel(rUnits);
        if nUnits < 1
            subplot(2,2,g); axis off; text(0.1,0.5,'No units passed thresholds','Units','normalized');
            subplot(2,2,g+2); axis off; text(0.1,0.5,'No units passed thresholds','Units','normalized');
            continue
        end

        isSig = isfinite(pUnits) & (pUnits < alpha_unit);

        % Fisher-mean r
        rOk = rUnits(isfinite(rUnits) & abs(rUnits) < 1);
        if ~isempty(rOk)
            rMeanFisher = tanh(mean(atanh(rOk)));
        else
            rMeanFisher = NaN;
        end

        % ------------------ TOP: DOT PLOT OF PER-UNIT r ------------------
        subplot(2,2,g); hold on;

        [rS, ord] = sort(rUnits, 'ascend');
        sigS = isSig(ord);

        xIdx = (1:nUnits)';
        jit  = (rand(nUnits,1)-0.5)*0.15;

        scatter(xIdx(~sigS)+jit(~sigS), rS(~sigS), 28, 'filled', ...
            'MarkerFaceColor',col_nonsig, 'MarkerFaceAlpha',0.95);
        scatter(xIdx(sigS)+jit(sigS), rS(sigS), 34, 'filled', ...
            'MarkerFaceColor',col_sig, 'MarkerFaceAlpha',1);

        yline(0,'k-','LineWidth',1);
        yline(rMeanFisher,'k--','LineWidth',1.8);

        xlim([0.5 nUnits+0.5]);
        ylabel('Per-unit Pearson r');
        xlabel('Units (sorted by r)');
        title(sprintf('%s: per-unit r (n=%d; sig=%d)', grpTitle, nUnits, sum(isSig)));

        txtTop = sprintf('Mean r (Fisher) = %.3f\nMedian r = %.3f', rMeanFisher, median(rUnits,'omitnan'));
        text(0.02,0.97,txtTop,'Units','normalized','VerticalAlignment','top', ...
            'FontSize',10,'BackgroundColor','w','Margin',4);

        box off

        % ------------------ BOTTOM: MEAN OF INDIVIDUAL REGRESSIONS ------------------
        subplot(2,2,g+2); hold on;

        yMean = mean(YhatUnits, 2, 'omitnan');

        % bootstrap CI over UNITS
        nU = size(YhatUnits,2);
        Yboot = nan(nGrid, nBoot);
        for b = 1:nBoot
            ii = randi(nU, nU, 1);
            Yboot(:,b) = mean(YhatUnits(:,ii), 2, 'omitnan');
        end
        lo = prctile(Yboot', 100*(alphaCI/2))';
        hi = prctile(Yboot', 100*(1-alphaCI/2))';

        okBand = isfinite(lo) & isfinite(hi);
        if any(okBand)
            fill([xGrid(okBand); flipud(xGrid(okBand))], ...
                 [lo(okBand);    flipud(hi(okBand))], ...
                 'k', 'FaceAlpha',0.12, 'EdgeColor','none');
        end

        % optional faint unit lines
        for i = 1:nU
            plot(xGrid, YhatUnits(:,i), '-', 'Color',[0.87 0.87 0.87], 'LineWidth',lw_nonsig);
        end

        plot(xGrid, yMean, 'k-', 'LineWidth',lw_mean);

        xlabel('CSM magnitude');
        ylabel('Neural response (Hz)');
        title(sprintf('%s: mean of unit regressions (unit-bootstrap CI)', grpTitle));
        box off

        % ------------------ STORE (SAFE FIELD NAME) ------------------
        Scale.unitSummary.(grpKey).rUnits      = rUnits;
        Scale.unitSummary.(grpKey).pUnits      = pUnits;
        Scale.unitSummary.(grpKey).slopes      = slopes;
        Scale.unitSummary.(grpKey).intercepts  = intercepts;
        Scale.unitSummary.(grpKey).xGrid       = xGrid;
        Scale.unitSummary.(grpKey).YhatUnits   = YhatUnits;
        Scale.unitSummary.(grpKey).rMeanFisher = rMeanFisher;
        Scale.unitSummary.(grpKey).nUnits      = nUnits;
        Scale.unitSummary.(grpKey).nSigUnits   = sum(isSig);

    end

    Scale.unitSummary.minTrialsPerUnit_line = minTrialsPerUnit_line;
    Scale.unitSummary.nBoot = nBoot;

else
    fprintf('\n[UNIT r + MEAN REGRESSIONS] Need Scale.table first.\n');
end

%% ===================== VIS: PER-UNIT r DOTS + MEAN OF INDIVIDUAL REGRESSIONS + MEAN-LINE STATS =====================
% TOP: each unit's Pearson r (dot), colored by p<0.05 + Fisher-mean r line
% BOTTOM: mean of individual unit regression lines + CI via bootstrap over UNITS
%         + stats for the mean line computed AGAINST DATA with equal-unit weighting:
%             - r_w: weighted corr(y, yhat_mean) with weights=1/nTrials(unit)
%             - R2_w: weighted R^2
%             - p_slope: unit-bootstrap p-value that mean slope != 0 (plus ttest backup)
%
% Requires Scale.table columns: unitIdx, neuralFR, csmMag, category

if do_scaling && isfield(Scale,'table') && ~isempty(Scale.table)

    T = Scale.table;

    % --- ensure category is string ---
    if iscell(T.category)
        T.category = string(T.category);
    elseif ischar(T.category)
        T.category = string(cellstr(T.category));
    elseif iscategorical(T.category)
        T.category = string(T.category);
    else
        T.category = string(T.category);
    end

    % groups (SAFE struct field names)
    grpKeys   = ["negLag_posCorr","negLag_negCorr"];
    grpTitles = ["posCorr units","negCorr units"];

    % settings
    minTrialsPerUnit_line = max(5, min_trials_for_unit_corr);
    nGrid   = 120;
    nBoot   = 2000;     % unit-bootstrap
    alphaCI = 0.05;     % 95% CI
    alpha_unit = 0.05;

    % styling
    lw_nonsig = 1;
    lw_mean   = 3;

    col_nonsig = [0.75 0.75 0.75];
    col_sig    = [0.10 0.10 0.10];

    figure('Color','w','Position',[200 120 980 720]);

    % init storage
    if ~isfield(Scale,'unitSummary') || ~isstruct(Scale.unitSummary)
        Scale.unitSummary = struct();
    end

    for g = 1:2

        grpKey   = grpKeys(g);
        grpTitle = grpTitles(g);

        Tg = T(T.category == grpKey,:);
        Tg = Tg(isfinite(Tg.neuralFR) & isfinite(Tg.csmMag),:);

        % ------------------ common x-grid ------------------
        if isempty(Tg)
            subplot(2,2,g); axis off; text(0.1,0.5,'Not enough data','Units','normalized');
            subplot(2,2,g+2); axis off; text(0.1,0.5,'Not enough data','Units','normalized');
            continue
        end

        xMin = min(Tg.csmMag,[],'omitnan');
        xMax = max(Tg.csmMag,[],'omitnan');
        if ~isfinite(xMin) || ~isfinite(xMax) || xMax==xMin
            subplot(2,2,g); axis off; text(0.1,0.5,'CSM magnitude has no range','Units','normalized');
            subplot(2,2,g+2); axis off; text(0.1,0.5,'CSM magnitude has no range','Units','normalized');
            continue
        end
        xGrid = linspace(xMin, xMax, nGrid)';

        % ------------------ per-unit computations ------------------
        uList = unique(Tg.unitIdx);

        rUnits     = [];
        pUnits     = [];
        slopes     = [];
        intercepts = [];
        nTrialsU   = [];

        YhatUnits  = [];   % nGrid x nUnitsKept

        % for mean-line stats (unit-weighted against real data)
        xAll = [];
        yAll = [];
        wAll = [];   % weights per trial = 1/nTrials(unit)

        for iU = 1:numel(uList)

            u  = uList(iU);
            ix = (Tg.unitIdx == u);

            x = Tg.csmMag(ix);
            y = Tg.neuralFR(ix);

            good = isfinite(x) & isfinite(y);
            if sum(good) < minTrialsPerUnit_line
                continue
            end

            xg = x(good);
            yg = y(good);

            if numel(unique(xg)) < 2
                continue
            end

            % per-unit Pearson r (for DOT plot)
            [rU, pU] = corr(xg, yg, 'type','Pearson', 'rows','complete');

            % per-unit regression line (for MEAN line)
            P = polyfit(xg, yg, 1);    % [slope intercept]
            yHat = polyval(P, xGrid);

            % store per-unit summaries
            rUnits(end+1,1)      = rU;       %#ok<SAGROW>
            pUnits(end+1,1)      = pU;       %#ok<SAGROW>
            slopes(end+1,1)      = P(1);     %#ok<SAGROW>
            intercepts(end+1,1)  = P(2);     %#ok<SAGROW>
            nTrialsU(end+1,1)    = numel(xg); %#ok<SAGROW>
            YhatUnits(:,end+1)   = yHat;     %#ok<SAGROW>

            % store trial-level data with equal-unit weights
            w = ones(numel(xg),1) ./ numel(xg);
            xAll = [xAll; xg(:)];            %#ok<SAGROW>
            yAll = [yAll; yg(:)];            %#ok<SAGROW>
            wAll = [wAll; w(:)];             %#ok<SAGROW>
        end

        nUnits = numel(rUnits);
        if nUnits < 1
            subplot(2,2,g); axis off; text(0.1,0.5,'No units passed thresholds','Units','normalized');
            subplot(2,2,g+2); axis off; text(0.1,0.5,'No units passed thresholds','Units','normalized');
            continue
        end

        isSig = isfinite(pUnits) & (pUnits < alpha_unit);

        % Fisher-mean r across units
        rOk = rUnits(isfinite(rUnits) & abs(rUnits) < 1);
        if ~isempty(rOk)
            rMeanFisher = tanh(mean(atanh(rOk)));
        else
            rMeanFisher = NaN;
        end

        % ------------------ TOP: DOT PLOT OF PER-UNIT r ------------------
        subplot(2,2,g); hold on;

        [rS, ord] = sort(rUnits, 'ascend');
        sigS = isSig(ord);

        xIdx = (1:nUnits)';
        jit  = (rand(nUnits,1)-0.5)*0.15;

        scatter(xIdx(~sigS)+jit(~sigS), rS(~sigS), 28, 'filled', ...
            'MarkerFaceColor',col_nonsig, 'MarkerFaceAlpha',0.95);
        scatter(xIdx(sigS)+jit(sigS), rS(sigS), 34, 'filled', ...
            'MarkerFaceColor',col_sig, 'MarkerFaceAlpha',1);

        yline(0,'k-','LineWidth',1);
        yline(rMeanFisher,'k--','LineWidth',1.8);

        xlim([0.5 nUnits+0.5]);
        ylabel('Per-unit Pearson r');
        xlabel('Units (sorted by r)');
        title(sprintf('%s: per-unit r (n=%d; sig=%d)', grpTitle, nUnits, sum(isSig)));

        txtTop = sprintf('Mean r (Fisher) = %.3f\nMedian r = %.3f', ...
            rMeanFisher, median(rUnits,'omitnan'));
        text(0.02,0.97,txtTop,'Units','normalized','VerticalAlignment','top', ...
            'FontSize',10,'BackgroundColor','w','Margin',4);

        box off

        % ------------------ BOTTOM: MEAN OF INDIVIDUAL REGRESSIONS + STATS ------------------
        subplot(2,2,g+2); hold on;

        % mean line across unit lines (definition of the plotted line)
        yMean = mean(YhatUnits, 2, 'omitnan');

        % bootstrap CI over UNITS for the mean line
        nU = size(YhatUnits,2);
        Yboot = nan(nGrid, nBoot);
        for b = 1:nBoot
            ii = randi(nU, nU, 1);
            Yboot(:,b) = mean(YhatUnits(:,ii), 2, 'omitnan');
        end
        lo = prctile(Yboot', 100*(alphaCI/2))';
        hi = prctile(Yboot', 100*(1-alphaCI/2))';

        okBand = isfinite(lo) & isfinite(hi);
        if any(okBand)
            fill([xGrid(okBand); flipud(xGrid(okBand))], ...
                 [lo(okBand);    flipud(hi(okBand))], ...
                 'k', 'FaceAlpha',0.12, 'EdgeColor','none');
        end

        % optional faint unit lines
        for i = 1:nU
            plot(xGrid, YhatUnits(:,i), '-', 'Color',[0.87 0.87 0.87], 'LineWidth',lw_nonsig);
        end

        % mean line on top
        plot(xGrid, yMean, 'k-', 'LineWidth',lw_mean);

        % ===== MEAN-LINE STATS computed AGAINST DATA (unit-weighted) =====
        % mean slope/intercept (equivalent to averaging regressions)
        mSlope = mean(slopes,'omitnan');
        mInt   = mean(intercepts,'omitnan');

        % predictions for each trial using the mean regression
        yHatAll_mean = mInt + mSlope .* xAll;

        % weighted correlation r_w
        r_w = weighted_corr(xAll, yAll, wAll);  % correlation between x and y (unit-weighted)
        % If you want correlation between y and yHat: use next line instead:
        r_y_yhat_w = weighted_corr(yAll, yHatAll_mean, wAll);

        % weighted R^2 of mean model: y ~ (mInt + mSlope*x)
        R2_w = weighted_R2(yAll, yHatAll_mean, wAll);

        % p-value for mean slope via unit-bootstrap (mean slope != 0)
        slopeBoot = nan(nBoot,1);
        for b = 1:nBoot
            ii = randi(nUnits, nUnits, 1);
            slopeBoot(b) = mean(slopes(ii),'omitnan');
        end
        p_slope_boot = 2 * min(mean(slopeBoot >= 0), mean(slopeBoot <= 0));

        % backup: one-sample t-test on slopes
        p_slope_t = NaN;
        try
            [~, p_slope_t] = ttest(slopes, 0);
        catch
            p_slope_t = NaN;
        end

        xlabel('CSM magnitude');
        ylabel('Neural response (Hz)');
        title(sprintf('%s: mean of unit regressions (unit-bootstrap CI)', grpTitle));
        box off

        txtBottom = sprintf(['Mean model: y = %.3f + %.3f*x\n' ...
                             'unit-wtd corr(x,y): r=%.3f\n' ...
                             'unit-wtd corr(y,yhat): r=%.3f\n' ...
                             'unit-wtd R^2(y~yhat)=%.3f\n' ...
                             'p(slope!=0): boot=%.3g | ttest=%.3g'], ...
                             mInt, mSlope, r_w, r_y_yhat_w, R2_w, p_slope_boot, p_slope_t);

        text(0.03,0.97,txtBottom,'Units','normalized','VerticalAlignment','top', ...
            'FontSize',9,'BackgroundColor','w','Margin',4);

        % ------------------ STORE (SAFE FIELD NAME) ------------------
        Scale.unitSummary.(grpKey).rUnits        = rUnits;
        Scale.unitSummary.(grpKey).pUnits        = pUnits;
        Scale.unitSummary.(grpKey).slopes        = slopes;
        Scale.unitSummary.(grpKey).intercepts    = intercepts;
        Scale.unitSummary.(grpKey).nTrialsU      = nTrialsU;
        Scale.unitSummary.(grpKey).xGrid         = xGrid;
        Scale.unitSummary.(grpKey).YhatUnits     = YhatUnits;
        Scale.unitSummary.(grpKey).yMean         = yMean;
        Scale.unitSummary.(grpKey).ci_lo         = lo;
        Scale.unitSummary.(grpKey).ci_hi         = hi;
        Scale.unitSummary.(grpKey).rMeanFisher   = rMeanFisher;
        Scale.unitSummary.(grpKey).nUnits        = nUnits;
        Scale.unitSummary.(grpKey).nSigUnits     = sum(isSig);

        Scale.unitSummary.(grpKey).meanSlope     = mSlope;
        Scale.unitSummary.(grpKey).meanIntercept = mInt;
        Scale.unitSummary.(grpKey).r_w_xy        = r_w;
        Scale.unitSummary.(grpKey).r_w_y_yhat    = r_y_yhat_w;
        Scale.unitSummary.(grpKey).R2_w          = R2_w;
        Scale.unitSummary.(grpKey).p_slope_boot  = p_slope_boot;
        Scale.unitSummary.(grpKey).p_slope_ttest = p_slope_t;

    end

    Scale.unitSummary.minTrialsPerUnit_line = minTrialsPerUnit_line;
    Scale.unitSummary.nBoot = nBoot;

else
    fprintf('\n[UNIT r + MEAN REGRESSIONS + STATS] Need Scale.table first.\n');
end

%% ===================== LOCAL HELPER FUNCTIONS (paste below in same script) =====================
function r = weighted_corr(a, b, w)
    a = a(:); b = b(:); w = w(:);
    good = isfinite(a) & isfinite(b) & isfinite(w) & (w>0);
    a = a(good); b = b(good); w = w(good);
    if numel(a) < 3
        r = NaN; return
    end
    w = w ./ sum(w);
    ma = sum(w .* a);
    mb = sum(w .* b);
    ca = a - ma;
    cb = b - mb;
    covab = sum(w .* ca .* cb);
    vara  = sum(w .* ca.^2);
    varb  = sum(w .* cb.^2);
    if vara<=0 || varb<=0
        r = NaN; return
    end
    r = covab ./ sqrt(vara*varb);
end

function R2 = weighted_R2(y, yhat, w)
    y = y(:); yhat = yhat(:); w = w(:);
    good = isfinite(y) & isfinite(yhat) & isfinite(w) & (w>0);
    y = y(good); yhat = yhat(good); w = w(good);
    if numel(y) < 3
        R2 = NaN; return
    end
    w = w ./ sum(w);
    ybar = sum(w .* y);
    SSres = sum(w .* (y - yhat).^2);
    SStot = sum(w .* (y - ybar).^2);
    if SStot<=0
        R2 = NaN; return
    end
    R2 = 1 - SSres./SStot;
end
%% ===================== WITHIN-UNIT NORMALIZED ANALYSIS =====================
% Outputs:
%   - Per-unit correlations (r, p)
%   - Mean of individual regression lines
%   - Bootstrap CI over units
%   - Histogram of per-unit r
%   - Results stored in Scale.unitCorr

if do_scaling && isfield(Scale,'table') && ~isempty(Scale.table)

    T = Scale.table;

    if ~isstring(T.category)
        T.category = string(T.category);
    end

    groups = ["negLag_posCorr","negLag_negCorr"];
    titles = ["posCorr units","negCorr units"];

    Scale.unitCorr = struct();

    figure('Color','w','Position',[200 200 1000 420]);

    for g = 1:2

        Tg = T(T.category == groups(g),:);
        Tg = Tg(isfinite(Tg.neuralFR) & isfinite(Tg.csmMag),:);

        subplot(1,2,g); hold on

        if isempty(Tg)
            text(0.1,0.5,'No data','Units','normalized');
            axis off; continue
        end

        uList = unique(Tg.unitIdx);

        x_all = [];
        y_all = [];

        slopes = [];
        intercepts = [];
        r_units = [];
        p_units = [];
        n_trials = [];

        % common grid for averaging regression lines
        xGrid = linspace(-3,3,120)';

        YhatUnits = [];

        for iU = 1:numel(uList)

            u = uList(iU);
            ix = Tg.unitIdx == u;

            x = Tg.csmMag(ix);
            y = Tg.neuralFR(ix);

            good = isfinite(x) & isfinite(y);
            x = x(good);
            y = y(good);

            if numel(x) < 5 || std(x)==0 || std(y)==0
                continue
            end

            % -------- normalize within unit --------
            xz = (x - mean(x)) ./ std(x);
            yz = (y - mean(y)) ./ std(y);

            x_all = [x_all; xz];
            y_all = [y_all; yz];

            % ---- per-unit correlation ----
            [rU, pU] = corr(xz, yz, 'type','Pearson');

            % ---- per-unit regression ----
            P = polyfit(xz, yz, 1);
            yHat = polyval(P, xGrid);

            % store
            slopes(end+1,1) = P(1);
            intercepts(end+1,1) = P(2);
            r_units(end+1,1) = rU;
            p_units(end+1,1) = pU;
            n_trials(end+1,1) = numel(xz);
            YhatUnits(:,end+1) = yHat;

        end

        nUnits = numel(r_units);

        % ===== pooled correlation on normalized data =====
        [r_pool, p_pool] = corr(x_all, y_all, 'type','Pearson');

        % ===== mean of individual regression lines =====
        yMean = mean(YhatUnits, 2, 'omitnan');

        % ===== bootstrap CI over units =====
        nBoot = 2000;
        Yboot = nan(numel(xGrid), nBoot);

        for b = 1:nBoot
            ii = randi(nUnits, nUnits, 1);
            Yboot(:,b) = mean(YhatUnits(:,ii), 2, 'omitnan');
        end

        lo = prctile(Yboot',2.5)';
        hi = prctile(Yboot',97.5)';

        % ===== scatter =====
        scatter(x_all, y_all, 10, 'k', 'filled', 'MarkerFaceAlpha',0.25);

        % CI band
        fill([xGrid; flipud(xGrid)], [lo; flipud(hi)], 'k', ...
            'FaceAlpha',0.12, 'EdgeColor','none');

        % mean regression line
        plot(xGrid, yMean, 'r', 'LineWidth',3);

        xlabel('CSM magnitude (z-scored within unit)');
        ylabel('Neural response (z-scored within unit)');
        title(titles(g));

        % ===== mean r (Fisher transform) =====
        z = atanh(r_units(abs(r_units)<1));
        r_mean = tanh(mean(z,'omitnan'));

        txt = sprintf(['nUnits=%d\n' ...
                       'Mean unit r (Fisher)=%.3f\n' ...
                       'Median r=%.3f\n' ...
                       'Pooled r=%.3f (p=%.3g)'], ...
                       nUnits, r_mean, median(r_units,'omitnan'), r_pool, p_pool);

        text(0.03,0.97,txt,'Units','normalized','VerticalAlignment','top',...
            'BackgroundColor','w','Margin',4);

        box off

        % ===== store results =====
        Scale.unitCorr.(groups(g)).r = r_units;
        Scale.unitCorr.(groups(g)).p = p_units;
        Scale.unitCorr.(groups(g)).slope = slopes;
        Scale.unitCorr.(groups(g)).intercept = intercepts;
        Scale.unitCorr.(groups(g)).nTrials = n_trials;
        Scale.unitCorr.(groups(g)).r_mean = r_mean;
        Scale.unitCorr.(groups(g)).r_median = median(r_units,'omitnan');

    end

    % ===== Histogram of per-unit correlations =====
    figure('Color','w','Position',[300 300 900 350]);

    for g = 1:2
        subplot(1,2,g); hold on
        r_units = Scale.unitCorr.(groups(g)).r;

        histogram(r_units, 20, 'FaceColor',[0.6 0.6 0.6]);
        xline(mean(r_units,'omitnan'),'r','LineWidth',2);
        xlabel('Per-unit Pearson r');
        ylabel('Count');
        title(titles(g));
        box off
    end

else
    fprintf('[UNIT NORMALIZED ANALYSIS] Need Scale.table first.\n');
end

%% ===================== UNIT-LEVEL r DOT PLOT + MEAN OF UNIT REGRESSIONS =====================
% NO trial-level scatter.
% Panel A (per group): each dot = one unit's Pearson r (computed within-unit on z-scored x/y)
% Panel B (per group): mean of individual unit regression lines + unit-bootstrap CI
%
% Requires Scale.table columns: unitIdx, neuralFR, csmMag, category

if do_scaling && isfield(Scale,'table') && ~isempty(Scale.table)

    T = Scale.table;

    % --- ensure category is string ---
    if iscell(T.category)
        T.category = string(T.category);
    elseif ischar(T.category)
        T.category = string(cellstr(T.category));
    elseif iscategorical(T.category)
        T.category = string(T.category);
    else
        T.category = string(T.category);
    end

    % groups
    grpNames    = ["negLag_posCorr","negLag_negCorr"];
    panelTitles = ["posCorr units","negCorr units"];

    % settings
    minTrialsPerUnit = max(5, min_trials_for_unit_corr);
    alpha_unit = 0.05;

    nGrid   = 120;   % for regression mean line
    nBoot   = 2000;  % unit-bootstrap
    alphaCI = 0.05;

    % plotting
    col_nonsig = [0.65 0.65 0.65];
    col_sig    = [0.10 0.10 0.10];

    figure('Color','w','Position',[180 120 1100 700]);

    Scale.unitLevel = struct();

    for g = 1:2

        Tg = T(T.category == grpNames(g),:);
        Tg = Tg(isfinite(Tg.neuralFR) & isfinite(Tg.csmMag),:);

        if isempty(Tg)
            % dot plot
            subplot(2,2,g); axis off
            text(0.1,0.5,'No data','Units','normalized');

            % mean-reg plot
            subplot(2,2,g+2); axis off
            text(0.1,0.5,'No data','Units','normalized');
            continue
        end

        uList = unique(Tg.unitIdx);

        % collectors
        unitID     = [];
        r_units    = [];
        p_units    = [];
        slopes     = [];
        intercepts = [];

        % fixed x-grid for normalized regressions (since xz, yz are z-scored)
        xGrid = linspace(-3, 3, nGrid)';

        YhatUnits = [];  % nGrid x nUnitsKept

        for iU = 1:numel(uList)

            u = uList(iU);
            ix = (Tg.unitIdx == u);

            x = Tg.csmMag(ix);
            y = Tg.neuralFR(ix);

            good = isfinite(x) & isfinite(y);
            x = x(good);
            y = y(good);

            if numel(x) < minTrialsPerUnit
                continue
            end
            if std(x)==0 || std(y)==0
                continue
            end

            % -------- within-unit normalize (z-score) --------
            xz = (x - mean(x)) ./ std(x);
            yz = (y - mean(y)) ./ std(y);

            % -------- per-unit correlation --------
            [rU, pU] = corr(xz, yz, 'type','Pearson', 'rows','complete');

            % -------- per-unit regression --------
            P = polyfit(xz, yz, 1);     % P(1)=slope, P(2)=intercept
            yHat = polyval(P, xGrid);   % evaluate on common grid

            % store
            unitID(end+1,1)     = u;          %#ok<SAGROW>
            r_units(end+1,1)    = rU;         %#ok<SAGROW>
            p_units(end+1,1)    = pU;         %#ok<SAGROW>
            slopes(end+1,1)     = P(1);       %#ok<SAGROW>
            intercepts(end+1,1) = P(2);       %#ok<SAGROW>
            YhatUnits(:,end+1)  = yHat;       %#ok<SAGROW>
        end

        nUnits = numel(r_units);

        % ------------------ (A) DOT PLOT: per-unit r only ------------------
        subplot(2,2,g); hold on;

        if nUnits < 1
            axis off
            text(0.1,0.5,'No units passed thresholds','Units','normalized');
        else
            % sort by r (nice visualization)
            [rS, ord] = sort(r_units, 'ascend');
            pS = p_units(ord);

            xIdx = (1:nUnits)';
            jit  = (rand(nUnits,1)-0.5)*0.15;

            isSig = isfinite(pS) & (pS < alpha_unit);

            scatter(xIdx(~isSig)+jit(~isSig), rS(~isSig), 28, 'filled', ...
                'MarkerFaceColor',col_nonsig, 'MarkerFaceAlpha',0.9);

            scatter(xIdx(isSig)+jit(isSig), rS(isSig), 34, 'filled', ...
                'MarkerFaceColor',col_sig, 'MarkerFaceAlpha',1);

            yline(0,'k-','LineWidth',1);

            % Fisher mean r
            rOk = r_units(isfinite(r_units) & abs(r_units)<1);
            if ~isempty(rOk)
                rMeanFisher = tanh(mean(atanh(rOk)));
            else
                rMeanFisher = NaN;
            end
            yline(rMeanFisher,'k--','LineWidth',1.6);

            xlabel('Units (sorted by r)');
            ylabel('Per-unit Pearson r (within-unit z-scored)');
            title(sprintf('%s: per-unit r (n=%d; sig=%d)', panelTitles(g), nUnits, sum(isfinite(p_units) & p_units<alpha_unit)));

            txt = sprintf('Mean r (Fisher) = %.3f\nMedian r = %.3f', rMeanFisher, median(r_units,'omitnan'));
            text(0.02,0.98,txt,'Units','normalized','VerticalAlignment','top', ...
                'BackgroundColor','w','Margin',4);

            box off
        end

        % ------------------ (B) MEAN OF INDIVIDUAL REGRESSIONS ------------------
        subplot(2,2,g+2); hold on;

        if size(YhatUnits,2) < 2
            axis off
            text(0.1,0.5,'Too few unit regressions to average','Units','normalized');
        else
            % (optional) show individual unit lines faint
            for i = 1:size(YhatUnits,2)
                plot(xGrid, YhatUnits(:,i), '-', 'Color',[0.85 0.85 0.85], 'LineWidth',1);
            end

            % mean line
            yMean = mean(YhatUnits, 2, 'omitnan');

            % unit-bootstrap CI
            nU = size(YhatUnits,2);
            Yboot = nan(nGrid, nBoot);
            for b = 1:nBoot
                ii = randi(nU, nU, 1);
                Yboot(:,b) = mean(YhatUnits(:,ii), 2, 'omitnan');
            end
            lo = prctile(Yboot', 100*(alphaCI/2))';
            hi = prctile(Yboot', 100*(1-alphaCI/2))';

            ok = isfinite(lo) & isfinite(hi);
            if any(ok)
                fill([xGrid(ok); flipud(xGrid(ok))], [lo(ok); flipud(hi(ok))], 'k', ...
                    'FaceAlpha',0.12, 'EdgeColor','none');
            end

            plot(xGrid, yMean, 'k-', 'LineWidth',3);

            xlabel('CSM magnitude (z-scored within unit)');
            ylabel('Neural response (z-scored within unit)');
            title(sprintf('%s: mean of unit regressions (n=%d)', panelTitles(g), nU));
            yline(0,'k:');
            xline(0,'k:');
            box off
        end

        % store
        Scale.unitLevel.(grpNames(g)).unitID     = unitID;
        Scale.unitLevel.(grpNames(g)).r          = r_units;
        Scale.unitLevel.(grpNames(g)).p          = p_units;
        Scale.unitLevel.(grpNames(g)).slope      = slopes;
        Scale.unitLevel.(grpNames(g)).intercept  = intercepts;
        Scale.unitLevel.(grpNames(g)).xGrid      = xGrid;
        Scale.unitLevel.(grpNames(g)).YhatUnits  = YhatUnits;

    end

else
    fprintf('\n[UNIT r + MEAN REGRESSIONS] Need Scale.table first.\n');
end
%% ===================== WITHIN-UNIT Z-SCORED CORRELATION =====================

if do_scaling && isfield(Scale,'table') && ~isempty(Scale.table)

    T = Scale.table;

    % ensure string category
    if ~isstring(T.category)
        T.category = string(T.category);
    end

    groups = ["negLag_posCorr","negLag_negCorr"];
    titles = ["posCorr units","negCorr units"];

    figure('Color','w','Position',[200 200 1000 420]);

    for g = 1:2

        Tg = T(T.category == groups(g),:);
        Tg = Tg(isfinite(Tg.neuralFR) & isfinite(Tg.csmMag),:);

        subplot(1,2,g); hold on;

        if isempty(Tg)
            text(0.1,0.5,'No data','Units','normalized');
            axis off; continue
        end

        uList = unique(Tg.unitIdx);

        x_all = [];
        y_all = [];
        r_units = [];

        for iU = 1:numel(uList)

            u = uList(iU);
            ix = Tg.unitIdx == u;

            x = Tg.csmMag(ix);
            y = Tg.neuralFR(ix);

            good = isfinite(x) & isfinite(y);
            x = x(good);
            y = y(good);

            if numel(x) < 5 || std(x)==0 || std(y)==0
                continue
            end

            % -------- within-unit z-score --------
            xz = (x - mean(x)) ./ std(x);
            yz = (y - mean(y)) ./ std(y);

            x_all = [x_all; xz];
            y_all = [y_all; yz];

            % per-unit correlation
            r_units(end+1,1) = corr(xz, yz, 'type','Pearson');
        end

        % pooled correlation on normalized data
        [r_pool, p_pool] = corr(x_all, y_all, 'type','Pearson');

        % regression line
        P = polyfit(x_all, y_all, 1);
        xGrid = linspace(min(x_all), max(x_all), 100);
        yFit = polyval(P, xGrid);

        % scatter (density-friendly)
        scatter(x_all, y_all, 10, 'k', 'filled', 'MarkerFaceAlpha',0.25);

        plot(xGrid, yFit, 'r', 'LineWidth',2.5);

        xlabel('CSM magnitude (z-scored within unit)');
        ylabel('Neural response (z-scored within unit)');
        title(sprintf('%s', titles(g)));

        txt = sprintf('Pooled r = %.3f (p=%.3g)\nMean per-unit r = %.3f\nnUnits=%d', ...
            r_pool, p_pool, mean(r_units,'omitnan'), numel(r_units));

        text(0.03,0.97,txt,'Units','normalized','VerticalAlignment','top',...
            'BackgroundColor','w','Margin',4);

        box off
    end
    % ---- Distribution of per-unit correlations ----
    figure('Color','w'); hold on

    for g = 1:2
        Tg = T(T.category == groups(g),:);
        uList = unique(Tg.unitIdx);

        r_units = [];

        for iU = 1:numel(uList)
            u = uList(iU);
            ix = Tg.unitIdx == u;

            x = Tg.csmMag(ix);
            y = Tg.neuralFR(ix);

            good = isfinite(x) & isfinite(y);
            x = x(good); y = y(good);

            if numel(x) < 5 || std(x)==0 || std(y)==0
                continue
            end

            xz = (x - mean(x)) ./ std(x);
            yz = (y - mean(y)) ./ std(y);

            r_units(end+1,1) = corr(xz, yz);
        end

        subplot(1,2,g); hold on
        histogram(r_units, 20, 'FaceColor',[0.5 0.5 0.5]);
        xline(mean(r_units,'omitnan'),'r','LineWidth',2);
        xlabel('Per-unit r (z-scored)');
        ylabel('Count');
        title(titles(g));
        box off
    end

else
    fprintf('[Z-SCORE CORR] Need Scale.table first.\n');
end

%%

function [vel_1ms, t_ms, accel_1ms] = getCSMAlignedKinematics_1ms(Reaches_use, win_pre_ms, win_post_ms, kinVar, kin_dim, posThresh, alignDir)
% Kinematics aligned to your CSM threshold-based event (1 ms)
% Uses AlignThreshCSMHF_1ms_edit (your function)

    C = AlignThreshCSMHF_1ms_edit(Reaches_use, posThresh, alignDir, win_pre_ms, win_post_ms);

    % If you want accel available, keep this (it will compute C.stimaccel if your stimvel_accel does)
    C = stimvel_accel(C);

    vel_1ms = squeeze(C.stimvel(:,kin_dim,:));
    t_ms    = C.stimvel(:,1,1);      % time is stored in column 1 (t_rel), should be same for all trials
    t_ms    = t_ms(:) * 1000;        % convert seconds -> ms (your script uses ms windows)

    if isfield(C,'stimaccel') && ~isempty(C.stimaccel)
        accel_1ms = squeeze(C.stimaccel(:,kin_dim,:));
    else
        accel_1ms = [];
    end
end

%% ===================== VIS: PER-UNIT CORRELATION SCATTER (POS vs NEG) =====================
% Computes per-unit Pearson r from Scale.table (trial-level rows)
% Then plots jittered scatter + big mean circle + bootstrap CI errorbar
%
% Requires Scale.table columns: unitIdx, neuralFR, csmMag, category
% Uses min_trials_for_unit_corr as the minimum trials per unit.

if do_scaling && isfield(Scale,'table') && ~isempty(Scale.table)

    T = Scale.table;

    % --- ensure category is string ---
    if iscell(T.category)
        T.category = string(T.category);
    elseif ischar(T.category)
        T.category = string(cellstr(T.category));
    elseif iscategorical(T.category)
        T.category = string(T.category);
    else
        T.category = string(T.category);
    end

    % groups
    isPos = (T.category == "negLag_posCorr");
    isNeg = (T.category == "negLag_negCorr");

    % Minimum trials per unit for correlation
    minTrials = min_trials_for_unit_corr;

    % ---------- compute per-unit Pearson r for POS ----------
    Tg = T(isPos,:);
    Tg = Tg(isfinite(Tg.neuralFR) & isfinite(Tg.csmMag),:);

    rPos = [];
    if ~isempty(Tg)
        uList = unique(Tg.unitIdx);
        rPos = nan(numel(uList),1);

        for iU = 1:numel(uList)
            u = uList(iU);
            ix = (Tg.unitIdx == u);

            x = Tg.csmMag(ix);
            y = Tg.neuralFR(ix);

            good = isfinite(x) & isfinite(y);
            if sum(good) < minTrials, continue; end
            if numel(unique(x(good))) < 2, continue; end

            rPos(iU) = corr(x(good), y(good), 'type','Pearson', 'rows','complete');
        end

        rPos = rPos(isfinite(rPos));
    end

    % ---------- compute per-unit Pearson r for NEG ----------
    Tg = T(isNeg,:);
    Tg = Tg(isfinite(Tg.neuralFR) & isfinite(Tg.csmMag),:);

    rNeg = [];
    if ~isempty(Tg)
        uList = unique(Tg.unitIdx);
        rNeg = nan(numel(uList),1);

        for iU = 1:numel(uList)
            u = uList(iU);
            ix = (Tg.unitIdx == u);

            x = Tg.csmMag(ix);
            y = Tg.neuralFR(ix);

            good = isfinite(x) & isfinite(y);
            if sum(good) < minTrials, continue; end
            if numel(unique(x(good))) < 2, continue; end

            rNeg(iU) = corr(x(good), y(good), 'type','Pearson', 'rows','complete');
        end

        rNeg = rNeg(isfinite(rNeg));
    end

    % ---------- plotting ----------
    figure('Color','w','Position',[300 300 430 420]); hold on;

    jit = 0.08;

    % scatter points
    if ~isempty(rPos)
        xPos = 1 + (rand(size(rPos))-0.5)*2*jit;
        scatter(xPos, rPos, 40, 'filled', 'MarkerFaceAlpha',0.6);
    end
    if ~isempty(rNeg)
        xNeg = 2 + (rand(size(rNeg))-0.5)*2*jit;
        scatter(xNeg, rNeg, 40, 'filled', 'MarkerFaceAlpha',0.6);
    end

    % ---------- big mean circles + bootstrap CI error bars ----------
    nBoot = 5000;
    alphaCI = 0.05;

    % POS
    if ~isempty(rPos)
        vals = rPos(:);
        vals = vals(isfinite(vals));
        if ~isempty(vals)
            bootMean = zeros(nBoot,1);
            n = numel(vals);
            for b = 1:nBoot
                samp = vals(randi(n,n,1));
                bootMean(b) = mean(samp,'omitnan');
            end
            mPos = mean(vals,'omitnan');
            ciPos = prctile(bootMean, [100*alphaCI/2, 100*(1-alphaCI/2)]);

            scatter(1, mPos, 160, 'k', 'filled');  % big circle
            errorbar(1, mPos, mPos-ciPos(1), ciPos(2)-mPos, 'k', ...
                'LineWidth',2, 'CapSize',12);
        end
    end

    % NEG
    if ~isempty(rNeg)
        vals = rNeg(:);
        vals = vals(isfinite(vals));
        if ~isempty(vals)
            bootMean = zeros(nBoot,1);
            n = numel(vals);
            for b = 1:nBoot
                samp = vals(randi(n,n,1));
                bootMean(b) = mean(samp,'omitnan');
            end
            mNeg = mean(vals,'omitnan');
            ciNeg = prctile(bootMean, [100*alphaCI/2, 100*(1-alphaCI/2)]);

            scatter(2, mNeg, 160, 'k', 'filled');  % big circle
            errorbar(2, mNeg, mNeg-ciNeg(1), ciNeg(2)-mNeg, 'k', ...
                'LineWidth',2, 'CapSize',12);
        end
    end

    yline(0,'--');
    xlim([0.5 2.5]);
    xticks([1 2]);
    xticklabels({'posCorr','negCorr'});
    ylabel('Per-unit Pearson r (neuralFR vs CSM magnitude)');
    title(sprintf('Per-unit scaling correlations (minTrials=%d)', minTrials));
    box off;

    % store
    Scale.rPos = rPos;
    Scale.rNeg = rNeg;

    fprintf('\n[UNIT r SCATTER] posCorr: nUnits=%d | negCorr: nUnits=%d\n', numel(rPos), numel(rNeg));

else
    fprintf('\n[UNIT r SCATTER] Need Scale.table first.\n');
end

%% SANITY CHECK 

T = Scale.table;
T.category = string(T.category);

cats = ["negLag_posCorr","negLag_negCorr"];
for c = 1:2
    Tg = T(T.category==cats(c),:);
    Tg = Tg(isfinite(Tg.neuralFR) & isfinite(Tg.csmMag),:);
    uList = unique(Tg.unitIdx);
    nTr = nan(numel(uList),1);
    for iU=1:numel(uList)
        nTr(iU) = sum(Tg.unitIdx==uList(iU));
    end
    fprintf('%s: units=%d | median trials/unit=%g | mean trials/unit=%g\n', ...
        cats(c), numel(uList), median(nTr,'omitnan'), mean(nTr,'omitnan'));
end

T = Scale.table;
T.category = string(T.category);

for cats = ["negLag_posCorr","negLag_negCorr"]
    Tg = T(T.category==cats,:);
    Tg = Tg(isfinite(Tg.neuralFR) & isfinite(Tg.csmMag),:);

    uList = unique(Tg.unitIdx);
    xC = []; yC = [];

    for iU=1:numel(uList)
        ix = Tg.unitIdx==uList(iU);
        x = Tg.csmMag(ix); y = Tg.neuralFR(ix);
        good = isfinite(x)&isfinite(y);
        if sum(good) < min_trials_for_unit_corr, continue; end

        % within-unit centering removes between-unit offsets
        xC = [xC; x(good)-mean(x(good),'omitnan')];
        yC = [yC; y(good)-mean(y(good),'omitnan')];
    end

    [rC,pC] = corr(xC,yC,'type','Pearson','rows','complete');
    fprintf('%s: within-unit-centered pooled r=%.3f (p=%.3g)\n', cats, rC, pC);
end

T = Scale.table; T.category = string(T.category);
for cats = ["negLag_posCorr","negLag_negCorr"]
    Tg = T(T.category==cats,:);
    Tg = Tg(isfinite(Tg.neuralFR) & isfinite(Tg.csmMag),:);

    uList = unique(Tg.unitIdx);
    xr = nan(numel(uList),1);
    for iU=1:numel(uList)
        ix = Tg.unitIdx==uList(iU);
        x = Tg.csmMag(ix);
        x = x(isfinite(x));
        if numel(x) < min_trials_for_unit_corr, continue; end
        xr(iU) = max(x)-min(x);
    end
    fprintf('%s: median csmMag range/unit = %.3f\n', cats, median(xr,'omitnan'));
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
time_idx     = find(full_timevec >= -250 & full_timevec <= 250);
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

    %xlabel('Time from max velocity (ms)');   % edit label if you want (csmT2 etc.)
    xlabel('Time from threshold crossing (ms)');   % edit label if you want (csmT2 etc.)
    ylabel('\Delta Z (Condition - Baseline)');
    xlim([-250 250]);
    title(sprintf('ALL-KEPT Δ: %s (N=%d)', cat, size(BL,1)), 'Interpreter','none');
    legend({'Stim-BL','', 'WO-BL',''}, 'Location','best','Box','off'); % blank entries hide fills
    box off;
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

title('ALL-KEPT Δ (AD - BL) : posCorr vs negCorr');

legend({'posCorr AD-BL','','negCorr AD-BL',''}, ...
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
