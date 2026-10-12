%% integrated_xcorr_and_psth_classification_NEG_LAG_ONLY.m
% EDITED: NEGATIVE LAG ONLY -> classify ONLY:
%   (1) negLag_posCorr
%   (2) negLag_negCorr
%
% Everything else left as close to your original as possible.

clear; clc;


cd("C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\ShamAnimals\ReachandCSMMod_allCSMStogther\CSMAll_wtihMaxV")
MyFolderInfo = dir('*CSMAll.mat');

%% ---------------------- ANALYSIS PARAMS ---------------------
rasterField   = 'max_vx';
kinVar        = 7;
gauss_sigma   = 10;
resample_fs   = 200;
win_pre_ms    = 100;
win_post_ms   = 150;
maxLag_ms     = 100;
lag_step      = 5;

%baseline_win_def = [-1000 -750];
baseline_win_def = [-500 -250];
%baseline_win_def = [-1000 -500];
%valid_peak_range_ms = [-100 100];
valid_peak_range_ms = [-250 250];

%%% ---------------- NEW CSM ANALYSIS PARAMS ----------------
csm_epoch_win      = [-100 100];   % use for CSM T0- larger window to see if more is captureed
csm_baseline_win   = [-500 -250];  % pre-CSM baseline (ms)

min_spikes_csm     = 5;           
nPerm              = 5000;
alpha_perm         = 0.05;

%%% time-resolved analysis
slide_win_ms       = 25;          
slide_step_ms      = 5;
cluster_alpha      = 0.05;        


XCorrData = struct;
AllGroups = struct();
NPdata = struct;

% ---------------- ONLY 2 NEG-LAG GROUPS ----------------
group_names = {'negLag_posCorr','negLag_negCorr'};
for g = 1:numel(group_names)
    gn = group_names{g};
    AllGroups.(gn).CSM_before = [];
    AllGroups.(gn).CSM_after  = [];
    AllGroups.(gn).NoCSM_before = [];
    AllGroups.(gn).NoCSM_after  = [];
end

%% ---------------------- MAIN SESSION LOOP ----------------------
for day = 1:length(MyFolderInfo)

    fprintf('\n========== SESSION %d / %d ==========\n', day, length(MyFolderInfo));
    load(MyFolderInfo(day).name, 'Reaches','data','filtered_modulated_indices','ExcitRR');

    % Save kinematics
    NPdata(day).CSM   = ExcitRR.ReReaches;
    NPdata(day).NoCSM = ExcitRR.NoCSM;

    session_units = filtered_modulated_indices;
    nUnits = length(session_units);

    % timevec for PSTHs (use your stored timebins)
    timevec = data.units(session_units(1)).max_vx.time_bins; % ms vector
    zero_idx = find(abs(timevec) == min(abs(timevec)), 1);
    baseline_win = find(timevec >= baseline_win_def(1) & timevec <= baseline_win_def(2));
    valid_peak_range = find(timevec >= valid_peak_range_ms(1) & timevec <= valid_peak_range_ms(2));

    % containers for session results
    before_idx_all = [];
    after_idx_all  = [];

    for uu = 1:nUnits
        unitIdx = session_units(uu);
        fprintf('  Analyzing UNIT %d...\n', unitIdx);

        % ----------------- LOAD SPIKE RASTER -----------------
        spikes_full = data.units(unitIdx).max_vx.Raster;
        [nTrials, nTime] = size(spikes_full);
        full_time_1ms = (-1000:1000)'; % original 1ms indexing used in your IFR calc

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
        C = AlignPeakV_1ms_edit(Reaches, kinVar, win_pre_ms, win_post_ms);
        C = stimvel_accel(C);
        vel_1ms_all = squeeze(C.stimvel(:,5,:)); % using velocity
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

                fr_z = (fr_seg - mean(fr_seg)) ./ std(fr_seg);
                v_z  = (v_seg - mean(v_seg)) ./ std(v_seg);

                if any(isnan(fr_z)) || any(isnan(v_z)), continue; end

                xcorr_trials(li,tr) = mean(fr_z .* v_z);
            end
        end

        % ---------------- MEDIAN + CONTINUOUS PEAK/MIN (NEGATIVE LAGS ONLY) ----------------
        medianCorr = nanmedian(xcorr_trials,2);

        neg_inds = find(lags_ms < 0);
        if isempty(neg_inds)
            peakLag_cont = NaN; peakVal_cont = NaN;
            minLag_cont  = NaN; minVal_cont  = NaN;
            chosenLag    = NaN; chosenVal    = NaN;
        else
            % Max within negative lags
            [~, iNmax] = max(medianCorr(neg_inds));
            idxMax = neg_inds(iNmax);

            if idxMax > 1 && idxMax < length(medianCorr)
                x = lags_ms(idxMax-1 : idxMax+1);
                y = medianCorr(idxMax-1 : idxMax+1);
                p = polyfit(x, y, 2);
                peakLag_cont = -p(2)/(2*p(1));
                peakVal_cont = polyval(p, peakLag_cont);

                % enforce negative lag (polyfit can drift slightly)
                if ~(peakLag_cont < 0)
                    peakLag_cont = lags_ms(idxMax);
                    peakVal_cont = medianCorr(idxMax);
                end
            else
                peakLag_cont = lags_ms(idxMax);
                peakVal_cont = medianCorr(idxMax);
            end

            % Min within negative lags
            [~, iNmin] = min(medianCorr(neg_inds));
            idxMin = neg_inds(iNmin);

            if idxMin > 1 && idxMin < length(medianCorr)
                x = lags_ms(idxMin-1 : idxMin+1);
                y = medianCorr(idxMin-1 : idxMin+1);
                p = polyfit(x, y, 2);
                minLag_cont = -p(2)/(2*p(1));
                minVal_cont = polyval(p, minLag_cont);

                % enforce negative lag
                if ~(minLag_cont < 0)
                    minLag_cont = lags_ms(idxMin);
                    minVal_cont = medianCorr(idxMin);
                end
            else
                minLag_cont = lags_ms(idxMin);
                minVal_cont = medianCorr(idxMin);
            end

            % SELECT BY LARGEST ABSOLUTE MAG (NEGATIVE LAGS ONLY)
            if abs(peakVal_cont) >= abs(minVal_cont)
                chosenVal = peakVal_cont;
                chosenLag = peakLag_cont;
            else
                chosenVal = minVal_cont;
                chosenLag = minLag_cont;
            end
        end

        % store xcorr results
        XCorrData(day).units(uu).unitIdx      = unitIdx;
        XCorrData(day).units(uu).lags_ms      = lags_ms;
        XCorrData(day).units(uu).xcorr_trials = xcorr_trials;
        XCorrData(day).units(uu).medianCorr   = medianCorr;

        % keep these fields for compatibility with your downstream code
        XCorrData(day).units(uu).peakVal      = peakVal_cont;
        XCorrData(day).units(uu).peakLag      = peakLag_cont;
        XCorrData(day).units(uu).minVal       = minVal_cont;
        XCorrData(day).units(uu).minLag       = minLag_cont;
        XCorrData(day).units(uu).chosenVal    = chosenVal;
        XCorrData(day).units(uu).chosenLag    = chosenLag;

        % ================== NEW (REVIEWER-FACING): EPOCH GLM + TIME-RESOLVED DELTAS (RAW + BASELINE-SUB) =================
        rng(0,'twister');

        % ================== LOAD RASTERS ==================
        BL_raster   = data.units(unitIdx).NoCSM_csmT2.Raster;   % Sham Control: NoCSM reaches
        Stim_raster = data.units(unitIdx).AllCSM_csmT2.Raster;  % Sham Control: AllCSM reaches

        % BL_raster   = data.units(unitIdx).NoCSM_csmT0.Raster;   % Sham Control: NoCSM reaches
        % Stim_raster = data.units(unitIdx).AllCSM_csmT0.Raster;  % Sham Control: AllCSM reaches

        % ================== TIME VECTOR ==================
        nTime_bl   = size(BL_raster,2);
        nTime_stim = size(Stim_raster,2);
        nTime2     = min(nTime_bl, nTime_stim);
        BL_raster   = BL_raster(:,1:nTime2);
        Stim_raster = Stim_raster(:,1:nTime2);

        full_time_1ms2 = (-1000:1000)';   % 2001 bins
        assert(nTime2 <= numel(full_time_1ms2), 'Raster longer than assumed -1000:1000 ms vector.');
        time_csm = full_time_1ms2(1:nTime2);

        % ================== EPOCH INDICES ==================
        csm_idx  = (time_csm >= csm_epoch_win(1))    & (time_csm <= csm_epoch_win(2));
        base_idx = (time_csm >= csm_baseline_win(1)) & (time_csm <= csm_baseline_win(2));

        assert(any(csm_idx),  'CSM epoch window is empty (check csm_epoch_win vs time_csm).');
        assert(any(base_idx), 'Baseline window is empty (check csm_baseline_win vs time_csm).');

        % ================== COUNTS PER TRIAL (EPOCH) ==================
        BL_base_spk   = sum(BL_raster(:,base_idx),2);
        Stim_base_spk = sum(Stim_raster(:,base_idx),2);
        BL_csm_spk    = sum(BL_raster(:,csm_idx),2);
        Stim_csm_spk  = sum(Stim_raster(:,csm_idx),2);

        % Inclusion gate
        if ( (sum(BL_csm_spk) + sum(Stim_csm_spk)) < min_spikes_csm )
            XCorrData(day).units(uu).CSM_epoch_glm = struct('skip',true,'reason','min_spikes_csm');
            XCorrData(day).units(uu).CSM_time_fast = struct('skip',true,'reason','min_spikes_csm');
        else

            % ================== EPOCH POISSON GLM + LR TEST (unchanged) ==================
            y    = [Stim_csm_spk; BL_csm_spk];
            cond = [ones(numel(Stim_csm_spk),1); zeros(numel(BL_csm_spk),1)];
            base = [Stim_base_spk; BL_base_spk];

            sd_base = std(base);
            if sd_base < 1e-6
                base_z = zeros(size(base));
            else
                base_z = (base - mean(base)) ./ sd_base;
            end

            Xfull = [cond, base_z];
            Xred  = [base_z];

            opts = statset('glmfit');
            opts.MaxIter = 200;

            if sum(y) == 0 || all(y==0)
                XCorrData(day).units(uu).CSM_epoch_glm = struct('skip',true,'reason','all_zero_counts');
            else
                try
                    [b_full, dev_full, stats_full] = glmfit(Xfull, y, 'poisson', 'link','log', 'options', opts);
                    [~,      dev_red]              = glmfit(Xred,  y, 'poisson', 'link','log', 'options', opts);

                    LR   = dev_red - dev_full;
                    p_LR = 1 - chi2cdf(LR, 1);

                    beta_cond = b_full(2);
                    RR        = exp(beta_cond);

                    se_cond = stats_full.se(2);
                    CI_beta = beta_cond + [-1.96 1.96]*se_cond;
                    CI_RR   = exp(CI_beta);

                    csm_win_s  = sum(csm_idx)  / 1000;
                    base_win_s = sum(base_idx) / 1000;

                    BL_netFR   = (BL_csm_spk./csm_win_s)    - (BL_base_spk./base_win_s);
                    Stim_netFR = (Stim_csm_spk./csm_win_s)  - (Stim_base_spk./base_win_s);
                    obs_delta_netFR = mean(Stim_netFR) - mean(BL_netFR);

                    XCorrData(day).units(uu).CSM_epoch_glm.skip        = false;
                    XCorrData(day).units(uu).CSM_epoch_glm.p_LR        = p_LR;
                    XCorrData(day).units(uu).CSM_epoch_glm.sig_LR      = (p_LR < alpha_perm);
                    XCorrData(day).units(uu).CSM_epoch_glm.LR          = LR;
                    XCorrData(day).units(uu).CSM_epoch_glm.beta_cond   = beta_cond;
                    XCorrData(day).units(uu).CSM_epoch_glm.RR          = RR;
                    XCorrData(day).units(uu).CSM_epoch_glm.CI_RR       = CI_RR;
                    XCorrData(day).units(uu).CSM_epoch_glm.delta_netFR = obs_delta_netFR;

                catch ME
                    XCorrData(day).units(uu).CSM_epoch_glm = struct('skip',true,'reason',['glmfit_failed: ' ME.message]);
                end
            end

            % ================== TIME-RESOLVED (FAST): RAW + NET DELTA + MAX-STAT PERM ==================
            bin_ms      = 10;
            bin_step_ms = 5;
            min_total_spikes_per_bin = 5;

            bin_starts = csm_epoch_win(1):bin_step_ms:(csm_epoch_win(2)-bin_ms);
            nBins      = numel(bin_starts);

            base_win_s = sum(base_idx) / 1000;
            bin_s      = bin_ms / 1000;

            BL_base_FR   = BL_base_spk   / base_win_s;
            Stim_base_FR = Stim_base_spk / base_win_s;

            nStim   = size(Stim_raster,1);
            nBL     = size(BL_raster,1);
            labels0 = [true(nStim,1); false(nBL,1)];

            rawFR_all = nan(nStim+nBL, nBins);
            netFR_all = nan(nStim+nBL, nBins);

            delta_raw = nan(nBins,1);
            delta_net = nan(nBins,1);

            for bbi = 1:nBins
                b0 = bin_starts(bbi);
                b_idx = (time_csm >= b0) & (time_csm < (b0 + bin_ms));
                if ~any(b_idx), continue; end

                Stim_bin_spk = sum(Stim_raster(:,b_idx),2);
                BL_bin_spk   = sum(BL_raster(:,b_idx),2);

                if (sum(Stim_bin_spk) + sum(BL_bin_spk)) < min_total_spikes_per_bin
                    continue;
                end

                stimFR = Stim_bin_spk / bin_s;
                blFR   = BL_bin_spk   / bin_s;

                rawFR_all(:,bbi) = [stimFR; blFR];
                netFR_all(:,bbi) = [stimFR - Stim_base_FR; blFR - BL_base_FR];

                delta_raw(bbi) = mean(rawFR_all(labels0,bbi), 'omitnan') - mean(rawFR_all(~labels0,bbi), 'omitnan');
                delta_net(bbi) = mean(netFR_all(labels0,bbi), 'omitnan') - mean(netFR_all(~labels0,bbi), 'omitnan');
            end

            max_null_raw = nan(nPerm,1);
            max_null_net = nan(nPerm,1);

            for pp = 1:nPerm
                permLab = labels0(randperm(numel(labels0)));

                perm_delta_raw = nan(nBins,1);
                perm_delta_net = nan(nBins,1);

                for bbi = 1:nBins
                    xr = rawFR_all(:,bbi);
                    xn = netFR_all(:,bbi);
                    if all(isnan(xr)), continue; end

                    perm_delta_raw(bbi) = mean(xr(permLab), 'omitnan') - mean(xr(~permLab), 'omitnan');
                    perm_delta_net(bbi) = mean(xn(permLab), 'omitnan') - mean(xn(~permLab), 'omitnan');
                end

                max_null_raw(pp) = max(abs(perm_delta_raw), [], 'omitnan');
                max_null_net(pp) = max(abs(perm_delta_net), [], 'omitnan');
            end

            p_corr_raw = nan(nBins,1);
            p_corr_net = nan(nBins,1);
            for bbi = 1:nBins
                if ~isnan(delta_raw(bbi))
                    p_corr_raw(bbi) = mean(max_null_raw >= abs(delta_raw(bbi)));
                end
                if ~isnan(delta_net(bbi))
                    p_corr_net(bbi) = mean(max_null_net >= abs(delta_net(bbi)));
                end
            end

            step_s = bin_step_ms/1000;
            AUC_raw = sum(delta_raw(isfinite(delta_raw))) * step_s;
            AUC_net = sum(delta_net(isfinite(delta_net))) * step_s;

            XCorrData(day).units(uu).CSM_time_fast.skip         = false;
            XCorrData(day).units(uu).CSM_time_fast.bin_ms       = bin_ms;
            XCorrData(day).units(uu).CSM_time_fast.step_ms      = bin_step_ms;
            XCorrData(day).units(uu).CSM_time_fast.starts_ms    = bin_starts(:);

            XCorrData(day).units(uu).CSM_time_fast.delta_rawFR  = delta_raw(:);
            XCorrData(day).units(uu).CSM_time_fast.delta_netFR  = delta_net(:);

            XCorrData(day).units(uu).CSM_time_fast.p_corr_raw   = p_corr_raw(:);
            XCorrData(day).units(uu).CSM_time_fast.p_corr_net   = p_corr_net(:);

            XCorrData(day).units(uu).CSM_time_fast.AUC_raw_Hzs  = AUC_raw;
            XCorrData(day).units(uu).CSM_time_fast.AUC_net_Hzs  = AUC_net;

            XCorrData(day).units(uu).CSM_time_fast.max_null_raw = max_null_raw;
            XCorrData(day).units(uu).CSM_time_fast.max_null_net = max_null_net;
        end
        % ================== END NEW BLOCK ==================

        % ---------------- PEAK-TIME BEFORE/AFTER & KEEP CHECK ----------------
        psth = data.units(unitIdx).max_vx.PSTH_smoothbin;
        [~, peak_idx] = max(psth);

        % only keep units whose peak falls in valid range
        if ~ismember(peak_idx, valid_peak_range)
            XCorrData(day).units(uu).keep = false;
            continue;
        else
            XCorrData(day).units(uu).keep = true;
        end

        % Before vs After based on psth peak relative to zero
        if timevec(peak_idx) < 0
            XCorrData(day).units(uu).before_after = 'before';
            before_idx_all(end+1) = unitIdx;
        else
            XCorrData(day).units(uu).before_after = 'after';
            after_idx_all(end+1) = unitIdx;
        end

        % ---------------- ASSIGN CATEGORY (NEGATIVE LAG ONLY; corr sign only) ----------------
        % If for some reason chosenLag drifted to >=0, force to NaN and skip grouping.
        if isnan(chosenLag) || ~(chosenLag < 0) || isnan(chosenVal)
            XCorrData(day).units(uu).category = 'negLag_posCorr'; % safe fallback label
            % NOTE: still allow downstream, but you can "continue" here if you'd rather drop:
            % continue;
        else
            if chosenVal >= 0
                cat_label = 'negLag_posCorr';
            else
                cat_label = 'negLag_negCorr';
            end
            XCorrData(day).units(uu).category = cat_label;
        end

        % ---------------- COMPUTE Z-SCORED PSTHs ----------------
        csm    = data.units(unitIdx).AllCSM_csmT2.PSTH_smoothbin;
        nocsm  = data.units(unitIdx).NoCSM_csmT2.PSTH_smoothbin;

        baseline_mean = mean(csm(baseline_win));
        baseline_std  = std(csm(baseline_win));
        if baseline_std == 0, baseline_std = 1; end

        csm_z   = (csm   - baseline_mean) / baseline_std;
        nocsm_z = (nocsm - baseline_mean) / baseline_std;

        % ---------------- AGGREGATE INTO GROUPS ----------------
        if exist('cat_label','var') && ismember(cat_label, group_names)
            if strcmp(XCorrData(day).units(uu).before_after,'before')
                AllGroups.(cat_label).CSM_before    = vertcat(AllGroups.(cat_label).CSM_before, csm_z);
                AllGroups.(cat_label).NoCSM_before  = vertcat(AllGroups.(cat_label).NoCSM_before, nocsm_z);
            else
                AllGroups.(cat_label).CSM_after     = vertcat(AllGroups.(cat_label).CSM_after, csm_z);
                AllGroups.(cat_label).NoCSM_after   = vertcat(AllGroups.(cat_label).NoCSM_after, nocsm_z);
            end
        end

        clear cat_label; % prevent carryover if a unit falls into fallback path

    end % unit loop

    % Save session-level lists
    XCorrData(day).before_peakv_units = before_idx_all;
    XCorrData(day).after_peakv_units  = after_idx_all;
    fprintf(' Session %d summary: %d before, %d after (kept units)\n', day, length(before_idx_all), length(after_idx_all));
end

fprintf('\n================ ALL SESSIONS COMPLETE ================\n');



%% ============================================================
% PUBLICATION PLOTTING: GLM RESULTS (CSM EPOCH) FOR TUNED+KEPT UNITS
% Split by NEG-LAG correlation sign:
%   - negLag_posCorr
%   - negLag_negCorr
%
% Plots (NO PSTHs):
%   Fig 1: Effect vs significance (beta_cond vs -log10(p_LR))
%   Fig 2: Effect size distribution (beta_cond and delta_netFR) by group
%   Fig 3: Fraction significant (sig_LR) by group + 95% binomial CI
%
% Requires XCorrData in workspace.
% Uses fields:
%   U.keep
%   U.category  ('negLag_posCorr' / 'negLag_negCorr')
%   U.CSM_epoch_glm.skip / sig_LR / p_LR / beta_cond / RR / CI_RR / delta_netFR / LR
% ============================================================

assert(exist('XCorrData','var')==1, 'XCorrData not found in workspace.');

cats = {'negLag_posCorr','negLag_negCorr'};

% colors (use rgb() if you have it; otherwise fixed RGB triplets)
if exist('rgb','file')==2 || exist('rgb','builtin')==5
    col_pos = rgb('SteelBlue');
    col_neg = rgb('DarkGreen');
else
    col_pos = [70 130 180]/255; % steelblue-ish
    col_neg = [0 100 0]/255;    % darkgreen-ish
end

% ---------------- COLLECT GLM METRICS ----------------
Cat      = {};
pLR      = [];
sigLR    = [];
betaCond = [];
RR       = [];
LR       = [];
dNetFR   = [];
dayIdx   = [];
unitIdxV = [];

for day = 1:numel(XCorrData)
    if ~isfield(XCorrData(day),'units') || isempty(XCorrData(day).units), continue; end
    for uu = 1:numel(XCorrData(day).units)
        U = XCorrData(day).units(uu);

        if ~isfield(U,'keep') || ~U.keep, continue; end
        if ~isfield(U,'category') || ~ismember(U.category, cats), continue; end
        if ~isfield(U,'CSM_epoch_glm') || ~isfield(U.CSM_epoch_glm,'skip') || U.CSM_epoch_glm.skip, continue; end
        if ~isfield(U.CSM_epoch_glm,'p_LR') || ~isfinite(U.CSM_epoch_glm.p_LR), continue; end

        Cat{end+1,1}      = U.category;
        pLR(end+1,1)      = U.CSM_epoch_glm.p_LR;
        sigLR(end+1,1)    = logical(U.CSM_epoch_glm.sig_LR);
        betaCond(end+1,1) = U.CSM_epoch_glm.beta_cond;
        RR(end+1,1)       = U.CSM_epoch_glm.RR;
        LR(end+1,1)       = U.CSM_epoch_glm.LR;

        if isfield(U.CSM_epoch_glm,'delta_netFR') && isfinite(U.CSM_epoch_glm.delta_netFR)
            dNetFR(end+1,1) = U.CSM_epoch_glm.delta_netFR;
        else
            dNetFR(end+1,1) = NaN;
        end

        if isfield(U,'unitIdx'), unitIdxV(end+1,1) = U.unitIdx; else, unitIdxV(end+1,1) = NaN; end
        dayIdx(end+1,1) = day;
    end
end

Cat = categorical(Cat, cats, cats);
nTot = numel(pLR);

fprintf('\nGLM summary (KEPT only): N=%d units\n', nTot);
for ci = 1:numel(cats)
    this = (Cat == cats{ci});
    fprintf('  %s: N=%d, tuned(sig)= %d (%.1f%%)\n', cats{ci}, sum(this), sum(sigLR(this)), 100*mean(sigLR(this)));
end

% transform for plotting
neglogp = -log10(max(pLR, realmin));
alpha_line = -log10(0.05);

isPos = (Cat == 'negLag_posCorr');
isNeg = (Cat == 'negLag_negCorr');

%% ---------------- FIG 1: EFFECT VS SIGNIFICANCE (beta vs -log10 p) ----------------
fig1 = figure('Color','w','Position',[200 120 900 520]);
ax = axes(fig1); hold(ax,'on'); box(ax,'off'); set(ax,'LineWidth',1,'FontSize',12);

% non-sig (faint) then sig (solid) so sig is on top
ns_pos = isPos & ~sigLR;
ns_neg = isNeg & ~sigLR;
s_pos  = isPos &  sigLR;
s_neg  = isNeg &  sigLR;

scatter(ax, betaCond(ns_pos), neglogp(ns_pos), 28, 'MarkerFaceColor', col_pos, 'MarkerEdgeColor', col_pos, 'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0.25);
scatter(ax, betaCond(ns_neg), neglogp(ns_neg), 28, 'MarkerFaceColor', col_neg, 'MarkerEdgeColor', col_neg, 'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0.25);

scatter(ax, betaCond(s_pos),  neglogp(s_pos),  34, 'MarkerFaceColor', col_pos, 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);
scatter(ax, betaCond(s_neg),  neglogp(s_neg),  34, 'MarkerFaceColor', col_neg, 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);

xline(ax, 0, '-', 'LineWidth', 1);
yline(ax, alpha_line, '--', 'LineWidth', 1);

xlabel(ax,'GLM \beta_{cond} (log rate ratio)');
ylabel(ax,'-log_{10}(p_{LR})');
title(ax,'CSM epoch GLM: effect vs significance (kept units)');

legend(ax, { ...
    sprintf('posCorr (ns)  n=%d', sum(ns_pos)), ...
    sprintf('negCorr (ns)  n=%d', sum(ns_neg)), ...
    sprintf('posCorr (sig) n=%d', sum(s_pos)), ...
    sprintf('negCorr (sig) n=%d', sum(s_neg))}, ...
    'Location','best', 'Box','off');

%% ---------------- FIG 2: DISTRIBUTIONS BY GROUP ----------------
% Panel A: beta_cond (log rate ratio)
% Panel B: delta_netFR (Hz) if present
fig2 = figure('Color','w','Position',[220 120 1100 420]);

% ---- A: beta_cond ----
ax1 = subplot(1,2,1); hold(ax1,'on'); box(ax1,'off'); set(ax1,'LineWidth',1,'FontSize',12);
boxchart(ax1, Cat, betaCond, 'MarkerStyle','.', 'JitterOutliers','on');
yline(ax1,0,'-','LineWidth',1);
ylabel(ax1,'\beta_{cond} (log rate ratio)');
title(ax1,'Effect size by correlation group');

% color boxes (manual)
bc = findobj(ax1,'Type','BoxChart');
% BoxChart objects are returned in reverse order sometimes; set by category name robustly:
for k = 1:numel(bc)
    if strcmp(string(bc(k).XData(1)), "negLag_posCorr")
        bc(k).BoxFaceColor = col_pos;
    else
        bc(k).BoxFaceColor = col_neg;
    end
end

% ---- B: delta_netFR ----
ax2 = subplot(1,2,2); hold(ax2,'on'); box(ax2,'off'); set(ax2,'LineWidth',1,'FontSize',12);
valid_d = isfinite(dNetFR);
if any(valid_d)
    boxchart(ax2, Cat(valid_d), dNetFR(valid_d), 'MarkerStyle','.', 'JitterOutliers','on');
    yline(ax2,0,'-','LineWidth',1);
    ylabel(ax2,'\Delta net FR (Hz): mean(Stim) - mean(BL)');
    title(ax2,'Net firing-rate change by group');

    bc2 = findobj(ax2,'Type','BoxChart');
    for k = 1:numel(bc2)
        if strcmp(string(bc2(k).XData(1)), "negLag_posCorr")
            bc2(k).BoxFaceColor = col_pos;
        else
            bc2(k).BoxFaceColor = col_neg;
        end
    end
else
    text(ax2, 0.5, 0.5, 'delta\_netFR not available for these units', 'HorizontalAlignment','center');
    set(ax2,'XTick',[],'YTick',[]);
end

%% ---------------- FIG 3: FRACTION SIGNIFICANT BY GROUP (+ 95% CI) ----------------
% Wilson CI (no toolboxes needed)
wilson_ci = @(k,n,z) deal( ...
    (k/n + z^2/(2*n)) / (1 + z^2/n) - (z/(1+z^2/n))*sqrt((k/n)*(1-k/n)/n + z^2/(4*n^2)), ...
    (k/n + z^2/(2*n)) / (1 + z^2/n) + (z/(1+z^2/n))*sqrt((k/n)*(1-k/n)/n + z^2/(4*n^2)) );

z = 1.96;
k_pos = sum(sigLR(isPos)); n_pos = sum(isPos);
k_neg = sum(sigLR(isNeg)); n_neg = sum(isNeg);

[lo_pos, hi_pos] = wilson_ci(k_pos, n_pos, z);
[lo_neg, hi_neg] = wilson_ci(k_neg, n_neg, z);

fig3 = figure('Color','w','Position',[240 120 650 420]);
ax = axes(fig3); hold(ax,'on'); box(ax,'off'); set(ax,'LineWidth',1,'FontSize',12);

x = [1 2];
y = [k_pos/n_pos, k_neg/n_neg];
err_lo = y - [lo_pos, lo_neg];
err_hi = [hi_pos, hi_neg] - y;

b = bar(ax, x, y, 0.6);
b.FaceColor = 'flat';
b.CData(1,:) = col_pos;
b.CData(2,:) = col_neg;

errorbar(ax, x, y, err_lo, err_hi, 'k.', 'LineWidth', 1.5, 'CapSize', 12);

set(ax,'XTick',x,'XTickLabel',{'negLag\_posCorr','negLag\_negCorr'});
ylim(ax,[0 1]);
ylabel(ax,'Fraction tuned (p_{LR}<0.05)');
title(ax,'CSM GLM: tuned fraction by correlation group (kept units)');

% annotate counts
text(1, y(1)+0.05, sprintf('%d/%d', k_pos, n_pos), 'HorizontalAlignment','center', 'FontSize',12);
text(2, y(2)+0.05, sprintf('%d/%d', k_neg, n_neg), 'HorizontalAlignment','center', 'FontSize',12);

%% ============================================================
% FIG: NEG-LAG ONLY — Epoch DELTA_NETFR scatter by correlation sign
%   - NOT significant: gray
%   - sig +corr (negLag_posCorr): SteelBlue
%   - sig -corr (negLag_negCorr): DarkGreen
%
% Uses (per unit):
%   U.keep == true
%   U.category in {'negLag_posCorr','negLag_negCorr'}
%   U.CSM_epoch_glm.skip == false
%   U.CSM_epoch_glm.sig_LR
%   U.CSM_epoch_glm.delta_netFR
%
% Drop AFTER main script (XCorrData exists; MyFolderInfo exists).
%% ============================================================

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

% -------- collect per-unit GLM results --------
Unit = struct('cat',{},'day',{},'unitIdx',{});
dEp  = [];   % delta_netFR
sigE = [];   % sig_LR

for day = 1:numel(XCorrData)
    if ~isfield(XCorrData(day),'units') || isempty(XCorrData(day).units), continue; end

    for uu = 1:numel(XCorrData(day).units)
        U = XCorrData(day).units(uu);

        if ~isfield(U,'keep') || ~U.keep, continue; end
        if ~isfield(U,'category') || ~ismember(U.category, groups), continue; end
        if ~isfield(U,'CSM_epoch_glm') || ~isfield(U.CSM_epoch_glm,'skip') || U.CSM_epoch_glm.skip, continue; end
        if ~isfield(U.CSM_epoch_glm,'delta_netFR') || ~isfinite(U.CSM_epoch_glm.delta_netFR), continue; end
        if ~isfield(U.CSM_epoch_glm,'sig_LR'), continue; end

        Unit(end+1).cat     = U.category;
        Unit(end).day       = day;
        Unit(end).unitIdx   = U.unitIdx;

        dEp(end+1,1)  = U.CSM_epoch_glm.delta_netFR;
        sigE(end+1,1) = logical(U.CSM_epoch_glm.sig_LR);
    end
end

if isempty(dEp)
    warning('No units found matching: kept + negLag_(pos/neg)Corr + valid GLM.');
else
    fprintf('Plotting %d kept units with GLM results (neg-lag only).\n', numel(dEp));
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
figure('Color','w','Position',[120 120 300 400]); hold on;

sz=55;
% non-sig (gray)
scatter(xj(idxNo), dEp(idxNo), 45, 'filled', ...
    'MarkerFaceColor', c_gray, 'MarkerFaceAlpha', 0.55, 'MarkerEdgeColor','none');

% sig posCorr (steel blue)
scatter(xj(idxSigPos), dEp(idxSigPos), sz, 'filled', ...
    'MarkerFaceColor', c_pos, 'MarkerFaceAlpha', 0.95, 'MarkerEdgeColor','k');

% sig negCorr (dark green)
scatter(xj(idxSigNeg), dEp(idxSigNeg), sz, 'filled', ...
    'MarkerFaceColor', c_neg, 'MarkerFaceAlpha', 0.95, 'MarkerEdgeColor','k');

xlim([0.5 numel(groups)+0.5]);
xticks(1:numel(groups));
xticklabels({'negLag posCorr','negLag negCorr'});
xtickangle(20);
%ylim([-110 60])
%ylim([-60 40])

yline(0,'k-','LineWidth',1);
ylabel('\Delta net FR (Hz): epoch mean (Stim - BL)');
title('Epoch effect size (GLM); significant neurons colored by correlation sign');

% % group means + counts
% for g = 1:numel(groups)
%     m = (x==g);
%     if any(m)
%         mu = mean(dEp(m),'omitnan');
%         plot(g, mu, 'kd', 'MarkerFaceColor','k', 'MarkerSize',8);
%         text(g, mu, sprintf('  %d/%d', sum(m & idxSig), sum(m)), ...
%             'FontWeight','bold','VerticalAlignment','bottom');
%     end
% end

for g = 1:numel(groups)
    m_all = (x==g);          % all units in this category
    m_sig = m_all & idxSig;  % only significant units in this category

    if any(m_sig)
        mu = mean(dEp(m_sig),'omitnan');
        plot(g, mu, 'kd','MarkerFaceColor','k','MarkerSize',8);

        text(g, mu, sprintf('  %d sig / %d total', sum(m_sig), sum(m_all)), ...
            'FontWeight','bold','VerticalAlignment','bottom');
    end
end


legend({'not sig','sig +corr','sig -corr','mean'}, 'Location','best', 'Box','off');

set(gca,'LineWidth',1,'FontSize',12);

%% -------- stats for significant cells only --------
SigModStats = struct();

for g = 1:numel(groups)

    m_all = (x == g);
    m_sig = m_all & idxSig;

    vals = dEp(m_sig);
    vals = vals(isfinite(vals));

    SigModStats(g).group = groups{g};
    SigModStats(g).n_sig = numel(vals);
    SigModStats(g).n_total = sum(m_all);

    if isempty(vals)
        SigModStats(g).mean = NaN;
        SigModStats(g).sd   = NaN;
        SigModStats(g).sem  = NaN;
        SigModStats(g).median = NaN;
        SigModStats(g).iqr    = NaN;
        SigModStats(g).cohens_d_vs0 = NaN;
        SigModStats(g).ttest_p = NaN;
        SigModStats(g).ttest_t = NaN;
        SigModStats(g).ttest_df = NaN;
    else
        SigModStats(g).mean   = mean(vals,'omitnan');
        SigModStats(g).sd     = std(vals,'omitnan');
        SigModStats(g).sem    = std(vals,'omitnan') ./ sqrt(numel(vals));
        SigModStats(g).median = median(vals,'omitnan');
        SigModStats(g).iqr    = iqr(vals);

        % one-sample test vs 0
        [~,p,~,stats] = ttest(vals,0);
        SigModStats(g).ttest_p  = p;
        SigModStats(g).ttest_t  = stats.tstat;
        SigModStats(g).ttest_df = stats.df;

        % Cohen's d for modulation relative to 0
        SigModStats(g).cohens_d_vs0 = mean(vals,'omitnan') ./ std(vals,'omitnan');
    end

    fprintf('\n[%s]\n', groups{g});
    fprintf('significant cells: %d / %d\n', SigModStats(g).n_sig, SigModStats(g).n_total);
    fprintf('mean modulation = %.2f Hz\n', SigModStats(g).mean);
    fprintf('SD = %.2f Hz | SEM = %.2f Hz\n', SigModStats(g).sd, SigModStats(g).sem);
    fprintf('median = %.2f Hz | IQR = %.2f Hz\n', SigModStats(g).median, SigModStats(g).iqr);
    fprintf('one-sample t-test vs 0: t(%d) = %.3f, p = %.5f\n', ...
        SigModStats(g).ttest_df, SigModStats(g).ttest_t, SigModStats(g).ttest_p);
    fprintf('Cohen''s d vs 0 = %.3f\n', SigModStats(g).cohens_d_vs0);
end

%% ===================== SIGNIFICANT MODULATED UNITS ONLY =====================
% Plots only GLM-significant units and reports:
%   - N sig units in each subset
%   - mean modulation
%   - SD / SEM
%   - one-sample t-test vs 0
%   - Cohen's d vs 0
%
% Assumes you already have:
%   dEp        = epoch mean delta FR (Stim - BL) for each unit
%   idxSig     = logical index of significantly modulated units
%   idxSigPos  = logical index of sig +Corr units
%   idxSigNeg  = logical index of sig -Corr units
%   c_pos, c_neg
%
% This code ignores non-significant units completely.

% -------- collect significant values only --------
vals_pos = dEp(idxSigPos);
vals_neg = dEp(idxSigNeg);

vals_pos = vals_pos(isfinite(vals_pos));
vals_neg = vals_neg(isfinite(vals_neg));

%% -------- descriptive stats --------
StatsSig = struct();

% +Corr
StatsSig.pos.group = 'negLag_posCorr';
StatsSig.pos.n     = numel(vals_pos);
StatsSig.pos.mean  = mean(vals_pos,'omitnan');
StatsSig.pos.sd    = std(vals_pos,'omitnan');
StatsSig.pos.sem   = StatsSig.pos.sd ./ sqrt(max(1,StatsSig.pos.n));
StatsSig.pos.median = median(vals_pos,'omitnan');
StatsSig.pos.iqr    = iqr(vals_pos);

if StatsSig.pos.n > 1 && StatsSig.pos.sd > 0
    [~,p_pos,~,st_pos] = ttest(vals_pos,0);
    StatsSig.pos.p_vs0   = p_pos;
    StatsSig.pos.t_vs0   = st_pos.tstat;
    StatsSig.pos.df_vs0  = st_pos.df;
    StatsSig.pos.d_vs0   = StatsSig.pos.mean / StatsSig.pos.sd;
else
    StatsSig.pos.p_vs0   = NaN;
    StatsSig.pos.t_vs0   = NaN;
    StatsSig.pos.df_vs0  = NaN;
    StatsSig.pos.d_vs0   = NaN;
end

% -Corr
StatsSig.neg.group = 'negLag_negCorr';
StatsSig.neg.n     = numel(vals_neg);
StatsSig.neg.mean  = mean(vals_neg,'omitnan');
StatsSig.neg.sd    = std(vals_neg,'omitnan');
StatsSig.neg.sem   = StatsSig.neg.sd ./ sqrt(max(1,StatsSig.neg.n));
StatsSig.neg.median = median(vals_neg,'omitnan');
StatsSig.neg.iqr    = iqr(vals_neg);

if StatsSig.neg.n > 1 && StatsSig.neg.sd > 0
    [~,p_neg,~,st_neg] = ttest(vals_neg,0);
    StatsSig.neg.p_vs0   = p_neg;
    StatsSig.neg.t_vs0   = st_neg.tstat;
    StatsSig.neg.df_vs0  = st_neg.df;
    StatsSig.neg.d_vs0   = StatsSig.neg.mean / StatsSig.neg.sd;
else
    StatsSig.neg.p_vs0   = NaN;
    StatsSig.neg.t_vs0   = NaN;
    StatsSig.neg.df_vs0  = NaN;
    StatsSig.neg.d_vs0   = NaN;
end

%% -------- between-group comparison: sig +Corr vs sig -Corr --------

BetweenSig = struct();

if ~isempty(vals_pos) && ~isempty(vals_neg)
    [~,p_between,ci_between,st_between] = ttest2(vals_pos, vals_neg, 'Vartype','unequal');
    p_rank = ranksum(vals_pos, vals_neg);

    n1 = numel(vals_pos);
    n2 = numel(vals_neg);
    s1 = std(vals_pos,'omitnan');
    s2 = std(vals_neg,'omitnan');
    m1 = mean(vals_pos,'omitnan');
    m2 = mean(vals_neg,'omitnan');

    spooled = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / max(1,(n1+n2-2)));
    d_between = (m1 - m2) / spooled;

    BetweenSig.group1 = 'negLag_posCorr';
    BetweenSig.group2 = 'negLag_negCorr';
    BetweenSig.n1 = n1;
    BetweenSig.n2 = n2;
    BetweenSig.mean1 = m1;
    BetweenSig.mean2 = m2;
    BetweenSig.meanDiff = m1 - m2;
    BetweenSig.ttest_p = p_between;
    BetweenSig.ttest_t = st_between.tstat;
    BetweenSig.ttest_df = st_between.df;
    BetweenSig.ttest_ci = ci_between;
    BetweenSig.ranksum_p = p_rank;
    BetweenSig.cohens_d = d_between;
else
    BetweenSig = [];
end

%% -------- print stats --------
fprintf('\n================ SIGNIFICANT MODULATED UNITS ONLY ================\n');

fprintf('\n[%s]\n', StatsSig.pos.group);
fprintf('N = %d\n', StatsSig.pos.n);
fprintf('Mean ± SD = %.2f ± %.2f Hz\n', StatsSig.pos.mean, StatsSig.pos.sd);
fprintf('Mean ± SEM = %.2f ± %.2f Hz\n', StatsSig.pos.mean, StatsSig.pos.sem);
fprintf('Median [IQR] = %.2f [%.2f] Hz\n', StatsSig.pos.median, StatsSig.pos.iqr);
fprintf('vs 0: t(%g) = %.3f, p = %.5f, Cohen''s d = %.3f\n', ...
    StatsSig.pos.df_vs0, StatsSig.pos.t_vs0, StatsSig.pos.p_vs0, StatsSig.pos.d_vs0);

fprintf('\n[%s]\n', StatsSig.neg.group);
fprintf('N = %d\n', StatsSig.neg.n);
fprintf('Mean ± SD = %.2f ± %.2f Hz\n', StatsSig.neg.mean, StatsSig.neg.sd);
fprintf('Mean ± SEM = %.2f ± %.2f Hz\n', StatsSig.neg.mean, StatsSig.neg.sem);
fprintf('Median [IQR] = %.2f [%.2f] Hz\n', StatsSig.neg.median, StatsSig.neg.iqr);
fprintf('vs 0: t(%g) = %.3f, p = %.5f, Cohen''s d = %.3f\n', ...
    StatsSig.neg.df_vs0, StatsSig.neg.t_vs0, StatsSig.neg.p_vs0, StatsSig.neg.d_vs0);

if ~isempty(BetweenSig)
    fprintf('\n---------------- BETWEEN-GROUP COMPARISON ----------------\n');
    fprintf('%s vs %s\n', BetweenSig.group1, BetweenSig.group2);
    fprintf('N = %d vs %d\n', BetweenSig.n1, BetweenSig.n2);
    fprintf('Mean = %.2f vs %.2f Hz\n', BetweenSig.mean1, BetweenSig.mean2);
    fprintf('Mean difference = %.2f Hz\n', BetweenSig.meanDiff);
    fprintf('Welch t-test: t(%g) = %.3f, p = %.5f\n', ...
        BetweenSig.ttest_df, BetweenSig.ttest_t, BetweenSig.ttest_p);
    fprintf('Rank-sum p = %.5f\n', BetweenSig.ranksum_p);
    fprintf('Cohen''s d = %.3f\n', BetweenSig.cohens_d);
end

%% -------- make plotting vectors: significant units only --------
plotVals = [vals_pos(:); vals_neg(:)];
plotGrp  = [ones(numel(vals_pos),1); 2*ones(numel(vals_neg),1)];

% jitter
rng(1)
xj = plotGrp + (rand(size(plotGrp)) - 0.5) * 0.18;

%% -------- plot --------
figure('Color','w','Position',[120 120 300 400]); hold on;

sz = 60;

% +Corr significant units
scatter(xj(plotGrp==1), plotVals(plotGrp==1), sz, 'filled', ...
    'MarkerFaceColor', c_pos, 'MarkerFaceAlpha', 0.95, ...
    'MarkerEdgeColor', 'k');

% -Corr significant units
scatter(xj(plotGrp==2), plotVals(plotGrp==2), sz, 'filled', ...
    'MarkerFaceColor', c_neg, 'MarkerFaceAlpha', 0.95, ...
    'MarkerEdgeColor', 'k');

% mean markers
if ~isempty(vals_pos)
    plot(1, StatsSig.pos.mean, 'kd', 'MarkerFaceColor','k', 'MarkerSize',8);
    text(1, StatsSig.pos.mean, sprintf('  %.1f \\pm %.1f Hz\n  n = %d', ...
        StatsSig.pos.mean, StatsSig.pos.sd, StatsSig.pos.n), ...
        'FontWeight','bold', 'VerticalAlignment','bottom');
end

if ~isempty(vals_neg)
    plot(2, StatsSig.neg.mean, 'kd', 'MarkerFaceColor','k', 'MarkerSize',8);
    text(2, StatsSig.neg.mean, sprintf('  %.1f \\pm %.1f Hz\n  n = %d', ...
        StatsSig.neg.mean, StatsSig.neg.sd, StatsSig.neg.n), ...
        'FontWeight','bold', 'VerticalAlignment','bottom');
end

xlim([0.5 2.5]);
xticks([1 2]);
xticklabels({'negLag posCorr','negLag negCorr'});
xtickangle(20);

yline(0,'k-','LineWidth',1);

ylabel('\Delta net FR (Hz): epoch mean (Stim - BL)');
title('Significantly modulated units only');

legend({'sig +corr','sig -corr','mean'}, 'Location','best', 'Box','off');
set(gca,'LineWidth',1,'FontSize',12);

%% -------- optional table for export --------
SigOnlyTable = table( ...
    [repmat("negLag_posCorr",numel(vals_pos),1); repmat("negLag_negCorr",numel(vals_neg),1)], ...
    plotVals, ...
    'VariableNames', {'group','deltaEpochFR_Hz'});