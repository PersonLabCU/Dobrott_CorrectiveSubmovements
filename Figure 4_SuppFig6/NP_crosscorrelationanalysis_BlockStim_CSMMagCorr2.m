% Just using negative lag to categorize the units
% NOTE: Keeps your ORIGINAL BL-based lag/corr classification EXACTLY.
% NEW ADD-ON: After computing BL neg-lag peakLag, use abs(peakLag) to define a
%            unit-specific window (relative to CSMT2=0) and extract per-trial
%            peak FR + peak outward velocity from ADCSM_csmT2 + ExcitRR_AD.
%
% IMPORTANT: Replace the ONE line marked "<<< REPLACE THIS >>>" with your actual
% CSMT2 alignment function (whatever you already use in your pipeline).

clear; clc;

% ---------------------- LOAD SESSIONS -----------------------

cd("C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\StimAnimals\BlockSeparated_CSMvsNoCSM_inducedalignment")
MyFolderInfo = dir('*Block_CSMvsNoCSM.mat');

%% ---------------------- ANALYSIS PARAMS ---------------------
rasterField   = 'max_vx';
kinVar        = 7;   % aligning to max outward velocity
gauss_sigma   = 10;
resample_fs   = 200;
win_pre_ms    = 100;
win_post_ms   = 150;
maxLag_ms     = 100;
lag_step      = 5;

% ---------------------- NEW PARAMS (CSMT2 PEAK WINDOW) ---------------------
peakWin_half_ms = 50;   % +/- window around the unit-specific center (ms)

% ---------------------- MASTER STRUCT ----------------------
XCorrData = struct;

for day = 1:length(MyFolderInfo)

    fprintf('\n========== SESSION %d / %d ==========\n', day, length(MyFolderInfo));
    load(MyFolderInfo(day).name, 'Reaches','data','filtered_modulated_indices','ExcitRR_BL','ExcitRR_AD','ExcitRR_WOE','ExcitRR_WO','Salt_index');

    if isempty(Salt_index)
        fprintf('Day %d: Salt_index is empty. Skipping...\n', day);
        continue;
    end

    session_units = Salt_index;
    nUnits = length(session_units);

    % Time window for indexing PSTHs
    timevec = data.units(session_units(1)).BLCSM_max_vx.time_bins;
    valid_peak_range = find(timevec >= -250 & timevec <= 250);

    before_idx_all = [];
    after_idx_all  = [];

    % ------------------ LOOP ACROSS UNITS -------------------
    for uu = 1:nUnits

        unitIdx = session_units(uu);
        fprintf('  Analyzing UNIT %d...\n', unitIdx);

        % ----------------- LOAD SPIKE RASTER -----------------
        spikes_full = data.units(unitIdx).BLCSM_max_vx.Raster;
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
        ReReaches = ExcitRR_BL.ReReaches;
        C = AlignPeakV_1ms_edit(ReReaches, kinVar, win_pre_ms, win_post_ms);
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

        % ---------------- MEDIAN + NEGATIVE-LAG ABS-MAX ONLY ----------------
        medianCorr = nanmedian(xcorr_trials,2);

        negMask = lags_ms < 0;
        peakLag_cont = NaN;
        peakVal_cont = NaN;

        if any(negMask) && any(~isnan(medianCorr(negMask)))
            idxNegAll = find(negMask);
            [~, iNeg] = max(abs(medianCorr(negMask)));      % <-- ONLY negative lags
            idxNeg = idxNegAll(iNeg);
            [peakLag_cont, peakVal_cont] = quadInterpExtremum(lags_ms, medianCorr, idxNeg);
        end

        % ---------------- SAVE UNIT-LEVEL RESULTS ----------------
        XCorrData(day).units(uu).unitIdx      = unitIdx;
        XCorrData(day).units(uu).lags_ms      = lags_ms;
        XCorrData(day).units(uu).xcorr_trials = xcorr_trials;
        XCorrData(day).units(uu).medianCorr   = medianCorr;

        % Keep original field names, but they now mean: NEGATIVE-LAG abs-max extremum
        XCorrData(day).units(uu).peakVal      = peakVal_cont;
        XCorrData(day).units(uu).peakLag      = peakLag_cont;

        % Not used anymore (kept for compatibility)
        XCorrData(day).units(uu).minVal       = NaN;
        XCorrData(day).units(uu).minLag       = NaN;

        % ================= PEAK-TIME GATE (DO THIS BEFORE PLOTTING) =================
        psth = data.units(unitIdx).BLCSM_max_vx.PSTH_smoothbin;
        [~, peak_idx] = max(psth);

        % Gate: only keep units whose PSTH peak occurs within your valid window
        if ~ismember(peak_idx, valid_peak_range)
            XCorrData(day).units(uu).keep = false;

            % still store placeholders for new analysis
            XCorrData(day).units(uu).unitLag_abs_ms   = NaN;
            XCorrData(day).units(uu).CSMT2_win_ms     = [NaN NaN];
            XCorrData(day).units(uu).trial_peakFR_AD  = [];
            XCorrData(day).units(uu).trial_peakVel_AD = [];
            XCorrData(day).units(uu).r_FR_vs_Vel_AD   = NaN;
            XCorrData(day).units(uu).p_FR_vs_Vel_AD   = NaN;
            continue;
        end
        XCorrData(day).units(uu).keep = true;

        % If in range, classify before/after and store IDs
        if timevec(peak_idx) < 0
            before_idx_all(end+1) = unitIdx;
        else
            after_idx_all(end+1)  = unitIdx;
        end

        % ==========================================================================================
        % NEW ADD-ON (AFTER KEEP GATE):
        % Use abs(negative-lag) from BL analysis to define unit-specific window relative to CSMT2=0
        % Then compute per-trial peak FR (from ADCSM_csmT2) and per-trial peak outward vel
        % (from ExcitRR_AD reaches aligned to CSMT2).
        % ==========================================================================================

        unitLag_abs_ms = abs(peakLag_cont); % <-- your requested positive lag
        XCorrData(day).units(uu).unitLag_abs_ms = unitLag_abs_ms;

        % ---- AD spike raster aligned to CSMT2 ----
        if ~isfield(data.units(unitIdx),'ADCSM_csmT2') || ~isfield(data.units(unitIdx).ADCSM_csmT2,'Raster') || ...
                isempty(data.units(unitIdx).ADCSM_csmT2.Raster)
            fprintf('    UNIT %d: missing ADCSM_csmT2 raster -> skipping peak-window analysis.\n', unitIdx);

            XCorrData(day).units(uu).CSMT2_win_ms     = [NaN NaN];
            XCorrData(day).units(uu).trial_peakFR_AD  = [];
            XCorrData(day).units(uu).trial_peakVel_AD = [];
            XCorrData(day).units(uu).r_FR_vs_Vel_AD   = NaN;
            XCorrData(day).units(uu).p_FR_vs_Vel_AD   = NaN;
            continue;
        end

        spikes_AD = data.units(unitIdx).ADCSM_csmT2.Raster; % trials x time(1ms)
        [nTrials_AD, nTime_AD] = size(spikes_AD);
     
        full_time_AD_1ms = (-1000:1000)';
       
        % ---- IFR for AD (same method as BL) ----
        IFR_AD_full = zeros(nTime_AD, nTrials_AD);

        kern_len = round(6*gauss_sigma);
        gk = fspecial('gaussian',[1 kern_len], gauss_sigma);
        gk = gk / sum(gk);

        for tr = 1:nTrials_AD
            sp_idx = find(spikes_AD(tr,:) > 0);
            if numel(sp_idx) < 2, continue; end

            spike_times = full_time_AD_1ms(sp_idx);
            isis = diff(spike_times)/1000;
            inst_rate = 1 ./ isis;

            rate_vec = zeros(nTime_AD,1);
            for s = 1:length(isis)
                interval_idx = sp_idx(s):sp_idx(s+1)-1;
                rate_vec(interval_idx) = inst_rate(s);
            end

            IFR_AD_full(:,tr) = conv(rate_vec, gk, 'same');
        end

        % ---- restrict to [-win_pre_ms, +win_post_ms] around CSMT2, then resample to 200Hz ----
        win_idx_AD = full_time_AD_1ms >= -win_pre_ms & full_time_AD_1ms <= win_post_ms;
        IFR_AD_1ms = IFR_AD_full(win_idx_AD, :);
        ifr_time_AD_1ms = full_time_AD_1ms(win_idx_AD);

        time_200Hz = (-win_pre_ms:1000/resample_fs:win_post_ms)';
        nSamples200 = length(time_200Hz);

        IFR_AD = nan(nSamples200, nTrials_AD);
        for tr = 1:nTrials_AD
            IFR_AD(:,tr) = interp1(ifr_time_AD_1ms, IFR_AD_1ms(:,tr), time_200Hz, 'linear', 0);
        end

        % ---- AD kinematics (ExcitRR_AD) aligned to CSMT2 ----
        % if ~isfield(ExcitRR_AD,'ReReaches') || isempty(ExcitRR_AD.ReReaches)
        %     fprintf('    UNIT %d: ExcitRR_AD.ReReaches missing/empty -> skipping peak-window analysis.\n', unitIdx);
        % 
        %     XCorrData(day).units(uu).CSMT2_win_ms     = [NaN NaN];
        %     XCorrData(day).units(uu).trial_peakFR_AD  = [];
        %     XCorrData(day).units(uu).trial_peakVel_AD = [];
        %     XCorrData(day).units(uu).r_FR_vs_Vel_AD   = NaN;
        %     XCorrData(day).units(uu).p_FR_vs_Vel_AD   = NaN;
        %     continue;
        % end

        ReReaches_AD = ExcitRR_AD.ReReaches;

        % <<< REPLACE THIS >>> with your actual CSMT2 alignment function
        % It must output a struct compatible with stimvel_accel and C.stimvel like your BL code.
      

        thresh = 0.1;   % <<< set to your positional threshold used for the CSM definition
        dir    = 1;      % 1 = outreach (matches your CSMT2-style alignment in the original function)

        C_AD = AlignThreshCSMHF_1ms_edit(ReReaches_AD, thresh, dir, win_pre_ms, win_post_ms);

        C_AD = stimvel_accel(C_AD);
        vel_1ms_all_AD = squeeze(C_AD.stimvel(:,5,:)); % outward velocity
       % vel_1ms_all_AD = squeeze(C_AD.stimvel(:,8,:)); % outward velocity
        vel_time_1ms_AD = linspace(-win_pre_ms, win_post_ms, size(vel_1ms_all_AD,1))';

        nTrials_vel_AD = size(vel_1ms_all_AD,2);

        % match trials between AD spikes and AD kinematics
        nTrials_use = min(nTrials_AD, nTrials_vel_AD);

        Vel_AD = nan(nSamples200, nTrials_use);
        for tr = 1:nTrials_use
            Vel_AD(:,tr) = interp1(vel_time_1ms_AD, vel_1ms_all_AD(:,tr), time_200Hz, 'linear', NaN);
        end

        IFR_AD = IFR_AD(:,1:nTrials_use);

        % ---- UNIT-SPECIFIC WINDOW centered at +abs(lag) ms relative to CSMT2 (0) ----
       % CSMT2_win_ms = [unitLag_abs_ms - peakWin_half_ms, unitLag_abs_ms + peakWin_half_ms];
       %CSMT2_win_ms = [-unitLag_abs_ms + peakWin_half_ms, unitLag_abs_ms + peakWin_half_ms];
       CSMT2_win_ms = [-50 100];
       % CSMT2_win_ms = [-unitLag_abs_ms, unitLag_abs_ms ];
        idxW = (time_200Hz >= CSMT2_win_ms(1)) & (time_200Hz <= CSMT2_win_ms(2));

        peakFR_AD  = nan(nTrials_use,1);
        peakVel_AD = nan(nTrials_use,1);

        for tr = 1:nTrials_use
            frw = IFR_AD(idxW,tr);
            vw  = Vel_AD(idxW,tr);

            if all(isnan(frw)) || all(isnan(vw)), continue; end
            peakFR_AD(tr)  = max(frw);
            peakVel_AD(tr) = max(vw);
        end

        ok = ~isnan(peakFR_AD) & ~isnan(peakVel_AD);
        if nnz(ok) >= 1
            [R,P] = corrcoef(peakVel_AD(ok), peakFR_AD(ok));
            rPV = R(1,2);
            pPV = P(1,2);
        else
            rPV = NaN; pPV = NaN;
        end

        % ---- SAVE RESULTS ----
        XCorrData(day).units(uu).CSMT2_win_ms     = CSMT2_win_ms;
        XCorrData(day).units(uu).trial_peakFR_AD  = peakFR_AD;
        XCorrData(day).units(uu).trial_peakVel_AD = peakVel_AD;
        XCorrData(day).units(uu).r_FR_vs_Vel_AD   = rPV;
        XCorrData(day).units(uu).p_FR_vs_Vel_AD   = pPV;


        % ===================== PER-UNIT PLOT + SAVE (PEAK FR vs PEAK VEL) =====================

        %saveFolder = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\CSMMagvsPeakFR_inidivudalunits";
        %saveFolder="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\CSMMagvsPeakFR_inidivudalunits\Shifted Window";
       saveFolder="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\StimAnimals_BlockStim\CSMMagvsPeakFR_inidivudalunits\Window -50 to 100";
        if ~exist(saveFolder,'dir'), mkdir(saveFolder); end

        ok = ~isnan(peakFR_AD) & ~isnan(peakVel_AD);

        if nnz(ok) >= 1

            % --- linear fit (for R^2) ---
            x = peakVel_AD(ok);
            y = peakFR_AD(ok);

            pfit = polyfit(x, y, 1);
            yhat = polyval(pfit, x);

            SSres = sum((y - yhat).^2);
            SStot = sum((y - mean(y)).^2);
            R2 = 1 - (SSres / SStot);

            if isnan(rPV)
                R = corrcoef(x, y);
                rPV = R(1,2);
            end

            fig = figure('Color','w','Position',[100 100 520 420]); hold on

            scatter(x, y, 45, 'k', 'filled', 'MarkerFaceAlpha', 0.7);

            xlinefit = linspace(min(x), max(x), 100);
            ylinefit = polyval(pfit, xlinefit);
            plot(xlinefit, ylinefit, 'k-', 'LineWidth', 2);

            xlabel('Peak outward velocity (CSMT2 unit window, AD)')
            ylabel('Peak firing rate (CSMT2 unit window, AD)')
            title(sprintf('Sess %d | Unit %d | abs(negLag)=%.1f ms', day, unitIdx, unitLag_abs_ms), ...
                'Interpreter','none');

            txt = sprintf('r = %.2f\nR^2 = %.2f\np = %.3g\nn = %d', rPV, R2, pPV, nnz(ok));
            xlim_curr = xlim; ylim_curr = ylim;
            text(xlim_curr(1) + 0.05*range(xlim_curr), ...
                ylim_curr(2) - 0.10*range(ylim_curr), ...
                txt, 'FontSize', 11, 'VerticalAlignment','top');

            box off

            % ---- SAVE AS SVG (robust) ----
            outBase = fullfile(saveFolder, sprintf('Sess%02d_Unit%03d_absLag%0.1fms', day, unitIdx, unitLag_abs_ms));

            outSvg = outBase + ".svg";                 % explicit extension
            set(fig, 'Renderer', 'painters');          % vector renderer

            print(fig, char(outSvg), '-dsvg');         % IMPORTANT: char(...) not string
            close(fig);

        end


       


    end % END UNIT LOOP

    % ---------------- SAVE SESSION GROUPING -------------------
    XCorrData(day).before_peakv_units = before_idx_all;
    XCorrData(day).after_peakv_units  = after_idx_all;

    fprintf('Session %d summary: %d before, %d after\n', ...
        day, length(before_idx_all), length(after_idx_all));
end

fprintf('\n================ ALL SESSIONS COMPLETE ================\n');

%% ========================= HELPER (put once at end of script) =========================
function [lag_cont, val_cont] = quadInterpExtremum(x, y, idx)
% Quadratic interpolation around idx (idx-1,idx,idx+1). Falls back if edge/bad fit.
    lag_cont = x(idx);
    val_cont = y(idx);

    if idx <= 1 || idx >= numel(y) || any(isnan(y(idx-1:idx+1)))
        return;
    end

    xx = x(idx-1:idx+1);
    yy = y(idx-1:idx+1);

    p = polyfit(xx, yy, 2);
    if abs(p(1)) < eps
        return;
    end

    lag0 = -p(2)/(2*p(1));
    if lag0 < min(xx) || lag0 > max(xx)
        return;
    end

    lag_cont = lag0;
    val_cont = polyval(p, lag_cont);
end
%%  Different ode 


%% ===================== PAIRED PER-UNIT (WORKING) + NEG-LAG POS/NEG CORR CATEGORIES =====================
% Keeps your paired analysis EXACTLY the same (ADCSM_pos_thresh vs eventopto; keep gate [-250,250]).
% NEW: For each kept unit, compute BL-based neg-lag abs-max cross-corr and classify:
%      posCorr (peakVal>0) vs negCorr (peakVal<0)
% Then: make paired plots split by category + paired t-tests per category.

clear; clc;

% ---------------------- LOAD SESSIONS -----------------------
cd("C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\StimAnimals\BlockSeparated_CSMvsNoCSM_inducedalignment")
MyFolderInfo = dir('*Block_CSMvsNoCSM.mat');

% ---------------------- PAIRED ANALYSIS SETTINGS (UNCHANGED) -----------------------
stimOn_ms    = 0;
%baseWin_ms   = [-50 0];
baseWin_ms   = [-500 -250];
stimWin_ms   = [stimOn_ms stimOn_ms+50];

sd_sigma_ms  = 5;
sd_kern_ms   = 6*sd_sigma_ms;

% Keep gate window for CONTROL PSTH peak
keepWin_ms   = [-250 250];

% ---------------------- NEG-LAG XCORR SETTINGS (FROM YOUR SCRIPT) -----------------------
kinVar        = 7;
gauss_sigma   = 10;
resample_fs   = 200;
win_pre_ms    = 100;
win_post_ms   = 150;
maxLag_ms     = 100;
lag_step      = 5;

%% ---------------------- OUTPUT ----------------------
behResp = [];
tagResp = [];
catSign = [];   % +1 = posCorr, -1 = negCorr
peakVal = [];
peakLag = [];
labels  = strings(0,1);

for day = 1:length(MyFolderInfo)

    fprintf('\n========== SESSION %d / %d ==========\n', day, length(MyFolderInfo));
    load(MyFolderInfo(day).name, 'data','ExcitRR_BL','Salt_index');

    if isempty(Salt_index)
        fprintf('Day %d: Salt_index empty. Skipping...\n', day);
        continue;
    end

    session_units = Salt_index(:);

    % === keep gate setup (uses BLCSM_max_vx, consistent with your pipelines) ===
    timevec = data.units(session_units(1)).BLCSM_max_vx.time_bins;
    valid_peak_range = find(timevec >= keepWin_ms(1) & timevec <= keepWin_ms(2));

    for uu = 1:numel(session_units)

        unitIdx = session_units(uu);
        u = data.units(unitIdx);

        % ===================== KEEP GATE (UNCHANGED) =====================
        if ~isfield(u,'BLCSM_max_vx') || isempty(u.BLCSM_max_vx.PSTH_smoothbin), continue; end
        psth_ctrl = u.BLCSM_max_vx.PSTH_smoothbin(:);
        [~, peak_idx] = max(psth_ctrl);
        if ~ismember(peak_idx, valid_peak_range), continue; end

        % ===================== REQUIRE THE TWO RASTERS (UNCHANGED) =====================
        if ~isfield(u,'ADCSM_pos_thresh') || isempty(u.ADCSM_pos_thresh.Raster), continue; end
        if ~isfield(u,'eventopto')        || isempty(u.eventopto.Raster),        continue; end

        R_beh = u.ADCSM_pos_thresh.Raster;
        R_tag = u.eventopto.Raster;

        % ===================== PAIRED METRICS (UNCHANGED) =====================
        b = computeResp(R_beh, baseWin_ms, stimWin_ms, sd_sigma_ms, sd_kern_ms);
        t = computeResp(R_tag, baseWin_ms, stimWin_ms, sd_sigma_ms, sd_kern_ms);
        if isnan(b) || isnan(t), continue; end

        % ===================== NEW: NEG-LAG XCORR CLASSIFICATION (BL) =====================
        % Uses your exact approach: BLCSM_max_vx raster + BL reaches aligned to peakV.
        [peakLag_cont, peakVal_cont] = negLagXCorr_oneUnit(u, ExcitRR_BL, kinVar, gauss_sigma, ...
            resample_fs, win_pre_ms, win_post_ms, maxLag_ms, lag_step);

        if isnan(peakVal_cont)
            continue;
        end

        behResp(end+1,1) = b; %#ok<AGROW>
        tagResp(end+1,1) = t; %#ok<AGROW>

        peakLag(end+1,1) = peakLag_cont; %#ok<AGROW>
        peakVal(end+1,1) = peakVal_cont; %#ok<AGROW>

        catSign(end+1,1) = sign(peakVal_cont); %#ok<AGROW>  % +1 posCorr, -1 negCorr, 0 if exactly 0

        labels(end+1,1)  = "S"+day+"_U"+unitIdx; %#ok<AGROW>
    end
end

fprintf('\nKept paired units with BL neg-lag category: %d\n', numel(behResp));
if isempty(behResp)
    error('No units survived keep + paired rasters + neg-lag xcorr.');
end

% Remove any catSign==0 (rare)
okCat = catSign ~= 0;
behResp = behResp(okCat);
tagResp = tagResp(okCat);
catSign = catSign(okCat);
peakLag = peakLag(okCat);
peakVal = peakVal(okCat);
labels  = labels(okCat);

%% ===================== CONNECTED PLOT SPLIT BY CATEGORY =====================
isPos = catSign > 0;
isNeg = catSign < 0;

figure('Color','w','Position',[200 200 900 420]);

subplot(1,2,1); hold on
plotConnected(behResp(isPos), tagResp(isPos));
title(sprintf('posCorr (BL neg-lag peakVal>0), N=%d', nnz(isPos)));
set(gca,'XTick',[1 2],'XTickLabel',{'Behavior','Optotag'}); xlim([0.7 2.3]);
ylabel('Peak FR (0–50) − baseline (Hz)'); box off

subplot(1,2,2); hold on
plotConnected(behResp(isNeg), tagResp(isNeg));
title(sprintf('negCorr (BL neg-lag peakVal<0), N=%d', nnz(isNeg)));
set(gca,'XTick',[1 2],'XTickLabel',{'Behavior','Optotag'}); xlim([0.7 2.3]);
ylabel('Peak FR (0–50) − baseline (Hz)'); box off


%% ===================== CONNECTED PLOT SPLIT BY CATEGORY (MEAN CIRCLES) =====================
isPos = catSign > 0;
isNeg = catSign < 0;

figure('Color','w','Position',[200 200 900 420]);

% -------- POSITIVE CORR --------
subplot(1,2,1); hold on
plotConnected(behResp(isPos), tagResp(isPos));

% --- means ---
mB = mean(behResp(isPos),'omitnan');
mT = mean(tagResp(isPos),'omitnan');

% --- large circles ---
scatter(1, mB, 140, 'k', 'filled', 'MarkerFaceAlpha',1);
scatter(2, mT, 140, 'k', 'filled', 'MarkerFaceAlpha',1);

% --- connecting line ---
plot([1 2], [mB mT], 'k-', 'LineWidth', 3);

title(sprintf('posCorr (BL neg-lag peakVal>0), N=%d', nnz(isPos)));
set(gca,'XTick',[1 2],'XTickLabel',{'Behavior','Optotag'});
xlim([0.7 2.3]);
ylabel('Peak FR (0–50) − baseline (Hz)');
box off


% -------- NEGATIVE CORR --------
subplot(1,2,2); hold on
plotConnected(behResp(isNeg), tagResp(isNeg));

mB = mean(behResp(isNeg),'omitnan');
mT = mean(tagResp(isNeg),'omitnan');

scatter(1, mB, 140, 'k', 'filled', 'MarkerFaceAlpha',1);
scatter(2, mT, 140, 'k', 'filled', 'MarkerFaceAlpha',1);

plot([1 2], [mB mT], 'k-', 'LineWidth', 3);

title(sprintf('negCorr (BL neg-lag peakVal<0), N=%d', nnz(isNeg)));
set(gca,'XTick',[1 2],'XTickLabel',{'Behavior','Optotag'});
xlim([0.7 2.3]);
ylabel('Peak FR (0–50) − baseline (Hz)');
box off

%% ===================== CONNECTED PLOT SPLIT BY CATEGORY (VERY COMPRESSED) =====================
% Opto on LEFT, Reach-related (Behavior) on RIGHT

isPos = catSign > 0;
isNeg = catSign < 0;

% <<< tighter spacing
xOpto  = 1.00;   % LEFT
xReach = 1.18;   % RIGHT

figure('Color','w','Position',[200 200 600 420]);

% -------- POSITIVE CORR --------
subplot(1,2,1); hold on

% swap order so opto plotted first
plotConnected_customX(tagResp(isPos), behResp(isPos), xOpto, xReach);

mOpto  = mean(tagResp(isPos),'omitnan');
mReach = mean(behResp(isPos),'omitnan');

scatter(xOpto,  mOpto,  140, 'k', 'filled');
scatter(xReach, mReach, 140, 'k', 'filled');
plot([xOpto xReach], [mOpto mReach], 'k-', 'LineWidth', 3);

title(sprintf('posCorr (BL neg-lag peakVal>0), N=%d', nnz(isPos)));
set(gca,'XTick',[xOpto xReach],'XTickLabel',{'Optotag','Reach-related'});
xlim([0.95 1.23]);
ylabel('Peak FR (0–50) − baseline (Hz)');
box off


% -------- NEGATIVE CORR --------
subplot(1,2,2); hold on

plotConnected_customX(tagResp(isNeg), behResp(isNeg), xOpto, xReach);

mOpto  = mean(tagResp(isNeg),'omitnan');
mReach = mean(behResp(isNeg),'omitnan');

scatter(xOpto,  mOpto,  140, 'k', 'filled');
scatter(xReach, mReach, 140, 'k', 'filled');
plot([xOpto xReach], [mOpto mReach], 'k-', 'LineWidth', 3);

title(sprintf('negCorr (BL neg-lag peakVal<0), N=%d', nnz(isNeg)));
set(gca,'XTick',[xOpto xReach],'XTickLabel',{'Optotag','Reach-related'});
xlim([0.95 1.23]);
ylabel('Peak FR (0–50) − baseline (Hz)');
box off

%% Box plot version 
%% ===================== BOX + CONNECTED LINES (BACKGROUND) =====================

isPos = catSign > 0;
isNeg = catSign < 0;

% spacing
xOpto  = 1.00;
xReach = 1.18;

figure('Color','w','Position',[200 200 600 420]);

% -------- POSITIVE CORR --------
subplot(1,2,1); hold on

% ---- background paired lines ----
plotConnected_customX(tagResp(isPos), behResp(isPos), xOpto, xReach);

% ---- boxplot data ----
data = [tagResp(isPos); behResp(isPos)];
group = [ones(sum(isPos),1); 2*ones(sum(isPos),1)];

boxplot(data, group, ...
    'Positions',[xOpto xReach], ...
    'Widths',0.05, ...
    'Colors','k', ...
    'Symbol','');

yline(0)

title(sprintf('posCorr (BL neg-lag peakVal>0), N=%d', nnz(isPos)));
set(gca,'XTick',[xOpto xReach],'XTickLabel',{'Optotag','Reach-related'});
xlim([0.95 1.23]);
ylabel('Peak FR (0–50) − baseline (Hz)');
box off


% -------- NEGATIVE CORR --------
subplot(1,2,2); hold on

% ---- background paired lines ----
plotConnected_customX(tagResp(isNeg), behResp(isNeg), xOpto, xReach);

% ---- boxplot data ----
data = [tagResp(isNeg); behResp(isNeg)];
group = [ones(sum(isNeg),1); 2*ones(sum(isNeg),1)];

boxplot(data, group, ...
    'Positions',[xOpto xReach], ...
    'Widths',0.05, ...
    'Colors','k', ...
    'Symbol','');

yline(0)

title(sprintf('negCorr (BL neg-lag peakVal<0), N=%d', nnz(isNeg)));
set(gca,'XTick',[xOpto xReach],'XTickLabel',{'Optotag','Reach-related'});
xlim([0.95 1.23]);
ylabel('Peak FR (0–50) − baseline (Hz)');
box off
%% ===================== CONNECTED PLOT SPLIT BY CATEGORY (VERY COMPRESSED) =====================
% Opto on LEFT, Reach-related on RIGHT
% Grand means shown as mean +/- SEM (no connecting mean line)

isPos = catSign > 0;
isNeg = catSign < 0;

% <<< tighter spacing
xOpto  = 1.00;   % LEFT
xReach = 1.18;   % RIGHT

figure('Color','w','Position',[200 200 600 420]);

% -------- POSITIVE CORR --------
subplot(1,2,1); hold on

% individual paired lines
plotConnected_customX(tagResp(isPos), behResp(isPos), xOpto, xReach);

% means
mOpto  = mean(tagResp(isPos), 'omitnan');
mReach = mean(behResp(isPos), 'omitnan');

% SEMs
semOpto  = std(tagResp(isPos),  'omitnan') ./ sqrt(sum(~isnan(tagResp(isPos))));
semReach = std(behResp(isPos), 'omitnan') ./ sqrt(sum(~isnan(behResp(isPos))));

% mean +/- SEM
errorbar(xOpto,  mOpto,  semOpto,  'ko', ...
    'MarkerFaceColor','k', 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
errorbar(xReach, mReach, semReach, 'ko', ...
    'MarkerFaceColor','k', 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);

yline(0)
title(sprintf('posCorr (BL neg-lag peakVal>0), N=%d', nnz(isPos)));
set(gca,'XTick',[xOpto xReach],'XTickLabel',{'Optotag','Reach-related'});
xlim([0.95 1.23]);
ylabel('Peak FR (0–50) − baseline (Hz)');
box off


% -------- NEGATIVE CORR --------
subplot(1,2,2); hold on

% individual paired lines
plotConnected_customX(tagResp(isNeg), behResp(isNeg), xOpto, xReach);

% means
mOpto  = mean(tagResp(isNeg), 'omitnan');
mReach = mean(behResp(isNeg), 'omitnan');

% SEMs
semOpto  = std(tagResp(isNeg),  'omitnan') ./ sqrt(sum(~isnan(tagResp(isNeg))));
semReach = std(behResp(isNeg), 'omitnan') ./ sqrt(sum(~isnan(behResp(isNeg))));

% mean +/- SEM
errorbar(xOpto,  mOpto,  semOpto,  'ko', ...
    'MarkerFaceColor','k', 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
errorbar(xReach, mReach, semReach, 'ko', ...
    'MarkerFaceColor','k', 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);

yline(0)
title(sprintf('negCorr (BL neg-lag peakVal<0), N=%d', nnz(isNeg)));
set(gca,'XTick',[xOpto xReach],'XTickLabel',{'Optotag','Reach-related'});
xlim([0.95 1.23]);
ylabel('Peak FR (0–50) − baseline (Hz)');
box off


%% FIgure edit code 
% ===================== CONNECTED PLOT SPLIT BY CATEGORY (VERY COMPRESSED) =====================
isPos = catSign > 0;
isNeg = catSign < 0;

% <<< tighter spacing
xB = 1.00;
xT = 1.18;

figure('Color','w','Position',[200 200 600 420]);

% -------- POSITIVE CORR --------
subplot(1,2,1); hold on
plotConnected_customX(behResp(isPos), tagResp(isPos), xB, xT);

mB = mean(behResp(isPos),'omitnan');
mT = mean(tagResp(isPos),'omitnan');

scatter(xB, mB, 140, 'k', 'filled');
scatter(xT, mT, 140, 'k', 'filled');
plot([xB xT], [mB mT], 'k-', 'LineWidth', 3);

title(sprintf('posCorr (BL neg-lag peakVal>0), N=%d', nnz(isPos)));
set(gca,'XTick',[xB xT],'XTickLabel',{'Behavior','Optotag'});
xlim([0.95 1.23]);   % <<< tight axis
ylabel('Peak FR (0–50) − baseline (Hz)');
box off


% -------- NEGATIVE CORR --------
subplot(1,2,2); hold on
plotConnected_customX(behResp(isNeg), tagResp(isNeg), xB, xT);

mB = mean(behResp(isNeg),'omitnan');
mT = mean(tagResp(isNeg),'omitnan');

scatter(xB, mB, 140, 'k', 'filled');
scatter(xT, mT, 140, 'k', 'filled');
plot([xB xT], [mB mT], 'k-', 'LineWidth', 3);

title(sprintf('negCorr (BL neg-lag peakVal<0), N=%d', nnz(isNeg)));
set(gca,'XTick',[xB xT],'XTickLabel',{'Behavior','Optotag'});
xlim([0.95 1.23]);
ylabel('Peak FR (0–50) − baseline (Hz)');
box off
%%
% %%
% function plotConnected_customX(a, b, x1, x2)
% 
% for i = 1:numel(a)
%     plot([x1 x2], [a(i) b(i)], '-', 'Color', [0 0 0 .2]);
% end
% 
% scatter(x1*ones(size(a)), a, 35, 'k', 'filled', 'MarkerFaceAlpha',0.65);
% scatter(x2*ones(size(b)), b, 35, 'k', 'filled', 'MarkerFaceAlpha',0.65);
% 
% end

function plotConnected_customX(y1, y2, x1, x2)

n = length(y1);

for i = 1:n

    if isnan(y1(i)) || isnan(y2(i))
        continue
    end

    % line color
    lineColor = [0.6 0.6 0.6];   % light gray (change if you want)

    % connecting line
    plot([x1 x2], [y1(i) y2(i)], '-', ...
        'Color', lineColor, ...
        'LineWidth', 1);

    % small circles (same color as line)
    scatter(x1, y1(i), 18, ...
        'MarkerEdgeColor', lineColor, ...
        'MarkerFaceColor', lineColor);

    scatter(x2, y2(i), 18, ...
        'MarkerEdgeColor', lineColor, ...
        'MarkerFaceColor', lineColor);

end
end

%% ===================== STATS: PAIRED t-test (overall + per category) =====================
fprintf('\n===== PAIRED t-TESTS =====\n');

[~,pAll,ciAll,stAll] = ttest(tagResp, behResp);
fprintf('ALL:    t(%d)=%.3f, p=%.4g, CI=[%.3f %.3f], meanΔ=%.3f\n', stAll.df, stAll.tstat, pAll, ciAll(1), ciAll(2), mean(tagResp-behResp));

if nnz(isPos) >= 3
    [~,pPos,ciPos,stPos] = ttest(tagResp(isPos), behResp(isPos));
    fprintf('posCorr t(%d)=%.3f, p=%.4g, CI=[%.3f %.3f], meanΔ=%.3f\n', stPos.df, stPos.tstat, pPos, ciPos(1), ciPos(2), mean(tagResp(isPos)-behResp(isPos)));
else
    fprintf('posCorr: too few units for t-test\n');
end

if nnz(isNeg) >= 3
    [~,pNeg,ciNeg,stNeg] = ttest(tagResp(isNeg), behResp(isNeg));
    fprintf('negCorr t(%d)=%.3f, p=%.4g, CI=[%.3f %.3f], meanΔ=%.3f\n', stNeg.df, stNeg.tstat, pNeg, ciNeg(1), ciNeg(2), mean(tagResp(isNeg)-behResp(isNeg)));
else
    fprintf('negCorr: too few units for t-test\n');
end

%% ===================== OPTIONAL: QUICK QC SCATTER peakLag vs peakVal =====================
figure('Color','w','Position',[1150 200 420 350]); hold on
scatter(peakLag(isPos), peakVal(isPos), 35, 'k', 'filled', 'MarkerFaceAlpha',0.55);
scatter(peakLag(isNeg), peakVal(isNeg), 35, 'k', 'filled', 'MarkerFaceAlpha',0.55);
xline(0,'--k'); yline(0,'--k');
xlabel('BL neg-lag peakLag (ms)'); ylabel('BL neg-lag peakVal');
title('BL classification inputs (neg-lag abs-max)');
box off

%% ===================== HELPERS =====================
function resp = computeResp(raster, baseWin_ms, stimWin_ms, sd_sigma_ms, sd_kern_ms)
% Peak FR in stim window minus mean baseline (mean PSTH spike density)
    [~, nT] = size(raster);

    if nT == 2001
        t = (-1000:1000)';
    else
        mid = ceil(nT/2);
        t = ((1:nT)' - mid);
    end

    idxBase = (t >= baseWin_ms(1)) & (t < baseWin_ms(2));
    idxStim = (t >= stimWin_ms(1)) & (t <= stimWin_ms(2));
    if nnz(idxBase) < 5 || nnz(idxStim) < 5
        resp = NaN; return;
    end

    g = gaussKernel1D(sd_sigma_ms, sd_kern_ms);
    m = mean(double(raster),1);            % spikes/ms
    fr = conv(m, g, 'same') * 1000;        % Hz

    resp = max(fr(idxStim)) - mean(fr(idxBase));
end

function [peakLag_cont, peakVal_cont] = negLagXCorr_oneUnit(u, ExcitRR_BL, kinVar, gauss_sigma, resample_fs, win_pre_ms, win_post_ms, maxLag_ms, lag_step)
% Your BL neg-lag abs-max cross-corr, wrapped as a function.
% Returns quad-interpolated extremum on medianCorr, restricted to negative lags.

    peakLag_cont = NaN; peakVal_cont = NaN;

    if ~isfield(u,'BLCSM_max_vx') || ~isfield(u.BLCSM_max_vx,'Raster') || isempty(u.BLCSM_max_vx.Raster)
        return;
    end

    spikes_full = u.BLCSM_max_vx.Raster; % [trial x time(1ms)]
    [nTrials, nTime] = size(spikes_full);
    full_time_1ms = (-1000:1000)';

    % IFR (same as your script)
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

    % restrict
    win_idx = full_time_1ms >= -win_pre_ms & full_time_1ms <= win_post_ms;
    IFR_1ms = IFR_full(win_idx, :);
    ifr_time_1ms = full_time_1ms(win_idx);

    % resample
    time_200Hz = (-win_pre_ms:1000/resample_fs:win_post_ms)'; % ms
    nSamples200 = length(time_200Hz);

    IFR = nan(nSamples200, nTrials);
    for tr = 1:nTrials
        IFR(:,tr) = interp1(ifr_time_1ms, IFR_1ms(:,tr), time_200Hz, 'linear', 0);
    end

    % BL kinematics aligned to peakV (same as you)
    ReReaches = ExcitRR_BL.ReReaches;
    C = AlignPeakV_1ms_edit(ReReaches, kinVar, win_pre_ms, win_post_ms);
    C = stimvel_accel(C);
    vel_1ms_all = squeeze(C.stimvel(:,5,:)); % outward velocity
    vel_time_1ms = linspace(-win_pre_ms, win_post_ms, size(vel_1ms_all,1))';

    Vel = nan(nSamples200, nTrials);
    for tr = 1:nTrials
        Vel(:,tr) = interp1(vel_time_1ms, vel_1ms_all(:,tr), time_200Hz, 'linear', NaN);
    end

    % sliding xcorr
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

            sd_fr = std(fr_seg); sd_v = std(v_seg);
            if sd_fr < 1e-9 || sd_v < 1e-9, continue; end

            fr_z = (fr_seg - mean(fr_seg)) ./ sd_fr;
            v_z  = (v_seg  - mean(v_seg )) ./ sd_v;

            if any(isnan(fr_z)) || any(isnan(v_z)), continue; end
            xcorr_trials(li,tr) = mean(fr_z .* v_z);
        end
    end

    medianCorr = nanmedian(xcorr_trials,2);

    negMask = lags_ms < 0;
    if ~any(negMask) || all(isnan(medianCorr(negMask)))
        return;
    end

    idxNegAll = find(negMask);
    [~, iNeg] = max(abs(medianCorr(negMask)));
    idxNeg = idxNegAll(iNeg);

    [peakLag_cont, peakVal_cont] = quadInterpExtremum(lags_ms, medianCorr, idxNeg);

    if isnan(peakLag_cont) || isnan(peakVal_cont)
        peakLag_cont = lags_ms(idxNeg);
        peakVal_cont = medianCorr(idxNeg);
    end
end

function plotConnected(a, b)
% connected paired plot helper
    if isempty(a) || isempty(b)
        text(0.5,0.5,'NO UNITS','Units','normalized','HorizontalAlignment','center'); return;
    end
    for i=1:numel(a)
        plot([1 2],[a(i) b(i)],'-','Color',[0 0 0 .2]);
    end
    scatter(ones(size(a)), a, 40, 'k', 'filled', 'MarkerFaceAlpha',0.65);
    scatter(2*ones(size(b)), b, 40, 'k', 'filled', 'MarkerFaceAlpha',0.65);
end



function g = gaussKernel1D(sigma_ms, kern_ms)
% Normalized 1D Gaussian kernel, odd length
    L = max(3, round(kern_ms));
    if mod(L,2)==0, L=L+1; end
    x = (-(L-1)/2:(L-1)/2);
    g = exp(-(x.^2)/(2*sigma_ms^2));
    g = g / sum(g);
end


%% DIFF METHOD - looking at the delta 

%% ===================== PAIRED PER-UNIT + NEG-LAG POS/NEG CORR CATEGORIES =====================
% Metric definitions:
%   tagResp = eventopto stim FR - eventopto prestim baseline FR
%   behResp = AD stim FR - BL stim FR
%
% Keeps:
%   - same keep gate
%   - same neg-lag xcorr classification
%   - same posCorr / negCorr split

clear; clc;

% ---------------------- LOAD SESSIONS -----------------------
cd("C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\StimAnimals\BlockSeparated_CSMvsNoCSM_inducedalignment")
MyFolderInfo = dir('*Block_CSMvsNoCSM.mat');

% ---------------------- ANALYSIS SETTINGS -----------------------
stimOn_ms    = 0;
baseWin_ms   = [-50 0];
stimWin_ms   = [stimOn_ms stimOn_ms+50];

sd_sigma_ms  = 5;
sd_kern_ms   = 6*sd_sigma_ms;

% Keep gate window for CONTROL PSTH peak
keepWin_ms   = [-250 250];

% ---------------------- NEG-LAG XCORR SETTINGS -----------------------
kinVar        = 7;
gauss_sigma   = 10;
resample_fs   = 200;
win_pre_ms    = 100;
win_post_ms   = 150;
maxLag_ms     = 100;
lag_step      = 5;

%% ---------------------- OUTPUT ----------------------
behResp = [];   % AD stim FR - BL stim FR
tagResp = [];   % eventopto stim FR - its own prestim baseline
catSign = [];   % +1 = posCorr, -1 = negCorr
peakVal = [];
peakLag = [];
labels  = strings(0,1);

for day = 1:length(MyFolderInfo)

    fprintf('\n========== SESSION %d / %d ==========\n', day, length(MyFolderInfo));
    load(MyFolderInfo(day).name, 'data','ExcitRR_BL','Salt_index');

    if isempty(Salt_index)
        fprintf('Day %d: Salt_index empty. Skipping...\n', day);
        continue;
    end

    session_units = Salt_index(:);

    % === keep gate setup ===
    timevec = data.units(session_units(1)).BLCSM_max_vx.time_bins;
    valid_peak_range = find(timevec >= keepWin_ms(1) & timevec <= keepWin_ms(2));

    for uu = 1:numel(session_units)

        unitIdx = session_units(uu);
        u = data.units(unitIdx);

        % ===================== KEEP GATE =====================
        if ~isfield(u,'BLCSM_max_vx') || isempty(u.BLCSM_max_vx.PSTH_smoothbin), continue; end
        psth_ctrl = u.BLCSM_max_vx.PSTH_smoothbin(:);
        [~, peak_idx] = max(psth_ctrl);
        if ~ismember(peak_idx, valid_peak_range), continue; end

        % ===================== REQUIRE RASTERS =====================
        if ~isfield(u,'ADCSM_pos_thresh') || isempty(u.ADCSM_pos_thresh.Raster), continue; end
        if ~isfield(u,'BLCSM_pos_thresh') || isempty(u.BLCSM_pos_thresh.Raster), continue; end
        if ~isfield(u,'eventopto')        || isempty(u.eventopto.Raster),        continue; end

        R_AD  = u.ADCSM_pos_thresh.Raster;
        R_BL  = u.BLCSM_pos_thresh.Raster;   % <-- change this field name if needed
        R_tag = u.eventopto.Raster;

        % ===================== METRICS =====================
        % eventopto: show increase relative to its own prestim baseline
        [tag_stimFR, tag_baseFR] = computeFRs(R_tag, baseWin_ms, stimWin_ms, sd_sigma_ms, sd_kern_ms);
        if isnan(tag_stimFR) || isnan(tag_baseFR), continue; end
        t = tag_stimFR - tag_baseFR;

        % behavioral comparison: AD stim FR - BL stim FR
        [ad_stimFR, ~] = computeFRs(R_AD, baseWin_ms, stimWin_ms, sd_sigma_ms, sd_kern_ms);
        [bl_stimFR, ~] = computeFRs(R_BL, baseWin_ms, stimWin_ms, sd_sigma_ms, sd_kern_ms);
        if isnan(ad_stimFR) || isnan(bl_stimFR), continue; end
        b = ad_stimFR - bl_stimFR;

        % ===================== NEG-LAG XCORR CLASSIFICATION (BL) =====================
        [peakLag_cont, peakVal_cont] = negLagXCorr_oneUnit(u, ExcitRR_BL, kinVar, gauss_sigma, ...
            resample_fs, win_pre_ms, win_post_ms, maxLag_ms, lag_step);

        if isnan(peakVal_cont)
            continue;
        end

        behResp(end+1,1) = b; %#ok<AGROW>
        tagResp(end+1,1) = t; %#ok<AGROW>

        peakLag(end+1,1) = peakLag_cont; %#ok<AGROW>
        peakVal(end+1,1) = peakVal_cont; %#ok<AGROW>

        catSign(end+1,1) = sign(peakVal_cont); %#ok<AGROW>
        labels(end+1,1)  = "S"+day+"_U"+unitIdx; %#ok<AGROW>
    end
end

fprintf('\nKept paired units with BL neg-lag category: %d\n', numel(behResp));
if isempty(behResp)
    error('No units survived keep + required rasters + neg-lag xcorr.');
end

% Remove any catSign==0
okCat = catSign ~= 0;
behResp = behResp(okCat);
tagResp = tagResp(okCat);
catSign = catSign(okCat);
peakLag = peakLag(okCat);
peakVal = peakVal(okCat);
labels  = labels(okCat);

%% ===================== OPTIONAL: QUICK CATEGORY COUNTS =====================
fprintf('posCorr units: %d\n', sum(catSign > 0));
fprintf('negCorr units: %d\n', sum(catSign < 0));

%% ===================== CONNECTED PLOT SPLIT BY CATEGORY =====================
% Behavior = AD stim FR - BL stim FR
% Event opto = eventopto stim FR - prestim baseline FR

isPos = catSign > 0;
isNeg = catSign < 0;

xB = 1.00;
xT = 1.18;

% shared y-limits across both panels
allVals = [behResp(:); tagResp(:)];
allVals = allVals(~isnan(allVals));

if isempty(allVals)
    error('No valid values to plot.');
end

yl = [min(allVals) max(allVals)];
if diff(yl) == 0
    yl = yl + [-1 1];
end
pad = 0.08 * diff(yl);
yl = yl + [-pad pad];

figure('Color','w','Position',[200 200 600 420]);

% -------- POSITIVE CORR --------
subplot(1,2,1); hold on

plotConnected_customX(behResp(isPos), tagResp(isPos), xB, xT);

mB = mean(behResp(isPos),'omitnan');
mT = mean(tagResp(isPos),'omitnan');

scatter(xB, mB, 140, 'k', 'filled');
scatter(xT, mT, 140, 'k', 'filled');
plot([xB xT], [mB mT], 'k-', 'LineWidth', 3);

title(sprintf('posCorr (N = %d)', nnz(isPos)));
set(gca,'XTick',[xB xT], ...
        'XTickLabel',{'AD - BL','Event opto'});
xlim([0.95 1.23]);
ylim(yl);
yline(0,'k:');
ylabel('\Delta firing rate (Hz)');
box off

% -------- NEGATIVE CORR --------
subplot(1,2,2); hold on

plotConnected_customX(behResp(isNeg), tagResp(isNeg), xB, xT);

mB = mean(behResp(isNeg),'omitnan');
mT = mean(tagResp(isNeg),'omitnan');

scatter(xB, mB, 140, 'k', 'filled');
scatter(xT, mT, 140, 'k', 'filled');
plot([xB xT], [mB mT], 'k-', 'LineWidth', 3);

title(sprintf('negCorr (N = %d)', nnz(isNeg)));
set(gca,'XTick',[xB xT], ...
        'XTickLabel',{'AD - BL','Event opto'});
xlim([0.95 1.23]);
ylim(yl);
yline(0,'k:');
ylabel('\Delta firing rate (Hz)');
box off

% ===================== OPTIONAL STATS =====================
[~,p_pos] = ttest(behResp(isPos), tagResp(isPos));
[~,p_neg] = ttest(behResp(isNeg), tagResp(isNeg));

fprintf('posCorr paired t-test: p = %.4g\n', p_pos);
fprintf('negCorr paired t-test: p = %.4g\n', p_neg);




%% ===================== HELPER =====================
function [stimFR, baseFR] = computeFRs(R, baseWin_ms, stimWin_ms, sd_sigma_ms, sd_kern_ms)

    if isempty(R)
        stimFR = NaN;
        baseFR = NaN;
        return;
    end

    % assumes columns are ms and time 0 is centered in raster
    nTime = size(R,2);
    tvec = ((1:nTime) - ceil(nTime/2));

    psth = mean(R,1) * 1000;  % spikes/ms -> Hz

    % Gaussian smoothing
    kern_x = -sd_kern_ms:sd_kern_ms;
    g = exp(-(kern_x.^2)/(2*sd_sigma_ms^2));
    g = g / sum(g);
    psth_sm = conv(psth, g, 'same');

    base_idx = tvec >= baseWin_ms(1) & tvec <  baseWin_ms(2);
    stim_idx = tvec >= stimWin_ms(1) & tvec <  stimWin_ms(2);

    if ~any(base_idx) || ~any(stim_idx)
        stimFR = NaN;
        baseFR = NaN;
        return;
    end

    baseFR = mean(psth_sm(base_idx), 'omitnan');
    stimFR = mean(psth_sm(stim_idx), 'omitnan');
end

%% New code edit- taking the delta 


%% ===================== PAIRED PER-UNIT: STIM EFFECT IN HZ =====================
% GOAL:
% Compare how much firing rate increases during the 50 ms stim epoch:
%
%   1) NON-REACH / EVENT OPTO:
%      tagResp = stim FR during 0-50 ms
%              - baseline FR from large windows outside stim
%
%   2) DURING BEHAVIOR:
%      behResp = AD stim FR during 0-50 ms
%              - BL stim FR during 0-50 ms
%
% Output units:
%   DELTA FIRING RATE (Hz)
%
% Keeps:
%   - same keep gate
%   - same BL neg-lag xcorr classification
%   - same posCorr / negCorr split

clear; clc;

% ---------------------- LOAD SESSIONS -----------------------
cd("C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\StimAnimals\BlockSeparated_CSMvsNoCSM_inducedalignment")
MyFolderInfo = dir('*Block_CSMvsNoCSM.mat');

%% ---------------------- SETTINGS -----------------------
stimWin_ms      = [0 50];        % 50 ms stim window

% Large baseline windows for eventopto, outside the stim epoch
tagBaseWin1_ms  = [-200 -50];
tagBaseWin2_ms  = [150 250];

% Keep gate window for BL control PSTH peak
keepWin_ms      = [-250 250];

% Neg-lag xcorr settings
kinVar        = 7;
gauss_sigma   = 10;
resample_fs   = 200;
win_pre_ms    = 100;
win_post_ms   = 150;
maxLag_ms     = 100;
lag_step      = 5;

% ---------------------- OUTPUT ARRAYS ----------------------
behResp = [];   % AD stim FR - BL stim FR
tagResp = [];   % eventopto stim FR - baseline FR
catSign = [];   % +1 = posCorr, -1 = negCorr
peakVal = [];
peakLag = [];
labels  = strings(0,1);

% Optional bookkeeping
tag_stimFR_all = [];
tag_baseFR_all = [];
ad_stimFR_all  = [];
bl_stimFR_all  = [];

%% ===================== MAIN LOOP =====================
for day = 1:length(MyFolderInfo)

    fprintf('\n========== SESSION %d / %d ==========\n', day, length(MyFolderInfo));
    load(MyFolderInfo(day).name, 'data', 'ExcitRR_BL', 'Salt_index');

    if isempty(Salt_index)
        fprintf('Day %d: Salt_index empty. Skipping...\n', day);
        continue;
    end

    session_units = Salt_index(:);

    % --- keep gate setup from BL control PSTH
    timevec = data.units(session_units(1)).BLCSM_max_vx.time_bins;
    valid_peak_range = find(timevec >= keepWin_ms(1) & timevec <= keepWin_ms(2));

    for uu = 1:numel(session_units)

        unitIdx = session_units(uu);
        u = data.units(unitIdx);

        % ---------------------- KEEP GATE ----------------------
        if ~isfield(u, 'BLCSM_max_vx') || isempty(u.BLCSM_max_vx.PSTH_smoothbin)
            continue
        end

        psth_ctrl = u.BLCSM_max_vx.PSTH_smoothbin(:);
        [~, peak_idx] = max(psth_ctrl);

        if ~ismember(peak_idx, valid_peak_range)
            continue
        end

        % ---------------------- REQUIRE RASTERS ----------------------
        if ~isfield(u, 'ADCSM_pos_thresh') || isempty(u.ADCSM_pos_thresh.Raster)
            continue
        end
        if ~isfield(u, 'BLCSM_pos_thresh') || isempty(u.BLCSM_pos_thresh.Raster)
            continue
        end
        if ~isfield(u, 'eventopto') || isempty(u.eventopto.Raster)
            continue
        end

        R_AD  = u.ADCSM_pos_thresh.Raster;
        R_BL  = u.BLCSM_pos_thresh.Raster;
        R_tag = u.eventopto.Raster;

        % ---------------------- EVENT OPTO METRIC (Hz) ----------------------
        % Stim FR during the 50 ms light pulse
        tag_stimFR = computeFR_inWindow(R_tag, stimWin_ms);

        % Baseline FR from large windows outside stim
        tag_baseFR = computeFR_multiWindow(R_tag, {tagBaseWin1_ms, tagBaseWin2_ms});

        if isnan(tag_stimFR) || isnan(tag_baseFR)
            continue
        end

        t = tag_stimFR - tag_baseFR;   % increase in FR during optotag

        % ---------------------- BEHAVIOR METRIC (Hz) ----------------------
        % Compare AD vs BL during the same 50 ms stim window
        ad_stimFR = computeFR_inWindow(R_AD, stimWin_ms);
        bl_stimFR = computeFR_inWindow(R_BL, stimWin_ms);

        if isnan(ad_stimFR) || isnan(bl_stimFR)
            continue
        end

        b = ad_stimFR - bl_stimFR;     % increase in FR during behavior due to stim

        % ---------------------- BL NEG-LAG XCORR CLASSIFICATION ----------------------
        [peakLag_cont, peakVal_cont] = negLagXCorr_oneUnit( ...
            u, ExcitRR_BL, kinVar, gauss_sigma, ...
            resample_fs, win_pre_ms, win_post_ms, maxLag_ms, lag_step);

        if isnan(peakVal_cont)
            continue
        end

        % ---------------------- SAVE ----------------------
        behResp(end+1,1) = b; %#ok<AGROW>
        tagResp(end+1,1) = t; %#ok<AGROW>

        tag_stimFR_all(end+1,1) = tag_stimFR; %#ok<AGROW>
        tag_baseFR_all(end+1,1) = tag_baseFR; %#ok<AGROW>
        ad_stimFR_all(end+1,1)  = ad_stimFR; %#ok<AGROW>
        bl_stimFR_all(end+1,1)  = bl_stimFR; %#ok<AGROW>

        peakLag(end+1,1) = peakLag_cont; %#ok<AGROW>
        peakVal(end+1,1) = peakVal_cont; %#ok<AGROW>

        catSign(end+1,1) = sign(peakVal_cont); %#ok<AGROW>
        labels(end+1,1)  = "S" + day + "_U" + unitIdx; %#ok<AGROW>
    end
end

fprintf('\nKept paired units with BL neg-lag category: %d\n', numel(behResp));
if isempty(behResp)
    error('No units survived keep gate + required rasters + neg-lag xcorr.');
end

%% ---------------------- REMOVE catSign == 0 ----------------------
okCat = catSign ~= 0;

behResp = behResp(okCat);
tagResp = tagResp(okCat);
catSign = catSign(okCat);
peakLag = peakLag(okCat);
peakVal = peakVal(okCat);
labels  = labels(okCat);

tag_stimFR_all = tag_stimFR_all(okCat);
tag_baseFR_all = tag_baseFR_all(okCat);
ad_stimFR_all  = ad_stimFR_all(okCat);
bl_stimFR_all  = bl_stimFR_all(okCat);

%% ---------------------- COUNTS / SUMMARY ----------------------
fprintf('posCorr units: %d\n', sum(catSign > 0));
fprintf('negCorr units: %d\n', sum(catSign < 0));

fprintf('Mean eventopto delta FR: %.4f Hz\n', mean(tagResp, 'omitnan'));
fprintf('Mean behavior delta FR (AD - BL): %.4f Hz\n', mean(behResp, 'omitnan'));

%% ===================== CONNECTED PLOT SPLIT BY CATEGORY =====================
% Behavior   = AD stim FR - BL stim FR
% Event opto = eventopto stim FR - baseline FR

isPos = catSign > 0;
isNeg = catSign < 0;

xB = 1.00;
xT = 1.18;

allVals = [behResp(:); tagResp(:)];
allVals = allVals(~isnan(allVals));

if isempty(allVals)
    error('No valid values to plot.');
end

yl = [min(allVals) max(allVals)];
if diff(yl) == 0
    yl = yl + [-1 1];
end
pad = 0.08 * diff(yl);
yl = yl + [-pad pad];

figure('Color', 'w', 'Position', [200 200 600 420]);

% -------- POSITIVE CORR --------
subplot(1,2,1); hold on
plotConnected_customX(behResp(isPos), tagResp(isPos), xB, xT);

mB = mean(behResp(isPos), 'omitnan');
mT = mean(tagResp(isPos), 'omitnan');

scatter(xB, mB, 140, 'k', 'filled');
scatter(xT, mT, 140, 'k', 'filled');
plot([xB xT], [mB mT], 'k-', 'LineWidth', 3);

title(sprintf('posCorr (N = %d)', nnz(isPos)));
set(gca, 'XTick', [xB xT], ...
         'XTickLabel', {'Behavior', 'Event opto'});
xlim([0.95 1.23]);
ylim(yl);
yline(0, 'k:');
ylabel('\Delta firing rate (Hz)');
box off

% -------- NEGATIVE CORR --------
subplot(1,2,2); hold on
plotConnected_customX(behResp(isNeg), tagResp(isNeg), xB, xT);

mB = mean(behResp(isNeg), 'omitnan');
mT = mean(tagResp(isNeg), 'omitnan');

scatter(xB, mB, 140, 'k', 'filled');
scatter(xT, mT, 140, 'k', 'filled');
plot([xB xT], [mB mT], 'k-', 'LineWidth', 3);

title(sprintf('negCorr (N = %d)', nnz(isNeg)));
set(gca, 'XTick', [xB xT], ...
         'XTickLabel', {'Behavior', 'Event opto'});
xlim([0.95 1.23]);
ylim(yl);
yline(0, 'k:');
ylabel('\Delta firing rate (Hz)');
box off

%% ===================== SUMMARY PLOT (MEAN ± SEM) =====================

isPos = catSign > 0;
isNeg = catSign < 0;

xB = 1.00;
xT = 1.15;   % tighter spacing

figure('Color','w','Position',[200 200 420 420]);  % narrower figure

% shared y-limits
allVals = [behResp(:); tagResp(:)];
allVals = allVals(~isnan(allVals));

yl = [min(allVals) max(allVals)];
if diff(yl) == 0
    yl = yl + [-1 1];
end
pad = 0.1 * diff(yl);
yl = yl + [-pad pad];

% -------- POSITIVE CORR --------
subplot(1,2,1); hold on

% means
mB = mean(behResp(isPos),'omitnan');
mT = mean(tagResp(isPos),'omitnan');

% SEMs
semB = std(behResp(isPos),'omitnan') ./ sqrt(sum(isPos));
semT = std(tagResp(isPos),'omitnan') ./ sqrt(sum(isPos));

% error bars
errorbar(xB, mB, semB, 'k', 'LineWidth', 2, 'CapSize', 0);
errorbar(xT, mT, semT, 'k', 'LineWidth', 2, 'CapSize', 0);

% large mean dots
scatter(xB, mB, 180, 'k', 'filled');
scatter(xT, mT, 180, 'k', 'filled');

% connecting mean line
plot([xB xT], [mB mT], 'k-', 'LineWidth', 2.5);

title(sprintf('posCorr (N = %d)', nnz(isPos)));
set(gca,'XTick',[xB xT],'XTickLabel',{'Behavior','Event opto'});
xlim([0.95 1.20]);
ylim(yl);
yline(0,'k:');
ylabel('\Delta firing rate (Hz)');
box off

% -------- NEGATIVE CORR --------
subplot(1,2,2); hold on

% means
mB = mean(behResp(isNeg),'omitnan');
mT = mean(tagResp(isNeg),'omitnan');

% SEMs
semB = std(behResp(isNeg),'omitnan') ./ sqrt(sum(isNeg));
semT = std(tagResp(isNeg),'omitnan') ./ sqrt(sum(isNeg));

% error bars
errorbar(xB, mB, semB, 'k', 'LineWidth', 2, 'CapSize', 0);
errorbar(xT, mT, semT, 'k', 'LineWidth', 2, 'CapSize', 0);

% large mean dots
scatter(xB, mB, 180, 'k', 'filled');
scatter(xT, mT, 180, 'k', 'filled');

% connecting mean line
plot([xB xT], [mB mT], 'k-', 'LineWidth', 2.5);

title(sprintf('negCorr (N = %d)', nnz(isNeg)));
set(gca,'XTick',[xB xT],'XTickLabel',{'Behavior','Event opto'});
xlim([0.95 1.20]);
ylim(yl);
yline(0,'k:');
ylabel('\Delta firing rate (Hz)');
box off
%% ---------------------- OPTIONAL STATS ----------------------
if sum(isPos) >= 2
    [~, p_pos] = ttest(behResp(isPos), tagResp(isPos));
    fprintf('posCorr paired t-test: p = %.4g\n', p_pos);
else
    fprintf('posCorr paired t-test: not enough units\n');
end

if sum(isNeg) >= 2
    [~, p_neg] = ttest(behResp(isNeg), tagResp(isNeg));
    fprintf('negCorr paired t-test: p = %.4g\n', p_neg);
else
    fprintf('negCorr paired t-test: not enough units\n');
end


%% ===================== CONNECTED PLOT SPLIT BY CATEGORY (VERY COMPRESSED) =====================
% Optotag on LEFT, Reach-related on RIGHT
% Grand means shown as mean +/- SEM (no connecting mean line)

isPos = catSign > 0;
isNeg = catSign < 0;

% <<< tighter spacing
xOpto  = 1.00;   % LEFT
xReach = 1.18;   % RIGHT

% shared y-limits across both panels
allVals = [tagResp(:); behResp(:)];
allVals = allVals(~isnan(allVals));

if isempty(allVals)
    error('No valid values to plot.');
end

yl = [min(allVals) max(allVals)];
if diff(yl) == 0
    yl = yl + [-1 1];
end
pad = 0.08 * diff(yl);
yl = yl + [-pad pad];

figure('Color','w','Position',[200 200 450 420]);


% -------- POSITIVE CORR --------
subplot(1,2,1); hold on

% individual paired lines
plotConnected_customX(tagResp(isPos), behResp(isPos), xOpto, xReach);

% means
mOpto  = mean(tagResp(isPos), 'omitnan');
mReach = mean(behResp(isPos), 'omitnan');

% SEMs
semOpto  = std(tagResp(isPos),  'omitnan') ./ sqrt(sum(~isnan(tagResp(isPos))));
semReach = std(behResp(isPos), 'omitnan') ./ sqrt(sum(~isnan(behResp(isPos))));

% mean +/- SEM
errorbar(xOpto,  mOpto,  semOpto,  'ko', ...
    'MarkerFaceColor','k', 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
errorbar(xReach, mReach, semReach, 'ko', ...
    'MarkerFaceColor','k', 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);

yline(0,'k:');
title(sprintf('posCorr (BL neg-lag peakVal > 0), N = %d', nnz(isPos)));
set(gca,'XTick',[xOpto xReach],'XTickLabel',{'Optotag','Reach-related'});
xlim([0.95 1.23]);
ylim([-100 200]);
ylabel('\Delta firing rate (Hz)');
box off

% -------- NEGATIVE CORR --------
subplot(1,2,2); hold on

% individual paired lines
plotConnected_customX(tagResp(isNeg), behResp(isNeg), xOpto, xReach);

% means
mOpto  = mean(tagResp(isNeg), 'omitnan');
mReach = mean(behResp(isNeg), 'omitnan');

% SEMs
semOpto  = std(tagResp(isNeg),  'omitnan') ./ sqrt(sum(~isnan(tagResp(isNeg))));
semReach = std(behResp(isNeg), 'omitnan') ./ sqrt(sum(~isnan(behResp(isNeg))));

% mean +/- SEM
errorbar(xOpto,  mOpto,  semOpto,  'ko', ...
    'MarkerFaceColor','k', 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
errorbar(xReach, mReach, semReach, 'ko', ...
    'MarkerFaceColor','k', 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);

yline(0,'k:');
title(sprintf('negCorr (BL neg-lag peakVal < 0), N = %d', nnz(isNeg)));
set(gca,'XTick',[xOpto xReach],'XTickLabel',{'Optotag','Reach-related'});
xlim([0.95 1.23]);
ylim([-100 200]);
ylabel('\Delta firing rate (Hz)');
box off

%% ===================== CONNECTED PLOT SPLIT BY CATEGORY (VERY COMPRESSED) =====================
% Optotag on LEFT, Reach-related on RIGHT
% Grand summary shown as median +/- IQR/2 (no connecting summary line)
% Also prints mean change in firing rate across all units

isPos = catSign > 0;
isNeg = catSign < 0;

% <<< tighter spacing
xOpto  = 1.00;   % LEFT
xReach = 1.18;   % RIGHT

% shared y-limits across both panels
allVals = [tagResp(:); behResp(:)];
allVals = allVals(~isnan(allVals));

if isempty(allVals)
    error('No valid values to plot.');
end

yl = [min(allVals) max(allVals)];
if diff(yl) == 0
    yl = yl + [-1 1];
end
pad = 0.08 * diff(yl);
yl = yl + [-pad pad];

figure('Color','w','Position',[200 200 450 420]);

% -------- POSITIVE CORR --------
subplot(1,2,1); hold on

% individual paired lines
plotConnected_customX(tagResp(isPos), behResp(isPos), xOpto, xReach);

% medians
medOpto  = median(tagResp(isPos), 'omitnan');
medReach = median(behResp(isPos), 'omitnan');

% IQR-based spread
iqrOpto  = iqr(tagResp(isPos));
iqrReach = iqr(behResp(isPos));

% median +/- IQR/2
errorbar(xOpto,  medOpto,  iqrOpto/2,  'ko', ...
    'MarkerFaceColor','k', 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
errorbar(xReach, medReach, iqrReach/2, 'ko', ...
    'MarkerFaceColor','k', 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);

yline(0,'k:');
title(sprintf('posCorr (BL neg-lag peakVal > 0), N = %d', nnz(isPos)));
set(gca,'XTick',[xOpto xReach],'XTickLabel',{'Optotag','Reach-related'});
xlim([0.95 1.23]);
ylim([-150 150]);
ylabel('\Delta firing rate (Hz)');
box off

% -------- NEGATIVE CORR --------
subplot(1,2,2); hold on

% individual paired lines
plotConnected_customX(tagResp(isNeg), behResp(isNeg), xOpto, xReach);

% medians
medOpto  = median(tagResp(isNeg), 'omitnan');
medReach = median(behResp(isNeg), 'omitnan');

% IQR-based spread
iqrOpto  = iqr(tagResp(isNeg));
iqrReach = iqr(behResp(isNeg));

% median +/- IQR/2
errorbar(xOpto,  medOpto,  iqrOpto/2,  'ko', ...
    'MarkerFaceColor','k', 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
errorbar(xReach, medReach, iqrReach/2, 'ko', ...
    'MarkerFaceColor','k', 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);

yline(0,'k:');
title(sprintf('negCorr (BL neg-lag peakVal < 0), N = %d', nnz(isNeg)));
set(gca,'XTick',[xOpto xReach],'XTickLabel',{'Optotag','Reach-related'});
xlim([0.95 1.23]);
ylim([-150 150]);
ylabel('\Delta firing rate (Hz)');
box off

% ===================== SUMMARY STATS =====================

% --- POSITIVE CORR ---
pos_tag   = tagResp(isPos);
pos_beh   = behResp(isPos);

mean_pos_tag = mean(pos_tag, 'omitnan');
sem_pos_tag  = std(pos_tag, 'omitnan') ./ sqrt(sum(~isnan(pos_tag)));

mean_pos_beh = mean(pos_beh, 'omitnan');
sem_pos_beh  = std(pos_beh, 'omitnan') ./ sqrt(sum(~isnan(pos_beh)));

% --- NEGATIVE CORR ---
neg_tag   = tagResp(isNeg);
neg_beh   = behResp(isNeg);

mean_neg_tag = mean(neg_tag, 'omitnan');
sem_neg_tag  = std(neg_tag, 'omitnan') ./ sqrt(sum(~isnan(neg_tag)));

mean_neg_beh = mean(neg_beh, 'omitnan');
sem_neg_beh  = std(neg_beh, 'omitnan') ./ sqrt(sum(~isnan(neg_beh)));

% --- PRINT ---
fprintf('\n==================== SUMMARY STATS ====================\n\n');

fprintf('POSITIVE CORR (N = %d)\n', sum(~isnan(pos_tag)));
fprintf('  Optotag:       %.3f ± %.3f Hz\n', mean_pos_tag, sem_pos_tag);
fprintf('  Reach-related: %.3f ± %.3f Hz\n\n', mean_pos_beh, sem_pos_beh);

fprintf('NEGATIVE CORR (N = %d)\n', sum(~isnan(neg_tag)));
fprintf('  Optotag:       %.3f ± %.3f Hz\n', mean_neg_tag, sem_neg_tag);
fprintf('  Reach-related: %.3f ± %.3f Hz\n', mean_neg_beh, sem_neg_beh);

%%
%%%% ===================== CONNECTED PLOT SPLIT BY CATEGORY (MATCHED PAIRS) =====================

isPos = catSign > 0;
isNeg = catSign < 0;

xOpto  = 1.00;
xReach = 1.18;

figure('Color','w','Position',[200 200 450 420]);

% -------- POSITIVE CORR --------
subplot(1,2,1); hold on

% only keep fully paired data
keepPos = isPos & ~isnan(tagResp) & ~isnan(behResp);

tagPos = tagResp(keepPos);
behPos = behResp(keepPos);

plotConnected_customX(tagPos, behPos, xOpto, xReach);

mOpto  = mean(tagPos);
mReach = mean(behPos);

semOpto  = std(tagPos) / sqrt(numel(tagPos));
semReach = std(behPos) / sqrt(numel(behPos));

errorbar(xOpto,  mOpto,  semOpto,  'ko', ...
    'MarkerFaceColor','k', 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
errorbar(xReach, mReach, semReach, 'ko', ...
    'MarkerFaceColor','k', 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);

yline(0)
title(sprintf('posCorr (N=%d)', numel(tagPos)));
set(gca,'XTick',[xOpto xReach],'XTickLabel',{'Optotag','Reach-related'});
xlim([0.95 1.23]);
ylabel('Peak FR (0–50) − baseline (Hz)');
box off

% -------- NEGATIVE CORR --------
subplot(1,2,2); hold on

keepNeg = isNeg & ~isnan(tagResp) & ~isnan(behResp);

tagNeg = tagResp(keepNeg);
behNeg = behResp(keepNeg);

plotConnected_customX(tagNeg, behNeg, xOpto, xReach);

mOpto  = mean(tagNeg);
mReach = mean(behNeg);

semOpto  = std(tagNeg) / sqrt(numel(tagNeg));
semReach = std(behNeg) / sqrt(numel(behNeg));

errorbar(xOpto,  mOpto,  semOpto,  'ko', ...
    'MarkerFaceColor','k', 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);
errorbar(xReach, mReach, semReach, 'ko', ...
    'MarkerFaceColor','k', 'MarkerSize',10, 'LineWidth',2, 'CapSize',10);

yline(0)
title(sprintf('negCorr (N=%d)', numel(tagNeg)));
set(gca,'XTick',[xOpto xReach],'XTickLabel',{'Optotag','Reach-related'});
xlim([0.95 1.23]);
ylabel('Peak FR (0–50) − baseline (Hz)');
box off

% ===================== PRINT MEAN ± SEM BY CATEGORY =====================

fprintf('\n==================== SUMMARY STATS ====================\n');

% helper function
calc_mean_sem = @(x) deal( ...
    mean(x,'omitnan'), ...
    std(x,'omitnan') ./ sqrt(sum(~isnan(x))) );

% -------- POSITIVE CORR --------
[mOpto_pos,  semOpto_pos]  = calc_mean_sem(tagResp(isPos));
[mReach_pos, semReach_pos] = calc_mean_sem(behResp(isPos));

fprintf('\nPOSITIVE CORR (N = %d)\n', nnz(isPos));
fprintf('  Optotag:      %.3f ± %.3f Hz\n', mOpto_pos,  semOpto_pos);
fprintf('  Reach-related: %.3f ± %.3f Hz\n', mReach_pos, semReach_pos);

% -------- NEGATIVE CORR --------
[mOpto_neg,  semOpto_neg]  = calc_mean_sem(tagResp(isNeg));
[mReach_neg, semReach_neg] = calc_mean_sem(behResp(isNeg));

fprintf('\nNEGATIVE CORR (N = %d)\n', nnz(isNeg));
fprintf('  Optotag:      %.3f ± %.3f Hz\n', mOpto_neg,  semOpto_neg);
fprintf('  Reach-related: %.3f ± %.3f Hz\n', mReach_neg, semReach_neg);

fprintf('======================================================\n\n');
%% ===================== HELPER FUNCTIONS =====================

function fr = computeFR_inWindow(R, win_ms)

    if isempty(R)
        fr = NaN;
        return
    end

    nTime = size(R,2);
    tvec = (1:nTime) - ceil(nTime/2);

    idx = tvec >= win_ms(1) & tvec < win_ms(2);
    if ~any(idx)
        fr = NaN;
        return
    end

    trialSpikes = sum(R(:,idx), 2);
    winDur_s = sum(idx) / 1000;    % assumes 1 ms bins

    if winDur_s <= 0
        fr = NaN;
        return
    end

    trialFR = trialSpikes ./ winDur_s;
    fr = mean(trialFR, 'omitnan');
end

function fr = computeFR_multiWindow(R, winCell)

    if isempty(R)
        fr = NaN;
        return
    end

    nTime = size(R,2);
    tvec = (1:nTime) - ceil(nTime/2);

    idx = false(size(tvec));

    for k = 1:numel(winCell)
        w = winCell{k};
        idx = idx | (tvec >= w(1) & tvec < w(2));
    end

    if ~any(idx)
        fr = NaN;
        return
    end

    trialSpikes = sum(R(:,idx), 2);
    winDur_s = sum(idx) / 1000;    % assumes 1 ms bins

    if winDur_s <= 0
        fr = NaN;
        return
    end

    trialFR = trialSpikes ./ winDur_s;
    fr = mean(trialFR, 'omitnan');
end

