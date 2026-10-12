%% integrated_xcorr_and_psth_classification.m
clear; clc;

cd("C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\StimAnimals\BlockSeparated_CSMvsNoCSM_inducedalignment")
MyFolderInfo = dir('*Block_CSMvsNoCSM.mat'); 

%% ---------------------- ANALYSIS PARAMS ---------------------
rasterField   = 'max_vx';
kinVar        = 7;
gauss_sigma   = 10;
resample_fs   = 200;
win_pre_ms    = 100;
win_post_ms   = 150;
maxLag_ms     = 100;

lag_step      = 5;

baseline_win_def = [-1000 -750];
%baseline_win_def = [-1000 -500];
valid_peak_range_ms = [-250 250];

%%% ---------------- NEW CSM ANALYSIS PARAMS ----------------
csm_epoch_win      = [-25 100];   % main CSM window (ms)
%csm_epoch_win      = [0 100];   % main CSM window (ms)
%csm_epoch_win      = [-100 0];   % main CSM window (ms)
%csm_baseline_win   = [-500 -250]; % pre-CSM baseline (ms)
csm_baseline_win   = [-1000 -750]; % pre-CSM baseline (ms) - using same baseline 
%csm_baseline_win   = [-1000 -500]; % pre-CSM baseline (ms) - using same baseline 

min_spikes_csm     = 5;           % unit inclusion threshold
nPerm              = 5000;
alpha_perm         = 0.05;

%%% time-resolved analysis
slide_win_ms       = 25;          % sliding window size
slide_step_ms      = 5;
cluster_alpha      = 0.05;        % cluster-forming threshold


XCorrData = struct;
AllGroups = struct();
NPdata = struct;

group_names = {'posLag_posCorr','posLag_negCorr','negLag_posCorr','negLag_negCorr'};
for g = 1:numel(group_names)
    gn = group_names{g};
    AllGroups.(gn).CSM_before = [];
    AllGroups.(gn).CSM_after  = [];
    AllGroups.(gn).NoCSM_before = [];
    AllGroups.(gn).NoCSM_after  = [];
end


%%
% ---------------------- MAIN SESSION LOOP ----------------------
for day = 1:length(MyFolderInfo)

    fprintf('\n========== SESSION %d / %d ==========\n', day, length(MyFolderInfo));
     load(MyFolderInfo(day).name, 'Reaches','data','filtered_modulated_indices','ExcitRR_BL','ExcitRR_AD','ExcitRR_WOE','ExcitRR_WO','Salt_index');

     if isempty(Salt_index)
         fprintf('Day %d: Salt_index is empty. Skipping...\n', day);
         continue;
     end



    %Save kinematics
    % NoCSM and ReReaches (for ExcitRR)
    NPdata(day).BLCSM = ExcitRR_BL.ReReaches;
    NPdata(day).BLNoCSM = ExcitRR_BL.NoReReach;

    NPdata(day).ADCSM = ExcitRR_AD.ReReaches;

    NPdata(day).WOCSM = ExcitRR_WO.ReReaches;
    NPdata(day).WONoCSM = ExcitRR_WO.NoReReach;
   
    NPdata(day).WOECSM = ExcitRR_WOE.ReReaches;
    NPdata(day).WOENoCSM = ExcitRR_WOE.NoReReach;

    session_units = Salt_index;
    nUnits = length(session_units);

    % timevec for PSTHs (use your stored timebins)
    timevec = data.units(session_units(1)).BLCSM_max_vx.time_bins; % ms vector
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
        %spikes_full = data.units(unitIdx).(rasterField).Raster; 
        spikes_full = data.units(unitIdx).BLCSM_max_vx.Raster;  
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
        %C = AlignPeakV_1ms_edit(Reaches, kinVar, win_pre_ms, win_post_ms); 
        Reaches=ExcitRR_BL.ReReaches;
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

        % ---------------- MEDIAN + CONTINUOUS PEAK/MIN ----------------
        medianCorr = nanmedian(xcorr_trials,2);

        % Max correlation interpolation
        [peakVal, idxMax] = max(medianCorr);
        if idxMax > 1 && idxMax < length(medianCorr)
            x = lags_ms(idxMax-1 : idxMax+1);
            y = medianCorr(idxMax-1 : idxMax+1);
            p = polyfit(x, y, 2);
            peakLag_cont = -p(2)/(2*p(1));
            peakVal_cont = polyval(p, peakLag_cont);
        else
            peakLag_cont = lags_ms(idxMax);
            peakVal_cont = peakVal;
        end

        % Min correlation interpolation
        [minVal, idxMin] = min(medianCorr);
        if idxMin > 1 && idxMin < length(medianCorr)
            x = lags_ms(idxMin-1 : idxMin+1);
            y = medianCorr(idxMin-1 : idxMin+1);
            p = polyfit(x, y, 2);
            minLag_cont = -p(2)/(2*p(1));
            minVal_cont = polyval(p, minLag_cont);
        else
            minLag_cont = lags_ms(idxMin);
            minVal_cont = minVal;
        end

        % ---------------- SELECT BY LARGEST ABSOLUTE MAG ----------------
        if abs(peakVal_cont) >= abs(minVal_cont)
            chosenVal = peakVal_cont;
            chosenLag = peakLag_cont;
        else
            chosenVal = minVal_cont;
            chosenLag = minLag_cont;
        end

        % store xcorr results
        XCorrData(day).units(uu).unitIdx      = unitIdx;
        XCorrData(day).units(uu).lags_ms      = lags_ms;
        XCorrData(day).units(uu).xcorr_trials = xcorr_trials;
        XCorrData(day).units(uu).medianCorr   = medianCorr;
        XCorrData(day).units(uu).peakVal      = peakVal_cont;
        XCorrData(day).units(uu).peakLag      = peakLag_cont;
        XCorrData(day).units(uu).minVal       = minVal_cont;
        XCorrData(day).units(uu).minLag       = minLag_cont;
        XCorrData(day).units(uu).chosenVal    = chosenVal;
        XCorrData(day).units(uu).chosenLag    = chosenLag;

        % ================== NEW (REVIEWER-FACING): PER-NEURON SIGNIFICANCE USING POISSON GLM (EPOCH) + FAST TIME-RESOLVED MAX-STAT =================
        % Goal:
        %  - Per-neuron significance using standard count model (Poisson GLM + LR test) at the EPOCH level
        %  - Time-resolved modulation in small bins using BASELINE-SUBTRACTED FR + MAX-STAT PERMUTATION (NO GLM INSIDE PERM LOOP)
        % Notes:
        %  - Epoch GLM uses baseline as covariate (controls trial-to-trial baseline variability)
        %  - Time-resolved uses baseline-subtracted FR per bin and max-stat correction (FWER) across bins
        %  - This avoids non-stop running + ill-conditioned GLM warnings from refitting GLMs inside permutations

        % --- reproducibility (optional) ---
        rng(0,'twister');

        % ================== LOAD RASTERS ==================
        BL_raster   = data.units(unitIdx).BLCSM_csmT2.Raster;   % [nTrials x nTime], 1 ms bins
        Stim_raster = data.units(unitIdx).ADCSM_csmT2.Raster;   % [nTrials x nTime], 1 ms bins

        % ================== TIME VECTOR ==================
        nTime_bl   = size(BL_raster,2);
        nTime_stim = size(Stim_raster,2);
        nTime      = min(nTime_bl, nTime_stim);  % enforce same time length (safest)
        BL_raster   = BL_raster(:,1:nTime);
        Stim_raster = Stim_raster(:,1:nTime);

        % Use stored time bins if available; otherwise assume -1000:1000 ms
        full_time_1ms = (-1000:1000)';   % 2001 bins
        assert(nTime <= numel(full_time_1ms), 'Raster longer than assumed -1000:1000 ms vector.');
        time_csm = full_time_1ms(1:nTime);

        % ================== EPOCH INDICES ==================
        csm_idx  = (time_csm >= csm_epoch_win(1))    & (time_csm <= csm_epoch_win(2));
        base_idx = (time_csm >= csm_baseline_win(1)) & (time_csm <= csm_baseline_win(2));

        assert(any(csm_idx),  'CSM epoch window is empty (check csm_epoch_win vs time_csm).');
        assert(any(base_idx), 'Baseline window is empty (check csm_baseline_win vs time_csm).');

        % ================== COUNTS PER TRIAL ==================
        BL_base_spk   = sum(BL_raster(:,base_idx),2);
        Stim_base_spk = sum(Stim_raster(:,base_idx),2);

        BL_csm_spk    = sum(BL_raster(:,csm_idx),2);
        Stim_csm_spk  = sum(Stim_raster(:,csm_idx),2);

        % Optional inclusion gate (your existing min_spikes_csm)
        % Make it depend on CSM window counts across BOTH conditions
        if ( (sum(BL_csm_spk) + sum(Stim_csm_spk)) < min_spikes_csm )
            XCorrData(day).units(uu).CSM_epoch_glm = struct('skip',true,'reason','min_spikes_csm');
            XCorrData(day).units(uu).CSM_time_fast = struct('skip',true,'reason','min_spikes_csm');
        else

            % ================== PER-NEURON POISSON GLM (EPOCH-LEVEL) + LR TEST ==================
            % Model:
            %   y = CSM_count
            %   predictors: condition (Stim vs BL), baseline_count (nuisance), intercept
            %
            % LR test compares:
            %   Reduced: intercept + baseline
            %   Full:    intercept + baseline + condition
            % and tests whether adding "condition" significantly improves fit.

            y    = [Stim_csm_spk; BL_csm_spk];
            cond = [ones(numel(Stim_csm_spk),1); zeros(numel(BL_csm_spk),1)];  % 1=Stim, 0=BL
            base = [Stim_base_spk; BL_base_spk];

            % Safe z-score for numerical stability
            sd_base = std(base);
            if sd_base < 1e-6
                base_z = zeros(size(base));
            else
                base_z = (base - mean(base)) ./ sd_base;
            end

            Xfull = [cond, base_z];  % glmfit adds intercept automatically
            Xred  = [base_z];

            opts = statset('glmfit');
            opts.MaxIter = 200;

            % Guard against all-zero outcomes
            if sum(y) == 0 || all(y==0)
                XCorrData(day).units(uu).CSM_epoch_glm = struct('skip',true,'reason','all_zero_counts');
            else
                try
                    [b_full, dev_full, stats_full] = glmfit(Xfull, y, 'poisson', 'link','log', 'options', opts);
                    [b_red,  dev_red]              = glmfit(Xred,  y, 'poisson', 'link','log', 'options', opts);

                    % Likelihood ratio test for "condition" effect:
                    LR   = dev_red - dev_full;
                    p_LR = 1 - chi2cdf(LR, 1);

                    % Effect size: rate ratio for Stim vs BL (exp(beta_condition))
                    beta_cond = b_full(2); % b_full = [intercept; cond; base]
                    RR        = exp(beta_cond);

                    % Wald CI for beta_cond -> CI on RR
                    se_cond = stats_full.se(2);
                    CI_beta = beta_cond + [-1.96 1.96]*se_cond;
                    CI_RR   = exp(CI_beta);

                    % Also store baseline-subtracted FR delta for plots
                    csm_win_s  = sum(csm_idx)  / 1000;
                    base_win_s = sum(base_idx) / 1000;

                    BL_netFR   = (BL_csm_spk./csm_win_s)   - (BL_base_spk./base_win_s);
                    Stim_netFR = (Stim_csm_spk./csm_win_s) - (Stim_base_spk./base_win_s);
                    obs_delta_netFR = mean(Stim_netFR) - mean(BL_netFR);

                    XCorrData(day).units(uu).CSM_epoch_glm.skip        = false;
                    XCorrData(day).units(uu).CSM_epoch_glm.p_LR        = p_LR;
                    XCorrData(day).units(uu).CSM_epoch_glm.sig_LR      = (p_LR < alpha_perm);
                    XCorrData(day).units(uu).CSM_epoch_glm.LR          = LR;
                    XCorrData(day).units(uu).CSM_epoch_glm.beta_cond   = beta_cond;
                    XCorrData(day).units(uu).CSM_epoch_glm.RR          = RR;
                    XCorrData(day).units(uu).CSM_epoch_glm.CI_RR       = CI_RR;
                    XCorrData(day).units(uu).CSM_epoch_glm.delta_netFR = obs_delta_netFR;
                    XCorrData(day).units(uu).CSM_epoch_glm.BL_netFR    = BL_netFR;
                    XCorrData(day).units(uu).CSM_epoch_glm.Stim_netFR  = Stim_netFR;

                catch ME
                    XCorrData(day).units(uu).CSM_epoch_glm = struct('skip',true,'reason',['glmfit_failed: ' ME.message]);
                end
            end

            % ================== TIME-RESOLVED (FAST + ROBUST): BASELINE-SUBTRACTED BIN FR + MAX-STAT PERMUTATION ==================
            % Avoids refitting GLMs inside permutation loops (which caused non-stop running).
            % Still provides a shuffled-control corrected p per bin across the CSM epoch.

            bin_ms      = 10;   % try 10 or 15; 5ms often too sparse
            bin_step_ms = 5;
            min_total_spikes_per_bin = 5;  % skip bins too sparse across ALL trials

            bin_starts = csm_epoch_win(1):bin_step_ms:(csm_epoch_win(2)-bin_ms);
            nBins      = numel(bin_starts);

            base_win_s = sum(base_idx) / 1000;
            bin_s      = bin_ms / 1000;

            % trial-level baseline FR (Hz)
            BL_base_FR   = BL_base_spk   / base_win_s;
            Stim_base_FR = Stim_base_spk / base_win_s;

            % stack baseline FR in same order [Stim; BL]
            baseFR_all = [Stim_base_FR; BL_base_FR];

            % condition labels in stacked order [Stim; BL]
            nStim   = size(Stim_raster,1);
            nBL     = size(BL_raster,1);
            labels0 = [true(nStim,1); false(nBL,1)];

            % Precompute baseline-subtracted FR values for ALL trials once
            % netFR_all(trial, bin) = (bin_count/bin_s) - baseline_FR(trial)
            netFR_all = nan(nStim+nBL, nBins);
            delta_t   = nan(nBins,1);

            for bbi = 1:nBins
                b0 = bin_starts(bbi);
                b_idx = (time_csm >= b0) & (time_csm < (b0 + bin_ms));
                if ~any(b_idx), continue; end

                Stim_bin_spk = sum(Stim_raster(:,b_idx),2);
                BL_bin_spk   = sum(BL_raster(:,b_idx),2);

                % Skip super sparse bins
                if (sum(Stim_bin_spk) + sum(BL_bin_spk)) < min_total_spikes_per_bin
                    continue;
                end

                binFR_all = [Stim_bin_spk; BL_bin_spk] / bin_s;
                netFR_all(:,bbi) = binFR_all - baseFR_all;

                delta_t(bbi) = mean(netFR_all(labels0,bbi), 'omitnan') - mean(netFR_all(~labels0,bbi), 'omitnan');
            end

            % Build max-stat null on |delta_t| across bins
            max_null = nan(nPerm,1);
            for pp = 1:nPerm
                permLab = labels0(randperm(numel(labels0)));

                perm_delta = nan(nBins,1);
                for bbi = 1:nBins
                    x = netFR_all(:,bbi);
                    if all(isnan(x)), continue; end
                    perm_delta(bbi) = mean(x(permLab), 'omitnan') - mean(x(~permLab), 'omitnan');
                end

                max_null(pp) = max(abs(perm_delta), [], 'omitnan');
            end

            % Corrected p per bin (FWER via max-stat)
            p_corr = nan(nBins,1);
            for bbi = 1:nBins
                if isnan(delta_t(bbi))
                    p_corr(bbi) = NaN;
                else
                    p_corr(bbi) = mean(max_null >= abs(delta_t(bbi)));
                end
            end

            % Store time results
            XCorrData(day).units(uu).CSM_time_fast.skip        = false;
            XCorrData(day).units(uu).CSM_time_fast.bin_ms      = bin_ms;
            XCorrData(day).units(uu).CSM_time_fast.step_ms     = bin_step_ms;
            XCorrData(day).units(uu).CSM_time_fast.starts_ms   = bin_starts(:);
            XCorrData(day).units(uu).CSM_time_fast.delta_netFR = delta_t(:);   % mean(Stim-BL) baseline-subtracted FR per bin
            XCorrData(day).units(uu).CSM_time_fast.p_corr      = p_corr(:);    % max-stat corrected p per bin
            XCorrData(day).units(uu).CSM_time_fast.max_null    = max_null;     % optional diagnostics
        end
        % ================== END NEW BLOCK ==================



        % ---------------- PEAK-TIME BEFORE/AFTER & KEEP CHECK ----------------
        psth = data.units(unitIdx).BLCSM_max_vx.PSTH_smoothbin;
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

        % ---------------- ASSIGN CATEGORY (lag sign & corr sign) ----------------
        lag_sign = sign(chosenLag);
        corr_sign = sign(chosenVal);
        if lag_sign >= 0, lag_label = 'posLag'; else lag_label = 'negLag'; end
        if corr_sign >= 0, corr_label = 'posCorr'; else corr_label = 'negCorr'; end
        cat_label = [lag_label '_' corr_label];
        XCorrData(day).units(uu).category = cat_label;

        % ---------------- COMPUTE Z-SCORED PSTHs (Wasout versus baseline ) ----------------
        % Pick the same PSTH types you used before
        csm    = data.units(unitIdx).BLCSM_csmT2.PSTH_smoothbin;
        nocsm  = data.units(unitIdx).WOECSM_csmT2.PSTH_smoothbin;

        % csm    = data.units(unitIdx).BLCSM_csmT2.PSTH_smoothbin;
        % nocsm  = data.units(unitIdx).ADCSM_csmT2.PSTH_smoothbin;

        %  csm    = data.units(unitIdx).BLCSM_pos_thresh.PSTH_smoothbin;
        % nocsm  = data.units(unitIdx).WOECSM_pos_thresh.PSTH_smoothbin;

        % csm    = data.units(unitIdx).BLCSM_pos_thresh.PSTH_smoothbin;
        % nocsm  = data.units(unitIdx).ADCSM_pos_thresh.PSTH_smoothbin;

        %  csm    = data.units(unitIdx).BLCSM_max_vx.PSTH_smoothbin;
        % nocsm  = data.units(unitIdx).WOECSM_max_vx.PSTH_smoothbin;

        %  csm    = data.units(unitIdx).BLCSM_max_vx.PSTH_smoothbin;
        % nocsm  = data.units(unitIdx).ADCSM_max_vx.PSTH_smoothbin;
        

        % csm    = data.units(unitIdx).AllCSM_csmT2.PSTH_smoothbin;
        % nocsm  = data.units(unitIdx).NoCSM_csmT2.PSTH_smoothbin;
        % 
        %   csm    = data.units(unitIdx).AllCSM_csmT0.PSTH_smoothbin;
        % nocsm  = data.units(unitIdx).NoCSM_csmT0.PSTH_smoothbin;

        % baseline from CSM as you did previously
        baseline_mean = mean(csm(baseline_win));
        baseline_std = std(csm(baseline_win));
        if baseline_std == 0, baseline_std = 1; end

        csm_z   = (csm - baseline_mean) / baseline_std;
        nocsm_z = (nocsm - baseline_mean) / baseline_std;
        % ---------------- AGGREGATE INTO GROUPS ----------------
        if ismember(cat_label, group_names)
            % store in AllGroups
            if strcmp(XCorrData(day).units(uu).before_after,'before')
                AllGroups.(cat_label).CSM_before = vertcat(AllGroups.(cat_label).CSM_before, csm_z);
                AllGroups.(cat_label).NoCSM_before = vertcat(AllGroups.(cat_label).NoCSM_before, nocsm_z);
            else
                AllGroups.(cat_label).CSM_after = vertcat(AllGroups.(cat_label).CSM_after, csm_z);
                AllGroups.(cat_label).NoCSM_after = vertcat(AllGroups.(cat_label).NoCSM_after, nocsm_z);
            end
        end

    end % unit loop

    % Save session-level lists
    XCorrData(day).before_peakv_units = before_idx_all;
    XCorrData(day).after_peakv_units  = after_idx_all;
    fprintf(' Session %d summary: %d before, %d after (kept units)\n', day, length(before_idx_all), length(after_idx_all));
end

fprintf('\n================ ALL SESSIONS COMPLETE ================\n');

%%
% ================== POPULATION SHUFFLED CONTROL (FIG 5D-STYLE, UPDATED) ==================
% Uses per-unit baseline-subtracted FR differences from CSM_time_fast
% Tests whether the POPULATION mean timecourse exceeds a shuffled null
% Controls family-wise error across time via max-stat

all_delta = [];
all_starts = [];

for day = 1:numel(XCorrData)
    if ~isfield(XCorrData(day),'units'), continue; end
    for uu = 1:numel(XCorrData(day).units)

        if isfield(XCorrData(day).units(uu),'CSM_time_fast') && ...
           ~XCorrData(day).units(uu).CSM_time_fast.skip

            d = XCorrData(day).units(uu).CSM_time_fast.delta_netFR(:)';
            s = XCorrData(day).units(uu).CSM_time_fast.starts_ms(:)';

            % initialize shared timebase
            if isempty(all_starts)
                all_starts = s;
            end

            % only include units with identical binning
            if numel(s)==numel(all_starts) && all(s==all_starts)
                all_delta = [all_delta; d]; %#ok<AGROW>
            end
        end
    end
end

if isempty(all_delta)
    warning('No units found for population shuffle (CSM_time_fast missing).');
else
    % ---------------- OBSERVED POPULATION EFFECT ----------------
    pop_delta = mean(all_delta,1,'omitnan');           % population mean FR difference
    obs_peak  = max(abs(pop_delta), [], 'omitnan');    % peak absolute effect

    % ---------------- SHUFFLED NULL (SIGN-FLIP) ----------------
    % Nonparametric test of whether population mean deviates from 0 consistently
    nPerm_pop = 5000;
    null_peak = nan(nPerm_pop,1);

    for p = 1:nPerm_pop
        flip = (rand(size(all_delta,1),1) > 0.5)*2 - 1;  % ±1 per neuron
        pop_perm = mean(all_delta .* flip, 1, 'omitnan');
        null_peak(p) = max(abs(pop_perm), [], 'omitnan');
    end

    p_pop = mean(null_peak >= obs_peak);

    % ---------------- STORE ----------------
    PopCSM = struct();
    PopCSM.starts_ms    = all_starts;
    PopCSM.pop_deltaFR  = pop_delta;
    PopCSM.obs_peak     = obs_peak;
    PopCSM.null_peak    = null_peak;
    PopCSM.p_pop_peak   = p_pop;
    PopCSM.n_units      = size(all_delta,1);

    XCorrData(1).PopulationCSM = PopCSM;

    fprintf('Population shuffled control (Fig 5D): N=%d units, p(peak|null)=%.4g\n', ...
            PopCSM.n_units, p_pop);
end
% ================== END POPULATION SHUFFLED CONTROL ==================

%% ===================== PLOT ALL RESULTS (RUN AFTER ANALYSIS FINISHES) =====================
% Assumes your main script has already filled XCorrData and (optionally) XCorrData(1).PopulationCSM
% Plots:
%   1) Epoch-level Poisson GLM (LR): delta distribution + volcano + sig counts (Early/Late)
%   2) Time-resolved: population mean timecourse + significant bins
%   3) Time-resolved: per-neuron heatmap of delta_netFR (sorted by peak)
%   4) Time-resolved: per-neuron heatmap of p_corr (or sig mask)
%   5) By-category population timecourses (xcorr category + before/after)
%
% NOTE: Uses fields created by the updated analysis:
%   XCorrData(day).units(uu).CSM_epoch_glm (p_LR, sig_LR, delta_netFR, RR, CI_RR)
%   XCorrData(day).units(uu).CSM_time_fast (starts_ms, delta_netFR, p_corr)
%   XCorrData(day).units(uu).before_after  ('before'/'after')
%   XCorrData(day).units(uu).category      ('posLag_posCorr', etc.)

clear P; close all;

% -------------------- COLLECT ALL UNITS --------------------
pLR      = [];
sigLR    = [];
deltaFR  = [];
RR       = [];
CI_RR_lo = [];
CI_RR_hi = [];

beforeAfter = {};
catLabel     = {};
unitKey      = {};  % "day-uu-unitIdx"

% time-resolved matrices
all_starts = [];
Dmat = [];   % [nUnits x nBins] delta_netFR
Pmat = [];   % [nUnits x nBins] p_corr

for day = 1:numel(XCorrData)
    if ~isfield(XCorrData(day),'units'), continue; end
    for uu = 1:numel(XCorrData(day).units)
        u = XCorrData(day).units(uu);

        % ---- epoch GLM ----
        if isfield(u,'CSM_epoch_glm') && isfield(u.CSM_epoch_glm,'skip') && ~u.CSM_epoch_glm.skip
            pLR(end+1,1)     = u.CSM_epoch_glm.p_LR;
            sigLR(end+1,1)   = u.CSM_epoch_glm.sig_LR;
            deltaFR(end+1,1) = u.CSM_epoch_glm.delta_netFR;

            if isfield(u.CSM_epoch_glm,'RR')
                RR(end+1,1) = u.CSM_epoch_glm.RR;
            else
                RR(end+1,1) = NaN;
            end

            if isfield(u.CSM_epoch_glm,'CI_RR')
                CI_RR_lo(end+1,1) = u.CSM_epoch_glm.CI_RR(1);
                CI_RR_hi(end+1,1) = u.CSM_epoch_glm.CI_RR(2);
            else
                CI_RR_lo(end+1,1) = NaN;
                CI_RR_hi(end+1,1) = NaN;
            end

            if isfield(u,'before_after'); beforeAfter{end+1,1} = u.before_after; else; beforeAfter{end+1,1} = 'unknown'; end
            if isfield(u,'category');     catLabel{end+1,1}    = u.category;     else; catLabel{end+1,1}    = 'unknown'; end
            if isfield(u,'unitIdx');      unitKey{end+1,1}     = sprintf('D%02d-U%03d-Idx%04d', day, uu, u.unitIdx);
            else;                         unitKey{end+1,1}     = sprintf('D%02d-U%03d', day, uu);
            end
        end

        % ---- time fast ----
        if isfield(u,'CSM_time_fast') && isfield(u.CSM_time_fast,'skip') && ~u.CSM_time_fast.skip
            s = u.CSM_time_fast.starts_ms(:)';
            d = u.CSM_time_fast.delta_netFR(:)';
            p = u.CSM_time_fast.p_corr(:)';

            if isempty(all_starts)
                all_starts = s;
            end

            % only keep units with identical binning
            if numel(s)==numel(all_starts) && all(s==all_starts)
                Dmat = [Dmat; d]; %#ok<AGROW>
                Pmat = [Pmat; p]; %#ok<AGROW>
            end
        end
    end
end

nEpochUnits = numel(pLR);
fprintf('Collected %d units with epoch GLM results.\n', nEpochUnits);
fprintf('Collected %d units with time-resolved results.\n', size(Dmat,1));

% -------------------- BASIC SANITY --------------------
if isempty(pLR)
    warning('No epoch GLM results found (CSM_epoch_glm.skip may be true for all units).');
end
if isempty(Dmat)
    warning('No time-resolved results found (CSM_time_fast missing/skip or binning mismatch).');
end

%% ===================== 1) EPOCH-LEVEL SUMMARY: COUNTS + DISTRIBUTIONS =====================
if ~isempty(pLR)
    alpha = 0.05;

    % counts (overall)
    nSig = sum(sigLR==1);
    fprintf('Epoch GLM significant (p<%.3f): %d / %d\n', alpha, nSig, nEpochUnits);

    % counts by before/after
    isBefore = strcmp(beforeAfter,'before');
    isAfter  = strcmp(beforeAfter,'after');

    nSig_before = sum(sigLR(isBefore)==1);
    nTot_before = sum(isBefore);
    nSig_after  = sum(sigLR(isAfter)==1);
    nTot_after  = sum(isAfter);

    fprintf('  BEFORE: %d / %d sig\n', nSig_before, nTot_before);
    fprintf('  AFTER : %d / %d sig\n', nSig_after,  nTot_after);

    % ---------- Figure A: delta_netFR distribution ----------
    fig = figure('Color','w','Position',[100 100 1100 500]);
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    nexttile; hold on;
    histogram(deltaFR, 30);
    xline(0,'k--','LineWidth',1);
    title('Epoch GLM: \Delta net FR (Stim - BL)');
    xlabel('\Delta net FR (Hz)');
    ylabel('# units');

    nexttile; hold on;
    % volcano: x=deltaFR, y=-log10(p)
    yv = -log10(max(pLR, realmin));
    scatter(deltaFR, yv, 25, 'filled');
    yline(-log10(alpha),'k--','LineWidth',1);
    xline(0,'k--','LineWidth',1);
    title('Epoch GLM: Volcano');
    xlabel('\Delta net FR (Hz)');
    ylabel('-log_{10}(p_{LR})');

    % annotate counts
    text(0.02,0.98, sprintf('N=%d\nsig=%d', nEpochUnits, nSig), 'Units','normalized', ...
        'VerticalAlignment','top','FontWeight','bold');

    % ---------- Figure B: RR (effect size) ----------
    fig = figure('Color','w','Position',[100 100 900 450]);
    hold on;
    validRR = ~isnan(RR);
    if any(validRR)
        histogram(log2(RR(validRR)), 30);
        xline(0,'k--','LineWidth',1);
        title('Epoch GLM Effect Size: log2(Rate Ratio)   (Stim vs BL)');
        xlabel('log2(RR)');
        ylabel('# units');
    else
        text(0.5,0.5,'No RR values found','HorizontalAlignment','center');
        axis off;
    end
end

%% ===================== 2) POPULATION TIMECOURSE (FAST DELTA) + SIG BINS =====================
if ~isempty(Dmat)
    t = all_starts(:);                      % bin start times (ms)
    pop = mean(Dmat,1,'omitnan')';          % mean across units
    sem = std(Dmat,0,1,'omitnan')' ./ sqrt(sum(~isnan(Dmat),1)');

    % If you ran the population shuffle block, use it; otherwise show mean±SEM only
    havePopShuffle = isfield(XCorrData(1),'PopulationCSM') && isfield(XCorrData(1).PopulationCSM,'pop_deltaFR');

    fig = figure('Color','w','Position',[100 100 1100 500]); hold on;

    % shaded SEM
    x = t;
    y = pop;
    y1 = pop - sem;
    y2 = pop + sem;
    fill([x; flipud(x)], [y1; flipud(y2)], 0.9*[1 1 1], 'EdgeColor','none');
    plot(x,y,'k','LineWidth',2);
    xline(0,'k--','LineWidth',1);
    yline(0,'k-','LineWidth',1);

    xlabel('Time (ms, bin start)');
    ylabel('Population \Delta net FR (Hz)');
    title('Population timecourse (mean \pm SEM): \Delta net FR (Stim - BL)');

    if havePopShuffle
        % Mark peak p-value from shuffle
        p_pop = XCorrData(1).PopulationCSM.p_pop_peak;
        txt = sprintf('Shuffle peak p = %.3g   (N=%d)', p_pop, XCorrData(1).PopulationCSM.n_units);
        text(0.02,0.98, txt, 'Units','normalized','VerticalAlignment','top','FontWeight','bold');
    end
end

%% ===================== 3) PER-NEURON HEATMAP: delta_netFR (SORTED BY PEAK TIME) =====================
if ~isempty(Dmat)
    % sort by time of peak abs effect
    [~, peakIdx] = max(abs(Dmat), [], 2, 'omitnan');
    [~, ord] = sort(peakIdx, 'ascend');

    fig = figure('Color','w','Position',[100 100 1100 650]);
    imagesc(all_starts, 1:size(Dmat,1), Dmat(ord,:));
    axis tight;
    xlabel('Time (ms, bin start)');
    ylabel('Units (sorted by peak time)');
    title('Per-unit \Delta net FR (Stim - BL) heatmap');
    colorbar;
    xline(0,'k--','LineWidth',1);
end

%% ===================== 4) PER-NEURON HEATMAP: corrected p-values (or sig mask) =====================
if ~isempty(Pmat)
    sigMask = Pmat < 0.05;

    fig = figure('Color','w','Position',[100 100 1100 650]);
    imagesc(all_starts, 1:size(sigMask,1), sigMask(ord,:));
    axis tight;
    xlabel('Time (ms, bin start)');
    ylabel('Units (same order as heatmap above)');
    title('Per-unit significance (p_{corr}<0.05) across time bins');
    colorbar;
    xline(0,'k--','LineWidth',1);
end

%% ===================== 5) POPULATION TIMECOURSE BY XCORR CATEGORY + BEFORE/AFTER =====================
groups = {'posLag_posCorr','posLag_negCorr','negLag_posCorr','negLag_negCorr'};
baList = {'before','after'};

% Initialize containers
GroupMat = struct();
for gi = 1:numel(groups)
    for bi = 1:numel(baList)
        GroupMat.(groups{gi}).(baList{bi}) = [];
    end
end

% ---- helper: safely convert a label to a char vector (or '' if not usable) ----
label2char = @(x) local_label2char(x);

% Collect timecourses
for day = 1:numel(XCorrData)
    if ~isfield(XCorrData(day),'units'), continue; end
    for uu = 1:numel(XCorrData(day).units)
        u = XCorrData(day).units(uu);

        % Must have time-resolved results
        if ~(isfield(u,'CSM_time_fast') && isfield(u.CSM_time_fast,'skip') && ~u.CSM_time_fast.skip)
            continue;
        end

        % Safe labels
        if isfield(u,'category')
            catC = label2char(u.category);
        else
            catC = '';
        end

        if isfield(u,'before_after')
            baC = label2char(u.before_after);
        else
            baC = '';
        end

        % Only keep valid labels
        if isempty(catC) || isempty(baC), continue; end
        if ~any(strcmp(catC, groups)), continue; end
        if ~any(strcmp(baC, baList)), continue; end

        % Grab time axis + data
        s = u.CSM_time_fast.starts_ms(:)';
        d = u.CSM_time_fast.delta_netFR(:)';

        % Initialize shared timebase if needed
        if isempty(all_starts)
            all_starts = s;
        end

        % Only include if binning matches
        if numel(s)==numel(all_starts) && all(s==all_starts)
            GroupMat.(catC).(baC) = [GroupMat.(catC).(baC); d]; 
        end
    end
end

% Plot groupwise timecourses
fig = figure('Color','w','Position',[100 100 1300 700]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

for gi = 1:numel(groups)
    nexttile; hold on;

    for bi = 1:numel(baList)
        D = GroupMat.(groups{gi}).(baList{bi});
        if isempty(D), continue; end

        m  = mean(D,1,'omitnan')';
        se = std(D,0,1,'omitnan')' ./ sqrt(sum(~isnan(D),1)');

        x  = all_starts(:);
        y1 = m - se;
        y2 = m + se;

        fill([x; flipud(x)], [y1; flipud(y2)], 0.92*[1 1 1], 'EdgeColor','none'); 
        plot(x, m, 'LineWidth', 2);
    end

    xline(0,'k--','LineWidth',1);
    yline(0,'k-','LineWidth',1);
    title(groups{gi}, 'Interpreter','none');
    xlabel('Time (ms, bin start)');
    ylabel('\Delta net FR (Hz)');

    nB = size(GroupMat.(groups{gi}).before,1);
    nA = size(GroupMat.(groups{gi}).after,1);
    legend({sprintf('before (N=%d)', nB), sprintf('after (N=%d)', nA)}, 'Location','best');
end

% --------------- local helper (put at END of your script) ---------------
function out = local_label2char(x)
% Convert label to char safely; return '' if not a usable text label.
    out = '';
    if isempty(x)
        return;
    end
    if ischar(x)
        out = x;
        return;
    end
    if isstring(x)
        if isscalar(x)
            out = char(x);
        else
            out = char(x(1));
        end
        return;
    end
    % anything else (double, logical, struct, etc.)
    out = '';
end



%% ===================== 6) SUMMARY BARPLOT: # SIG NEURONS (EPOCH GLM) BY GROUP + BEFORE/AFTER =====================
if ~isempty(pLR)

    groups = {'posLag_posCorr','posLag_negCorr','negLag_posCorr','negLag_negCorr'};
    baList = {'before','after'};

    % Ensure catLabel and beforeAfter are strings
    catS = strings(nEpochUnits,1);
    baS  = strings(nEpochUnits,1);
    for i = 1:nEpochUnits
        if ischar(catLabel{i}) || (isstring(catLabel{i}) && isscalar(catLabel{i}))
            catS(i) = string(catLabel{i});
        else
            catS(i) = "";
        end

        if ischar(beforeAfter{i}) || (isstring(beforeAfter{i}) && isscalar(beforeAfter{i}))
            baS(i) = string(beforeAfter{i});
        else
            baS(i) = "";
        end
    end

    sigCounts = zeros(numel(groups), numel(baList));
    totCounts = zeros(numel(groups), numel(baList));

    for i = 1:nEpochUnits
        if catS(i)=="" || baS(i)=="", continue; end
        gi = find(strcmp(catS(i), groups), 1);
        bi = find(strcmp(baS(i),  baList), 1);
        if isempty(gi) || isempty(bi), continue; end

        totCounts(gi,bi) = totCounts(gi,bi) + 1;
        sigCounts(gi,bi) = sigCounts(gi,bi) + (sigLR(i)==1);
    end

    fig = figure('Color','w','Position',[100 100 1200 450]); hold on;
    bar(sigCounts);
    xticks(1:numel(groups)); xticklabels(groups); xtickangle(20);
    ylabel('# significant units (epoch GLM p<0.05)');
    legend({'before','after'},'Location','best');
    title('Epoch GLM significant neurons by xcorr category and peak timing');

    % annotate with sig/total
    for gi = 1:numel(groups)
        for bi = 1:numel(baList)
            text(gi, sigCounts(gi,bi) + 0.25, sprintf('%d/%d', sigCounts(gi,bi), totCounts(gi,bi)), ...
                'HorizontalAlignment','center', 'FontSize',10, 'FontWeight','bold');
        end
    end
end


%% ===================== SIG-ONLY BY CATEGORY (Epoch GLM sig_LR) =====================
% Uses:
%   - XCorrData(day).units(uu).CSM_epoch_glm.sig_LR  (epoch-level significance)
%   - XCorrData(day).units(uu).CSM_time_fast.delta_netFR + starts_ms (time-resolved effect)
%   - XCorrData(day).units(uu).category  (posLag_negCorr, etc.)
%   - XCorrData(day).units(uu).before_after (optional)

groups = {'posLag_posCorr','posLag_negCorr','negLag_posCorr','negLag_negCorr'};
baList = {'before','after'};   % optional split

alpha = 0.05;  % significance threshold already used in sig_LR

% ---------- local helper: safe label -> char ----------
label2char = @(x) local_label2char(x);

% ---------- containers ----------
Counts = struct();
for gi = 1:numel(groups)
    Counts.(groups{gi}).total = 0;
    Counts.(groups{gi}).sig   = 0;

    for bi = 1:numel(baList)
        Counts.(groups{gi}).(baList{bi}).total = 0;
        Counts.(groups{gi}).(baList{bi}).sig   = 0;
    end
end

SigMat = struct();   % store sig-only timecourses per category
AllMat = struct();   % store all timecourses per category (optional for comparison)
for gi = 1:numel(groups)
    SigMat.(groups{gi}).all = [];
    AllMat.(groups{gi}).all = [];

    for bi = 1:numel(baList)
        SigMat.(groups{gi}).(baList{bi}) = [];
        AllMat.(groups{gi}).(baList{bi}) = [];
    end
end

all_starts = [];  % common time base

% ---------- walk through all units ----------
for day = 1:numel(XCorrData)
    if ~isfield(XCorrData(day),'units'), continue; end

    for uu = 1:numel(XCorrData(day).units)
        u = XCorrData(day).units(uu);

        % need category + epoch significance info
        if ~isfield(u,'category'), continue; end
        if ~isfield(u,'CSM_epoch_glm') || ~isfield(u.CSM_epoch_glm,'skip') || u.CSM_epoch_glm.skip
            continue;
        end

        catC = label2char(u.category);
        if isempty(catC) || ~any(strcmp(catC, groups)), continue; end

        % before/after (optional)
        if isfield(u,'before_after')
            baC = label2char(u.before_after);
        else
            baC = '';
        end
        if ~any(strcmp(baC, baList))
            baC = ''; % treat as "unknown" -> only count overall
        end

        % significance
        isSig = false;
        if isfield(u.CSM_epoch_glm,'sig_LR')
            isSig = logical(u.CSM_epoch_glm.sig_LR);
        elseif isfield(u.CSM_epoch_glm,'p_LR')
            isSig = (u.CSM_epoch_glm.p_LR < alpha);
        end

        % update counts
        Counts.(catC).total = Counts.(catC).total + 1;
        Counts.(catC).sig   = Counts.(catC).sig   + double(isSig);

        if ~isempty(baC)
            Counts.(catC).(baC).total = Counts.(catC).(baC).total + 1;
            Counts.(catC).(baC).sig   = Counts.(catC).(baC).sig   + double(isSig);
        end

        % timecourse (optional, but needed for plots below)
        if ~(isfield(u,'CSM_time_fast') && isfield(u.CSM_time_fast,'skip') && ~u.CSM_time_fast.skip)
            continue;
        end

        s = u.CSM_time_fast.starts_ms(:)';
        d = u.CSM_time_fast.delta_netFR(:)';

        if isempty(all_starts)
            all_starts = s;
        end
        if ~(numel(s)==numel(all_starts) && all(s==all_starts))
            continue; % skip mismatched binning
        end

        % store ALL timecourses for that category (optional)
        AllMat.(catC).all = [AllMat.(catC).all; d]; 
        if ~isempty(baC)
            AllMat.(catC).(baC) = [AllMat.(catC).(baC); d]; 
        end

        % store SIG-only timecourses
        if isSig
            SigMat.(catC).all = [SigMat.(catC).all; d]; 
            if ~isempty(baC)
                SigMat.(catC).(baC) = [SigMat.(catC).(baC); d];
            end
        end
    end
end

% ===================== A) PRINT COUNTS (sig / total) =====================
fprintf('\n=== Significant units by category (epoch GLM sig_LR) ===\n');
for gi = 1:numel(groups)
    g = groups{gi};
    fprintf('%s: %d / %d sig\n', g, Counts.(g).sig, Counts.(g).total);
end

fprintf('\n=== Significant units by category AND before/after ===\n');
for gi = 1:numel(groups)
    g = groups{gi};
    for bi = 1:numel(baList)
        b = baList{bi};
        fprintf('%s - %s: %d / %d sig\n', g, b, Counts.(g).(b).sig, Counts.(g).(b).total);
    end
end

%% ===================== B) BARPLOT: SIG COUNTS BY CATEGORY =====================
sigCounts = zeros(numel(groups),1);
totCounts = zeros(numel(groups),1);
for gi = 1:numel(groups)
    g = groups{gi};
    sigCounts(gi) = Counts.(g).sig;
    totCounts(gi) = Counts.(g).total;
end

figure('Color','w','Position',[100 100 900 450]); hold on;
bar(sigCounts);
xticks(1:numel(groups)); xticklabels(groups); xtickangle(20);
ylabel('# significant units (epoch GLM p<0.05)');
title('Significant CSM-modulated units by xcorr category');

for gi = 1:numel(groups)
    text(gi, sigCounts(gi)+0.25, sprintf('%d/%d', sigCounts(gi), totCounts(gi)), ...
        'HorizontalAlignment','center','FontWeight','bold');
end

%% ===================== C) POPULATION TIMECOURSE (SIG-ONLY) BY CATEGORY =====================
% Plots mean ± SEM using ONLY epoch-significant units
figure('Color','w','Position',[100 100 1300 700]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

for gi = 1:numel(groups)
    g = groups{gi};
    nexttile; hold on;

    Dsig = SigMat.(g).all;
    if isempty(Dsig)
        title(sprintf('%s (NO SIG UNITS)', g), 'Interpreter','none');
        axis off;
        continue;
    end

    m  = mean(Dsig,1,'omitnan')';
    se = std(Dsig,0,1,'omitnan')' ./ sqrt(sum(~isnan(Dsig),1)');

    x  = all_starts(:);
    y1 = m - se;
    y2 = m + se;

    fill([x; flipud(x)], [y1; flipud(y2)], 0.92*[1 1 1], 'EdgeColor','none');
    plot(x, m, 'k', 'LineWidth', 2);
    xline(0,'k--','LineWidth',1);
    yline(0,'k-','LineWidth',1);

    title(sprintf('%s (sig N=%d)', g, size(Dsig,1)), 'Interpreter','none');
    xlabel('Time (ms, bin start)');
    ylabel('\Delta net FR (Hz)');
    xlim([-25 80])
end

%% ===================== D) OPTIONAL: SIG-ONLY HEATMAPS BY CATEGORY =====================
% Comment out if you don’t want heatmaps
for gi = 1:numel(groups)
    g = groups{gi};
    Dsig = SigMat.(g).all;
    if isempty(Dsig), continue; end

    % sort units by peak abs time
    [~, pk] = max(abs(Dsig), [], 2, 'omitnan');
    [~, ord] = sort(pk,'ascend');

    figure('Color','w','Position',[100 100 1100 550]);
    imagesc(all_starts, 1:size(Dsig,1), Dsig(ord,:));
    axis tight;
    colorbar;
    xline(0,'k--','LineWidth',1);
    xlabel('Time (ms, bin start)');
    ylabel('Sig units (sorted by peak time)');
    title(sprintf('SIG-ONLY heatmap: %s', g), 'Interpreter','none');
end


%% ===================== SCATTER: DELTA FR PER UNIT BY CATEGORY (SIG IN RED) =====================
groups = {'posLag_posCorr','posLag_negCorr','negLag_posCorr','negLag_negCorr'};
alpha  = 0.05;

cat_all   = {};
delta_all = [];
sig_all   = [];   % will store as logical

for day = 1:numel(XCorrData)
    if ~isfield(XCorrData(day),'units'), continue; end
    for uu = 1:numel(XCorrData(day).units)
        u = XCorrData(day).units(uu);

        % must have epoch glm
        if ~isfield(u,'CSM_epoch_glm') || ~isfield(u.CSM_epoch_glm,'skip') || u.CSM_epoch_glm.skip
            continue;
        end
        if ~isfield(u,'category'), continue; end

        catC = local_label2char(u.category);
        if isempty(catC) || ~any(strcmp(catC, groups))
            continue;
        end

        % delta FR
        if ~isfield(u.CSM_epoch_glm,'delta_netFR') || isempty(u.CSM_epoch_glm.delta_netFR) || ~isnumeric(u.CSM_epoch_glm.delta_netFR)
            continue;
        end
        d = double(u.CSM_epoch_glm.delta_netFR);
        if ~isfinite(d)
            continue;
        end

        % significance (FORCE LOGICAL)
        if isfield(u.CSM_epoch_glm,'sig_LR') && ~isempty(u.CSM_epoch_glm.sig_LR)
            isSig = logical(u.CSM_epoch_glm.sig_LR);
        elseif isfield(u.CSM_epoch_glm,'p_LR') && ~isempty(u.CSM_epoch_glm.p_LR)
            isSig = (double(u.CSM_epoch_glm.p_LR) < alpha);
        else
            isSig = false;
        end
        if ~isscalar(isSig), isSig = false; end  % safety

        % store
        cat_all{end+1,1}   = catC; 
        delta_all(end+1,1) = d;   
        sig_all(end+1,1)   = isSig;
    end
end

if isempty(delta_all)
    warning('No units found with CSM_epoch_glm.delta_netFR + valid category.');
else
    % map categories to x positions
    x_cat = nan(numel(cat_all),1);
    for i = 1:numel(cat_all)
        gi = find(strcmp(cat_all{i}, groups), 1);
        if isempty(gi), gi = NaN; end
        x_cat(i) = gi;
    end

    % keep only valid mapped points
    keep = isfinite(x_cat) & isfinite(delta_all) & ~isnan(sig_all);
    x_cat = x_cat(keep);
    delta_all = delta_all(keep);
    sig_all = sig_all(keep);

    % FORCE logical indices
    idx1 = (sig_all == true);
    idx0 = ~idx1;

    % jitter
    jitter = (rand(size(x_cat)) - 0.5) * 0.25;
    xj = x_cat + jitter;

    figure('Color','w','Position',[100 100 1200 520]); hold on;

    % non-sig points
    scatter(xj(idx0), delta_all(idx0), 28, 'filled', ...
        'MarkerFaceAlpha', 0.65, 'MarkerEdgeColor','none');

    % sig points (red)
    scatter(xj(idx1), delta_all(idx1), 40, 'filled', ...
        'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor','none');

    xlim([0.5 numel(groups)+0.5]);
    xticks(1:numel(groups));
    xticklabels(groups);
    xtickangle(20);

    yline(0,'k-','LineWidth',1);
    ylabel('\Delta net FR (Stim - BL), Hz');
    title(sprintf('Per-unit epoch effect by category (sig in red; p<%.2f)', alpha));

    % overlay means and annotate sig/total
    for gi = 1:numel(groups)
        mask = (x_cat == gi);
        if any(mask)
            m = mean(delta_all(mask),'omitnan');
            plot(gi, m, 'kd', 'MarkerFaceColor','k', 'MarkerSize',8);

            nTot = sum(mask);
            nSig = sum(mask & idx1);
            text(gi, m, sprintf('  %d/%d', nSig, nTot), ...
                'VerticalAlignment','bottom', 'HorizontalAlignment','left', 'FontWeight','bold');
        end
    end

    legend({'not sig','sig','category mean'}, 'Location','best');
end


%% ===================== TIME-RESOLVED SIGNIFICANCE (SUMMARY, NOT FULL TIMECOURSE) =====================
groups = {'posLag_posCorr','posLag_negCorr','negLag_posCorr','negLag_negCorr'};
alpha  = 0.05;

cat_all            = {};
peak_delta_abs     = [];
peak_delta_signed  = [];
anySig             = false(0,1);
peakSig            = false(0,1);
peakTime_ms        = [];

for day = 1:numel(XCorrData)
    if ~isfield(XCorrData(day),'units'), continue; end
    for uu = 1:numel(XCorrData(day).units)
        u = XCorrData(day).units(uu);

        if ~isfield(u,'CSM_time_fast') || ~isfield(u.CSM_time_fast,'skip') || u.CSM_time_fast.skip
            continue;
        end
        if ~isfield(u,'category'), continue; end

        catC = local_label2char(u.category);
        if isempty(catC) || ~any(strcmp(catC, groups)), continue; end

        d = u.CSM_time_fast.delta_netFR(:);
        p = u.CSM_time_fast.p_corr(:);
        t = u.CSM_time_fast.starts_ms(:);

        if isempty(d) || isempty(p) || numel(d)~=numel(p)
            continue;
        end

        % Find peak bin: max absolute delta, ignoring NaNs
        dabs = abs(d);
        dabs(isnan(dabs)) = -Inf;
        [~, pkIdx] = max(dabs);
        if isempty(pkIdx) || ~isfinite(dabs(pkIdx))
            continue;
        end

        pkAbs    = abs(d(pkIdx));
        pkSigned = d(pkIdx);

        % Any-bin significance (NaN-safe; no 'omitnan' needed)
        isAnySig = any((p < alpha) & ~isnan(p));

        % Is the peak bin significant?
        isPeakSig = false;
        if pkIdx <= numel(p) && ~isnan(p(pkIdx))
            isPeakSig = (p(pkIdx) < alpha);
        end

        % Store
        cat_all{end+1,1}           = catC; 
        peak_delta_abs(end+1,1)    = pkAbs; 
        peak_delta_signed(end+1,1) = pkSigned; 
        anySig(end+1,1)            = logical(isAnySig); 
        peakSig(end+1,1)           = logical(isPeakSig);
        if pkIdx <= numel(t)
            peakTime_ms(end+1,1)   = t(pkIdx); 
        else
            peakTime_ms(end+1,1)   = NaN; 
        end
    end
end

if isempty(peak_delta_signed)
    warning('No units found with CSM_time_fast results + valid category.');
else
    % Map categories -> x positions
    x_cat = nan(numel(cat_all),1);
    for i = 1:numel(cat_all)
        x_cat(i) = find(strcmp(cat_all{i}, groups), 1);
    end

    % Keep valid entries
    keep = isfinite(x_cat) & isfinite(peak_delta_signed) & ~isnan(anySig);
    x_cat = x_cat(keep);
    peak_delta_abs = peak_delta_abs(keep);
    peak_delta_signed = peak_delta_signed(keep);
    anySig = anySig(keep);
    peakSig = peakSig(keep);
    peakTime_ms = peakTime_ms(keep);

    % Jitter
    jitter = (rand(size(x_cat)) - 0.5) * 0.25;
    xj = x_cat + jitter;

    % --------- PLOT 1: signed peak delta (sig if ANY bin significant) ----------
    figure('Color','w','Position',[100 100 1200 520]); hold on;
    idxSig = (anySig == true);
    idxNo  = ~idxSig;

    scatter(xj(idxNo), peak_delta_signed(idxNo), 28, 'filled', ...
        'MarkerFaceAlpha', 0.60, 'MarkerEdgeColor','none');

    scatter(xj(idxSig), peak_delta_signed(idxSig), 42, 'filled', ...
        'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor','none');

    xlim([0.5 numel(groups)+0.5]);
    xticks(1:numel(groups)); xticklabels(groups); xtickangle(20);
    yline(0,'k-','LineWidth',1);
    ylabel('Peak \Delta net FR (Stim - BL) at max |effect|, Hz');
    title(sprintf('Time-resolved summary: per-unit peak effect (sig = ANY bin p_{corr}<%.2f)', alpha));

    for gi = 1:numel(groups)
        mask = (x_cat == gi);
        if any(mask)
            m = mean(peak_delta_signed(mask));
            plot(gi, m, 'kd', 'MarkerFaceColor','k', 'MarkerSize',8);

            nTot = sum(mask);
            nSig = sum(mask & idxSig);
            text(gi, m, sprintf('  %d/%d', nSig, nTot), ...
                'VerticalAlignment','bottom','HorizontalAlignment','left','FontWeight','bold');
        end
    end
    legend({'not sig','sig (any bin)','category mean'}, 'Location','best');

    % --------- PLOT 2: peak ABS delta (sig if PEAK bin significant) ----------
    figure('Color','w','Position',[100 100 1200 520]); hold on;
    idxSig2 = (peakSig == true);
    idxNo2  = ~idxSig2;

    scatter(xj(idxNo2), peak_delta_abs(idxNo2), 28, 'filled', ...
        'MarkerFaceAlpha', 0.60, 'MarkerEdgeColor','none');

    scatter(xj(idxSig2), peak_delta_abs(idxSig2), 42, 'filled', ...
        'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor','none');

    xlim([0.5 numel(groups)+0.5]);
    xticks(1:numel(groups)); xticklabels(groups); xtickangle(20);
    yline(0,'k-','LineWidth',1);
    ylabel('Peak | \Delta net FR |, Hz');
    title(sprintf('Time-resolved summary: per-unit peak magnitude (sig = PEAK bin p_{corr}<%.2f)', alpha));

    for gi = 1:numel(groups)
        mask = (x_cat == gi);
        if any(mask)
            m = mean(peak_delta_abs(mask));
            plot(gi, m, 'kd', 'MarkerFaceColor','k', 'MarkerSize',8);

            nTot = sum(mask);
            nSig = sum(mask & idxSig2);
            text(gi, m, sprintf('  %d/%d', nSig, nTot), ...
                'VerticalAlignment','bottom','HorizontalAlignment','left','FontWeight','bold');
        end
    end
    legend({'not sig','sig (peak bin)','category mean'}, 'Location','best');

    % --------- OPTIONAL: peak-time histogram (sig vs not) ----------
    figure('Color','w','Position',[100 100 900 450]); hold on;
    histogram(peakTime_ms(idxNo), 20);
    histogram(peakTime_ms(idxSig), 20);
    xline(0,'k--','LineWidth',1);
    xlabel('Time of peak |effect| (bin start ms)');
    ylabel('# units');
    title('Peak-time distribution (sig any-bin vs not)');
    legend({'not sig','sig (any bin)'}, 'Location','best');
end


%% ===================== TIME-RESOLVED SIGNED AUC (CSM epoch) + SIG HIGHLIGHT BY CATEGORY =====================
% For each unit with CSM_time_fast:
%   - Signed AUC over the CSM epoch using delta_netFR(t) across bins
%     AUC_units: (Hz * s) ~ spikes (baseline-subtracted extra spikes)
%   - Time-resolved significance: any bin with p_corr < alpha
% Then plot per category with sig units in red.

groups = {'posLag_posCorr','posLag_negCorr','negLag_posCorr','negLag_negCorr'};
alpha  = 0.05;

cat_all   = {};
auc_signed = [];
sig_anybin = false(0,1);   % logical
n_bins_used = [];
unit_id = {};              % optional identifiers

for day = 1:numel(XCorrData)
    if ~isfield(XCorrData(day),'units'), continue; end
    for uu = 1:numel(XCorrData(day).units)
        u = XCorrData(day).units(uu);

        % must have time-resolved results
        if ~isfield(u,'CSM_time_fast') || ~isfield(u.CSM_time_fast,'skip') || u.CSM_time_fast.skip
            continue;
        end
        if ~isfield(u,'category'), continue; end

        catC = local_label2char(u.category);
        if isempty(catC) || ~any(strcmp(catC, groups)), continue; end

        d = u.CSM_time_fast.delta_netFR(:);   % Hz per bin (Stim-BL net FR)
        p = u.CSM_time_fast.p_corr(:);        % corrected p per bin
        s = u.CSM_time_fast.starts_ms(:);     % bin starts (ms)

        if isempty(d) || isempty(s) || numel(d)~=numel(s)
            continue;
        end

        % Determine effective step size in seconds from starts_ms (robust)
        if numel(s) >= 2
            step_ms = median(diff(s(~isnan(s))));
        else
            % fallback: use stored step if available
            if isfield(u.CSM_time_fast,'step_ms')
                step_ms = u.CSM_time_fast.step_ms;
            else
                step_ms = NaN;
            end
        end
        if ~isfinite(step_ms) || step_ms <= 0
            continue;
        end
        step_s = step_ms / 1000;

        % Signed AUC over available (non-NaN) bins
        good = isfinite(d);
        if ~any(good)
            continue;
        end

        % Signed AUC (Hz * s); interpretable as baseline-subtracted "extra spikes" in epoch
        auc = sum(d(good)) * step_s;

        % Time-resolved significance: ANY corrected bin significant (NaN-safe)
        if ~isempty(p) && numel(p)==numel(d)
            isSig = any((p < alpha) & ~isnan(p));
        else
            isSig = false;
        end

        % store
        cat_all{end+1,1}   = catC; 
        auc_signed(end+1,1)= auc;  
        sig_anybin(end+1,1)= logical(isSig);
        n_bins_used(end+1,1)= sum(good);

        if isfield(u,'unitIdx')
            unit_id{end+1,1} = sprintf('D%02d-U%03d-Idx%04d', day, uu, u.unitIdx);
        else
            unit_id{end+1,1} = sprintf('D%02d-U%03d', day, uu);
        end
    end
end

if isempty(auc_signed)
    warning('No units found with CSM_time_fast.delta_netFR for AUC computation.');
else
    % map categories -> x positions
    x_cat = nan(numel(cat_all),1);
    for i = 1:numel(cat_all)
        x_cat(i) = find(strcmp(cat_all{i}, groups), 1);
    end

    % keep valid
    keep = isfinite(x_cat) & isfinite(auc_signed);
    x_cat = x_cat(keep);
    auc_signed = auc_signed(keep);
    sig_anybin = sig_anybin(keep);

    % jitter
    jitter = (rand(size(x_cat)) - 0.5) * 0.25;
    xj = x_cat + jitter;

    % plot
    figure('Color','w','Position',[100 100 1200 520]); hold on;

    idxSig = (sig_anybin == true);
    idxNo  = ~idxSig;

    % non-sig
    scatter(xj(idxNo), auc_signed(idxNo), 28, 'filled', ...
        'MarkerFaceAlpha', 0.60, 'MarkerEdgeColor','none');

    % sig (any bin)
    scatter(xj(idxSig), auc_signed(idxSig), 42, 'filled', ...
        'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor','none');

    xlim([0.5 numel(groups)+0.5]);
    xticks(1:numel(groups));
    xticklabels(groups);
    xtickangle(20);

    yline(0,'k-','LineWidth',1);
    ylabel('Signed AUC of \Delta net FR over CSM epoch (Hz·s ≈ spikes)');
    title(sprintf('Time-resolved summary: signed AUC (sig = ANY bin p_{corr}<%.2f)', alpha));

    % overlay category means + annotate sig/total
    for gi = 1:numel(groups)
        mask = (x_cat == gi);
        if any(mask)
            m = mean(auc_signed(mask));
            plot(gi, m, 'kd', 'MarkerFaceColor','k', 'MarkerSize',8);

            nTot = sum(mask);
            nSig = sum(mask & idxSig);
            text(gi, m, sprintf('  %d/%d', nSig, nTot), ...
                'VerticalAlignment','bottom','HorizontalAlignment','left','FontWeight','bold');
        end
    end

    legend({'not sig','sig (any bin)','category mean'}, 'Location','best');
end


%% NEW CODE TO PLOT MORE OF THE AUC

% ===================== TIME-RESOLVED SUMMARY PLOTS (AUC, absAUC, peak, peak time, RR/GLM) =====================
% Makes:
%  1) Signed AUC by category (sig ANY-bin in red)
%  2) Abs AUC by category (sig ANY-bin in red)
%  3) Peak signed delta by category (sig ANY-bin in red)
%  4) Peak |delta| by category (sig PEAK-bin in red)
%  5) Peak-time histogram (sig vs not)
%  6) Scatter: signed AUC vs epoch GLM RR (color: epoch sig, marker: time sig)
%  7) Volcano: epoch delta_netFR vs -log10(p_LR) (color: epoch sig)
%
% Requires per unit fields:
%   u.category
%   u.CSM_time_fast.delta_netFR
%   u.CSM_time_fast.p_corr
%   u.CSM_time_fast.starts_ms
% Optional:
%   u.CSM_epoch_glm.RR, u.CSM_epoch_glm.p_LR, u.CSM_epoch_glm.sig_LR, u.CSM_epoch_glm.delta_netFR

groups = {'posLag_posCorr','posLag_negCorr','negLag_posCorr','negLag_negCorr'};
alpha_time  = 0.05;
alpha_epoch = 0.05;

% ---------- collect per-unit summaries ----------
cat_all = {};
x_cat   = [];

auc_signed = [];
auc_abs    = [];
peak_signed = [];
peak_abs    = [];
peak_time_ms = [];

sig_anybin_time = false(0,1);
sig_peakbin_time = false(0,1);

% epoch GLM (optional)
has_epoch = false(0,1);
epoch_RR  = nan(0,1);
epoch_pLR = nan(0,1);
epoch_sig = false(0,1);
epoch_delta = nan(0,1);

for day = 1:numel(XCorrData)
    if ~isfield(XCorrData(day),'units'), continue; end
    for uu = 1:numel(XCorrData(day).units)
        u = XCorrData(day).units(uu);

        % ---------- need time-resolved data ----------
        if ~isfield(u,'CSM_time_fast') || ~isfield(u.CSM_time_fast,'skip') || u.CSM_time_fast.skip
            continue;
        end
        if ~isfield(u,'category'), continue; end

        catC = local_label2char(u.category);
        if isempty(catC) || ~any(strcmp(catC, groups)), continue; end

        d = u.CSM_time_fast.delta_netFR(:);   % Hz
        p = u.CSM_time_fast.p_corr(:);        % corrected p
        s = u.CSM_time_fast.starts_ms(:);     % ms

        if isempty(d) || isempty(s) || numel(d)~=numel(s)
            continue;
        end

        % step (s) for AUC integration
        if numel(s) >= 2
            step_ms = median(diff(s(~isnan(s))));
        elseif isfield(u.CSM_time_fast,'step_ms')
            step_ms = u.CSM_time_fast.step_ms;
        else
            step_ms = NaN;
        end
        if ~isfinite(step_ms) || step_ms <= 0, continue; end
        step_s = step_ms/1000;

        good = isfinite(d);
        if ~any(good), continue; end

        % peak indices ignoring NaN
        dabs = abs(d); dabs(~good) = -Inf;
        [~, pkIdx] = max(dabs);
        if isempty(pkIdx) || ~isfinite(dabs(pkIdx)), continue; end

        % --- compute summaries ---
        AUCs = sum(d(good)) * step_s;          % Hz*s ~ spikes
        AUCa = sum(abs(d(good))) * step_s;     % magnitude (no sign)
        Pks  = d(pkIdx);
        Pka  = abs(d(pkIdx));
        Pkt  = (pkIdx <= numel(s)) * s(pkIdx) + (pkIdx > numel(s)) * NaN;

        % time significance
        if ~isempty(p) && numel(p)==numel(d)
            anySig = any((p < alpha_time) & ~isnan(p));
            peakSig = ~isnan(p(pkIdx)) && (p(pkIdx) < alpha_time);
        else
            anySig = false;
            peakSig = false;
        end

        % epoch info (optional)
        e_has = false; e_RR = NaN; e_p = NaN; e_sig = false; e_del = NaN;
        if isfield(u,'CSM_epoch_glm') && isfield(u.CSM_epoch_glm,'skip') && ~u.CSM_epoch_glm.skip
            e_has = true;
            if isfield(u.CSM_epoch_glm,'RR'), e_RR = double(u.CSM_epoch_glm.RR); end
            if isfield(u.CSM_epoch_glm,'p_LR'), e_p = double(u.CSM_epoch_glm.p_LR); end
            if isfield(u.CSM_epoch_glm,'sig_LR'), e_sig = logical(u.CSM_epoch_glm.sig_LR); end
            if isfield(u.CSM_epoch_glm,'delta_netFR'), e_del = double(u.CSM_epoch_glm.delta_netFR); end
            if isfinite(e_p), e_sig = (e_p < alpha_epoch); end
        end

        % store
        cat_all{end+1,1} = catC; 
        x_cat(end+1,1)   = find(strcmp(catC, groups),1); 

        auc_signed(end+1,1) = AUCs; 
        auc_abs(end+1,1)    = AUCa; 
        peak_signed(end+1,1)= Pks; 
        peak_abs(end+1,1)   = Pka;  
        peak_time_ms(end+1,1)= Pkt; 

        sig_anybin_time(end+1,1)  = logical(anySig);  
        sig_peakbin_time(end+1,1) = logical(peakSig); 

        has_epoch(end+1,1)  = logical(e_has); 
        epoch_RR(end+1,1)   = e_RR; 
        epoch_pLR(end+1,1)  = e_p;  
        epoch_sig(end+1,1)  = logical(e_sig); 
        epoch_delta(end+1,1)= e_del; 
    end
end

if isempty(x_cat)
    warning('No units found with CSM_time_fast + category.');
    return;
end

% jitter
jitter = (rand(size(x_cat)) - 0.5) * 0.25;
xj = x_cat + jitter;

%% ===================== 1) SIGNED AUC BY CATEGORY =====================
local_cat_scatter(groups, x_cat, xj, auc_signed, sig_anybin_time, ...
    sprintf('Signed AUC of \\Delta net FR (Hz·s ≈ spikes); sig = ANY bin p_{corr}<%.2f', alpha_time), ...
    'Signed AUC (Hz·s ≈ spikes)');

%% ===================== 2) ABS AUC BY CATEGORY =====================
local_cat_scatter(groups, x_cat, xj, auc_abs, sig_anybin_time, ...
    sprintf('Abs AUC of |\\Delta net FR| (Hz·s); sig = ANY bin p_{corr}<%.2f', alpha_time), ...
    'Abs AUC (Hz·s)');

%% ===================== 3) PEAK SIGNED DELTA BY CATEGORY =====================
local_cat_scatter(groups, x_cat, xj, peak_signed, sig_anybin_time, ...
    sprintf('Peak \\Delta net FR at max |effect|; sig = ANY bin p_{corr}<%.2f', alpha_time), ...
    'Peak \Delta net FR (Hz)');

%% ===================== 4) PEAK |DELTA| BY CATEGORY (peak-bin sig) =====================
local_cat_scatter(groups, x_cat, xj, peak_abs, sig_peakbin_time, ...
    sprintf('Peak |\\Delta net FR|; sig = PEAK bin p_{corr}<%.2f', alpha_time), ...
    'Peak |\Delta net FR| (Hz)');

%% ===================== 5) PEAK TIME HISTOGRAM (sig vs not) =====================
figure('Color','w','Position',[100 100 950 420]); hold on;
idxSig = sig_anybin_time;
idxNo  = ~idxSig;
histogram(peak_time_ms(idxNo & isfinite(peak_time_ms)), 20);
histogram(peak_time_ms(idxSig & isfinite(peak_time_ms)), 20);
xline(0,'k--','LineWidth',1);
xlabel('Time of peak |effect| (bin start, ms)');
ylabel('# units');
title('Peak-time distribution (time-sig any-bin vs not)');
legend({'not sig','sig (any bin)'}, 'Location','best');

%% ===================== 6) AUC VS EPOCH GLM RR (epoch sig color, time sig marker) =====================
% Use only units with epoch info and finite RR
use = has_epoch & isfinite(epoch_RR) & epoch_RR>0 & isfinite(auc_signed);
if any(use)
    figure('Color','w','Position',[100 100 900 520]); hold on;

    % marker style for time sig vs not
    mkSig = 'o';
    mkNo  = 'o';

    % epoch non-sig (gray)
    idx0 = use & ~epoch_sig;
    scatter(log2(epoch_RR(idx0)), auc_signed(idx0), 36, 'filled', ...
        'MarkerFaceAlpha',0.6, 'MarkerEdgeColor','none');

    % epoch sig (red)
    idx1 = use & epoch_sig;
    scatter(log2(epoch_RR(idx1)), auc_signed(idx1), 44, 'filled', ...
        'MarkerFaceColor',[1 0 0], 'MarkerEdgeColor','none');

    % overlay time-sig ring
    ring = use & sig_anybin_time;
    scatter(log2(epoch_RR(ring)), auc_signed(ring), 70, ...
        'MarkerEdgeColor','k', 'LineWidth',1);

    xline(0,'k--'); yline(0,'k-');
    xlabel('log2(RR) from epoch Poisson GLM (Stim vs BL)');
    ylabel('Signed AUC (Hz·s ≈ spikes)');
    title('Consistency check: time-resolved AUC vs epoch GLM effect');
    legend({'epoch not sig','epoch sig','time-sig outline'}, 'Location','best');
else
    warning('No usable epoch GLM RR values found for AUC vs RR plot.');
end

%% ===================== 7) VOLCANO: epoch delta_netFR vs -log10(p_LR) =====================
useV = has_epoch & isfinite(epoch_delta) & isfinite(epoch_pLR) & epoch_pLR>0;
if any(useV)
    figure('Color','w','Position',[100 100 900 520]); hold on;

    yv = -log10(epoch_pLR(useV));
    xv = epoch_delta(useV);
    sigv = epoch_sig(useV);

    scatter(xv(~sigv), yv(~sigv), 30, 'filled', 'MarkerFaceAlpha',0.6, 'MarkerEdgeColor','none');
    scatter(xv(sigv),  yv(sigv),  42, 'filled', 'MarkerFaceColor',[1 0 0], 'MarkerEdgeColor','none');

    xline(0,'k-');
    yline(-log10(alpha_epoch),'k--');
    xlabel('\Delta net FR (epoch) (Stim-BL), Hz');
    ylabel('-log10(p_{LR})');
    title('Epoch GLM volcano: effect size vs significance');
    legend({'epoch not sig','epoch sig','p=0.05'}, 'Location','best');
else
    warning('No usable epoch GLM delta/p values found for volcano plot.');
end

%% ===================== LOCAL HELPERS (PUT AT END OF SCRIPT) =====================
function local_cat_scatter(groups, x_cat, xj, y, sigMask, figTitle, ylab)
    figure('Color','w','Position',[100 100 1200 520]); hold on;

    sigMask = (sigMask == true);
    idxNo = ~sigMask;
    idxSig = sigMask;

    scatter(xj(idxNo), y(idxNo), 28, 'filled', 'MarkerFaceAlpha',0.60, 'MarkerEdgeColor','none');
    scatter(xj(idxSig), y(idxSig), 42, 'filled', 'MarkerFaceColor',[1 0 0], 'MarkerEdgeColor','none');

    xlim([0.5 numel(groups)+0.5]);
    xticks(1:numel(groups));
    xticklabels(groups);
    xtickangle(20);
    yline(0,'k-','LineWidth',1);

    ylabel(ylab);
    title(figTitle);

    for gi = 1:numel(groups)
        mask = (x_cat == gi);
        if any(mask)
            m = mean(y(mask));
            plot(gi, m, 'kd', 'MarkerFaceColor','k', 'MarkerSize',8);
            nTot = sum(mask);
            nSig = sum(mask & sigMask);
            text(gi, m, sprintf('  %d/%d', nSig, nTot), ...
                'VerticalAlignment','bottom','HorizontalAlignment','left','FontWeight','bold');
        end
    end

    legend({'not sig','sig','category mean'}, 'Location','best');
end


%% PLOT to axtually use 
%% ===================== MAIN REVIEWER-FACING PLOT: SIGNED AUC BY CATEGORY =====================
groups = {'posLag_posCorr','posLag_negCorr','negLag_posCorr','negLag_negCorr'};
alpha  = 0.05;

cat_all = {};
auc_signed = [];
sig_anybin = false(0,1);

for day = 1:numel(XCorrData)
    if ~isfield(XCorrData(day),'units'), continue; end
    for uu = 1:numel(XCorrData(day).units)
        u = XCorrData(day).units(uu);

        % Must have time-resolved results
        if ~isfield(u,'CSM_time_fast') || u.CSM_time_fast.skip
            continue;
        end
        if ~isfield(u,'category'), continue; end

        catC = local_label2char(u.category);
        if isempty(catC) || ~any(strcmp(catC, groups)), continue; end

        d = u.CSM_time_fast.delta_netFR(:);   % Hz
        p = u.CSM_time_fast.p_corr(:);        % corrected p
        s = u.CSM_time_fast.starts_ms(:);     % bin starts (ms)

        if isempty(d) || numel(d) < 2, continue; end

        % Bin step (seconds)
        step_ms = median(diff(s));
        step_s  = step_ms / 1000;

        % Signed AUC (Hz*s ≈ extra spikes)
        good = isfinite(d);
        auc  = sum(d(good)) * step_s;

        % Time-resolved significance: ANY bin survives correction
        isSig = any((p < alpha) & ~isnan(p));

        % Store
        cat_all{end+1,1}   = catC; %#ok<AGROW>
        auc_signed(end+1,1)= auc;  %#ok<AGROW>
        sig_anybin(end+1,1)= logical(isSig); %#ok<AGROW>
    end
end

% Map categories to x
x_cat = nan(numel(cat_all),1);
for i = 1:numel(cat_all)
    x_cat(i) = find(strcmp(cat_all{i}, groups),1);
end

% Jitter for visibility
xj = x_cat + (rand(size(x_cat))-0.5)*0.25;

%% ---------- PLOT ----------
figure('Color','w','Position',[100 100 1200 520]); hold on;

idxSig = sig_anybin;
idxNo  = ~idxSig;

scatter(xj(idxNo), auc_signed(idxNo), 28, 'filled', ...
    'MarkerFaceAlpha',0.6,'MarkerEdgeColor','none');

scatter(xj(idxSig), auc_signed(idxSig), 42, 'filled', ...
    'MarkerFaceColor',[1 0 0],'MarkerEdgeColor','none');

xlim([0.5 numel(groups)+0.5]);
xticks(1:numel(groups));
xticklabels(groups);
xtickangle(20);
yline(0,'k-','LineWidth',1);

ylabel('Signed AUC of \Delta net FR (Hz·s ≈ extra spikes)');
title(sprintf('CSM modulation by category (red = time-resolved sig, p_{corr}<%.2f)', alpha));

% Category means + counts
for gi = 1:numel(groups)
    mask = (x_cat == gi);
    if any(mask)
        m = mean(auc_signed(mask));
        plot(gi, m, 'kd','MarkerFaceColor','k','MarkerSize',8);

        nTot = sum(mask);
        nSig = sum(mask & idxSig);
        text(gi, m, sprintf('  %d/%d', nSig, nTot), ...
            'FontWeight','bold','VerticalAlignment','bottom');
    end
end

legend({'not sig','sig','category mean'},'Location','best');

