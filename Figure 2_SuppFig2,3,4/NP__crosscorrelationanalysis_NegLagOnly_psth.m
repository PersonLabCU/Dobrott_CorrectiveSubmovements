


%% ============================================================
% integrated_xcorr_and_psth_classification_NEG_LAG_2CATS_CSMMOD_FULL.m
%
% FULLY INTEGRATED:
%   - XCorr: NEGATIVE LAG ONLY -> 2 categories:
%       1) negLag_posCorr
%       2) negLag_negCorr
%   - CSM modulation: epoch Poisson GLM LR-test (sig_LR)
%   - PSTH collection (Z-scored) for:
%       (A) ALL units (no keep filter)
%       (B) KEPT units only (max_vx peak within valid_peak_range_ms)
%     ...and for both:
%       - ALL units pooled + per-category
%       - CSM-tuned pooled + per-category
%   - PLOTTING:
%       1) Pooled PSTHs: ALL vs TUNED, and KEPT vs TUNED-KEPT
%       2) Category PSTHs + Category DELTAs for ALL units
%       3) Category PSTHs + Category DELTAs for KEPT units
%       4) Category scatter of delta_netFR (epoch) with labels Day+UnitID:
%            - ALL units
%            - KEPT units
%
% NOTES:
%   - No before/after designation anywhere.
%   - Uses AllCSM_csmT2 / NoCSM_csmT2 by default (edit fields below).
%   - delta_netFR is computed and stored in CSM_epoch_glm.delta_netFR.
% ============================================================

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

baseline_win_def     = [-1000 -500];  % zscore baseline window (ms)
%valid_peak_range_ms  = [-150 150];    % keep if max_vx PSTH peak falls here (ms)
valid_peak_range_ms  = [-250 250];    % keep if max_vx PSTH peak falls here (ms)

%%% ---------------- CSM MODULATION PARAMS ----------------
csm_epoch_win      = [-100 100];   % epoch for CSM modulation test
csm_baseline_win   = [-500 -250];  % baseline for epoch test
min_spikes_csm     = 5;
alpha_perm         = 0.05;

%% ---------------- PSTH FIELDS (EDIT AS NEEDED) ----------------
% CSM_RASTER_FIELD   = 'AllCSM_csmT2';
% NoCSM_RASTER_FIELD = 'NoCSM_csmT2';
% 


% CSM_RASTER_FIELD   = 'AllCSM_csmT0';
% NoCSM_RASTER_FIELD = 'NoCSM_csmT0';
% 

% Alternative:
CSM_RASTER_FIELD   = 'AllCSM_pos_thresh';
NoCSM_RASTER_FIELD = 'NoCSM_pos_thresh';


% ---------------- OUTPUT STRUCTS ----------------
XCorrData = struct;
NPdata    = struct;

% 2 categories only
group_names = {'negLag_posCorr','negLag_negCorr'};

% ---- pooled containers ----
AllAll   = struct('CSM',[],'NoCSM',[]);
AllKept  = struct('CSM',[],'NoCSM',[]);
TunedAll = struct('CSM',[],'NoCSM',[]);
TunedKept= struct('CSM',[],'NoCSM',[]);

% ---- per-category containers ----
AllGroupsAll    = struct();
AllGroupsKept   = struct();
TunedGroupsAll  = struct();
TunedGroupsKept = struct();

for g = 1:numel(group_names)
    gn = group_names{g};

    AllGroupsAll.(gn).CSM     = [];
    AllGroupsAll.(gn).NoCSM   = [];

    AllGroupsKept.(gn).CSM    = [];
    AllGroupsKept.(gn).NoCSM  = [];

    TunedGroupsAll.(gn).CSM   = [];
    TunedGroupsAll.(gn).NoCSM = [];

    TunedGroupsKept.(gn).CSM   = [];
    TunedGroupsKept.(gn).NoCSM = [];
end

timevec_global = [];
baseline_win_global = [];
valid_peak_range_global = [];

%% ---------------------- MAIN SESSION LOOP ----------------------
for day = 1:length(MyFolderInfo)

    fprintf('\n========== SESSION %d / %d ==========\n', day, length(MyFolderInfo));
    load(MyFolderInfo(day).name, 'Reaches','data','filtered_modulated_indices','ExcitRR');

    % Save kinematics (optional)
    NPdata(day).CSM   = ExcitRR.ReReaches;
    NPdata(day).NoCSM = ExcitRR.NoCSM;

    session_units = filtered_modulated_indices;
    if isempty(session_units)
        fprintf(' Day %d: no modulated units. Skipping.\n', day);
        continue;
    end

    nUnits = length(session_units);

    % timevec for PSTHs
    timevec = data.units(session_units(1)).max_vx.time_bins; % ms
    baseline_win     = find(timevec >= baseline_win_def(1)    & timevec <= baseline_win_def(2));
    valid_peak_range = find(timevec >= valid_peak_range_ms(1) & timevec <= valid_peak_range_ms(2));

    if isempty(timevec_global)
        timevec_global = timevec;
        baseline_win_global = baseline_win;
        valid_peak_range_global = valid_peak_range;
    end

    kept_units_session  = [];
    tuned_units_session = [];

    for uu = 1:nUnits

        unitIdx = session_units(uu);
        fprintf('  Analyzing UNIT %d...\n', unitIdx);

        % ----------------- LOAD SPIKE RASTER (for IFR/xCorr) -----------------
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

        % Resample IFR -> 200 Hz
        time_200Hz = (-win_pre_ms:1000/resample_fs:win_post_ms)';
        nSamples200 = length(time_200Hz);

        IFR = nan(nSamples200, nTrials);
        for tr = 1:nTrials
            IFR(:,tr) = interp1(ifr_time_1ms, IFR_1ms(:,tr), time_200Hz, 'linear', 0);
        end

        % ----------------- LOAD KINEMATICS ----------------------
        C = AlignPeakV_1ms_edit(Reaches, kinVar, win_pre_ms, win_post_ms);
        C = stimvel_accel(C);
        vel_1ms_all  = squeeze(C.stimvel(:,5,:));
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

                fr_sd = std(fr_seg); v_sd = std(v_seg);
                if fr_sd==0 || v_sd==0 || isnan(fr_sd) || isnan(v_sd), continue; end

                fr_z = (fr_seg - mean(fr_seg)) ./ fr_sd;
                v_z  = (v_seg  - mean(v_seg))  ./ v_sd;

                if any(isnan(fr_z)) || any(isnan(v_z)), continue; end
                xcorr_trials(li,tr) = mean(fr_z .* v_z);
            end
        end

        medianCorr = nanmedian(xcorr_trials,2);

        % ---------------- NEGATIVE-LAG ONLY: choose max-|corr| within lags<0 ------------
        negMask = (lags_ms < 0);
        negIdx  = find(negMask);

        negLag_peakLag = NaN;
        negLag_peakVal = NaN;
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

            if negLag_peakVal >= 0
                cat_label = 'negLag_posCorr';
            else
                cat_label = 'negLag_negCorr';
            end
        end

        % ---------------- KEEP FILTER (PSTH PEAK IN WINDOW) ----------------
        psth_ref = data.units(unitIdx).max_vx.PSTH_smoothbin;
        [~, peak_idx] = max(psth_ref);
        keep = ismember(peak_idx, valid_peak_range);

        % ---------------- CSM MODULATION TEST (EPOCH GLM) ----------------
        rng(0,'twister');

        isTuned = false;
        delta_netFR = NaN;

        try
            BL_raster   = data.units(unitIdx).(NoCSM_RASTER_FIELD).Raster;
            Stim_raster = data.units(unitIdx).(CSM_RASTER_FIELD).Raster;

            nTime_bl   = size(BL_raster,2);
            nTime_stim = size(Stim_raster,2);
            nTime_min  = min(nTime_bl, nTime_stim);
            BL_raster   = BL_raster(:,1:nTime_min);
            Stim_raster = Stim_raster(:,1:nTime_min);

            full_time_1ms = (-1000:1000)';
            if nTime_min > numel(full_time_1ms)
                error('Raster longer than assumed -1000:1000 ms vector.');
            end
            time_csm = full_time_1ms(1:nTime_min);

            csm_idx  = (time_csm >= csm_epoch_win(1))    & (time_csm <= csm_epoch_win(2));
            base_idx = (time_csm >= csm_baseline_win(1)) & (time_csm <= csm_baseline_win(2));

            BL_base_spk   = sum(BL_raster(:,base_idx),2);
            Stim_base_spk = sum(Stim_raster(:,base_idx),2);
            BL_csm_spk    = sum(BL_raster(:,csm_idx),2);
            Stim_csm_spk  = sum(Stim_raster(:,csm_idx),2);

            % compute delta_netFR regardless of GLM success (if windows valid)
            csm_win_s  = sum(csm_idx)  / 1000;
            base_win_s = sum(base_idx) / 1000;

            BL_netFR   = (BL_csm_spk./csm_win_s)   - (BL_base_spk./base_win_s);
            Stim_netFR = (Stim_csm_spk./csm_win_s) - (Stim_base_spk./base_win_s);
            delta_netFR = mean(Stim_netFR) - mean(BL_netFR);

            if ( (sum(BL_csm_spk) + sum(Stim_csm_spk)) < min_spikes_csm )
                CSM_epoch_glm = struct('skip',true,'reason','min_spikes_csm');
            else
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

                opts = statset('glmfit'); opts.MaxIter = 200;

                if sum(y)==0 || all(y==0)
                    CSM_epoch_glm = struct('skip',true,'reason','all_zero_counts');
                else
                    [~, dev_full] = glmfit(Xfull, y, 'poisson', 'link','log', 'options', opts);
                    [~, dev_red ] = glmfit(Xred,  y, 'poisson', 'link','log', 'options', opts);

                    LR   = dev_red - dev_full;
                    p_LR = 1 - chi2cdf(LR, 1);

                    CSM_epoch_glm.skip   = false;
                    CSM_epoch_glm.p_LR   = p_LR;
                    CSM_epoch_glm.sig_LR = (p_LR < alpha_perm);
                    CSM_epoch_glm.LR     = LR;

                    isTuned = CSM_epoch_glm.sig_LR;
                end
            end

        catch ME
            CSM_epoch_glm = struct('skip',true,'reason',['epoch_block_failed: ' ME.message]);
            isTuned = false;
        end

        % always attach delta_netFR if we computed it
        if ~isfield(CSM_epoch_glm,'delta_netFR') || isempty(CSM_epoch_glm.delta_netFR)
            CSM_epoch_glm.delta_netFR = delta_netFR;
        end

        % ---------------- PSTHs (Z-SCORED) ----------------
        % PSTH extraction can fail for some units if fields missing
        % got_psth = true;
        % try
        %     % csm_psth   = data.units(unitIdx).(CSM_PSTH_FIELD).PSTH_smoothbin;
        %     % nocsm_psth = data.units(unitIdx).(NoCSM_PSTH_FIELD).PSTH_smoothbin;
        % 
        %     csm_psth   = data.units(unitIdx).AllCSM_pos_thresh.PSTH_smoothbin;
        %     nocsm_psth = data.units(unitIdx).NoCSM_pos_thresh.PSTH_smoothbin;
        % catch
        %     got_psth = false;
        % end

         %Recomputing PSTH

         csm_raster   = data.units(unitIdx).AllCSM_pos_thresh.Raster;
         nocsm_raster = data.units(unitIdx).NoCSM_pos_thresh.Raster;

         % csm_raster   = data.units(unitIdx).AllCSM_csmT2.Raster;
         % nocsm_raster = data.units(unitIdx).NoCSM_csmT2.Raster;

         % csm_raster   = data.units(unitIdx).AllCSM_csmT0.Raster;
         % nocsm_raster = data.units(unitIdx).NoCSM_csmT0.Raster;

        bin_size = 10;
        smoothing_bins = 3; % originally 3 

        % --- compute PSTHs ---
        [psth_raw_csm,    csm_psth,    time_bins] = computePSTH(csm_raster,    bin_size, smoothing_bins);
        [psth_raw_nocsm,  nocsm_psth,    ~]         = computePSTH(nocsm_raster,    bin_size, smoothing_bins);
    
        % % --- row vectors ---
        % psth_norm    = psth_smooth_norm(:)';
        % psth_base    = psth_smooth_base(:)';
        % 
        % --- use returned time bins directly ---
        psth_timevec = time_bins(:)';

       
        baseline_mean = mean(csm_psth(baseline_win));
        baseline_std  = std(csm_psth(baseline_win));
        % 
        % csm_z   = (csm_psth   - baseline_mean) / baseline_std;
        % nocsm_z = (nocsm_psth - baseline_mean) / baseline_std;

        csm_z   = csm_psth ;
        nocsm_z = nocsm_psth;

            % --- pooled ALL units ---
            AllAll.CSM   = vertcat(AllAll.CSM,   csm_z);
            AllAll.NoCSM = vertcat(AllAll.NoCSM, nocsm_z);

            if ismember(cat_label, group_names)
                AllGroupsAll.(cat_label).CSM   = vertcat(AllGroupsAll.(cat_label).CSM,   csm_z);
                AllGroupsAll.(cat_label).NoCSM = vertcat(AllGroupsAll.(cat_label).NoCSM, nocsm_z);
            end

            % --- tuned ALL units ---
            if isTuned
                TunedAll.CSM   = vertcat(TunedAll.CSM,   csm_z);
                TunedAll.NoCSM = vertcat(TunedAll.NoCSM, nocsm_z);

                if ismember(cat_label, group_names)
                    TunedGroupsAll.(cat_label).CSM   = vertcat(TunedGroupsAll.(cat_label).CSM,   csm_z);
                    TunedGroupsAll.(cat_label).NoCSM = vertcat(TunedGroupsAll.(cat_label).NoCSM, nocsm_z);
                end
            end

            % --- kept subsets (only when keep==true) ---
            if keep
                AllKept.CSM   = vertcat(AllKept.CSM,   csm_z);
                AllKept.NoCSM = vertcat(AllKept.NoCSM, nocsm_z);

                if ismember(cat_label, group_names)
                    AllGroupsKept.(cat_label).CSM   = vertcat(AllGroupsKept.(cat_label).CSM,   csm_z);
                    AllGroupsKept.(cat_label).NoCSM = vertcat(AllGroupsKept.(cat_label).NoCSM, nocsm_z);
                end

                if isTuned
                    TunedKept.CSM   = vertcat(TunedKept.CSM,   csm_z);
                    TunedKept.NoCSM = vertcat(TunedKept.NoCSM, nocsm_z);

                    if ismember(cat_label, group_names)
                        TunedGroupsKept.(cat_label).CSM   = vertcat(TunedGroupsKept.(cat_label).CSM,   csm_z);
                        TunedGroupsKept.(cat_label).NoCSM = vertcat(TunedGroupsKept.(cat_label).NoCSM, nocsm_z);
                    end
                end
            end


        % ---------------- SAVE UNIT-LEVEL RESULTS ----------------
        XCorrData(day).units(uu).unitIdx        = unitIdx;
        XCorrData(day).units(uu).day            = day;

        XCorrData(day).units(uu).lags_ms        = lags_ms;
        XCorrData(day).units(uu).xcorr_trials   = xcorr_trials;
        XCorrData(day).units(uu).medianCorr     = medianCorr;

        XCorrData(day).units(uu).negLag_peakLag = negLag_peakLag;
        XCorrData(day).units(uu).negLag_peakVal = negLag_peakVal;
        XCorrData(day).units(uu).category       = cat_label;

        XCorrData(day).units(uu).keep           = keep;
        XCorrData(day).units(uu).psth_peak_time = timevec(peak_idx);

        XCorrData(day).units(uu).CSM_epoch_glm  = CSM_epoch_glm;
        XCorrData(day).units(uu).isCSMtuned     = isTuned;

        if keep
            kept_units_session(end+1) = unitIdx;
            if isTuned
                tuned_units_session(end+1) = unitIdx;
            end
        end

    end % unit loop

    XCorrData(day).kept_units  = kept_units_session;
    XCorrData(day).tuned_units = tuned_units_session;

    fprintf(' Session %d kept=%d | tuned_kept=%d\n', day, numel(kept_units_session), numel(tuned_units_session));
end

fprintf('\n================ ALL SESSIONS COMPLETE ================\n');


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
    %figure('Color','w','Position',[200 100 360 650]); hold on
    figure('Color','w','Position',[200 100 650 400]); hold on

    % faint row guides
    for i = 1:nU
        plot([valid_peak_range_ms(1) valid_peak_range_ms(2)], [i i], '-', ...
            'Color',[0.92 0.92 0.92], 'LineWidth',0.5);
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
            'Color', col, 'LineWidth', 3.0);
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

%% ===================== PEAK TIME SUMMARY STATS BY SUBSET =====================
% Requires:
%   PeakPlotTable with variables:
%       - peakTime_ms
%       - category
%
% Computes:
%   - N, mean, SD, SEM, median, IQR for each subset
%   - two-sample stats comparing:
%         negLag_posCorr  vs  negLag_negCorr
%   - effect size (Cohen's d)
%
% Outputs:
%   PeakTimeStatsTable
%   PeakTimeCompare

if ~exist('PeakPlotTable','var') || isempty(PeakPlotTable)
    error('PeakPlotTable not found. Run the peak timing collection block first.');
end

% ---- define subsets of interest ----
subset_names = {'negLag_posCorr','negLag_negCorr','other/unknown'};

% initialize
subset_label = {};
N_sub        = [];
mean_sub     = [];
sd_sub       = [];
sem_sub      = [];
median_sub   = [];
iqr_sub      = [];
min_sub      = [];
max_sub      = [];

for si = 1:numel(subset_names)

    this_subset = subset_names{si};

    if strcmp(this_subset,'other/unknown')
        idx = ~strcmp(PeakPlotTable.category,'negLag_posCorr') & ...
              ~strcmp(PeakPlotTable.category,'negLag_negCorr');
    else
        idx = strcmp(PeakPlotTable.category,this_subset);
    end

    vals = PeakPlotTable.peakTime_ms(idx);
    vals = vals(isfinite(vals));

    subset_label{end+1,1} = this_subset;
    N_sub(end+1,1)        = numel(vals);

    if isempty(vals)
        mean_sub(end+1,1)   = NaN;
        sd_sub(end+1,1)     = NaN;
        sem_sub(end+1,1)    = NaN;
        median_sub(end+1,1) = NaN;
        iqr_sub(end+1,1)    = NaN;
        min_sub(end+1,1)    = NaN;
        max_sub(end+1,1)    = NaN;
    else
        mean_sub(end+1,1)   = mean(vals,'omitnan');
        sd_sub(end+1,1)     = std(vals,'omitnan');
        sem_sub(end+1,1)    = std(vals,'omitnan') ./ sqrt(sum(~isnan(vals)));
        median_sub(end+1,1) = median(vals,'omitnan');
        iqr_sub(end+1,1)    = iqr(vals);
        min_sub(end+1,1)    = min(vals);
        max_sub(end+1,1)    = max(vals);
    end
end

% ---- table of descriptive stats ----
PeakTimeStatsTable = table( ...
    string(subset_label), N_sub, mean_sub, sd_sub, sem_sub, ...
    median_sub, iqr_sub, min_sub, max_sub, ...
    'VariableNames', {'Subset','N','Mean_ms','SD_ms','SEM_ms','Median_ms','IQR_ms','Min_ms','Max_ms'});

disp(' ');
disp('================ PEAK TIME SUMMARY STATS ================');
disp(PeakTimeStatsTable);

% ---- comparison: posCorr vs negCorr ----
pos_vals = PeakPlotTable.peakTime_ms(strcmp(PeakPlotTable.category,'negLag_posCorr'));
neg_vals = PeakPlotTable.peakTime_ms(strcmp(PeakPlotTable.category,'negLag_negCorr'));

pos_vals = pos_vals(isfinite(pos_vals));
neg_vals = neg_vals(isfinite(neg_vals));

PeakTimeCompare = struct();

if isempty(pos_vals) || isempty(neg_vals)
    warning('One or both comparison groups are empty. No inferential stats computed.');
else
    % Welch two-sample t-test (does not assume equal variance)
    [h_ttest,p_ttest,ci_ttest,stats_ttest] = ttest2(pos_vals, neg_vals, 'Vartype','unequal');

    % Nonparametric comparison
    p_ranksum = ranksum(pos_vals, neg_vals);

    % effect size: Cohen's d
    n1 = numel(pos_vals);
    n2 = numel(neg_vals);
    s1 = std(pos_vals,'omitnan');
    s2 = std(neg_vals,'omitnan');
    m1 = mean(pos_vals,'omitnan');
    m2 = mean(neg_vals,'omitnan');

    spooled = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / max(1,(n1+n2-2)));
    cohens_d = (m1 - m2) / spooled;

    % store results
    PeakTimeCompare.group1 = 'negLag_posCorr';
    PeakTimeCompare.group2 = 'negLag_negCorr';
    PeakTimeCompare.N1 = n1;
    PeakTimeCompare.N2 = n2;
    PeakTimeCompare.mean1_ms = m1;
    PeakTimeCompare.mean2_ms = m2;
    PeakTimeCompare.meanDiff_ms = m1 - m2;

    PeakTimeCompare.ttest_h = h_ttest;
    PeakTimeCompare.ttest_p = p_ttest;
    PeakTimeCompare.ttest_tstat = stats_ttest.tstat;
    PeakTimeCompare.ttest_df = stats_ttest.df;
    PeakTimeCompare.ttest_ci = ci_ttest;

    PeakTimeCompare.ranksum_p = p_ranksum;
    PeakTimeCompare.cohens_d = cohens_d;

    % print nicely
    fprintf('\n================ PEAK TIME COMPARISON ================\n');
    fprintf('Comparison: %s vs %s\n', PeakTimeCompare.group1, PeakTimeCompare.group2);
    fprintf('N = %d vs %d\n', PeakTimeCompare.N1, PeakTimeCompare.N2);
    fprintf('Mean peak time = %.2f ms vs %.2f ms\n', PeakTimeCompare.mean1_ms, PeakTimeCompare.mean2_ms);
    fprintf('Mean difference = %.2f ms\n', PeakTimeCompare.meanDiff_ms);
    fprintf('Welch t-test: t(%.2f) = %.3f, p = %.5f\n', ...
        PeakTimeCompare.ttest_df, PeakTimeCompare.ttest_tstat, PeakTimeCompare.ttest_p);
    fprintf('Rank-sum p = %.5f\n', PeakTimeCompare.ranksum_p);
    fprintf('Cohen''s d = %.3f\n', PeakTimeCompare.cohens_d);
end
%% ============================================================
% ============================ PLOTS ===========================
% ============================================================

% Use global timevec (from first session encountered)
timevec = timevec_global;

% -------------------- PSTH plotting helpers --------------------
Color_CSM   = rgb('DarkMagenta');
Color_NoCSM = rgb('Black');

full_timevec = timevec;
time_idx     = find(full_timevec >= -250 & full_timevec <= 250);
time_vector  = full_timevec(time_idx);

mean_sem = @(mat) deal( ...
    nanmean(mat(:,time_idx),1), ...
    nanstd(mat(:,time_idx),0,1) ./ sqrt(max(1,size(mat,1))) );

delta_sem = @(csm,no) deal( ...
    nanmean(csm(:,time_idx) - no(:,time_idx),1), ...
    nanstd(csm(:,time_idx) - no(:,time_idx),0,1) ./ sqrt(max(1,size(csm,1))) );

%% ============================================================
% FIG SET 0: POOLED PSTHs
%   - ALL units vs Tuned ALL
%   - KEPT units vs Tuned KEPT
% ============================================================

% --- pooled: ALL units ---
if ~isempty(AllAll.CSM) && ~isempty(AllAll.NoCSM)
    [mCSM, semCSM] = mean_sem(AllAll.CSM);
    [mNo,  semNo ] = mean_sem(AllAll.NoCSM);

    figure('Name','POOLED ALL UNITS','Position',[80,80,650,320]); hold on;
    plot(time_vector, mCSM,'Color',Color_CSM,'LineWidth',2);
    fill([time_vector fliplr(time_vector)], [mCSM+semCSM fliplr(mCSM-semCSM)], Color_CSM,'FaceAlpha',0.2,'EdgeColor','none');
    plot(time_vector, mNo,'Color',Color_NoCSM,'LineWidth',2);
    fill([time_vector fliplr(time_vector)], [mNo+semNo fliplr(mNo-semNo)], Color_NoCSM,'FaceAlpha',0.15,'EdgeColor','none');
    xline(0,'--k'); xlabel('Time from CSM T2 (ms)'); ylabel('Z-scored FR');
    xlim([-250 250]);
    title(sprintf('Pooled ALL units (N=%d)', size(AllAll.CSM,1)));
    box off;
end

% --- pooled: TUNED ALL units ---
if ~isempty(TunedAll.CSM) && ~isempty(TunedAll.NoCSM)
    [mCSM, semCSM] = mean_sem(TunedAll.CSM);
    [mNo,  semNo ] = mean_sem(TunedAll.NoCSM);

    figure('Name','POOLED TUNED (ALL UNITS)','Position',[80,420,650,320]); hold on;
    plot(time_vector, mCSM,'Color',Color_CSM,'LineWidth',2);
    fill([time_vector fliplr(time_vector)], [mCSM+semCSM fliplr(mCSM-semCSM)], Color_CSM,'FaceAlpha',0.2,'EdgeColor','none');
    plot(time_vector, mNo,'Color',Color_NoCSM,'LineWidth',2);
    fill([time_vector fliplr(time_vector)], [mNo+semNo fliplr(mNo-semNo)], Color_NoCSM,'FaceAlpha',0.15,'EdgeColor','none');
    xline(0,'--k'); xlabel('Time from CSM T2 (ms)'); ylabel('Z-scored FR');
    xlim([-250 250]);
    title(sprintf('Pooled TUNED (ALL units) (N=%d; p_{LR}<%.3f)', size(TunedAll.CSM,1), alpha_perm));
    box off;
end

% --- pooled: KEPT units ---
if ~isempty(AllKept.CSM) && ~isempty(AllKept.NoCSM)
    [mCSM, semCSM] = mean_sem(AllKept.CSM);
    [mNo,  semNo ] = mean_sem(AllKept.NoCSM);

    figure('Name','POOLED KEPT UNITS','Position',[780,80,650,320]); hold on;
    plot(time_vector, mCSM,'Color',Color_CSM,'LineWidth',2);
    fill([time_vector fliplr(time_vector)], [mCSM+semCSM fliplr(mCSM-semCSM)], Color_CSM,'FaceAlpha',0.2,'EdgeColor','none');
    plot(time_vector, mNo,'Color',Color_NoCSM,'LineWidth',2);
    fill([time_vector fliplr(time_vector)], [mNo+semNo fliplr(mNo-semNo)], Color_NoCSM,'FaceAlpha',0.15,'EdgeColor','none');
    xline(0,'--k'); xlabel('Time from CSM T2 (ms)'); ylabel('Z-scored FR');
    xlim([-250 250]);
    title(sprintf('Pooled KEPT units (N=%d)', size(AllKept.CSM,1)));
    box off;
end

% --- pooled: TUNED KEPT units ---
if ~isempty(TunedKept.CSM) && ~isempty(TunedKept.NoCSM)
    [mCSM, semCSM] = mean_sem(TunedKept.CSM);
    [mNo,  semNo ] = mean_sem(TunedKept.NoCSM);

    figure('Name','POOLED TUNED (KEPT UNITS)','Position',[780,420,650,320]); hold on;
    plot(time_vector, mCSM,'Color',Color_CSM,'LineWidth',2);
    fill([time_vector fliplr(time_vector)], [mCSM+semCSM fliplr(mCSM-semCSM)], Color_CSM,'FaceAlpha',0.2,'EdgeColor','none');
    plot(time_vector, mNo,'Color',Color_NoCSM,'LineWidth',2);
    fill([time_vector fliplr(time_vector)], [mNo+semNo fliplr(mNo-semNo)], Color_NoCSM,'FaceAlpha',0.15,'EdgeColor','none');
    xline(0,'--k'); xlabel('Time from CSM T2 (ms)'); ylabel('Z-scored FR');
    xlim([-250 250]);
    title(sprintf('Pooled TUNED (KEPT units) (N=%d; p_{LR}<%.3f)', size(TunedKept.CSM,1), alpha_perm));
    box off;
end

%% ============================================================
% FIG SET 1: CATEGORY PSTHs + DELTAS (ALL units)
%% ============================================================

DeltaColor_pos = rgb('SteelBlue'); % negLag_posCorr
DeltaColor_neg = rgb('DarkGreen'); % negLag_negCorr

for ci = 1:numel(group_names)
    cat = group_names{ci};

    if strcmp(cat,'negLag_posCorr')
        dcol = DeltaColor_pos;

    else
        dcol = DeltaColor_neg;
    end

    CSM_all = AllGroupsAll.(cat).CSM;
    No_all  = AllGroupsAll.(cat).NoCSM;

    if ~isempty(CSM_all) && ~isempty(No_all)

        [mCSM, semCSM] = mean_sem(CSM_all);
        [mNo,  semNo ] = mean_sem(No_all);

        figure('Name',['ALL PSTH ' cat]); hold on;
        plot(time_vector, mCSM,'Color',Color_CSM,'LineWidth',2);
        fill([time_vector fliplr(time_vector)], [mCSM+semCSM fliplr(mCSM-semCSM)], Color_CSM,'FaceAlpha',0.2,'EdgeColor','none');
        plot(time_vector, mNo,'Color',Color_NoCSM,'LineWidth',2);
        fill([time_vector fliplr(time_vector)], [mNo+semNo fliplr(mNo-semNo)], Color_NoCSM,'FaceAlpha',0.15,'EdgeColor','none');
        xline(0,'--k'); xlabel('Time from CSM T2 (ms)'); ylabel('Z-scored FR');
        xlim([-250 250]);
        title(sprintf('ALL units: %s (N=%d)', cat, size(CSM_all,1)),'Interpreter','none');
        box off;

        [mD, semD] = delta_sem(CSM_all, No_all);
        figure('Name',['ALL DELTA ' cat]); hold on;
        plot(time_vector, mD,'Color',dcol,'LineWidth',2);
        fill([time_vector fliplr(time_vector)], [mD+semD fliplr(mD-semD)], dcol,'FaceAlpha',0.2,'EdgeColor','none');
        xline(0,'--k'); yline(0,'--k'); yline(-2,':'); yline(2,':');
        xlabel('Time from CSM T2 (ms)'); ylabel('\Delta Z (CSM - NoCSM)');
        xlim([-250 250]);
        title(sprintf('ALL units Δ: %s (N=%d)', cat, size(CSM_all,1)),'Interpreter','none');
        box off;
    end
end

%% ============================================================
% FIG SET 2: CATEGORY PSTHs + DELTAS (KEPT units)
%% ============================================================

for ci = 1:numel(group_names)
    cat = group_names{ci};

    if strcmp(cat,'negLag_posCorr')
        dcol = DeltaColor_pos;
    else
        dcol = DeltaColor_neg;
    end

    CSM_k = AllGroupsKept.(cat).CSM;
    No_k  = AllGroupsKept.(cat).NoCSM;

    if ~isempty(CSM_k) && ~isempty(No_k)

        [mCSM, semCSM] = mean_sem(CSM_k);
        [mNo,  semNo ] = mean_sem(No_k);

        figure('Name',['KEPT PSTH ' cat],'Position',[100,100,600,300]); hold on;
        plot(time_vector, mCSM,'Color',Color_CSM,'LineWidth',2);
        fill([time_vector fliplr(time_vector)], [mCSM+semCSM fliplr(mCSM-semCSM)], Color_CSM,'FaceAlpha',0.2,'EdgeColor','none');
        plot(time_vector, mNo,'Color',Color_NoCSM,'LineWidth',2);
        fill([time_vector fliplr(time_vector)], [mNo+semNo fliplr(mNo-semNo)], Color_NoCSM,'FaceAlpha',0.15,'EdgeColor','none');
        xline(0,'--k'); 
       % xlabel('Time from CSM T2 (ms)');
        %xlabel('Time from CSM T0 (ms)');
        xlabel('Time from threshold crossing(ms)');
        ylabel('Z-scored FR');
        xlim([-250 250]);
        ylim([0 12])
        title(sprintf('KEPT units: %s (N=%d)', cat, size(CSM_k,1)),'Interpreter','none');
        box off;

        [mD, semD] = delta_sem(CSM_k, No_k);
        figure('Name',['KEPT DELTA ' cat],'Position',[100,100,600,300]); hold on;
        plot(time_vector, mD,'Color',dcol,'LineWidth',2);
        fill([time_vector fliplr(time_vector)], [mD+semD fliplr(mD-semD)], dcol,'FaceAlpha',0.2,'EdgeColor','none');
        xline(0,'--k'); yline(0,'--k'); yline(-2,':'); yline(2,':');
        %xlabel('Time from CSM T2 (ms)'); 
        %xlabel('Time from CSM T0 (ms)'); 
        xlabel('Time from threshold crossing (ms)'); 
        ylabel('\Delta Z (CSM - NoCSM)');
        xlim([-250 250]);
        ylim([-3 3])
        title(sprintf('KEPT units Δ: %s (N=%d)', cat, size(CSM_k,1)),'Interpreter','none');
        box off;
    end
end
%% ============================================================
% COMBINED DELTA PLOT (KEPT units — BOTH categories)
%% ============================================================

figure('Name','KEPT DELTA (both categories)','Position',[100,100,650,320]); 
hold on;

legend_entries = {};

for ci = 1:numel(group_names)

    cat = group_names{ci};

    % --- color selection ---
    if strcmp(cat,'negLag_posCorr')
        dcol = DeltaColor_pos;
    else
        dcol = DeltaColor_neg;
    end

    CSM_k = AllGroupsKept.(cat).CSM;
    No_k  = AllGroupsKept.(cat).NoCSM;

    if isempty(CSM_k) || isempty(No_k)
        fprintf('Skipping %s (empty)\n', cat);
        continue;
    end

    % --- delta ---
    [mD, semD] = delta_sem(CSM_k, No_k);

    % --- plot ---
    plot(time_vector, mD, 'Color', dcol, 'LineWidth', 2);

    fill([time_vector fliplr(time_vector)], ...
         [mD+semD fliplr(mD-semD)], ...
         dcol, 'FaceAlpha', 0.2, 'EdgeColor','none');

    legend_entries{end+1} = sprintf('%s (N=%d)', cat, size(CSM_k,1));

end

xline(0,'--k');
yline(0,'--k');
yline(-2,':');
yline(2,':');

xlabel('Time from threshold crossing (ms)');
ylabel('\Delta Z (CSM - NoCSM)');
xlim([-250 250]);
%ylim([-3 3]);

legend(legend_entries,'Location','best','Interpreter','none');
title('KEPT units Δ (both categories)');
box off;


%% ============================================================
% KEPT DELTA (both categories) with smoothing
% ============================================================

figure('Name','KEPT DELTA (both categories)','Position',[100,100,650,320]); 
hold on;

legend_entries = {};

% ---- smoothing settings ----
smooth_bins = 4;   % try 3, 5, or 7
smooth_plot = @(x) smoothdata(x, 'gaussian', smooth_bins);

for ci = 1:numel(group_names)

    cat = group_names{ci};

    % --- color selection ---
    if strcmp(cat,'negLag_posCorr')
        dcol = DeltaColor_pos;
    else
        dcol = DeltaColor_neg;
    end

    CSM_k = AllGroupsKept.(cat).CSM;
    No_k  = AllGroupsKept.(cat).NoCSM;

    if isempty(CSM_k) || isempty(No_k)
        fprintf('Skipping %s (empty)\n', cat);
        continue;
    end

    % --- delta ---
    [mD, semD] = delta_sem(CSM_k, No_k);

    % --- smooth mean and SEM for display ---
    mD_s   = smooth_plot(mD);
    semD_s = smooth_plot(semD);

    % --- plot ---
    fill([time_vector fliplr(time_vector)], ...
         [mD_s + semD_s, fliplr(mD_s - semD_s)], ...
         dcol, 'FaceAlpha', 0.2, 'EdgeColor','none');

    plot(time_vector, mD_s, 'Color', dcol, 'LineWidth', 2);

    legend_entries{end+1} = sprintf('%s (N=%d)', cat, size(CSM_k,1));

end

xline(0,'--k');
yline(0,'--k');
yline(-2,':');
yline(2,':');

xlabel('Time from threshold crossing (ms)');
ylabel('\Delta Z (CSM - NoCSM)');
xlim([-250 250]);
% ylim([-3 3]);

legend(legend_entries,'Location','best','Interpreter','none');
title('KEPT units \Delta (both categories)');
box off;

%% ============================================================
% KEPT DELTA (both categories)
%   - restricts data to same time window as time_vector
%   - smooths for display only
%   - finds max Z
%   - reports p-value at peak timepoint
% ============================================================

figure('Name','KEPT DELTA (both categories)','Position',[100,100,700,340]);
hold on;

legend_entries = {};
PeakStats = struct;

% ---------------- settings ----------------
smooth_bins = 4;   % display smoothing only
plot_win    = [-250 250];

toRow       = @(x) x(:)';
smooth_plot = @(x) smoothdata(toRow(x), 'gaussian', smooth_bins);

% -----------------------------------------------------------------
% IMPORTANT:
% full_timevec must correspond to the full number of columns in CSM/NoCSM
% Example:
%   full_timevec = timevec_global;
% -----------------------------------------------------------------
tv_full = toRow(full_timevec);

% time indices for plotting window
time_idx = find(tv_full >= plot_win(1) & tv_full <= plot_win(2));
tv       = tv_full(time_idx);

for ci = 1:numel(group_names)

    cat = group_names{ci};

    % ---------------- color ----------------
    if strcmp(cat,'negLag_posCorr')
        dcol = DeltaColor_pos;
    else
        dcol = DeltaColor_neg;
    end

    % ---------------- data ----------------
    CSM_k = AllGroupsKept.(cat).CSM;
    No_k  = AllGroupsKept.(cat).NoCSM;

    if isempty(CSM_k) || isempty(No_k)
        fprintf('Skipping %s (empty)\n', cat);
        continue;
    end

    % ---------------- safety checks ----------------
    if size(CSM_k,2) ~= size(No_k,2)
        error('CSM and NoCSM do not have the same number of time bins for %s.', cat);
    end

    if size(CSM_k,2) ~= numel(tv_full)
        error(['Number of columns in %s data (%d) does not match length of full_timevec (%d). ' ...
               'Make sure full_timevec is the full trace time vector.'], ...
               cat, size(CSM_k,2), numel(tv_full));
    end

    % ---------------- restrict to same plotted time window ----------------
    CSM_plot = CSM_k(:, time_idx);
    No_plot  = No_k(:,  time_idx);

    % rows = units, cols = time
    delta_units = CSM_plot - No_plot;

    % ---------------- mean + sem ----------------
    mD = nanmean(delta_units, 1);

    nPerT = sum(~isnan(delta_units), 1);
    semD  = nanstd(delta_units, 0, 1) ./ sqrt(max(nPerT,1));

    % force row vectors
    mD   = toRow(mD);
    semD = toRow(semD);

    % smooth for display only
    mD_s   = smooth_plot(mD);
    semD_s = smooth_plot(semD);

    % ---------------- peak ----------------
    [~, peak_idx] = max(abs(mD_s));
    peak_t = tv(peak_idx);
    peak_z = mD_s(peak_idx);

    % p-value at peak timepoint
    peak_vals = delta_units(:, peak_idx);
    peak_vals = peak_vals(~isnan(peak_vals));

    if numel(peak_vals) >= 2
        p_peak = signrank(peak_vals, 0);
    else
        p_peak = NaN;
    end

    % store
    PeakStats.(cat).peak_idx = peak_idx;
    PeakStats.(cat).peak_t   = peak_t;
    PeakStats.(cat).peak_z   = peak_z;
    PeakStats.(cat).p_peak   = p_peak;
    PeakStats.(cat).n        = size(delta_units,1);

    % ---------------- final size check before plotting ----------------
    if numel(tv) ~= numel(mD_s) || numel(tv) ~= numel(semD_s)
        error('Length mismatch before plotting for %s: tv=%d, mD_s=%d, semD_s=%d', ...
            cat, numel(tv), numel(mD_s), numel(semD_s));
    end

    % ---------------- plot ----------------
    fill([tv, fliplr(tv)], ...
         [mD_s + semD_s, fliplr(mD_s - semD_s)], ...
         dcol, ...
         'FaceAlpha', 0.2, ...
         'EdgeColor', 'none');

    plot(tv, mD_s, 'Color', dcol, 'LineWidth', 2);

    scatter(peak_t, peak_z, 55, 'filled', ...
        'MarkerFaceColor', dcol, ...
        'MarkerEdgeColor', 'k');

    peak_txt = sprintf('max Z = %.2f\np = %.3g', peak_z, p_peak);
    text(peak_t + 8, peak_z, peak_txt, ...
        'Color', dcol, ...
        'FontSize', 9, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle');

    legend_entries{end+1} = sprintf('%s (N=%d)', cat, size(delta_units,1));

end

% ---------------- reference lines ----------------
xline(0,'--k');
yline(0,'--k');
yline(-2,':');
yline(2,':');

xlabel('Time from threshold crossing (ms)');
ylabel('\Delta Z (CSM - NoCSM)');
xlim(plot_win);
title('KEPT units \Delta (both categories)');
legend(legend_entries, 'Location', 'best', 'Interpreter', 'none');
box off;

% ============================================================
% Print summary
% ============================================================
fprintf('\n===== Peak stats =====\n');
for ci = 1:numel(group_names)
    cat = group_names{ci};
    if isfield(PeakStats, cat)
        fprintf('%s: max Z = %.3f at %d ms, p = %.4g, N = %d\n', ...
            cat, ...
            PeakStats.(cat).peak_z, ...
            round(PeakStats.(cat).peak_t), ...
            PeakStats.(cat).p_peak, ...
            PeakStats.(cat).n);
    end
end

%% ============================================================
% KEPT DELTA (both categories)
%   - restricts data to same time window as full_timevec
%   - smooths mean/SEM for display only
%   - finds ONLY the single max peak per category
%   - reports peak Z ± SEM and p-value at that peak timepoint
% ============================================================

figure('Name','KEPT DELTA (both categories)','Position',[100,100,780,360]);
hold on;

legend_entries = {};
PeakStats = struct;

% ---------------- settings ----------------
smooth_bins = 5;      % display smoothing only
plot_win    = [-250 250];
minPeakDist = 15;     % only used if local peaks are found first

toRow       = @(x) x(:)';
smooth_plot = @(x) smoothdata(toRow(x), 'gaussian', smooth_bins);

% ------------------------------------------------------------
% full_timevec must match the full number of columns in CSM/NoCSM
% Example:
% full_timevec = timevec_global;
% ------------------------------------------------------------
tv_full = toRow(full_timevec);

% time indices for the plotting window
time_idx = find(tv_full >= plot_win(1) & tv_full <= plot_win(2));
tv       = tv_full(time_idx);

for ci = 1:numel(group_names)

    cat = group_names{ci};

    % ---------------- color ----------------
    if strcmp(cat,'negLag_posCorr')
        dcol = DeltaColor_pos;
    else
        dcol = DeltaColor_neg;
    end

    % ---------------- data ----------------
    CSM_k = AllGroupsKept.(cat).CSM;
    No_k  = AllGroupsKept.(cat).NoCSM;

    if isempty(CSM_k) || isempty(No_k)
        fprintf('Skipping %s (empty)\n', cat);
        continue;
    end

    % ---------------- safety checks ----------------
    if size(CSM_k,2) ~= size(No_k,2)
        error('CSM and NoCSM do not have the same number of time bins for %s.', cat);
    end

    if size(CSM_k,2) ~= numel(tv_full)
        error(['Number of columns in %s data (%d) does not match length of full_timevec (%d). ' ...
               'Make sure full_timevec is the full trace time vector.'], ...
               cat, size(CSM_k,2), numel(tv_full));
    end

    % ---------------- restrict to plot window ----------------
    CSM_plot = CSM_k(:, time_idx);
    No_plot  = No_k(:,  time_idx);

    % rows = units, cols = time
    delta_units = CSM_plot - No_plot;

    % ---------------- mean + sem from raw data ----------------
    mD    = nanmean(delta_units, 1);
    nPerT = sum(~isnan(delta_units), 1);
    semD  = nanstd(delta_units, 0, 1) ./ sqrt(max(nPerT,1));

    % force row vectors
    mD   = toRow(mD);
    semD = toRow(semD);

    % ---------------- smooth for display only ----------------
    mD_s   = smooth_plot(mD);
    semD_s = smooth_plot(semD);

    % upper/lower SEM envelopes
    upper_s = mD_s + semD_s;
    lower_s = mD_s - semD_s;

    % ---------------- final size check ----------------
    if numel(tv) ~= numel(mD_s) || numel(tv) ~= numel(semD_s)
        error('Length mismatch before plotting for %s: tv=%d, mD_s=%d, semD_s=%d', ...
            cat, numel(tv), numel(mD_s), numel(semD_s));
    end

    % ---------------- plot SEM and mean ----------------
    fill([tv fliplr(tv)], ...
         [upper_s fliplr(lower_s)], ...
         dcol, ...
         'FaceAlpha', 0.2, ...
         'EdgeColor', 'none');

    plot(tv, mD_s, 'Color', dcol, 'LineWidth', 2);

    % ========================================================
    % SINGLE MAX PEAK ONLY
    %   use the absolute max of the smoothed mean trace
    % ========================================================
    [~, peak_idx] = max(abs(mD_s));

    peak_t   = tv(peak_idx);
    peak_z   = mD_s(peak_idx);
    peak_sem = semD_s(peak_idx);
    peak_up  = upper_s(peak_idx);
    peak_lo  = lower_s(peak_idx);

    % p-value at this timepoint across units
    peak_vals = delta_units(:, peak_idx);
    peak_vals = peak_vals(~isnan(peak_vals));

    if numel(peak_vals) >= 2
        p_peak = signrank(peak_vals, 0);
    else
        p_peak = NaN;
    end

    % store
    PeakStats.(cat).peak_idx   = peak_idx;
    PeakStats.(cat).peak_t     = peak_t;
    PeakStats.(cat).peak_z     = peak_z;
    PeakStats.(cat).peak_sem   = peak_sem;
    PeakStats.(cat).peak_upper = peak_up;
    PeakStats.(cat).peak_lower = peak_lo;
    PeakStats.(cat).p_peak     = p_peak;
    PeakStats.(cat).n          = size(delta_units,1);

    % ---------------- plot marker at mean peak ----------------
    scatter(peak_t, peak_z, 60, 'filled', ...
        'MarkerFaceColor', dcol, ...
        'MarkerEdgeColor', 'k');

    % vertical line showing mean ± SEM at the peak
    plot([peak_t peak_t], [peak_lo peak_up], 'Color', dcol, 'LineWidth', 1.3);

    % label
    peak_txt = sprintf('%.2f \\pm %.2f\np=%.2g', peak_z, peak_sem, p_peak);
    text(peak_t + 6, peak_z, peak_txt, ...
        'Color', dcol, ...
        'FontSize', 8, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle');

    legend_entries{end+1} = sprintf('%s (N=%d)', cat, size(delta_units,1));

end

% ---------------- reference lines ----------------
xline(0,'--k');
yline(0,'--k');
yline(-2,':');
yline(2,':');

% ---------------- labels ----------------
xlabel('Time from threshold crossing (ms)');
ylabel('\Delta Z (CSM - NoCSM)');
xlim(plot_win);
title('KEPT units \Delta (both categories)');
legend(legend_entries, 'Location', 'best', 'Interpreter', 'none');
box off;

% ============================================================
% Print summary to command window
% ============================================================
fprintf('\n===== Peak stats =====\n');
for ci = 1:numel(group_names)

    cat = group_names{ci};

    if ~isfield(PeakStats, cat)
        continue;
    end

    fprintf('%s (N = %d): Z = %.3f +/- %.3f at %d ms (upper = %.3f, lower = %.3f), p = %.4g\n', ...
        cat, ...
        PeakStats.(cat).n, ...
        PeakStats.(cat).peak_z, ...
        PeakStats.(cat).peak_sem, ...
        round(PeakStats.(cat).peak_t), ...
        PeakStats.(cat).peak_upper, ...
        PeakStats.(cat).peak_lower, ...
        PeakStats.(cat).p_peak);
end
%% ============================================================
% KEPT DELTA (both categories)
%   - restricts data to same time window as full_timevec
%   - smooths mean/SEM for display only
%   - shows multiple top peak values on the plot
%   - reports p-value at each peak timepoint
% ============================================================

figure('Name','KEPT DELTA (both categories)','Position',[100,100,760,360]);
hold on;

legend_entries = {};
PeakStats = struct;

% ---------------- settings ----------------
smooth_bins = 5;      % display smoothing only
plot_win    = [-250 250];
nPeaksShow  = 3;      % number of top peaks to display per category
minPeakDist = 15;     % in samples, to avoid peaks too close together

toRow       = @(x) x(:)';
smooth_plot = @(x) smoothdata(toRow(x), 'gaussian', smooth_bins);

% ------------------------------------------------------------
% full_timevec must match the full number of columns in CSM/NoCSM
% Example:
% full_timevec = timevec_global;
% ------------------------------------------------------------
tv_full = toRow(full_timevec);

% time indices for the plotting window
time_idx = find(tv_full >= plot_win(1) & tv_full <= plot_win(2));
tv       = tv_full(time_idx);

for ci = 1:numel(group_names)

    cat = group_names{ci};

    % ---------------- color ----------------
    if strcmp(cat,'negLag_posCorr')
        dcol = DeltaColor_pos;
    else
        dcol = DeltaColor_neg;
    end

    % ---------------- data ----------------
    CSM_k = AllGroupsKept.(cat).CSM;
    No_k  = AllGroupsKept.(cat).NoCSM;

    if isempty(CSM_k) || isempty(No_k)
        fprintf('Skipping %s (empty)\n', cat);
        continue;
    end

    % ---------------- safety checks ----------------
    if size(CSM_k,2) ~= size(No_k,2)
        error('CSM and NoCSM do not have the same number of time bins for %s.', cat);
    end

    if size(CSM_k,2) ~= numel(tv_full)
        error(['Number of columns in %s data (%d) does not match length of full_timevec (%d). ' ...
               'Make sure full_timevec is the full trace time vector.'], ...
               cat, size(CSM_k,2), numel(tv_full));
    end

    % ---------------- restrict to plot window ----------------
    CSM_plot = CSM_k(:, time_idx);
    No_plot  = No_k(:,  time_idx);

    % rows = units, cols = time
    delta_units = CSM_plot - No_plot;

    % ---------------- mean + sem from raw data ----------------
    mD    = nanmean(delta_units, 1);
    nPerT = sum(~isnan(delta_units), 1);
    semD  = nanstd(delta_units, 0, 1) ./ sqrt(max(nPerT,1));

    % force row vectors
    mD   = toRow(mD);
    semD = toRow(semD);

    % ---------------- smooth for display only ----------------
    mD_s   = smooth_plot(mD);
    semD_s = smooth_plot(semD);

    % ---------------- final size check ----------------
    if numel(tv) ~= numel(mD_s) || numel(tv) ~= numel(semD_s)
        error('Length mismatch before plotting for %s: tv=%d, mD_s=%d, semD_s=%d', ...
            cat, numel(tv), numel(mD_s), numel(semD_s));
    end

    % ---------------- plot SEM and mean ----------------
    fill([tv fliplr(tv)], ...
         [mD_s + semD_s, fliplr(mD_s - semD_s)], ...
         dcol, ...
         'FaceAlpha', 0.2, ...
         'EdgeColor', 'none');

    plot(tv, mD_s, 'Color', dcol, 'LineWidth', 2);

    % ========================================================
    % MULTIPLE TOP PEAKS
    %   find peaks on absolute smoothed trace so both positive
    %   and negative deflections can be captured
    % ========================================================
    [pks_abs, locs] = findpeaks(abs(mD_s), 'MinPeakDistance', minPeakDist);

    % fallback in case no local peaks are found
    if isempty(locs)
        [~, locs] = max(abs(mD_s));
        pks_abs   = abs(mD_s(locs));
    end

    % sort by peak magnitude
    [~, sort_idx] = sort(pks_abs, 'descend');

    % keep top N
    nUse = min(nPeaksShow, numel(sort_idx));
    top_locs = locs(sort_idx(1:nUse));

    % sort selected peaks in time order for nicer labeling/storage
    [top_locs, order_idx] = sort(top_locs, 'ascend');

    % initialize storage
    PeakStats.(cat).peak_idx = nan(1,nUse);
    PeakStats.(cat).peak_t   = nan(1,nUse);
    PeakStats.(cat).peak_z   = nan(1,nUse);
    PeakStats.(cat).p_peak   = nan(1,nUse);
    PeakStats.(cat).n        = size(delta_units,1);

    for pi = 1:nUse

        peak_idx = top_locs(pi);
        peak_t   = tv(peak_idx);
        peak_z   = mD_s(peak_idx);

        % p-value at this timepoint across units
        peak_vals = delta_units(:, peak_idx);
        peak_vals = peak_vals(~isnan(peak_vals));

        if numel(peak_vals) >= 2
            p_peak = signrank(peak_vals, 0);
        else
            p_peak = NaN;
        end

        % store
        PeakStats.(cat).peak_idx(pi) = peak_idx;
        PeakStats.(cat).peak_t(pi)   = peak_t;
        PeakStats.(cat).peak_z(pi)   = peak_z;
        PeakStats.(cat).p_peak(pi)   = p_peak;

        % plot marker
        scatter(peak_t, peak_z, 55, 'filled', ...
            'MarkerFaceColor', dcol, ...
            'MarkerEdgeColor', 'k');

        % stagger text labels a little to reduce overlap
        y_offset = 0.18 * (-1)^(pi+1);

        peak_txt = sprintf('%.2f\np=%.2g', peak_z, p_peak);
        text(peak_t + 6, peak_z + y_offset, peak_txt, ...
            'Color', dcol, ...
            'FontSize', 8, ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle');
    end

    legend_entries{end+1} = sprintf('%s (N=%d)', cat, size(delta_units,1));

end

% ---------------- reference lines ----------------
xline(0,'--k');
yline(0,'--k');
yline(-2,':');
yline(2,':');

% ---------------- labels ----------------
xlabel('Time from threshold crossing (ms)');
ylabel('\Delta Z (CSM - NoCSM)');
xlim(plot_win);
title('KEPT units \Delta (both categories)');
legend(legend_entries, 'Location', 'best', 'Interpreter', 'none');
box off;
% ============================================================
% Print summary to command window
% ============================================================
fprintf('\n===== Peak stats =====\n');
for ci = 1:numel(group_names)

    cat = group_names{ci};

    if ~isfield(PeakStats, cat)
        continue;
    end

    fprintf('\n%s (N = %d)\n', cat, PeakStats.(cat).n);

    for pi = 1:numel(PeakStats.(cat).peak_t)
        fprintf('  Peak %d: Z = %.3f at %d ms, p = %.4g\n', ...
            pi, ...
            PeakStats.(cat).peak_z(pi), ...
            round(PeakStats.(cat).peak_t(pi)), ...
            PeakStats.(cat).p_peak(pi));
    end
end
%% ============================================================
% COMBINED PLOT: Δ Neural (left y) + Δ Behavior (right y)
% Illustrator-editable export (vector PDF): ax2 (behavior) FIRST, ax1 (neural) SECOND
%% ============================================================

% ---------- SETTINGS ----------
xlim_ms = [-250 250];

% Behavior xlines originally in seconds:
xlines_behavior_s = [-0.055, 0.087, 0];           % sec
xlines_behavior_ms = xlines_behavior_s * 1000;    % ms

% If you want the event line at t=0 and horizontal refs on neural:
do_neural_ref_lines = true;

% ---------- PREP BEHAVIOR DELTA (CSM - NoCSM) ----------
% Column indices you were using:
colY = 8;  % e.g., y-velocity is column 8 in Average/Std (your code)

diff_avg_y = RR_big.Average(:,colY) - NRR_big.Average(:,colY);
diff_avg_y = diff_avg_y * 10; % your scaling

SE_diff_y  = sqrt((RR_big.Std(:,colY).^2 ./ RR_big.N) + (NRR_big.Std(:,colY).^2 ./ NRR_big.N));

time_s  = RR_big.Average(:,1);          % seconds
time_ms = (time_s(:) * 1000)';          % row vector in ms

% Make sure behavior vectors are row vectors for fill()
diff_avg_y = diff_avg_y(:)'; 
SE_diff_y  = SE_diff_y(:)';

% Optional: window behavior to match x-lims (keeps plots tidy)
beh_keep = (time_ms >= xlim_ms(1)) & (time_ms <= xlim_ms(2));
time_ms_b = time_ms(beh_keep);
diff_b    = diff_avg_y(beh_keep);
SE_b      = SE_diff_y(beh_keep);

% ---------- FIGURE ----------
%fig = figure('Name','Δ Neural + Δ Behavior','Color','w','Position',[100,100,700,340]);
fig = figure('Name','Δ Neural + Δ Behavior','Color','w','Position',[100,100,650,320]);

% ===================== ax2: BEHAVIOR (RIGHT Y) FIRST =====================
ax2 = axes('Parent',fig);
hold(ax2,'on');

% Shaded SEM for behavior (vector-friendly fill)
fill(ax2, ...
    [time_ms_b fliplr(time_ms_b)], ...
    [diff_b+SE_b fliplr(diff_b-SE_b)], ...
    'k', 'FaceAlpha',0.18, 'EdgeColor','none');

% Mean line
plot(ax2, time_ms_b, diff_b, 'k', 'LineWidth',2);

% ax2 formatting
ax2.YAxisLocation = 'right';
ax2.Box = 'off';
ax2.XLim = xlim_ms;
ax2.XColor = 'k';
ax2.YColor = 'k';
ylabel(ax2, 'Δ Behavior (CSM - NoCSM)');

% Behavior reference lines (convert your seconds -> ms)
yline(ax2, 0, '--k');
for xl = xlines_behavior_ms
    xline(ax2, xl, '--k');
end

% Hide ax2 x-labels/ticks (ax1 will own x-axis visually)
ax2.XTickLabel = [];
ax2.XLabel.String = '';

% ===================== ax1: NEURAL (LEFT Y) SECOND (OVERLAY) =====================
ax1 = axes('Parent',fig);
hold(ax1,'on');

legend_entries = {};

for ci = 1:numel(group_names)

    cat = group_names{ci};

    % --- color selection ---
    if strcmp(cat,'negLag_posCorr')
        dcol = DeltaColor_pos;
    else
        dcol = DeltaColor_neg;
    end

    CSM_k = AllGroupsKept.(cat).CSM;
    No_k  = AllGroupsKept.(cat).NoCSM;

    if isempty(CSM_k) || isempty(No_k)
        fprintf('Skipping %s (empty)\n', cat);
        continue;
    end

    % --- delta ---
    [mD, semD] = delta_sem(CSM_k, No_k);

    % --- plot neural mean + SEM ---
    plot(ax1, time_vector, mD, 'Color', dcol, 'LineWidth',2);

    fill(ax1, [time_vector fliplr(time_vector)], ...
              [mD+semD      fliplr(mD-semD)], ...
              dcol, 'FaceAlpha',0.20, 'EdgeColor','none');

    legend_entries{end+1} = sprintf('%s (N=%d)', cat, size(CSM_k,1));
end

% ax1 formatting
ax1.Color = 'none';          % critical so ax1 doesn't cover ax2 (Illustrator-friendly)
ax1.Box   = 'off';
ax1.XLim  = xlim_ms;

xlabel(ax1, 'Time from threshold crossing (ms)');
ylabel(ax1, '\Delta Z (CSM - NoCSM)');

% Neural reference lines
if do_neural_ref_lines
    xline(ax1, 0, '--k');
    yline(ax1, 0, '--k');
    yline(ax1, -2, ':k');
    yline(ax1,  2, ':k');
end

% Put both axes in same position + link x
ax1.Position = ax2.Position;
linkaxes([ax1 ax2],'x');

% Legend on ax1 (neural categories)
if ~isempty(legend_entries)
    legend(ax1, legend_entries, 'Location','best', 'Interpreter','none');
end

ylim([-3 4])
title(ax1, 'Δ Neural (left) + Δ Behavior (right)');

% ---------- EXPORT (Illustrator-editable vector PDF) ----------
set(fig,'Renderer','painters');

% If you want to save, uncomment and set output_folder/filename:
% if ~exist(output_folder,'dir'); mkdir(output_folder); end
% outFile = fullfile(output_folder, 'CombinedDelta_NeuralBehavior.pdf');
% print(fig, outFile, '-dpdf', '-painters');

%% ============================================================
% COMBINED DELTA PLOT: NEURAL (left y) + BEHAVIOR (right y)
% FULL REWRITE WITH ZERO ALIGNED ACROSS BOTH Y-AXES
% Illustrator-editable (painters renderer)
% ============================================================

% -------------------------- SETTINGS --------------------------
xlim_ms = [-250 250];

% Behavior event lines originally in seconds
xlines_behavior_s  = [-0.055, 0.087, 0];
xlines_behavior_ms = xlines_behavior_s * 1000;

do_neural_ref_lines = true;
colY = 8;   % behavior column to use

% ------------------ PREP BEHAVIOR DELTA ----------------------
% Δ behavior = CSM - NoCSM
diff_avg_y = RR_big.Average(:,colY) - NRR_big.Average(:,colY);
diff_avg_y = diff_avg_y * 10;   % keep your scaling

SE_diff_y = sqrt( ...
    (RR_big.Std(:,colY).^2  ./ RR_big.N) + ...
    (NRR_big.Std(:,colY).^2 ./ NRR_big.N) );

time_s  = RR_big.Average(:,1);
time_ms = time_s(:)' * 1000;

diff_avg_y = diff_avg_y(:)';
SE_diff_y  = SE_diff_y(:)';

% Window behavior to x-lims
beh_keep  = (time_ms >= xlim_ms(1)) & (time_ms <= xlim_ms(2));
time_ms_b = time_ms(beh_keep);
diff_b    = diff_avg_y(beh_keep);
SE_b      = SE_diff_y(beh_keep);

% -------------------------- FIGURE ---------------------------
fig = figure('Name','Δ Neural + Δ Behavior', ...
             'Color','w', ...
             'Position',[100,100,650,320]);
set(fig,'Renderer','painters');

% ============================================================
% AX2: BEHAVIOR (RIGHT Y-AXIS) FIRST
% ============================================================
ax2 = axes('Parent',fig);
hold(ax2,'on');
box(ax2,'off');

% Behavior SEM
fill(ax2, ...
    [time_ms_b fliplr(time_ms_b)], ...
    [diff_b + SE_b, fliplr(diff_b - SE_b)], ...
    'k', ...
    'FaceAlpha',0.18, ...
    'EdgeColor','none');

% Behavior mean
plot(ax2, time_ms_b, diff_b, 'k', 'LineWidth', 2);

% Behavior reference lines
yline(ax2, 0, '--k', 'LineWidth', 1);
for xl = xlines_behavior_ms
    xline(ax2, xl, '--k', 'LineWidth', 1);
end

% Behavior axis formatting
ax2.YAxisLocation = 'right';
ax2.XLim          = xlim_ms;
ax2.XColor        = 'k';
ax2.YColor        = 'k';
ax2.LineWidth     = 1;
ax2.FontSize      = 11;
ax2.TickDir       = 'out';
ax2.XTickLabel    = [];
ax2.XLabel.String = '';
ylabel(ax2, 'Δ Behavior (CSM - NoCSM)');

% ============================================================
% AX1: NEURAL (LEFT Y-AXIS) SECOND
% ============================================================
ax1 = axes('Parent',fig);
hold(ax1,'on');
box(ax1,'off');

legend_entries = {};

for ci = 1:numel(group_names)

    cat = group_names{ci};

    % Color by category
    if strcmp(cat, 'negLag_posCorr')
        dcol = DeltaColor_pos;
    else
        dcol = DeltaColor_neg;
    end

    CSM_k = AllGroupsKept.(cat).CSM;
    No_k  = AllGroupsKept.(cat).NoCSM;

    if isempty(CSM_k) || isempty(No_k)
        fprintf('Skipping %s (empty)\n', cat);
        continue;
    end

    % Neural delta
    [mD, semD] = delta_sem(CSM_k, No_k);

    % SEM patch first
    fill(ax1, ...
        [time_vector fliplr(time_vector)], ...
        [mD + semD, fliplr(mD - semD)], ...
        dcol, ...
        'FaceAlpha',0.20, ...
        'EdgeColor','none');

    % Mean trace
    plot(ax1, time_vector, mD, 'Color', dcol, 'LineWidth', 2);

    legend_entries{end+1} = sprintf('%s (N=%d)', cat, size(CSM_k,1));
end

% Neural reference lines
if do_neural_ref_lines
    xline(ax1, 0, '--k', 'LineWidth', 1);
    yline(ax1, 0, '--k', 'LineWidth', 1);
    yline(ax1, -2, ':k', 'LineWidth', 1);
    yline(ax1,  2, ':k', 'LineWidth', 1);
end

% Neural axis formatting
ax1.Color     = 'none';
ax1.XLim      = xlim_ms;
ax1.XColor    = 'k';
ax1.YColor    = 'k';
ax1.LineWidth = 1;
ax1.FontSize  = 11;
ax1.TickDir   = 'out';
xlabel(ax1, 'Time from threshold crossing (ms)');
ylabel(ax1, '\Delta Z (CSM - NoCSM)');

% Match positions exactly and link x
ax1.Position = ax2.Position;
linkaxes([ax1 ax2], 'x');

% ============================================================
% FORCE ZERO TO BE AT THE SAME HEIGHT ON BOTH AXES
% ============================================================
% Make each y-axis symmetric around zero so the 0-lines align.

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

% Re-apply exact matching position after limit updates
ax1.Position = ax2.Position;

% Legend
if ~isempty(legend_entries)
    legend(ax1, legend_entries, ...
        'Location','best', ...
        'Interpreter','none', ...
        'Box','off');
end

title(ax1, 'Δ Neural (left) + Δ Behavior (right)');

% ---------------- OPTIONAL EXPORT ----------------
% if ~exist(output_folder,'dir'); mkdir(output_folder); end
% outFile = fullfile(output_folder, 'CombinedDelta_NeuralBehavior_zeroAligned.pdf');
% print(fig, outFile, '-dpdf', '-painters');
%% With smoothing added

%% ============================================================
% COMBINED DELTA PLOT: NEURAL (left y) + BEHAVIOR (right y)
% WITH LIGHT DISPLAY-ONLY SMOOTHING
% FULL REWRITE WITH ZERO ALIGNED ACROSS BOTH Y-AXES
% Illustrator-editable (painters renderer)
% ============================================================

% -------------------------- SETTINGS --------------------------
xlim_ms = [-250 250];

% Behavior event lines originally in seconds
xlines_behavior_s  = [-0.055, 0.087, 0];
xlines_behavior_ms = xlines_behavior_s * 1000;

do_neural_ref_lines = true;
colY = 8;   % behavior column to use

% ---- light smoothing for DISPLAY ONLY ----
smooth_bins = 4;   % same style as your previous plots
smooth_plot = @(x) smoothdata(x, 'gaussian', smooth_bins);

% ------------------ PREP BEHAVIOR DELTA ----------------------
% Δ behavior = CSM - NoCSM
diff_avg_y = RR_big.Average(:,colY) - NRR_big.Average(:,colY);
diff_avg_y = diff_avg_y * 10;   % keep your scaling

SE_diff_y = sqrt( ...
    (RR_big.Std(:,colY).^2  ./ RR_big.N) + ...
    (NRR_big.Std(:,colY).^2 ./ NRR_big.N) );

time_s  = RR_big.Average(:,1);
time_ms = time_s(:)' * 1000;

diff_avg_y = diff_avg_y(:)';
SE_diff_y  = SE_diff_y(:)';

% Window behavior to x-lims
beh_keep  = (time_ms >= xlim_ms(1)) & (time_ms <= xlim_ms(2));
time_ms_b = time_ms(beh_keep);
diff_b    = diff_avg_y(beh_keep);
SE_b      = SE_diff_y(beh_keep);

% ---- smooth behavior for plotting only ----
diff_b_plot = smooth_plot(diff_b);
SE_b_plot   = smooth_plot(SE_b);

% -------------------------- FIGURE ---------------------------
fig = figure('Name','Δ Neural + Δ Behavior', ...
             'Color','w', ...
             'Position',[100,100,650,320]);
set(fig,'Renderer','painters');

% ============================================================
% AX2: BEHAVIOR (RIGHT Y-AXIS) FIRST
% ============================================================
ax2 = axes('Parent',fig);
hold(ax2,'on');
box(ax2,'off');

% Behavior SEM
fill(ax2, ...
    [time_ms_b fliplr(time_ms_b)], ...
    [diff_b_plot + SE_b_plot, fliplr(diff_b_plot - SE_b_plot)], ...
    'k', ...
    'FaceAlpha',0.18, ...
    'EdgeColor','none');

% Behavior mean
plot(ax2, time_ms_b, diff_b_plot, 'k', 'LineWidth', 2);

% Behavior reference lines
yline(ax2, 0, '--k', 'LineWidth', 1);
for xl = xlines_behavior_ms
    xline(ax2, xl, '--k', 'LineWidth', 1);
end

% Behavior axis formatting
ax2.YAxisLocation = 'right';
ax2.XLim          = xlim_ms;
ax2.XColor        = 'k';
ax2.YColor        = 'k';
ax2.LineWidth     = 1;
ax2.FontSize      = 11;
ax2.TickDir       = 'out';
ax2.XTickLabel    = [];
ax2.XLabel.String = '';
ylabel(ax2, 'Δ Behavior (CSM - NoCSM)');

% ============================================================
% AX1: NEURAL (LEFT Y-AXIS) SECOND
% ============================================================
ax1 = axes('Parent',fig);
hold(ax1,'on');
box(ax1,'off');

legend_handles = [];
legend_entries = {};

for ci = 1:numel(group_names)

    cat = group_names{ci};

    % Color by category
    if strcmp(cat, 'negLag_posCorr')
        dcol = DeltaColor_pos;
    else
        dcol = DeltaColor_neg;
    end

    CSM_k = AllGroupsKept.(cat).CSM;
    No_k  = AllGroupsKept.(cat).NoCSM;

    if isempty(CSM_k) || isempty(No_k)
        fprintf('Skipping %s (empty)\n', cat);
        continue;
    end

    % Neural delta
    [mD, semD] = delta_sem(CSM_k, No_k);

    % ---- smooth neural for plotting only ----
    mD_plot   = smooth_plot(mD);
    semD_plot = smooth_plot(semD);

    % SEM patch first
    fill(ax1, ...
        [time_vector fliplr(time_vector)], ...
        [mD_plot + semD_plot, fliplr(mD_plot - semD_plot)], ...
        dcol, ...
        'FaceAlpha',0.20, ...
        'EdgeColor','none');

    % Mean trace
    h = plot(ax1, time_vector, mD_plot, 'Color', dcol, 'LineWidth', 2);

    legend_handles(end+1) = h; %#ok<SAGROW>
    legend_entries{end+1} = sprintf('%s (N=%d)', cat, size(CSM_k,1)); %#ok<SAGROW>
end

% Neural reference lines
if do_neural_ref_lines
    xline(ax1, 0, '--k', 'LineWidth', 1);
    yline(ax1, 0, '--k', 'LineWidth', 1);
    yline(ax1, -2, ':k', 'LineWidth', 1);
    yline(ax1,  2, ':k', 'LineWidth', 1);
end

% Neural axis formatting
ax1.Color     = 'none';
ax1.XLim      = xlim_ms;
ax1.XColor    = 'k';
ax1.YColor    = 'k';
ax1.LineWidth = 1;
ax1.FontSize  = 11;
ax1.TickDir   = 'out';
xlabel(ax1, 'Time from threshold crossing (ms)');
ylabel(ax1, '\Delta Z (CSM - NoCSM)');

% Match positions exactly and link x
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

% Re-apply exact matching position after limit updates
ax1.Position = ax2.Position;

% Legend
if ~isempty(legend_handles)
    legend(ax1, legend_handles, legend_entries, ...
        'Location','best', ...
        'Interpreter','none', ...
        'Box','off');
end

title(ax1, 'Δ Neural (left) + Δ Behavior (right)');

% ---------------- OPTIONAL EXPORT ----------------
% if ~exist(output_folder,'dir'); mkdir(output_folder); end
% outFile = fullfile(output_folder, 'CombinedDelta_NeuralBehavior_zeroAligned.pdf');
% print(fig, outFile, '-dpdf', '-painters');
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
figure('Position', [100, 100, 600, 200]);  
shadedErrorBar(time, diff_avg_y, SE_diff_y, 'lineProps', {'k', 'LineWidth', 2});
xlabel('Time from threshold crossing (s)');
xlim([-0.25 0.25])
%xlim('auto')
ylabel('Y-Velocity Difference (CSM - NoCSM)');
title('CSM vs NoCSM Kinematic Difference (Y-Direction)');
yline(0, '--');
xline(-0.055)
xline(0.087)
xline(0,'--')
%% ============================================================
% DELTA DIFFERENCE: (posCorr Δ) - (negCorr Δ)
% where Δ = CSM - NoCSM
% ============================================================

catPos = 'negLag_posCorr';
catNeg = 'negLag_negCorr';

CSM_pos = AllGroupsKept.(catPos).CSM;
No_pos  = AllGroupsKept.(catPos).NoCSM;

CSM_neg = AllGroupsKept.(catNeg).CSM;
No_neg  = AllGroupsKept.(catNeg).NoCSM;

% --- checks ---
if isempty(CSM_pos) || isempty(No_pos) || isempty(CSM_neg) || isempty(No_neg)
    error('One of the matrices is empty — cannot compute delta difference.');
end

% --- compute deltas ---
[mD_pos, semD_pos] = delta_sem(CSM_pos, No_pos);
[mD_neg, semD_neg] = delta_sem(CSM_neg, No_neg);

% --- difference of deltas ---
deltaDiff = mD_pos - mD_neg;

% --- SEM propagation ---
sem_deltaDiff = sqrt(semD_pos.^2 + semD_neg.^2);

% ------------------ PLOT ------------------
figure('Name','Δ(posCorr) - Δ(negCorr)','Position',[120 120 650 320]); 
hold on;

plot(time_vector, deltaDiff, 'k', 'LineWidth', 2);

fill([time_vector fliplr(time_vector)], ...
     [deltaDiff+sem_deltaDiff fliplr(deltaDiff-sem_deltaDiff)], ...
     [0.2 0.2 0.2], 'FaceAlpha',0.25, 'EdgeColor','none');

xline(0,'--k');
yline(0,'--k');

xlabel('Time from threshold crossing (ms)');
ylabel('\Delta difference  (posCorr - negCorr)');
xlim([-250 250]);

title('Difference of Δ traces between categories');
box off;
%% ============================================================
% NET (posCorr - negCorr) FOR CSM vs NoCSM  (same plot)
% Uses KEPT units from AllGroupsKept.(cat).CSM / .NoCSM
% Requires: mean_sem() that returns [meanTrace, semTrace]
%% ============================================================

catPos = 'negLag_posCorr';
catNeg = 'negLag_negCorr';

% ---------- pull matrices ----------
CSM_pos = AllGroupsKept.(catPos).CSM;     % [units x time]
CSM_neg = AllGroupsKept.(catNeg).CSM;

No_pos  = AllGroupsKept.(catPos).NoCSM;
No_neg  = AllGroupsKept.(catNeg).NoCSM;

% ---------- basic checks ----------
if isempty(CSM_pos) || isempty(CSM_neg)
    error('CSM matrices empty for pos and/or neg categories.');
end
if isempty(No_pos) || isempty(No_neg)
    error('NoCSM matrices empty for pos and/or neg categories.');
end

% ---------- means + SEMs ----------
[mCSM_pos, semCSM_pos] = mean_sem(CSM_pos);
[mCSM_neg, semCSM_neg] = mean_sem(CSM_neg);

[mNo_pos,  semNo_pos ] = mean_sem(No_pos);
[mNo_neg,  semNo_neg ] = mean_sem(No_neg);

% ---------- NETS (pos - neg) ----------
net_CSM = mCSM_pos - mCSM_neg;
net_No  = mNo_pos  - mNo_neg;

% SEM propagation (assumes independence; good default for group compare)
sem_net_CSM = sqrt(semCSM_pos.^2 + semCSM_neg.^2);
sem_net_No  = sqrt(semNo_pos.^2  + semNo_neg.^2);

% ---------- plot ----------
figure('Name','NET pos-neg (CSM vs NoCSM)','Position',[120 120 650 320]); hold on;


% CSM net
plot(time_vector, net_CSM, 'Color', Color_CSM, 'LineWidth', 2);
fill([time_vector fliplr(time_vector)], ...
     [net_CSM+sem_net_CSM fliplr(net_CSM-sem_net_CSM)], ...
     Color_CSM, 'FaceAlpha', 0.20, 'EdgeColor','none');

% NoCSM net
plot(time_vector, net_No, 'Color', Color_NoCSM, 'LineWidth', 2);
fill([time_vector fliplr(time_vector)], ...
     [net_No+sem_net_No fliplr(net_No-sem_net_No)], ...
     Color_NoCSM, 'FaceAlpha', 0.15, 'EdgeColor','none');

xline(0,'--k');
yline(0,'--k');

xlabel('Time from threshold crossing (ms)');
ylabel('Net Z (posCorr - negCorr)');
xlim([-250 250]);
% ylim([-3 3]);  % uncomment if you want fixed limits

legend({sprintf('CSM net (pos-neg): Npos=%d, Nneg=%d', size(CSM_pos,1), size(CSM_neg,1)), ...
        sprintf('NoCSM net (pos-neg): Npos=%d, Nneg=%d', size(No_pos,1),  size(No_neg,1))}, ...
       'Location','best','Interpreter','none');

title('NET difference between categories (posCorr - negCorr)');
box off;
%% ============================================================
% ADD-ON PLOTTING: TUNED (CSM-modulated) KEPT units BY CATEGORY
% (uses TunedGroupsKept.(negLag_posCorr / negLag_negCorr))
% Also makes category DELTA plots for tuned kept.
% Drop this AFTER your main script finishes (after TunedGroupsKept exists).
%% ============================================================

Color_CSM   = rgb('DarkMagenta');
Color_NoCSM = rgb('Black');

DeltaColor_pos = rgb('SteelBlue');  % tuned negLag_posCorr
DeltaColor_neg = rgb('DarkGreen');  % tuned negLag_negCorr

% from main script
full_timevec = timevec_global;
time_idx     = find(full_timevec >= -250 & full_timevec <= 250);
time_vector  = full_timevec(time_idx);

mean_sem = @(mat) deal( ...
    nanmean(mat(:,time_idx),1), ...
    nanstd(mat(:,time_idx),0,1) ./ sqrt(max(1,size(mat,1))) );

delta_sem = @(csm,no) deal( ...
    nanmean(csm(:,time_idx) - no(:,time_idx),1), ...
    nanstd(csm(:,time_idx) - no(:,time_idx),0,1) ./ sqrt(max(1,size(csm,1))) );

group_names = {'negLag_posCorr','negLag_negCorr'};

for ci = 1:numel(group_names)
    cat = group_names{ci};

    if strcmp(cat,'negLag_posCorr')
        dcol = DeltaColor_pos;
    else
        dcol = DeltaColor_neg;
    end

    CSM_tk = TunedGroupsKept.(cat).CSM;
    No_tk  = TunedGroupsKept.(cat).NoCSM;

    if isempty(CSM_tk) || isempty(No_tk)
        fprintf('TunedGroupsKept.%s empty -> skipping.\n', cat);
        continue;
    end

    % -------- PSTH (tuned kept) --------
    [mCSM, semCSM] = mean_sem(CSM_tk);
    [mNo,  semNo ] = mean_sem(No_tk);

    figure('Name',['TUNED-KEPT PSTH ' cat],'Position',[100,100,600,300]); hold on;
    plot(time_vector, mCSM,'Color',Color_CSM,'LineWidth',2);
    fill([time_vector fliplr(time_vector)], ...
         [mCSM+semCSM fliplr(mCSM-semCSM)], ...
         Color_CSM,'FaceAlpha',0.2,'EdgeColor','none');

    plot(time_vector, mNo,'Color',Color_NoCSM,'LineWidth',2);
    fill([time_vector fliplr(time_vector)], ...
         [mNo+semNo fliplr(mNo-semNo)], ...
         Color_NoCSM,'FaceAlpha',0.15,'EdgeColor','none');

    xline(0,'--k');
    xlabel('Time from CSM T2 (ms)'); ylabel('Z-scored FR');
    xlim([-250 250]);
    title(sprintf('TUNED-KEPT: %s (N=%d; p_{LR}<%.3f)', cat, size(CSM_tk,1), alpha_perm), 'Interpreter','none');
    box off;

    % -------- DELTA (tuned kept) --------
    [mD, semD] = delta_sem(CSM_tk, No_tk);

    figure('Name',['TUNED-KEPT DELTA ' cat],'Position',[100,100,600,300]); hold on;
    plot(time_vector, mD, 'Color', dcol, 'LineWidth', 2);
    fill([time_vector fliplr(time_vector)], ...
         [mD+semD fliplr(mD-semD)], ...
         dcol,'FaceAlpha',0.2,'EdgeColor','none');

    xline(0,'--k'); yline(0,'--k');
    yline(-2,':'); yline(2,':');
    %xlabel('Time from CSM T2 (ms)');
    xlabel('Time from threshold crossing (ms)');
    ylabel('\Delta Z (CSM - NoCSM)');
    xlim([-250 250]);
    ylim([-8 8])
    %ylim('auto')
    title(sprintf('TUNED-KEPT Δ: %s (N=%d)', cat, size(CSM_tk,1)), 'Interpreter','none');
    box off;
end

%%
%% ============================================================
% OVERLAY POS + NEG CORR (TUNED KEPT)
% ============================================================

Color_CSM   = rgb('DarkMagenta');
Color_NoCSM = rgb('Black');

DeltaColor_pos = rgb('SteelBlue');  
DeltaColor_neg = rgb('DarkGreen');  

% ---- time window ----
full_timevec = timevec_global;
time_idx     = find(full_timevec >= -250 & full_timevec <= 250);
time_vector  = full_timevec(time_idx);

% ---- helpers ----
mean_sem = @(mat) deal( ...
    nanmean(mat(:,time_idx),1), ...
    nanstd(mat(:,time_idx),0,1) ./ sqrt(max(1,size(mat,1))) );

delta_sem = @(csm,no) deal( ...
    nanmean(csm(:,time_idx) - no(:,time_idx),1), ...
    nanstd(csm(:,time_idx) - no(:,time_idx),0,1) ./ sqrt(max(1,size(csm,1))) );

% ---- groups ----
groups = {'negLag_posCorr','negLag_negCorr'};
colors = {DeltaColor_pos, DeltaColor_neg};

% ============================================================
% -------------------- PSTH (OVERLAY) ------------------------
% ============================================================

figure('Name','TUNED-KEPT PSTH (POS vs NEG)','Position',[100,100,650,350]); 
hold on;

for ci = 1:numel(groups)

    cat = groups{ci};
    dcol = colors{ci};

    CSM_tk = TunedGroupsKept.(cat).CSM;
    No_tk  = TunedGroupsKept.(cat).NoCSM;

    if isempty(CSM_tk) || isempty(No_tk)
        fprintf('Skipping %s (empty)\n', cat);
        continue;
    end

    % ---- means ----
    [mCSM, semCSM] = mean_sem(CSM_tk);
    [mNo,  semNo ] = mean_sem(No_tk);

    % ---- plot CSM ----
    plot(time_vector, mCSM, 'Color', dcol, 'LineWidth', 2.5);
    fill([time_vector fliplr(time_vector)], ...
         [mCSM+semCSM fliplr(mCSM-semCSM)], ...
         dcol, 'FaceAlpha', 0.2, 'EdgeColor','none');

    % ---- plot NoCSM (lighter same color) ----
    plot(time_vector, mNo, '--', 'Color', dcol, 'LineWidth', 2);
end

xline(0,'--k');
xlabel('Time from CSM T2 (ms)');
ylabel('Z-scored FR');
xlim([-250 250]);
title('TUNED-KEPT PSTH (posCorr vs negCorr)');
box off;

legend({'posCorr CSM','posCorr SEM','posCorr NoCSM', ...
        'negCorr CSM','negCorr SEM','negCorr NoCSM'}, ...
        'Location','best');


% ============================================================
% -------------------- DELTA (OVERLAY) ------------------------
% ============================================================

figure('Name','TUNED-KEPT DELTA (POS vs NEG)','Position',[100,100,650,350]); 
hold on;

for ci = 1:numel(groups)

    cat = groups{ci};
    dcol = colors{ci};

    CSM_tk = TunedGroupsKept.(cat).CSM;
    No_tk  = TunedGroupsKept.(cat).NoCSM;

    if isempty(CSM_tk) || isempty(No_tk)
        continue;
    end

    % ---- delta ----
    [mD, semD] = delta_sem(CSM_tk, No_tk);

    plot(time_vector, mD, 'Color', dcol, 'LineWidth', 2.5);
    fill([time_vector fliplr(time_vector)], ...
         [mD+semD fliplr(mD-semD)], ...
         dcol, 'FaceAlpha', 0.2, 'EdgeColor','none');
end

xline(0,'--k'); 
yline(0,'--k');
yline(2)
yline(-2)

xlabel('Time from CSMT2 (ms)');
ylabel('\Delta Z (CSM - NoCSM)');
xlim([-250 250]);
%ylim([-8 8]);
title('TUNED-KEPT Δ (posCorr vs negCorr)');
box off;

legend({'posCorr','posCorr SEM','negCorr','negCorr SEM'}, ...
        'Location','best');
%% ============================================================
% FIG SET 3: CATEGORY SCATTER of delta_netFR (epoch)
%   - one for ALL units
%   - one for KEPT units
% Labels: D#:U#
% Sig in red: CSM_epoch_glm.sig_LR
%% ============================================================

plot_epoch_scatter_from_XCorrData(XCorrData, group_names, false); % ALL units
plot_epoch_scatter_from_XCorrData(XCorrData, group_names, true ); % KEPT only

%% ===================== LOCAL FUNCTION =====================
function plot_epoch_scatter_from_XCorrData(XCorrData, groups, useKeepOnly)

    jitter_amp = 0.25;
    dx_label   = 0.02;
    label_fs   = 8;

    Unit = struct('day',{},'cat',{},'unitIdx',{},'dEp',{},'sigE',{});

    for day = 1:numel(XCorrData)
        if ~isfield(XCorrData(day),'units') || isempty(XCorrData(day).units)
            continue;
        end

        for uu = 1:numel(XCorrData(day).units)
            u = XCorrData(day).units(uu);

            if useKeepOnly
                if ~isfield(u,'keep') || ~u.keep
                    continue;
                end
            end

            if ~isfield(u,'category') || isempty(u.category) || ~ismember(u.category, groups)
                continue;
            end

            if ~isfield(u,'CSM_epoch_glm') || isempty(u.CSM_epoch_glm)
                continue;
            end
            if isfield(u.CSM_epoch_glm,'skip') && u.CSM_epoch_glm.skip
                continue;
            end

            if ~isfield(u.CSM_epoch_glm,'delta_netFR') || ~isfinite(u.CSM_epoch_glm.delta_netFR)
                continue;
            end

            sigE = false;
            if isfield(u.CSM_epoch_glm,'sig_LR')
                sigE = logical(u.CSM_epoch_glm.sig_LR);
            end

            Unit(end+1).day    = day;
            Unit(end).cat      = u.category;
            Unit(end).unitIdx  = u.unitIdx;
            Unit(end).dEp      = u.CSM_epoch_glm.delta_netFR;
            Unit(end).sigE     = sigE;
        end
    end

    if isempty(Unit)
        warning('No units for epoch scatter (useKeepOnly=%d).', useKeepOnly);
        return;
    end

    nU  = numel(Unit);
    dEp = nan(nU,1);
    sigE= false(nU,1);
    x   = nan(nU,1);

    for i = 1:nU
        dEp(i)  = Unit(i).dEp;
        sigE(i) = logical(Unit(i).sigE);
        x(i)    = find(strcmp(Unit(i).cat, groups), 1, 'first');
    end

    valid = isfinite(x) & isfinite(dEp);
    Unit = Unit(valid);
    x    = x(valid);
    dEp  = dEp(valid);
    sigE = sigE(valid);

    rng(0,'twister');
    jitter = (rand(numel(x),1)-0.5) * jitter_amp;
    xj = x + jitter;

    idxSig = logical(sigE);
    idxNo  = ~idxSig;

    if useKeepOnly
        figName = 'EPOCH DELTA_NETFR SCATTER (KEPT ONLY)';
        ttl     = 'Epoch effect size by category (KEPT only); labels Day+UnitID; sig in red';
    else
        figName = 'EPOCH DELTA_NETFR SCATTER (ALL UNITS)';
        ttl     = 'Epoch effect size by category (ALL units); labels Day+UnitID; sig in red';
    end

    figure('Color','w','Name',figName,'Position',[100 100 1200 520]); hold on;

    scatter(xj(idxNo), dEp(idxNo), 28, 'filled', 'MarkerFaceAlpha', 0.55, 'MarkerEdgeColor', 'none');
    scatter(xj(idxSig), dEp(idxSig), 46, 'filled', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none');

    xlim([0.5 numel(groups)+0.5]);
    xticks(1:numel(groups));
    xticklabels(groups);
    xtickangle(20);
    yline(0,'k-','LineWidth',1);

    ylabel('\Delta net FR (Hz): epoch mean (CSM - NoCSM)');
    title(ttl,'Interpreter','none');
    box off; set(gca,'TickDir','out');

    % group means + counts
    for g = 1:numel(groups)
        m = (x == g);
        if any(m)
            mu = mean(dEp(m), 'omitnan');
            plot(g, mu, 'kd', 'MarkerFaceColor','k', 'MarkerSize',8);
            text(g, mu, sprintf('  %d/%d', sum(m & idxSig), sum(m)), ...
                'FontWeight','bold', 'VerticalAlignment','bottom');
        end
    end

    dy_label = 0.02 * range(dEp);
    if ~isfinite(dy_label) || dy_label == 0
        dy_label = 0.05;
    end

    for i = 1:numel(Unit)

        uid = Unit(i).unitIdx;
        dd  = Unit(i).day;

        if idxSig(i)
            txtColor = [0.6 0 0];
            fontW    = 'bold';
        else
            txtColor = [0.2 0.2 0.2];
            fontW    = 'normal';
        end

        lab = sprintf('D%d:U%d', dd, uid);

        text(xj(i) + dx_label, dEp(i) + dy_label, lab, ...
            'FontSize', label_fs, ...
            'Color', txtColor, ...
            'FontWeight', fontW, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'bottom');
    end

    legend({'not sig','sig','mean'},'Location','best');
end
