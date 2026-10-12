%% Master Kinematic analysis looking at per session effects 

%% ============================================================
% WITHIN-SESSION behavioral effects
% Separate each day into BL, AD, WO, WOE
% Align to threshold crossing
% Measure deceleration from 0 to 100 ms after threshold crossing
% ============================================================

clear; clc;

ConFon = struct;

fname = uigetdir;
direc = dir(fullfile(fname, '*.mat'));
%%
% ---------------- settings ----------------
filt_state = 1;
pct   = 0.8;
leng  = 5;
thresh = 0.1;   
align_ms = 200;  % amount of alignment window
win_ms   = [0 100];   % deceleration window after threshold crossing

for i = 1:length(direc)

    load(fullfile(fname, direc(i).name));   % loads anim_data

    fullname   = direc(i).name;
    mouse_name = regexprep(fullname, '.mat', '');

    ConFon.(mouse_name).anim_data = anim_data;

    dnum = fieldnames(anim_data);

    for d = 1:length(dnum)

        day_name = dnum{d};

        % ---------------- raw reaches for this day ----------------
        Reaches = anim_data.(day_name).Reaches;
        rn = fieldnames(Reaches);

        blk_array = nan(length(rn),1);
        for r = 1:length(rn)
            blk_array(r) = Reaches.(rn{r}).block;
        end

        % ---------------- split blocks ----------------
        bl_inx = find(blk_array == 1);
        ad_inx = find(blk_array == 2);
        %wo_inx = find(blk_array == 3);
        wo_inx = find(ismember(blk_array, [3 4]));

        BASE  = Extract(Reaches, bl_inx);
        ADAPT = Extract(Reaches, ad_inx);
        WASH  = Extract(Reaches, wo_inx);

       

        % =========================================================
        % BASELINE
        % BL = all control reaches from baseline
        % =========================================================
        [BASE_S, BASE_C] = StimFilter(BASE); 

        % use ALL control reaches from baseline
        BL_block = BASE_C;

        if filt_state == 1 && ~isempty(fieldnames(BL_block))
            [BL_block, ~] = rfilts_standard(BL_block, thresh, leng, pct);
        end

        if ~isempty(fieldnames(BL_block))
            BL_block = accel_toreaches2(BL_block);
        end

        % =========================================================
        % ADAPTATION
        % AD = all stim reaches from adaptation
        % =========================================================
        %[ADAPT_S, ADAPT_C] = StimFilter(ADAPT); 
        [ADAPT_S, ADAPT_C] = StimFilter(Reaches); 

        % use ALL stim reaches from adaptation
        AD_block = ADAPT_S;

        if filt_state == 1 && ~isempty(fieldnames(AD_block))
            [AD_block, ~] = rfilts_standard(AD_block, thresh, leng, pct);
        end

        if ~isempty(fieldnames(AD_block))
            AD_block = accel_toreaches2(AD_block);
        end

        % =========================================================
        % WASHOUT
        % WO  = all washout control reaches
        % WOE = first 5 washout reaches
        % =========================================================
        [WASH_S, WASH_C] = StimFilter(WASH);

        % all washout control reaches
        WO_block = WASH_C;

        if filt_state == 1 && ~isempty(fieldnames(WO_block))
            [WO_block, ~] = rfilts_standard(WO_block, thresh, leng, pct);
        end

        if ~isempty(fieldnames(WO_block))
            WO_block = accel_toreaches2(WO_block);
        end

        % first 5 washout reaches = WOE
        rn_wash = fieldnames(WO_block);
        if isempty(rn_wash)
            WOE_block = struct;
        else
            nWOE = min(5, numel(rn_wash));
            WOE_block = Extract(WO_block, 1:nWOE);
        end

        % =========================================================
        % SAVE RAW BLOCKS PER DAY
        % =========================================================

        ConFon.(mouse_name).Session.(day_name).BL   = BL_block;
        ConFon.(mouse_name).Session.(day_name).AD   = AD_block;
        ConFon.(mouse_name).Session.(day_name).WO   = WO_block;
        ConFon.(mouse_name).Session.(day_name).WOE  = WOE_block;
        
        %Alignment section 

        if ~isempty(fieldnames(BL_block))
            BL_aligned = AlignThresh(BL_block, thresh, 1, 250);
            BL_aligned=stimvel_accel(BL_aligned);
        else
            BL_aligned = struct;
        end

        if ~isempty(fieldnames(AD_block))
            AD_aligned = AlignStim(AD_block,250);
            AD_aligned=stimvel_accel(AD_aligned);
        else
            AD_aligned = struct;
        end

        if ~isempty(fieldnames(WO_block))
            WO_aligned = AlignThresh(WO_block, thresh, 1, 250);
            WO_aligned=stimvel_accel(WO_aligned);
        else
            WO_aligned = struct;
        end

        if ~isempty(fieldnames(WOE_block))
            WOE_aligned = AlignThresh(WOE_block, thresh, 1, 250);
            WOE_aligned=stimvel_accel(WOE_aligned);
        else
            WOE_aligned = struct;
        end

 

        ConFon.(mouse_name).Session.(day_name).BL_align  = BL_aligned;
        ConFon.(mouse_name).Session.(day_name).AD_align  = AD_aligned;
        ConFon.(mouse_name).Session.(day_name).WO_align  = WO_aligned;
        ConFon.(mouse_name).Session.(day_name).WOE_align = WOE_aligned;

        % =========================================================
        % THRESHOLD-ALIGNED DECELERATION: 0 to 100 ms
        % =========================================================
        ConFon.(mouse_name).Session.(day_name).BL_decel0_100  = ...
            aavg_decel_70(BL_block);

        ConFon.(mouse_name).Session.(day_name).AD_decel0_100  = ...
            aavg_decel_70(AD_block);

        ConFon.(mouse_name).Session.(day_name).WO_decel0_100  = ...
            aavg_decel_70(WO_block);

        ConFon.(mouse_name).Session.(day_name).WOE_decel0_100 = ...
            aavg_decel_70(WOE_block);

        % Optional summary means
        ConFon.(mouse_name).Session.(day_name).BL_decel0_100_mean  = ...
            mean(ConFon.(mouse_name).Session.(day_name).BL_decel0_100(:,2), 'omitnan');

        ConFon.(mouse_name).Session.(day_name).AD_decel0_100_mean  = ...
            mean(ConFon.(mouse_name).Session.(day_name).AD_decel0_100(:,2), 'omitnan');

        ConFon.(mouse_name).Session.(day_name).WO_decel0_100_mean  = ...
            mean(ConFon.(mouse_name).Session.(day_name).WO_decel0_100(:,2), 'omitnan');

        ConFon.(mouse_name).Session.(day_name).WOE_decel0_100_mean = ...
            mean(ConFon.(mouse_name).Session.(day_name).WOE_decel0_100(:,2), 'omitnan');

    end

    clearvars -except ConFon fname direc i filt_state pct leng thresh align_ms win_ms
end

clearvars -except ConFon

%% ============================================================
% PLOT DECELERATION FOR BL, AD, WOE, WO
% Uses session means already saved in:
% ConFon.(mouse).Session.(day).BL_decel0_100
% ConFon.(mouse).Session.(day).AD_decel0_100
% ConFon.(mouse).Session.(day).WOE_decel0_100
% ConFon.(mouse).Session.(day).WO_decel0_100
%
% decel_col = which column of the decel output to plot:
%   1 = col 10 all
%   2 = col 11 out 
%   3 = col 12 up 
%   4 = col 13 lat
% ============================================================

decel_col = 2;   % <-- change if needed

allRows = {};

mice = fieldnames(ConFon);

for m = 1:numel(mice)
    mouse = mice{m};

    if ~isfield(ConFon.(mouse),'Session')
        continue
    end

    days = fieldnames(ConFon.(mouse).Session);

    for d = 1:numel(days)
        day = days{d};
        S = ConFon.(mouse).Session.(day);

        % --- session means ---
        BLm  = get_decel_mean_safe(S, 'BL_decel0_100',  decel_col);
        ADm  = get_decel_mean_safe(S, 'AD_decel0_100',  decel_col);
        WOEm = get_decel_mean_safe(S, 'WOE_decel0_100', decel_col);
        WOm  = get_decel_mean_safe(S, 'WO_decel0_100',  decel_col);

        % --- n reaches ---
        nBL  = get_n_safe(S, 'BL_decel0_100');
        nAD  = get_n_safe(S, 'AD_decel0_100');
        nWOE = get_n_safe(S, 'WOE_decel0_100');
        nWO  = get_n_safe(S, 'WO_decel0_100');

        allRows(end+1,:) = {mouse, day, BLm, ADm, WOEm, WOm, nBL, nAD, nWOE, nWO}; %#ok<SAGROW>
    end
end

DecelTable = cell2table(allRows, ...
    'VariableNames', {'Mouse','Session','BL','AD','WOE','WO','nBL','nAD','nWOE','nWO'});

%% ============================================================
% CONNECTED SESSION PLOT
% Each line = one session
% ============================================================

X = [1 2 3 4];
Y = [DecelTable.BL, DecelTable.AD, DecelTable.WOE, DecelTable.WO];

figure('Color','w','Position',[100 100 520 360]); hold on;

for i = 1:size(Y,1)
    plot(X, Y(i,:), '-o', ...
        'Color', [0.6 0.6 0.6], ...
        'MarkerFaceColor', [0.6 0.6 0.6], ...
        'LineWidth', 1.0, ...
        'MarkerSize', 5);
end

mY = mean(Y,1,'omitnan');
sY = std(Y,0,1,'omitnan') ./ sqrt(sum(~isnan(Y),1));

errorbar(X, mY, sY, 'k-o', ...
    'LineWidth', 2.2, ...
    'MarkerFaceColor', 'k', ...
    'MarkerSize', 7);

xlim([0.7 4.3]);
xticks(X);
xticklabels({'BL','AD','WOE','WO'});
ylabel(sprintf('Mean decel (col %d)', decel_col));
title('Session mean deceleration (0-100 ms after threshold)');
box off;
set(gca,'FontSize',12,'LineWidth',1.2);

%% ============================================================
% BOXPLOT + SESSION POINTS
% ============================================================

vals   = [DecelTable.BL; DecelTable.AD; DecelTable.WOE; DecelTable.WO];
group  = [repmat({'BL'},height(DecelTable),1); ...
          repmat({'AD'},height(DecelTable),1); ...
          repmat({'WOE'},height(DecelTable),1); ...
          repmat({'WO'},height(DecelTable),1)];

figure('Color','w','Position',[100 100 500 360]); hold on;

boxplot(vals, group, 'Symbol','');

jit = 0.08;
xBL  = 1 + (rand(height(DecelTable),1)-0.5)*2*jit;
xAD  = 2 + (rand(height(DecelTable),1)-0.5)*2*jit;
xWOE = 3 + (rand(height(DecelTable),1)-0.5)*2*jit;
xWO  = 4 + (rand(height(DecelTable),1)-0.5)*2*jit;

scatter(xBL,  DecelTable.BL,  28, 'filled', 'MarkerFaceAlpha',0.7);
scatter(xAD,  DecelTable.AD,  28, 'filled', 'MarkerFaceAlpha',0.7);
scatter(xWOE, DecelTable.WOE, 28, 'filled', 'MarkerFaceAlpha',0.7);
scatter(xWO,  DecelTable.WO,  28, 'filled', 'MarkerFaceAlpha',0.7);

ylabel(sprintf('Mean decel (col %d)', decel_col));
title('Deceleration across sessions');
box off;
set(gca,'FontSize',12,'LineWidth',1.2);

%% ============================================================
% OPTIONAL: BAR PLOT OF GRAND MEAN +/- SEM
% ============================================================

figure('Color','w','Position',[100 100 450 340]); hold on;

grandMean = mean(Y,1,'omitnan');
grandSEM  = std(Y,0,1,'omitnan') ./ sqrt(sum(~isnan(Y),1));

bar(X, grandMean, 0.65, 'FaceColor',[0.75 0.75 0.75], 'EdgeColor','k');
errorbar(X, grandMean, grandSEM, 'k.', 'LineWidth',1.8, 'CapSize',0);

scatter(repmat(X(1),height(DecelTable),1), DecelTable.BL,  24, 'filled', 'MarkerFaceAlpha',0.65);
scatter(repmat(X(2),height(DecelTable),1), DecelTable.AD,  24, 'filled', 'MarkerFaceAlpha',0.65);
scatter(repmat(X(3),height(DecelTable),1), DecelTable.WOE, 24, 'filled', 'MarkerFaceAlpha',0.65);
scatter(repmat(X(4),height(DecelTable),1), DecelTable.WO,  24, 'filled', 'MarkerFaceAlpha',0.65);

xlim([0.4 4.6]);
xticks(X);
xticklabels({'BL','AD','WOE','WO'});
ylabel(sprintf('Mean decel (col %d)', decel_col));
title('Grand mean deceleration');
box off;
set(gca,'FontSize',12,'LineWidth',1.2);

%% ============================================================
% HELPERS
% ============================================================

function out = get_decel_mean_safe(S, fieldname, decel_col)
if isfield(S, fieldname) && ~isempty(S.(fieldname)) && size(S.(fieldname),2) >= decel_col
    out = mean(S.(fieldname)(:,decel_col), 'omitnan');
else
    out = NaN;
end
end

function n = get_n_safe(S, fieldname)
if isfield(S, fieldname) && ~isempty(S.(fieldname))
    n = size(S.(fieldname),1);
else
    n = 0;
end
end

%% Same analysis but looking at per anmial 

%% ============================================================
% PLOT DECELERATION FOR BL, AD, WOE, WO
% AVERAGED WITHIN ANIMAL ACROSS SESSIONS
% ============================================================

decel_col = 2;   % choose which decel output column to use

mice = fieldnames(ConFon);

AnimalRows = {};

for m = 1:numel(mice)
    mouse = mice{m};

    if ~isfield(ConFon.(mouse),'Session')
        continue
    end

    days = fieldnames(ConFon.(mouse).Session);

    BL_all  = [];
    AD_all  = [];
    WOE_all = [];
    WO_all  = [];

    for d = 1:numel(days)
        day = days{d};
        S = ConFon.(mouse).Session.(day);

        BLm  = get_decel_mean_safe(S, 'BL_decel0_100',  decel_col);
        ADm  = get_decel_mean_safe(S, 'AD_decel0_100',  decel_col);
        WOEm = get_decel_mean_safe(S, 'WOE_decel0_100', decel_col);
        WOm  = get_decel_mean_safe(S, 'WO_decel0_100',  decel_col);

        BL_all(end+1,1)  = BLm;   %#ok<SAGROW>
        AD_all(end+1,1)  = ADm;   %#ok<SAGROW>
        WOE_all(end+1,1) = WOEm;  %#ok<SAGROW>
        WO_all(end+1,1)  = WOm;   %#ok<SAGROW>
    end

    AnimalRows(end+1,:) = { ...
        mouse, ...
        mean(BL_all,  'omitnan'), ...
        mean(AD_all,  'omitnan'), ...
        mean(WOE_all, 'omitnan'), ...
        mean(WO_all,  'omitnan'), ...
        sum(~isnan(BL_all)), ...
        sum(~isnan(AD_all)), ...
        sum(~isnan(WOE_all)), ...
        sum(~isnan(WO_all))}; %#ok<SAGROW>
end

AnimalTable = cell2table(AnimalRows, ...
    'VariableNames', {'Mouse','BL','AD','WOE','WO','nBLsess','nADsess','nWOEsess','nWOsess'});

%% ============================================================
% CONNECTED ANIMAL PLOT
% Each line = one animal
% ============================================================

X = [1 2 3 4];
Y = [AnimalTable.BL, AnimalTable.AD, AnimalTable.WOE, AnimalTable.WO];

figure('Color','w','Position',[100 100 520 360]); hold on;

for i = 1:size(Y,1)
    plot(X, Y(i,:), '-o', ...
        'Color', [0.65 0.65 0.65], ...
        'MarkerFaceColor', [0.65 0.65 0.65], ...
        'LineWidth', 1.2, ...
        'MarkerSize', 6);
end

mY = mean(Y,1,'omitnan');
sY = std(Y,0,1,'omitnan') ./ sqrt(sum(~isnan(Y),1));

errorbar(X, mY, sY, 'k-o', ...
    'LineWidth', 2.5, ...
    'MarkerFaceColor', 'k', ...
    'MarkerSize', 7);

xlim([0.7 4.3]);
xticks(X);
xticklabels({'BL','AD','WOE','WO'});
ylabel(sprintf('Mean decel (col %d)', decel_col));
title('Animal-averaged deceleration');
box off;
set(gca,'FontSize',12,'LineWidth',1.2);

%% ============================================================
% BOXPLOT / DOT PLOT ACROSS ANIMALS
% ============================================================

vals  = [AnimalTable.BL; AnimalTable.AD; AnimalTable.WOE; AnimalTable.WO];
group = [repmat({'BL'}, height(AnimalTable),1); ...
         repmat({'AD'}, height(AnimalTable),1); ...
         repmat({'WOE'},height(AnimalTable),1); ...
         repmat({'WO'}, height(AnimalTable),1)];

figure('Color','w','Position',[100 100 500 360]); hold on;
boxplot(vals, group, 'Symbol','');

jit = 0.08;
scatter(1 + (rand(height(AnimalTable),1)-0.5)*2*jit, AnimalTable.BL,  32, 'filled', 'MarkerFaceAlpha',0.75);
scatter(2 + (rand(height(AnimalTable),1)-0.5)*2*jit, AnimalTable.AD,  32, 'filled', 'MarkerFaceAlpha',0.75);
scatter(3 + (rand(height(AnimalTable),1)-0.5)*2*jit, AnimalTable.WOE, 32, 'filled', 'MarkerFaceAlpha',0.75);
scatter(4 + (rand(height(AnimalTable),1)-0.5)*2*jit, AnimalTable.WO,  32, 'filled', 'MarkerFaceAlpha',0.75);

ylabel(sprintf('Mean decel (col %d)', decel_col));
title('Animal means');
box off;
set(gca,'FontSize',12,'LineWidth',1.2);

%% ============================================================
% BAR PLOT OF GRAND MEAN +/- SEM ACROSS ANIMALS
% ============================================================

figure('Color','w','Position',[100 100 450 340]); hold on;

grandMean = mean(Y,1,'omitnan');
grandSEM  = std(Y,0,1,'omitnan') ./ sqrt(sum(~isnan(Y),1));

bar(X, grandMean, 0.65, 'FaceColor',[0.75 0.75 0.75], 'EdgeColor','k');
errorbar(X, grandMean, grandSEM, 'k.', 'LineWidth',1.8, 'CapSize',0);

scatter(repmat(1,height(AnimalTable),1), AnimalTable.BL,  28, 'filled', 'MarkerFaceAlpha',0.7);
scatter(repmat(2,height(AnimalTable),1), AnimalTable.AD,  28, 'filled', 'MarkerFaceAlpha',0.7);
scatter(repmat(3,height(AnimalTable),1), AnimalTable.WOE, 28, 'filled', 'MarkerFaceAlpha',0.7);
scatter(repmat(4,height(AnimalTable),1), AnimalTable.WO,  28, 'filled', 'MarkerFaceAlpha',0.7);

xlim([0.4 4.6]);
xticks(X);
xticklabels({'BL','AD','WOE','WO'});
ylabel(sprintf('Mean decel (col %d)', decel_col));
title('Grand mean deceleration across animals');
box off;
set(gca,'FontSize',12,'LineWidth',1.2);

%% ============================================================
% THRESHOLD-ALIGNED DECELERATION AUC FROM SESSION AVERAGE TRACES
% Uses:
%   ConFon.(mouse).Session.(day).BL_align.Average
%   ConFon.(mouse).Session.(day).AD_align.Average
%   ConFon.(mouse).Session.(day).WOE_align.Average
%   ConFon.(mouse).Session.(day).WO_align.Average
%
% Computes AUC of NEGATIVE acceleration only from 0 to 100 ms
% ============================================================

acc_col = 8;          % <-- choose acceleration column from .Average
time_col = 1;          % assumes col 1 is time in seconds
win_s = [0 0.100];     % 0 to 100 ms after threshold crossing

mice = fieldnames(ConFon);

AnimalRows = {};

for m = 1:numel(mice)
    mouse = mice{m};

    if ~isfield(ConFon.(mouse), 'Session')
        continue
    end

    days = fieldnames(ConFon.(mouse).Session);

    BL_auc_all  = [];
    AD_auc_all  = [];
    WOE_auc_all = [];
    WO_auc_all  = [];

    for d = 1:numel(days)
        day = days{d};
        S = ConFon.(mouse).Session.(day);

        BL_auc  = get_auc_from_align_safe(S, 'BL_align',  acc_col, time_col, win_s);
        AD_auc  = get_auc_from_align_safe(S, 'AD_align',  acc_col, time_col, win_s);
        WOE_auc = get_auc_from_align_safe(S, 'WOE_align', acc_col, time_col, win_s);
        WO_auc  = get_auc_from_align_safe(S, 'WO_align',  acc_col, time_col, win_s);

        BL_auc_all(end+1,1)  = BL_auc;   %#ok<SAGROW>
        AD_auc_all(end+1,1)  = AD_auc;   %#ok<SAGROW>
        WOE_auc_all(end+1,1) = WOE_auc;  %#ok<SAGROW>
        WO_auc_all(end+1,1)  = WO_auc;   %#ok<SAGROW>
    end

    AnimalRows(end+1,:) = { ...
        mouse, ...
        mean(BL_auc_all,  'omitnan'), ...
        mean(AD_auc_all,  'omitnan'), ...
        mean(WOE_auc_all, 'omitnan'), ...
        mean(WO_auc_all,  'omitnan'), ...
        sum(~isnan(BL_auc_all)), ...
        sum(~isnan(AD_auc_all)), ...
        sum(~isnan(WOE_auc_all)), ...
        sum(~isnan(WO_auc_all))}; %#ok<SAGROW>
end

AnimalTable = cell2table(AnimalRows, ...
    'VariableNames', {'Mouse','BL','AD','WOE','WO', ...
                      'nBLsess','nADsess','nWOEsess','nWOsess'});


function auc_val = get_auc_from_align_safe(S, align_field, acc_col, time_col, win_s)

auc_val = NaN;

if ~isfield(S, align_field)
    return
end

A = S.(align_field);

if ~isstruct(A) || ~isfield(A, 'Average') || isempty(A.Average)
    return
end

M = A.Average;

if size(M,2) < max(acc_col, time_col)
    return
end

t = M(:,time_col);
a = M(:,acc_col);

% restrict to desired window
idx = t >= win_s(1) & t <= win_s(2);

if ~any(idx)
    return
end

t_win = t(idx);
a_win = a(idx);

% keep only negative acceleration (deceleration)
a_neg = a_win;
a_neg(a_neg > 0) = 0;

% AUC of deceleration magnitude
auc_val = trapz(t_win, abs(a_neg));
end
%% ============================================================
% CONNECTED PLOT: one line per animal
% ============================================================

X = [1 2 3 4];
Y = [AnimalTable.BL, AnimalTable.AD, AnimalTable.WOE, AnimalTable.WO];

figure('Color','w','Position',[100 100 520 360]); hold on;

for i = 1:size(Y,1)
    plot(X, Y(i,:), '-o', ...
        'Color', [0.7 0.7 0.7], ...
        'MarkerFaceColor', [0.7 0.7 0.7], ...
        'LineWidth', 1.2, ...
        'MarkerSize', 6);
end

mY = mean(Y,1,'omitnan');
sY = std(Y,0,1,'omitnan') ./ sqrt(sum(~isnan(Y),1));

errorbar(X, mY, sY, 'k-o', ...
    'LineWidth', 2.5, ...
    'MarkerFaceColor', 'k', ...
    'MarkerSize', 7);

xlim([0.7 4.3]);
xticks(X);
xticklabels({'BL','AD','WOE','WO'});
ylabel('Negative acceleration AUC (0-100 ms)');
title('Threshold-aligned deceleration AUC');
box off;
set(gca,'FontSize',12,'LineWidth',1.2);

%% ============================================================
% EXTRACT THRESHOLD-ALIGNED MEAN DECELERATION (0-70 ms)
% FROM SESSION AVERAGE TRACES, THEN PLOT ANIMAL AVERAGES
%
% Assumes data live here:
%   ConFon.(mouse).Session.(day).BL_align.Average
%   ConFon.(mouse).Session.(day).AD_align.Average
%   ConFon.(mouse).Session.(day).WOE_align.Average
%   ConFon.(mouse).Session.(day).WO_align.Average
%
% Assumes:
%   - rows 26:33 correspond to 0 to 70 ms
%   - acc_col is the column containing the deceleration signal
%
% Output:
%   AnimalTable = one row per animal
%   SessionTable = one row per session
%   figure with one line per animal + grand mean +/- SEM
% ============================================================



% ---------------- SETTINGS ----------------
acc_col   = 8;      % out accel
row_win   = 26:33;   % 0 to 70 ms
%row_stim=26:31;
row_stim=25:31;
%row_stim=26:36;
row_nostim=26:36;
neg_only  = false;    % true = keep only negative values
flip_sign = false;    % true = larger positive number = stronger deceleration

% ---------------- EXTRACT ----------------
mice = fieldnames(ConFon);

AnimalRows  = {};
SessionRows = {};

for m = 1:numel(mice)
    mouse = mice{m};

    if ~isfield(ConFon.(mouse), 'Session')
        continue
    end

    days = fieldnames(ConFon.(mouse).Session);

    BL_all  = [];
    AD_all  = [];
    WOE_all = [];
    WO_all  = [];

    for d = 1:numel(days)
        day = days{d};
        S   = ConFon.(mouse).Session.(day);
        % 
        % BL_val  = get_mean_from_align_idx(S, 'BL_align',  acc_col, row_win, neg_only, flip_sign);
        % AD_val  = get_mean_from_align_idx(S, 'AD_align',  acc_col, row_win, neg_only, flip_sign);
        % WOE_val = get_mean_from_align_idx(S, 'WOE_align', acc_col, row_win, neg_only, flip_sign);
        % WO_val  = get_mean_from_align_idx(S, 'WO_align',  acc_col, row_win, neg_only, flip_sign);

         BL_val  = get_mean_from_align_idx(S, 'BL_align',  acc_col, row_nostim, neg_only, flip_sign);
        AD_val  = get_mean_from_align_idx(S, 'AD_align',  acc_col, row_stim, neg_only, flip_sign);
        WOE_val = get_mean_from_align_idx(S, 'WOE_align', acc_col, row_nostim, neg_only, flip_sign);
        WO_val  = get_mean_from_align_idx(S, 'WO_align',  acc_col, row_nostim, neg_only, flip_sign);

        BL_all(end+1,1)  = BL_val;   %#ok<SAGROW>
        AD_all(end+1,1)  = AD_val;   %#ok<SAGROW>
        WOE_all(end+1,1) = WOE_val;  %#ok<SAGROW>
        WO_all(end+1,1)  = WO_val;   %#ok<SAGROW>

        SessionRows(end+1,:) = { ...
            mouse, day, BL_val, AD_val, WOE_val, WO_val}; %#ok<SAGROW>
    end

    AnimalRows(end+1,:) = { ...
        mouse, ...
        mean(BL_all,  'omitnan'), ...
        mean(AD_all,  'omitnan'), ...
        mean(WOE_all, 'omitnan'), ...
        mean(WO_all,  'omitnan'), ...
        sum(~isnan(BL_all)), ...
        sum(~isnan(AD_all)), ...
        sum(~isnan(WOE_all)), ...
        sum(~isnan(WO_all))}; %#ok<SAGROW>
end

AnimalTable = cell2table(AnimalRows, ...
    'VariableNames', {'Mouse','BL','AD','WOE','WO', ...
                      'nBLsess','nADsess','nWOEsess','nWOsess'});

SessionTable = cell2table(SessionRows, ...
    'VariableNames', {'Mouse','Session','BL','AD','WOE','WO'});

disp(AnimalTable);

%% ============================================================
% PLOT 1: CONNECTED ANIMAL MEANS
% ============================================================

X = [1 2 3 4];
Y = [AnimalTable.BL, AnimalTable.AD, AnimalTable.WOE, AnimalTable.WO];

fig1 = figure('Color','w','Position',[100 100 560 380]); 
hold on;

for i = 1:size(Y,1)
    plot(X, Y(i,:), '-o', ...
        'Color', [0.72 0.72 0.72], ...
        'MarkerFaceColor', [0.72 0.72 0.72], ...
        'MarkerEdgeColor', [0.72 0.72 0.72], ...
        'LineWidth', 1.2, ...
        'MarkerSize', 6);
end

mY = mean(Y, 1, 'omitnan');
sY = std(Y, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(Y),1));

errorbar(X, mY, sY, 'k-o', ...
    'LineWidth', 2.2, ...
    'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'k', ...
    'MarkerSize', 7, ...
    'CapSize', 0);

xticks(X);
xticklabels({'BL','AD','WOE','WO'});
xlim([0.7 4.3]);

if flip_sign
    ylabel('Mean deceleration (0-70 ms)');
else
    ylabel('Mean negative acceleration (0-70 ms)');
end

title('Animal-averaged threshold-aligned deceleration');
box off;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);

%% ============================================================
% PLOT 2: DOT PLOT + GRAND MEAN +/- SEM
% ============================================================

fig2 = figure('Color','w','Position',[100 100 500 360]);
hold on;

grandMean = mean(Y,1,'omitnan');
grandSEM  = std(Y,0,1,'omitnan') ./ sqrt(sum(~isnan(Y),1));

bar(X, grandMean, 0.68, ...
    'FaceColor', [0.85 0.85 0.85], ...
    'EdgeColor', 'k', ...
    'LineWidth', 1.2);

errorbar(X, grandMean, grandSEM, 'k.', ...
    'LineWidth', 1.8, ...
    'CapSize', 0);

jit = 0.08;
scatter(1 + (rand(height(AnimalTable),1)-0.5)*2*jit, AnimalTable.BL,  34, 'filled', ...
    'MarkerFaceColor',[0.25 0.25 0.25], 'MarkerFaceAlpha',0.8);
scatter(2 + (rand(height(AnimalTable),1)-0.5)*2*jit, AnimalTable.AD,  34, 'filled', ...
    'MarkerFaceColor',[0.25 0.25 0.25], 'MarkerFaceAlpha',0.8);
scatter(3 + (rand(height(AnimalTable),1)-0.5)*2*jit, AnimalTable.WOE, 34, 'filled', ...
    'MarkerFaceColor',[0.25 0.25 0.25], 'MarkerFaceAlpha',0.8);
scatter(4 + (rand(height(AnimalTable),1)-0.5)*2*jit, AnimalTable.WO,  34, 'filled', ...
    'MarkerFaceColor',[0.25 0.25 0.25], 'MarkerFaceAlpha',0.8);

xticks(X);
xticklabels({'BL','AD','WOE','WO'});
xlim([0.5 4.5]);

if flip_sign
    ylabel('Mean deceleration (0-70 ms)');
else
    ylabel('Mean negative acceleration (0-70 ms)');
end

title('Animal means');
box off;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);


%% ============================================================
% DELTA FROM BASELINE ANALYSIS
%
% Computes mean deceleration from threshold-aligned average traces
% using rows 26:33 (0 to 70 ms), then computes:
%   AD_delta  = AD  - BL
%   WOE_delta = WOE - BL
%   WO_delta  = WO  - BL
%
% Outputs:
%   SessionTable : one row per session
%   AnimalTable  : one row per animal (animal mean of session deltas)
%
% Plots:
%   1) connected animal plot of delta from BL
%   2) bar + dots plot of animal mean deltas
%
% Optional:
%   linear mixed effects model on session-level deltas
%
% Assumes data are in:
%   ConFon.(mouse).Session.(day).BL_align.Average
%   ConFon.(mouse).Session.(day).AD_align.Average
%   ConFon.(mouse).Session.(day).WOE_align.Average
%   ConFon.(mouse).Session.(day).WO_align.Average
% ============================================================


% ---------------- SETTINGS ----------------
acc_col   = 8;      % <-- CHANGE to your decel/accel column
row_win   = 26:31;   % rows corresponding to 0 to 70 ms
row_nostim=26:33;
neg_only  = true;    % keep only negative values
flip_sign = false;    % flip sign so larger positive = stronger deceleration

%% ============================================================
% EXTRACT SESSION-LEVEL VALUES AND DELTAS
% ============================================================

mice = fieldnames(ConFon);

SessionRows = {};

for m = 1:numel(mice)
    mouse = mice{m};

    if ~isfield(ConFon.(mouse), 'Session')
        continue
    end

    days = fieldnames(ConFon.(mouse).Session);

    for d = 1:numel(days)
        day = days{d};
        S   = ConFon.(mouse).Session.(day);

        BL_val  = get_mean_from_align_idx(S, 'BL_align',  acc_col, row_win, neg_only, flip_sign);
        AD_val  = get_mean_from_align_idx(S, 'AD_align',  acc_col, row_win, neg_only, flip_sign);
        WOE_val = get_mean_from_align_idx(S, 'WOE_align', acc_col, row_win, neg_only, flip_sign);
        WO_val  = get_mean_from_align_idx(S, 'WO_align',  acc_col, row_win, neg_only, flip_sign);

        % BL_val  = get_mean_from_align_idx(S, 'BL_align',  acc_col, row_nostim, neg_only, flip_sign);
        % AD_val  = get_mean_from_align_idx(S, 'AD_align',  acc_col, row_win, neg_only, flip_sign);
        % WOE_val = get_mean_from_align_idx(S, 'WOE_align', acc_col, row_nostim, neg_only, flip_sign);
        % WO_val  = get_mean_from_align_idx(S, 'WO_align',  acc_col, row_nostim, neg_only, flip_sign);

        AD_delta  = AD_val  - BL_val;
        WOE_delta = WOE_val - BL_val;
        WO_delta  = WO_val  - BL_val;
        
        %converting values 
        AD_delta=AD_delta*10; 
        WOE_delta=WOE_delta*10; 
        WO_delta=WO_delta*10; 

        SessionRows(end+1,:) = { ...
            mouse, day, d, ...
            BL_val, AD_val, WOE_val, WO_val, ...
            AD_delta, WOE_delta, WO_delta}; %#ok<SAGROW>
    end
end

SessionTable = cell2table(SessionRows, ...
    'VariableNames', {'Mouse','Session','SessionNum', ...
                      'BL','AD','WOE','WO', ...
                      'AD_delta','WOE_delta','WO_delta'});

disp(SessionTable)

%% ============================================================
% COMPUTE ANIMAL-LEVEL MEAN DELTAS
% ============================================================

AnimalRows = {};

for m = 1:numel(mice)
    mouse = mice{m};

    idx = strcmp(SessionTable.Mouse, mouse);
    if ~any(idx)
        continue
    end

    Tm = SessionTable(idx,:);

    AnimalRows(end+1,:) = { ...
        mouse, ...
        mean(Tm.AD_delta,  'omitnan'), ...
        mean(Tm.WOE_delta, 'omitnan'), ...
        mean(Tm.WO_delta,  'omitnan'), ...
        sum(~isnan(Tm.AD_delta)), ...
        sum(~isnan(Tm.WOE_delta)), ...
        sum(~isnan(Tm.WO_delta))}; %#ok<SAGROW>
end

AnimalTable = cell2table(AnimalRows, ...
    'VariableNames', {'Mouse','AD_delta','WOE_delta','WO_delta', ...
                      'nADsess','nWOEsess','nWOsess'});

disp(AnimalTable)

%% ============================================================
% PLOT 1: CONNECTED ANIMAL DELTA PLOT
% ============================================================

X = [1 2 3];
Y = [AnimalTable.AD_delta, AnimalTable.WOE_delta, AnimalTable.WO_delta];

figure('Color','w','Position',[100 100 500 360]); hold on;

for i = 1:size(Y,1)
    plot(X, Y(i,:), '-o', ...
        'Color', [0.72 0.72 0.72], ...
        'MarkerFaceColor', [0.72 0.72 0.72], ...
        'MarkerEdgeColor', [0.72 0.72 0.72], ...
        'LineWidth', 1.2, ...
        'MarkerSize', 6);
end

mY = mean(Y,1,'omitnan');
sY = std(Y,0,1,'omitnan') ./ sqrt(sum(~isnan(Y),1));

errorbar(X, mY, sY, 'k-o', ...
    'LineWidth', 2.2, ...
    'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'k', ...
    'MarkerSize', 7, ...
    'CapSize', 0);

yline(0,'k--','LineWidth',1);

xticks(X);
xticklabels({'AD-BL','WOE-BL','WO-BL'});
xlim([0.7 3.3]);
ylabel('\Delta deceleration from baseline');
title('Animal-averaged delta from baseline');
box off;
set(gca,'FontSize',12,'LineWidth',1.2);

%% ============================================================
% STATS: ARE DELTAS DIFFERENT FROM ZERO?
% ============================================================

% --- extract vectors ---
AD  = AnimalTable.AD_delta;
WOE = AnimalTable.WOE_delta;
WO  = AnimalTable.WO_delta;

% --- t-tests ---
[~,p_AD,~,stats_AD]   = ttest(AD);
[~,p_WOE,~,stats_WOE] = ttest(WOE);
[~,p_WO,~,stats_WO]   = ttest(WO);

fprintf('\n=== One-sample t-tests vs 0 ===\n');
fprintf('AD-BL:  p = %.4g, t(%d) = %.3f\n',  p_AD,  stats_AD.df,  stats_AD.tstat);
fprintf('WOE-BL: p = %.4g, t(%d) = %.3f\n', p_WOE, stats_WOE.df, stats_WOE.tstat);
fprintf('WO-BL:  p = %.4g, t(%d) = %.3f\n',  p_WO,  stats_WO.df,  stats_WO.tstat);

% --- optional: nonparametric ---
p_AD_np  = signrank(AD);
p_WOE_np = signrank(WOE);
p_WO_np  = signrank(WO);

fprintf('\n=== Wilcoxon signed-rank ===\n');
fprintf('AD-BL:  p = %.4g\n', p_AD_np);
fprintf('WOE-BL: p = %.4g\n', p_WOE_np);
fprintf('WO-BL:  p = %.4g\n', p_WO_np);

%% ============================================================
% PLOT 2: BAR + DOTS OF ANIMAL DELTAS
% ============================================================

figure('Color','w','Position',[100 100 460 340]); hold on;

grandMean = mean(Y,1,'omitnan');
grandSEM  = std(Y,0,1,'omitnan') ./ sqrt(sum(~isnan(Y),1));

bar(X, grandMean, 0.65, ...
    'FaceColor', [0.86 0.86 0.86], ...
    'EdgeColor', 'k', ...
    'LineWidth', 1.2);

errorbar(X, grandMean, grandSEM, 'k.', ...
    'LineWidth', 1.8, ...
    'CapSize', 0);

jit = 0.08;
scatter(1 + (rand(height(AnimalTable),1)-0.5)*2*jit, AnimalTable.AD_delta,  34, ...
    'filled', 'MarkerFaceColor',[0.2 0.2 0.2], 'MarkerFaceAlpha',0.8);
scatter(2 + (rand(height(AnimalTable),1)-0.5)*2*jit, AnimalTable.WOE_delta, 34, ...
    'filled', 'MarkerFaceColor',[0.2 0.2 0.2], 'MarkerFaceAlpha',0.8);
scatter(3 + (rand(height(AnimalTable),1)-0.5)*2*jit, AnimalTable.WO_delta,  34, ...
    'filled', 'MarkerFaceColor',[0.2 0.2 0.2], 'MarkerFaceAlpha',0.8);

yline(0,'k--','LineWidth',1);

xticks(X);
xticklabels({'AD-BL','WOE-BL','WO-BL'});
xlim([0.5 3.5]);
ylabel('\Delta deceleration from baseline');
title('Animal mean deltas');
box off;
set(gca,'FontSize',12,'LineWidth',1.2);

%% ============================================================
% OPTIONAL: SESSION-LEVEL DELTA PLOT ACROSS SESSIONS
% ============================================================

figure('Color','w','Position',[100 100 700 300]);

subplot(1,3,1); hold on;
plot_session_delta(SessionTable, 'AD_delta');
title('AD - BL');
xlabel('Session #'); ylabel('\Delta deceleration');
yline(0,'k--');

subplot(1,3,2); hold on;
plot_session_delta(SessionTable, 'WOE_delta');
title('WOE - BL');
xlabel('Session #'); ylabel('\Delta deceleration');
yline(0,'k--');

subplot(1,3,3); hold on;
plot_session_delta(SessionTable, 'WO_delta');
title('WO - BL');
xlabel('Session #'); ylabel('\Delta deceleration');
yline(0,'k--');

%% ============================================================
% OPTIONAL: LINEAR MIXED EFFECTS MODEL
% Session-level deltas from baseline
% ============================================================

DeltaRows = {};

for i = 1:height(SessionTable)
    DeltaRows(end+1,:) = {SessionTable.Mouse{i}, SessionTable.SessionNum(i), 'AD',  SessionTable.AD_delta(i)}; %#ok<SAGROW>
    DeltaRows(end+1,:) = {SessionTable.Mouse{i}, SessionTable.SessionNum(i), 'WOE', SessionTable.WOE_delta(i)}; %#ok<SAGROW>
    DeltaRows(end+1,:) = {SessionTable.Mouse{i}, SessionTable.SessionNum(i), 'WO',  SessionTable.WO_delta(i)}; %#ok<SAGROW>
end

DeltaTable = cell2table(DeltaRows, ...
    'VariableNames', {'Mouse','SessionNum','Type','Delta'});

DeltaTable.Mouse = categorical(DeltaTable.Mouse);
DeltaTable.Type  = categorical(DeltaTable.Type);

DeltaTable = DeltaTable(~isnan(DeltaTable.Delta), :);

% Full model: does delta differ by type and change across session?
lme_delta = fitlme(DeltaTable, 'Delta ~ Type * SessionNum + (1|Mouse)');
disp(lme_delta)
disp(anova(lme_delta))

% Individual tests by type
for typeName = {'AD','WOE','WO'}
    tt = DeltaTable(DeltaTable.Type == categorical(typeName), :);
    if height(tt) > 0
        fprintf('\n================ %s delta vs session ================\n', typeName{1});
        lme_one = fitlme(tt, 'Delta ~ SessionNum + (1|Mouse)');
        disp(lme_one)
        disp(anova(lme_one))
    end
end

% One-sample tests across animal means (simple summary)
fprintf('\n================ Animal mean delta vs zero ================\n');
[~,p_ad,~,stats_ad]   = ttest(AnimalTable.AD_delta);
[~,p_woe,~,stats_woe] = ttest(AnimalTable.WOE_delta);
[~,p_wo,~,stats_wo]   = ttest(AnimalTable.WO_delta);

fprintf('AD-BL:  p = %.4g, t(%d) = %.3f\n',  p_ad,  stats_ad.df,  stats_ad.tstat);
fprintf('WOE-BL: p = %.4g, t(%d) = %.3f\n', p_woe, stats_woe.df, stats_woe.tstat);
fprintf('WO-BL:  p = %.4g, t(%d) = %.3f\n',  p_wo,  stats_wo.df,  stats_wo.tstat);

%% ============================================================
% OPTIONAL: SAVE TABLES
% ============================================================

% writetable(SessionTable, 'SessionTable_deltaFromBL.csv');
% writetable(AnimalTable,  'AnimalTable_deltaFromBL.csv');
% writetable(DeltaTable,   'DeltaTable_LME.csv');

%% ============================================================
% HELPER FUNCTIONS
% ============================================================


function plot_session_delta(SessionTable, fieldName)

mice = unique(SessionTable.Mouse);
hold on

for i = 1:numel(mice)
    idx = strcmp(SessionTable.Mouse, mice{i});
    Tm  = SessionTable(idx,:);

    [sess_sorted, ord] = sort(Tm.SessionNum);
    vals = Tm.(fieldName)(ord);

    plot(sess_sorted, vals, '-o', ...
        'Color', [0.65 0.65 0.65], ...
        'MarkerFaceColor', [0.65 0.65 0.65], ...
        'LineWidth', 1.1, ...
        'MarkerSize', 5);
end

sess_all = unique(SessionTable.SessionNum);
m = nan(size(sess_all));
s = nan(size(sess_all));

for k = 1:numel(sess_all)
    idx = SessionTable.SessionNum == sess_all(k);
    vals = SessionTable.(fieldName)(idx);
    m(k) = mean(vals, 'omitnan');
    s(k) = std(vals, 0, 'omitnan') ./ sqrt(sum(~isnan(vals)));
end

errorbar(sess_all, m, s, 'k-o', ...
    'LineWidth', 2, ...
    'MarkerFaceColor', 'k', ...
    'MarkerSize', 6, ...
    'CapSize', 0);

box off
set(gca,'FontSize',11,'LineWidth',1.2)
end

%% ============================================================
% OPTIONAL: SAVE TABLES
% ============================================================

% writetable(AnimalTable, 'AnimalTable_deceleration_0_70ms.csv');
% writetable(SessionTable, 'SessionTable_deceleration_0_70ms.csv');

%% ============================================================
% HELPER FUNCTION
% ============================================================

function out = get_mean_from_align_idx(S, align_field, acc_col, row_win, neg_only, flip_sign)

out = NaN;

if ~isfield(S, align_field)
    return
end

A = S.(align_field);

if ~isstruct(A) || ~isfield(A, 'Average') || isempty(A.Average)
    return
end

M = A.Average;

if size(M,2) < acc_col || size(M,1) < max(row_win)
    return
end

a_win = M(row_win, acc_col);

if neg_only
    a_win(a_win > 0) = NaN;
end

out = mean(a_win, 'omitnan');

if flip_sign
    out = -out;
end
end

%%


%% ============================================================
% SESSION-LEVEL DELTA STATS FROM BASELINE
%
% Uses existing SessionTable with columns:
%   Mouse
%   Session
%   SessionNum
%   BL
%   AD
%   WOE
%   WO
%   AD_delta
%   WOE_delta
%   WO_delta
%
% Runs:
%   1) LME for AD-BL vs 0
%   2) LME for WOE-BL vs 0
%   3) LME for WO-BL vs 0
%   4) Combined LME comparing AD / WOE / WO
%   5) Combined LME testing change across sessions
%   6) One figure with session-level delta plots
%
% Notes:
%   - Random intercept for Mouse accounts for repeated sessions
%   - Intercept test in each one-condition model tests whether mean delta ~= 0
% ============================================================

%% ============================================================
% CHECK REQUIRED TABLE
% ============================================================

reqVars = {'Mouse','SessionNum','AD_delta','WOE_delta','WO_delta'};
for i = 1:numel(reqVars)
    if ~ismember(reqVars{i}, SessionTable.Properties.VariableNames)
        error('SessionTable is missing required variable: %s', reqVars{i});
    end
end

%% ============================================================
% BUILD LONG DELTA TABLE
% ============================================================

rows = {};

for i = 1:height(SessionTable)
    rows(end+1,:) = {SessionTable.Mouse{i}, SessionTable.SessionNum(i), 'AD',  SessionTable.AD_delta(i)};  %#ok<SAGROW>
    rows(end+1,:) = {SessionTable.Mouse{i}, SessionTable.SessionNum(i), 'WOE', SessionTable.WOE_delta(i)}; %#ok<SAGROW>
    rows(end+1,:) = {SessionTable.Mouse{i}, SessionTable.SessionNum(i), 'WO',  SessionTable.WO_delta(i)};  %#ok<SAGROW>
end

DeltaTable = cell2table(rows, ...
    'VariableNames', {'Mouse','SessionNum','Type','Delta'});

DeltaTable.Mouse = categorical(DeltaTable.Mouse);
DeltaTable.Type  = categorical(DeltaTable.Type);

DeltaTable = DeltaTable(~isnan(DeltaTable.Delta), :);

disp('================ DELTA TABLE ================')
disp(DeltaTable)

%% ============================================================
% INDIVIDUAL MODELS: IS EACH DELTA DIFFERENT FROM 0?
% Model:
%   Delta ~ 1 + (1|Mouse)
% The intercept tests whether the mean delta differs from zero
% ============================================================

% ---------------- AD - BL ----------------
tbl_AD = DeltaTable(DeltaTable.Type == 'AD', :);
lme_AD = fitlme(tbl_AD, 'Delta ~ 1 + (1|Mouse)');

disp(' ')
disp('================ AD - BL : LME vs 0 ================')
disp(lme_AD)
disp(anova(lme_AD))

coef_AD = lme_AD.Coefficients;
ad_row = strcmp(coef_AD.Name,'(Intercept)');
AD_est = coef_AD.Estimate(ad_row);
AD_se  = coef_AD.SE(ad_row);
AD_t   = coef_AD.tStat(ad_row);
AD_p   = coef_AD.pValue(ad_row);

fprintf('AD-BL intercept: estimate = %.6f, SE = %.6f, t = %.3f, p = %.4g\n', ...
    AD_est, AD_se, AD_t, AD_p);

% ---------------- WOE - BL ----------------
tbl_WOE = DeltaTable(DeltaTable.Type == 'WOE', :);
lme_WOE = fitlme(tbl_WOE, 'Delta ~ 1 + (1|Mouse)');

disp(' ')
disp('================ WOE - BL : LME vs 0 ================')
disp(lme_WOE)
disp(anova(lme_WOE))

coef_WOE = lme_WOE.Coefficients;
woe_row = strcmp(coef_WOE.Name,'(Intercept)');
WOE_est = coef_WOE.Estimate(woe_row);
WOE_se  = coef_WOE.SE(woe_row);
WOE_t   = coef_WOE.tStat(woe_row);
WOE_p   = coef_WOE.pValue(woe_row);

fprintf('WOE-BL intercept: estimate = %.6f, SE = %.6f, t = %.3f, p = %.4g\n', ...
    WOE_est, WOE_se, WOE_t, WOE_p);

% ---------------- WO - BL ----------------
tbl_WO = DeltaTable(DeltaTable.Type == 'WO', :);
lme_WO = fitlme(tbl_WO, 'Delta ~ 1 + (1|Mouse)');

disp(' ')
disp('================ WO - BL : LME vs 0 ================')
disp(lme_WO)
disp(anova(lme_WO))

coef_WO = lme_WO.Coefficients;
wo_row = strcmp(coef_WO.Name,'(Intercept)');
WO_est = coef_WO.Estimate(wo_row);
WO_se  = coef_WO.SE(wo_row);
WO_t   = coef_WO.tStat(wo_row);
WO_p   = coef_WO.pValue(wo_row);

fprintf('WO-BL intercept: estimate = %.6f, SE = %.6f, t = %.3f, p = %.4g\n', ...
    WO_est, WO_se, WO_t, WO_p);

%% ============================================================
% MULTIPLE COMPARISON CORRECTION ACROSS THE 3 PRIMARY TESTS
% (Benjamini-Hochberg FDR)
% ============================================================

p_primary = [AD_p; WOE_p; WO_p];
[p_fdr, crit_p, adj_p] = fdr_bh_local(p_primary, 0.05);

disp(' ')
disp('================ FDR-CORRECTED PRIMARY TESTS ================')
PrimaryStats = table( ...
    ["AD-BL"; "WOE-BL"; "WO-BL"], ...
    p_primary, ...
    adj_p, ...
    p_fdr, ...
    'VariableNames', {'Comparison','p_raw','p_FDR','sig_FDR'});
disp(PrimaryStats)

%% ============================================================
% COMBINED MODEL: DO DELTAS DIFFER BY TYPE?
% ============================================================

lme_type = fitlme(DeltaTable, 'Delta ~ Type + (1|Mouse)');

disp(' ')
disp('================ COMBINED MODEL: TYPE EFFECT ================')
disp(lme_type)
disp(anova(lme_type))

%% ============================================================
% COMBINED MODEL: DO DELTAS CHANGE ACROSS SESSIONS?
% ============================================================

lme_type_session = fitlme(DeltaTable, 'Delta ~ Type * SessionNum + (1|Mouse)');

disp(' ')
disp('================ COMBINED MODEL: TYPE * SESSION ================')
disp(lme_type_session)
disp(anova(lme_type_session))

%% ============================================================
% OPTIONAL SIMPLE PAIRED ANIMAL-LEVEL TESTS BETWEEN DELTA TYPES
% ============================================================

if exist('AnimalTable','var') && all(ismember({'AD_delta','WOE_delta','WO_delta'}, AnimalTable.Properties.VariableNames))
    [~,p_AD_vs_WOE,~,stats_AD_vs_WOE] = ttest(AnimalTable.AD_delta, AnimalTable.WOE_delta);
    [~,p_AD_vs_WO,~,stats_AD_vs_WO]   = ttest(AnimalTable.AD_delta, AnimalTable.WO_delta);
    [~,p_WOE_vs_WO,~,stats_WOE_vs_WO] = ttest(AnimalTable.WOE_delta, AnimalTable.WO_delta);

    disp(' ')
    disp('================ ANIMAL-LEVEL PAIRED TESTS BETWEEN DELTA TYPES ================')
    fprintf('AD vs WOE: p = %.4g, t(%d) = %.3f\n', p_AD_vs_WOE, stats_AD_vs_WOE.df, stats_AD_vs_WOE.tstat);
    fprintf('AD vs WO:  p = %.4g, t(%d) = %.3f\n', p_AD_vs_WO,  stats_AD_vs_WO.df,  stats_AD_vs_WO.tstat);
    fprintf('WOE vs WO: p = %.4g, t(%d) = %.3f\n', p_WOE_vs_WO, stats_WOE_vs_WO.df, stats_WOE_vs_WO.tstat);
end

%% ============================================================
% PLOT: SESSION-LEVEL DELTAS
% ============================================================

fig = figure('Color','w','Position',[100 100 900 300]);

subplot(1,3,1); hold on;
plot_session_delta_local(SessionTable, 'AD_delta');
title(sprintf('AD - BL\nLME p = %.3g', AD_p));
xlabel('Session #');
ylabel('\Delta from baseline');
yline(0,'k--','LineWidth',1);
box off;
set(gca,'FontSize',11,'LineWidth',1.2);

subplot(1,3,2); hold on;
plot_session_delta_local(SessionTable, 'WOE_delta');
title(sprintf('WOE - BL\nLME p = %.3g', WOE_p));
xlabel('Session #');
ylabel('\Delta from baseline');
yline(0,'k--','LineWidth',1);
box off;
set(gca,'FontSize',11,'LineWidth',1.2);

subplot(1,3,3); hold on;
plot_session_delta_local(SessionTable, 'WO_delta');
title(sprintf('WO - BL\nLME p = %.3g', WO_p));
xlabel('Session #');
ylabel('\Delta from baseline');
yline(0,'k--','LineWidth',1);
box off;
set(gca,'FontSize',11,'LineWidth',1.2);

%% ============================================================
% OPTIONAL SUMMARY TEXT FOR MANUSCRIPT / NOTES
% ============================================================

disp(' ')
disp('================ SUMMARY ================')
fprintf('AD-BL vs 0:  estimate = %.6f, p = %.4g\n', AD_est, AD_p);
fprintf('WOE-BL vs 0: estimate = %.6f, p = %.4g\n', WOE_est, WOE_p);
fprintf('WO-BL vs 0:  estimate = %.6f, p = %.4g\n', WO_est, WO_p);

%% ============================================================
% HELPER FUNCTIONS
% ============================================================

function plot_session_delta_local(SessionTable, fieldName)

mice = unique(SessionTable.Mouse);
hold on

for i = 1:numel(mice)
    idx = strcmp(SessionTable.Mouse, mice{i});
    Tm  = SessionTable(idx,:);

    [sess_sorted, ord] = sort(Tm.SessionNum);
    vals = Tm.(fieldName)(ord);

    plot(sess_sorted, vals, '-o', ...
        'Color', [0.75 0.75 0.75], ...
        'MarkerFaceColor', [0.75 0.75 0.75], ...
        'MarkerEdgeColor', [0.75 0.75 0.75], ...
        'LineWidth', 1.0, ...
        'MarkerSize', 5);
end

sess_all = unique(SessionTable.SessionNum);
m = nan(size(sess_all));
s = nan(size(sess_all));

for k = 1:numel(sess_all)
    idx = SessionTable.SessionNum == sess_all(k);
    vals = SessionTable.(fieldName)(idx);
    m(k) = mean(vals, 'omitnan');
    s(k) = std(vals, 0, 'omitnan') ./ sqrt(sum(~isnan(vals)));
end

errorbar(sess_all, m, s, 'k-o', ...
    'LineWidth', 2, ...
    'MarkerFaceColor', 'k', ...
    'MarkerEdgeColor', 'k', ...
    'MarkerSize', 6, ...
    'CapSize', 0);
end

function [h, crit_p, adj_p] = fdr_bh_local(pvals, q)
% Simple Benjamini-Hochberg FDR correction
% pvals: column vector
% q: desired FDR, e.g. 0.05

p = pvals(:);
m = numel(p);

[ps, sort_ids] = sort(p);
thresh = (1:m)'/m * q;

w = find(ps <= thresh);
h = false(size(p));

if isempty(w)
    crit_p = 0;
else
    crit_p = ps(max(w));
    h = p <= crit_p;
end

adj_p = nan(size(p));
for i = 1:m
    adj_p(sort_ids(i)) = min(1, min(m * ps(i:end) ./ (i:m)'));
end
end

%% 
% Now looking at just 4 sessions across each mouse to make it fair 

%% ============================================================
% SESSION-LEVEL DELTA STATS FROM BASELINE
% ONLY FIRST 4 SESSIONS PER MOUSE
% ========================================================================================================================
% CHECK REQUIRED TABLE
% ============================================================

reqVars = {'Mouse','SessionNum','AD_delta','WOE_delta','WO_delta'};
for i = 1:numel(reqVars)
    if ~ismember(reqVars{i}, SessionTable.Properties.VariableNames)
        error('SessionTable is missing required variable: %s', reqVars{i});
    end
end

%% ============================================================
% KEEP ONLY FIRST 4 SESSIONS PER MOUSE
% ============================================================

SessionTable_4 = table();

allMice = unique(SessionTable.Mouse);

for i = 1:numel(allMice)
    thisMouse = allMice{i};
    idx_mouse = strcmp(SessionTable.Mouse, thisMouse);
    Tm = SessionTable(idx_mouse,:);

    [~, ord] = sort(Tm.SessionNum);
    Tm = Tm(ord,:);

    nKeep = min(4, height(Tm));
    SessionTable_4 = [SessionTable_4; Tm(1:nKeep,:)]; %#ok<AGROW>
end

SessionTable_4 = sortrows(SessionTable_4, {'Mouse','SessionNum'});

disp('================ SESSION TABLE (FIRST 4 SESSIONS / MOUSE) ================')
disp(SessionTable_4)

%% ============================================================
% BUILD LONG DELTA TABLE
% ============================================================

rows = {};

for i = 1:height(SessionTable_4)
    rows(end+1,:) = {SessionTable_4.Mouse{i}, SessionTable_4.SessionNum(i), 'AD',  SessionTable_4.AD_delta(i)};  %#ok<SAGROW>
    rows(end+1,:) = {SessionTable_4.Mouse{i}, SessionTable_4.SessionNum(i), 'WOE', SessionTable_4.WOE_delta(i)}; %#ok<SAGROW>
    rows(end+1,:) = {SessionTable_4.Mouse{i}, SessionTable_4.SessionNum(i), 'WO',  SessionTable_4.WO_delta(i)};  %#ok<SAGROW>
end

DeltaTable = cell2table(rows, ...
    'VariableNames', {'Mouse','SessionNum','Type','Delta'});

DeltaTable.Mouse = categorical(DeltaTable.Mouse);
DeltaTable.Type  = categorical(DeltaTable.Type);
DeltaTable = DeltaTable(~isnan(DeltaTable.Delta), :);

disp('================ DELTA TABLE ================')
disp(DeltaTable)

%% ============================================================
% INDIVIDUAL MODELS: IS EACH DELTA DIFFERENT FROM 0?
% ============================================================

% ---------------- AD - BL ----------------
tbl_AD = DeltaTable(DeltaTable.Type == 'AD', :);
lme_AD = fitlme(tbl_AD, 'Delta ~ 1 + (1|Mouse)');

disp(' ')
disp('================ AD - BL : LME vs 0 ================')
disp(lme_AD)
disp(anova(lme_AD))

coef_AD = lme_AD.Coefficients;
ad_row = strcmp(coef_AD.Name,'(Intercept)');
AD_est = coef_AD.Estimate(ad_row);
AD_se  = coef_AD.SE(ad_row);
AD_t   = coef_AD.tStat(ad_row);
AD_p   = coef_AD.pValue(ad_row);

fprintf('AD-BL intercept: estimate = %.6f, SE = %.6f, t = %.3f, p = %.4g\n', ...
    AD_est, AD_se, AD_t, AD_p);

% ---------------- WOE - BL ----------------
tbl_WOE = DeltaTable(DeltaTable.Type == 'WOE', :);
lme_WOE = fitlme(tbl_WOE, 'Delta ~ 1 + (1|Mouse)');

disp(' ')
disp('================ WOE - BL : LME vs 0 ================')
disp(lme_WOE)
disp(anova(lme_WOE))

coef_WOE = lme_WOE.Coefficients;
woe_row = strcmp(coef_WOE.Name,'(Intercept)');
WOE_est = coef_WOE.Estimate(woe_row);
WOE_se  = coef_WOE.SE(woe_row);
WOE_t   = coef_WOE.tStat(woe_row);
WOE_p   = coef_WOE.pValue(woe_row);

fprintf('WOE-BL intercept: estimate = %.6f, SE = %.6f, t = %.3f, p = %.4g\n', ...
    WOE_est, WOE_se, WOE_t, WOE_p);

% ---------------- WO - BL ----------------
tbl_WO = DeltaTable(DeltaTable.Type == 'WO', :);
lme_WO = fitlme(tbl_WO, 'Delta ~ 1 + (1|Mouse)');

disp(' ')
disp('================ WO - BL : LME vs 0 ================')
disp(lme_WO)
disp(anova(lme_WO))

coef_WO = lme_WO.Coefficients;
wo_row = strcmp(coef_WO.Name,'(Intercept)');
WO_est = coef_WO.Estimate(wo_row);
WO_se  = coef_WO.SE(wo_row);
WO_t   = coef_WO.tStat(wo_row);
WO_p   = coef_WO.pValue(wo_row);

fprintf('WO-BL intercept: estimate = %.6f, SE = %.6f, t = %.3f, p = %.4g\n', ...
    WO_est, WO_se, WO_t, WO_p);

%% ============================================================
% FDR CORRECTION
% ============================================================

p_primary = [AD_p; WOE_p; WO_p];
[p_fdr, crit_p, adj_p] = fdr_bh_local(p_primary, 0.05);

disp(' ')
disp('================ FDR-CORRECTED PRIMARY TESTS ================')
PrimaryStats = table( ...
    ["AD-BL"; "WOE-BL"; "WO-BL"], ...
    p_primary, ...
    adj_p, ...
    p_fdr, ...
    'VariableNames', {'Comparison','p_raw','p_FDR','sig_FDR'});
disp(PrimaryStats)

%% ============================================================
% COMBINED MODEL: TYPE EFFECT
% ============================================================

lme_type = fitlme(DeltaTable, 'Delta ~ Type + (1|Mouse)');

disp(' ')
disp('================ COMBINED MODEL: TYPE EFFECT ================')
disp(lme_type)
disp(anova(lme_type))

%% ============================================================
% COMBINED MODEL: TYPE * SESSION
% ============================================================

lme_type_session = fitlme(DeltaTable, 'Delta ~ Type * SessionNum + (1|Mouse)');

disp(' ')
disp('================ COMBINED MODEL: TYPE * SESSION ================')
disp(lme_type_session)
disp(anova(lme_type_session))

%% ============================================================
% ANIMAL TABLE RECOMPUTED USING ONLY FIRST 4 SESSIONS
% ============================================================

AnimalRows_4 = {};

mice4 = unique(SessionTable_4.Mouse);

for i = 1:numel(mice4)
    thisMouse = mice4{i};
    idx = strcmp(SessionTable_4.Mouse, thisMouse);
    Tm = SessionTable_4(idx,:);

    AnimalRows_4(end+1,:) = { ...
        thisMouse, ...
        mean(Tm.AD_delta,  'omitnan'), ...
        mean(Tm.WOE_delta, 'omitnan'), ...
        mean(Tm.WO_delta,  'omitnan')}; %#ok<SAGROW>
end

AnimalTable_4 = cell2table(AnimalRows_4, ...
    'VariableNames', {'Mouse','AD_delta','WOE_delta','WO_delta'});

disp('================ ANIMAL TABLE (FIRST 4 SESSIONS / MOUSE) ================')
disp(AnimalTable_4)

%% ============================================================
% OPTIONAL PAIRED TESTS ON ANIMAL MEANS
% ============================================================

[~,p_AD_vs_WOE,~,stats_AD_vs_WOE] = ttest(AnimalTable_4.AD_delta, AnimalTable_4.WOE_delta);
[~,p_AD_vs_WO,~,stats_AD_vs_WO]   = ttest(AnimalTable_4.AD_delta, AnimalTable_4.WO_delta);
[~,p_WOE_vs_WO,~,stats_WOE_vs_WO] = ttest(AnimalTable_4.WOE_delta, AnimalTable_4.WO_delta);

disp(' ')
disp('================ ANIMAL-LEVEL PAIRED TESTS BETWEEN DELTA TYPES ================')
fprintf('AD vs WOE: p = %.4g, t(%d) = %.3f\n', p_AD_vs_WOE, stats_AD_vs_WOE.df, stats_AD_vs_WOE.tstat);
fprintf('AD vs WO:  p = %.4g, t(%d) = %.3f\n', p_AD_vs_WO,  stats_AD_vs_WO.df,  stats_AD_vs_WO.tstat);
fprintf('WOE vs WO: p = %.4g, t(%d) = %.3f\n', p_WOE_vs_WO, stats_WOE_vs_WO.df, stats_WOE_vs_WO.tstat);

%% ============================================================
% PLOT: SESSION-LEVEL DELTAS
% ============================================================

figure('Color','w','Position',[100 100 900 300]);

subplot(1,3,1); hold on;
plot_session_delta_local(SessionTable_4, 'AD_delta');
title(sprintf('AD - BL\nLME p = %.3g', AD_p));
xlabel('Session #');
ylabel('\Delta from baseline');
yline(0,'k--','LineWidth',1);
box off;
set(gca,'FontSize',11,'LineWidth',1.2);

subplot(1,3,2); hold on;
plot_session_delta_local(SessionTable_4, 'WOE_delta');
title(sprintf('WOE - BL\nLME p = %.3g', WOE_p));
xlabel('Session #');
ylabel('\Delta from baseline');
yline(0,'k--','LineWidth',1);
box off;
set(gca,'FontSize',11,'LineWidth',1.2);

subplot(1,3,3); hold on;
plot_session_delta_local(SessionTable_4, 'WO_delta');
title(sprintf('WO - BL\nLME p = %.3g', WO_p));
xlabel('Session #');
ylabel('\Delta from baseline');
yline(0,'k--','LineWidth',1);
box off;
set(gca,'FontSize',11,'LineWidth',1.2);


%% ============================================================
% SUMMARY
% ============================================================

disp(' ')
disp('================ SUMMARY ================')
fprintf('AD-BL vs 0:  estimate = %.6f, p = %.4g\n', AD_est, AD_p);
fprintf('WOE-BL vs 0: estimate = %.6f, p = %.4g\n', WOE_est, WOE_p);
fprintf('WO-BL vs 0:  estimate = %.6f, p = %.4g\n', WO_est, WO_p);

%%
% ============================================================
% BOX PLOT: AD vs WOE (ANIMAL-LEVEL, FIRST 4 SESSIONS)
% WITH PAIRED LINES + p-VALUE
% ============================================================

AD  = AnimalTable_4.AD_delta;
WOE = AnimalTable_4.WOE_delta;

% paired t-test
[~,p_AD_vs_WOE,~,stats] = ttest(AD, WOE);

fprintf('AD vs WOE: p = %.4g, t(%d) = %.3f\n', ...
    p_AD_vs_WOE, stats.df, stats.tstat);

% ---------------- plot ----------------
figure('Color','w','Position',[100 100 320 380]); hold on;

x1 = 1;
x2 = 2;

% --- paired lines (background) ---
for i = 1:numel(AD)
    plot([x1 x2], [AD(i) WOE(i)], '-', ...
        'Color', [0.75 0.75 0.75], ...
        'LineWidth', 1.2);
end

% --- boxplot ---
boxplot([AD WOE], ...
    'Positions',[x1 x2], ...
    'Widths',0.5, ...
    'Colors','k', ...
    'Symbol','');

% --- overlay points ---
scatter(repmat(x1,size(AD)), AD, 55, 'filled', ...
    'MarkerFaceColor',[0.3 0.3 0.3], ...
    'MarkerEdgeColor','k');

scatter(repmat(x2,size(WOE)), WOE, 55, 'filled', ...
    'MarkerFaceColor',[0.3 0.3 0.3], ...
    'MarkerEdgeColor','k');

% --- mean markers (optional but nice) ---
plot(x1, mean(AD,'omitnan'), 'ko', 'MarkerFaceColor','b', 'MarkerSize',7);
plot(x2, mean(WOE,'omitnan'), 'ko', 'MarkerFaceColor','m', 'MarkerSize',7);

% --- axes ---
xlim([0.5 2.5])
xticks([x1 x2])
xticklabels({'AD-BL','WOE-BL'})
ylabel('\Delta deceleration from baseline')
yline(0,'k--')

box off
set(gca,'FontSize',12,'LineWidth',1.2)

% ---------------- p-value annotation ----------------
yl = ylim;
yr = yl(2) - yl(1);

y_bracket = yl(2) + 0.05*yr;
h = 0.04*yr;

ylim([yl(1), yl(2) + 0.2*yr]);

plot([x1 x1 x2 x2], [y_bracket y_bracket+h y_bracket+h y_bracket], ...
    'k-', 'LineWidth', 1.5);

text(mean([x1 x2]), y_bracket + h + 0.01*yr, ...
    sprintf('p = %.3f', p_AD_vs_WOE), ...
    'HorizontalAlignment','center', ...
    'FontSize',11);

%% ============================================================
% STATS: AD and WOE vs BASELINE (0)
% (ANIMAL-LEVEL, FIRST 4 SESSIONS)
% ============================================================

AD  = AnimalTable_4.AD_delta;
WOE = AnimalTable_4.WOE_delta;

% --- one-sample t-tests vs 0 ---
[~,p_AD,~,stats_AD]   = ttest(AD);
[~,p_WOE,~,stats_WOE] = ttest(WOE);

fprintf('\n=== VS BASELINE (0) ===\n')
fprintf('AD-BL:  p = %.4g, t(%d) = %.3f\n',  p_AD,  stats_AD.df,  stats_AD.tstat);
fprintf('WOE-BL: p = %.4g, t(%d) = %.3f\n', p_WOE, stats_WOE.df, stats_WOE.tstat);

% ============================================================
% BOX PLOT WITH p-VALUES VS 0
% ============================================================

figure('Color','w','Position',[100 100 320 380]); hold on;

x1 = 1;
x2 = 2;

% paired lines (optional but nice)
for i = 1:numel(AD)
    plot([x1 x2], [AD(i) WOE(i)], '-', ...
        'Color', [0.75 0.75 0.75], ...
        'LineWidth', 1.2);
end

% boxplot
boxplot([AD WOE], ...
    'Positions',[x1 x2], ...
    'Widths',0.5, ...
    'Colors','k', ...
    'Symbol','');

% points
scatter(repmat(x1,size(AD)), AD, 55, 'filled', ...
    'MarkerFaceColor',[0.3 0.3 0.3], 'MarkerEdgeColor','k');

scatter(repmat(x2,size(WOE)), WOE, 55, 'filled', ...
    'MarkerFaceColor',[0.3 0.3 0.3], 'MarkerEdgeColor','k');

% means
plot(x1, mean(AD,'omitnan'), 'ko', 'MarkerFaceColor','k', 'MarkerSize',7);
plot(x2, mean(WOE,'omitnan'), 'ko', 'MarkerFaceColor','k', 'MarkerSize',7);

% axis
xlim([0.5 2.5])
xticks([x1 x2])
xticklabels({'AD-BL','WOE-BL'})
ylabel('\Delta deceleration from baseline')
yline(0,'k--','LineWidth',1)

box off
set(gca,'FontSize',12,'LineWidth',1.2)

% ============================================================
% ADD p-values vs 0 (above each box)
% ============================================================

yl = ylim;
yr = yl(2) - yl(1);

y_text = yl(2) + 0.05*yr;
ylim([yl(1), yl(2) + 0.2*yr]);

text(x1, y_text, sprintf('p = %.3f', p_AD), ...
    'HorizontalAlignment','center', ...
    'FontSize',11);

text(x2, y_text, sprintf('p = %.3f', p_WOE), ...
    'HorizontalAlignment','center', ...
    'FontSize',11);

%% TRying out a hierarchial mixed effects model: 

%% ============================================================
% HIERARCHICAL MIXED-EFFECTS MODEL
% WASHOUT vs BASELINE USING PER-SESSION DATA
%
% Uses session-level means already stored in SessionTable:
%   Mouse, Session, SessionNum, BL, AD, WOE, WO
%
% Model structure:
%   Decel ~ Condition + SessionNum + (1|Mouse) + (1|Mouse:SessionID)
%
% This tests:
%   - whether WOE differs from BL
%   - whether WO differs from BL
%   - optional change across session number
%
% BL is the reference level
% ============================================================

%% ---------------- CHECK REQUIRED TABLE ----------------
reqVars = {'Mouse','Session','SessionNum','BL','WOE','WO'};
for i = 1:numel(reqVars)
    if ~ismember(reqVars{i}, SessionTable.Properties.VariableNames)
        error('SessionTable is missing required variable: %s', reqVars{i});
    end
end

% ---------------- BUILD LONG TABLE ----------------
LongRows = {};

for i = 1:height(SessionTable)

    mouse_i   = SessionTable.Mouse{i};
    sess_i    = SessionTable.Session{i};
    sessnum_i = SessionTable.SessionNum(i);

    % BL
    LongRows(end+1,:) = {mouse_i, sess_i, sessnum_i, 'BL',  SessionTable.BL(i)};   %#ok<SAGROW>

    % WOE
    LongRows(end+1,:) = {mouse_i, sess_i, sessnum_i, 'WOE', SessionTable.WOE(i)};  %#ok<SAGROW>

    % WO
    LongRows(end+1,:) = {mouse_i, sess_i, sessnum_i, 'WO',  SessionTable.WO(i)};   %#ok<SAGROW>
end

LongTable = cell2table(LongRows, ...
    'VariableNames', {'Mouse','Session','SessionNum','Condition','Decel'});

% remove missing
LongTable = LongTable(~isnan(LongTable.Decel), :);

% categorical variables
LongTable.Mouse     = categorical(LongTable.Mouse);
LongTable.Session   = categorical(LongTable.Session);
LongTable.Condition = categorical(LongTable.Condition);

% make BL the reference
LongTable.Condition = reordercats(LongTable.Condition, {'BL','WOE','WO'});

% unique session ID nested within mouse
LongTable.SessionID = categorical(strcat(string(LongTable.Mouse), "__", string(LongTable.Session)));

disp('================ LONG TABLE FOR HIERARCHICAL MODEL ================')
disp(LongTable)

%% ============================================================
% MODEL 1: MAIN EFFECT OF CONDITION
% Tests whether WOE and WO differ from BL
% ============================================================

lme_wash = fitlme(LongTable, ...
    'Decel ~ Condition + (1|Mouse) + (1|SessionID)');

disp(' ')
disp('================ HIERARCHICAL MODEL: CONDITION ONLY ================')
disp(lme_wash)
disp(anova(lme_wash))

coefTbl = lme_wash.Coefficients;
disp(' ')
disp('================ FIXED EFFECT COEFFICIENTS ================')
disp(coefTbl)

% Extract comparisons vs BL
row_WOE = strcmp(coefTbl.Name, 'Condition_WOE');
row_WO  = strcmp(coefTbl.Name, 'Condition_WO');

if any(row_WOE)
    fprintf('\nWOE vs BL: estimate = %.6f, SE = %.6f, t = %.3f, p = %.4g\n', ...
        coefTbl.Estimate(row_WOE), coefTbl.SE(row_WOE), ...
        coefTbl.tStat(row_WOE), coefTbl.pValue(row_WOE));
end

if any(row_WO)
    fprintf('WO vs BL:  estimate = %.6f, SE = %.6f, t = %.3f, p = %.4g\n', ...
        coefTbl.Estimate(row_WO), coefTbl.SE(row_WO), ...
        coefTbl.tStat(row_WO), coefTbl.pValue(row_WO));
end

%% ============================================================
% MODEL 2: ADD SESSION NUMBER
% Tests whether washout effects remain after accounting for session number
% ============================================================

lme_wash_sess = fitlme(LongTable, ...
    'Decel ~ Condition + SessionNum + (1|Mouse) + (1|SessionID)');

disp(' ')
disp('================ HIERARCHICAL MODEL: CONDITION + SESSION ================')
disp(lme_wash_sess)
disp(anova(lme_wash_sess))

coefTbl2 = lme_wash_sess.Coefficients;
disp(' ')
disp('================ FIXED EFFECT COEFFICIENTS (WITH SESSION) ================')
disp(coefTbl2)

row_WOE2 = strcmp(coefTbl2.Name, 'Condition_WOE');
row_WO2  = strcmp(coefTbl2.Name, 'Condition_WO');
row_Sess = strcmp(coefTbl2.Name, 'SessionNum');

if any(row_WOE2)
    fprintf('\nWOE vs BL (controlling for SessionNum): estimate = %.6f, SE = %.6f, t = %.3f, p = %.4g\n', ...
        coefTbl2.Estimate(row_WOE2), coefTbl2.SE(row_WOE2), ...
        coefTbl2.tStat(row_WOE2), coefTbl2.pValue(row_WOE2));
end

if any(row_WO2)
    fprintf('WO vs BL  (controlling for SessionNum): estimate = %.6f, SE = %.6f, t = %.3f, p = %.4g\n', ...
        coefTbl2.Estimate(row_WO2), coefTbl2.SE(row_WO2), ...
        coefTbl2.tStat(row_WO2), coefTbl2.pValue(row_WO2));
end

if any(row_Sess)
    fprintf('SessionNum effect: estimate = %.6f, SE = %.6f, t = %.3f, p = %.4g\n', ...
        coefTbl2.Estimate(row_Sess), coefTbl2.SE(row_Sess), ...
        coefTbl2.tStat(row_Sess), coefTbl2.pValue(row_Sess));
end
%%
% ============================================================
% MODEL 3: CONDITION x SESSION INTERACTION
% Tests whether WOE or WO change across sessions differently from BL
% ============================================================

lme_wash_int = fitlme(LongTable, ...
    'Decel ~ Condition * SessionNum + (1|Mouse) + (1|SessionID)');

disp(' ')
disp('================ HIERARCHICAL MODEL: CONDITION x SESSION ================')
disp(lme_wash_int)
disp(anova(lme_wash_int))

% ============================================================
% OPTIONAL: RANDOM SLOPE MODEL
% Use only if your dataset is large/stable enough
% This lets condition effects vary by mouse
% ============================================================

try
    lme_wash_rs = fitlme(LongTable, ...
        'Decel ~ Condition + SessionNum + (Condition|Mouse) + (1|SessionID)');

    disp(' ')
    disp('================ HIERARCHICAL MODEL: RANDOM SLOPE ================')
    disp(lme_wash_rs)
    disp(anova(lme_wash_rs))
catch ME
    warning('Random slope model failed to converge: %s', ME.message);
end

% ============================================================
% OPTIONAL: PREDICTED CONDITION MEANS FROM MODEL
% ============================================================

PredTable = table( ...
    categorical({'dummy';'dummy';'dummy'}), ...
    categorical({'dummy__dummy';'dummy__dummy';'dummy__dummy'}), ...
    [mean(LongTable.SessionNum,'omitnan'); mean(LongTable.SessionNum,'omitnan'); mean(LongTable.SessionNum,'omitnan')], ...
    categorical({'BL';'WOE';'WO'}, {'BL','WOE','WO'}), ...
    'VariableNames', {'Mouse','SessionID','SessionNum','Condition'});

PredVals = predict(lme_wash_sess, PredTable);

figure('Color','w','Position',[100 100 360 320]); hold on;
bar(1:3, PredVals, 0.65, 'FaceColor',[0.8 0.8 0.8], 'EdgeColor','k');
plot(1:3, PredVals, 'ko', 'MarkerFaceColor','k', 'MarkerSize',6);
xticks(1:3);
xticklabels({'BL','WOE','WO'});
ylabel('Model-predicted deceleration');
title('Hierarchical mixed-effects model');
box off;
set(gca,'FontSize',12,'LineWidth',1.2);

%% WITH STIMULATD TRIALS

%% ============================================================
% HIERARCHICAL MIXED-EFFECTS MODEL
% BL vs AD vs WOE vs WO (PER-SESSION DATA)
%
% Tests:
%   AD vs BL   (stimulation effect)
%   WOE vs BL  (early washout)
%   WO vs BL   (washout)
%
% Model:
%   Decel ~ Condition + SessionNum + (1|Mouse) + (1|SessionID)
%
% Key:
%   BL = reference level
% ============================================================

% ---------------- BUILD LONG TABLE ----------------
LongRows = {};

for i = 1:height(SessionTable)

    mouse_i   = SessionTable.Mouse{i};
    sess_i    = SessionTable.Session{i};
    sessnum_i = SessionTable.SessionNum(i);

    % --- BL ---
    LongRows(end+1,:) = {mouse_i, sess_i, sessnum_i, 'BL',  SessionTable.BL(i)}; %#ok<SAGROW>

    % --- AD (stimulated) ---
    LongRows(end+1,:) = {mouse_i, sess_i, sessnum_i, 'AD',  SessionTable.AD(i)}; %#ok<SAGROW>

    % --- WOE ---
    LongRows(end+1,:) = {mouse_i, sess_i, sessnum_i, 'WOE', SessionTable.WOE(i)}; %#ok<SAGROW>

    % --- WO ---
    LongRows(end+1,:) = {mouse_i, sess_i, sessnum_i, 'WO',  SessionTable.WO(i)}; %#ok<SAGROW>
end

LongTable = cell2table(LongRows, ...
    'VariableNames', {'Mouse','Session','SessionNum','Condition','Decel'});

% remove NaNs
LongTable = LongTable(~isnan(LongTable.Decel), :);

% categorical variables
LongTable.Mouse     = categorical(LongTable.Mouse);
LongTable.Session   = categorical(LongTable.Session);
LongTable.Condition = categorical(LongTable.Condition);

% set BL as reference
LongTable.Condition = reordercats(LongTable.Condition, {'BL','AD','WOE','WO'});

% nested session ID
LongTable.SessionID = categorical(strcat(string(LongTable.Mouse), "__", string(LongTable.Session)));

disp('================ LONG TABLE ================')
disp(LongTable)

%% ============================================================
% MODEL 1: CONDITION EFFECT (PRIMARY MODEL)
% ============================================================

lme_all = fitlme(LongTable, ...
    'Decel ~ Condition + (1|Mouse) + (1|SessionID)');

disp(' ')
disp('================ MODEL: CONDITION ONLY ================')
disp(lme_all)
disp(anova(lme_all))

coefTbl = lme_all.Coefficients;
disp(' ')
disp('================ FIXED EFFECTS ================')
disp(coefTbl)

% --- extract comparisons vs BL ---
row_AD  = strcmp(coefTbl.Name, 'Condition_AD');
row_WOE = strcmp(coefTbl.Name, 'Condition_WOE');
row_WO  = strcmp(coefTbl.Name, 'Condition_WO');

fprintf('\n===== EFFECTS VS BASELINE =====\n')

if any(row_AD)
    fprintf('AD vs BL:  estimate = %.6f, SE = %.6f, t = %.3f, p = %.4g\n', ...
        coefTbl.Estimate(row_AD), coefTbl.SE(row_AD), ...
        coefTbl.tStat(row_AD), coefTbl.pValue(row_AD));
end

if any(row_WOE)
    fprintf('WOE vs BL: estimate = %.6f, SE = %.6f, t = %.3f, p = %.4g\n', ...
        coefTbl.Estimate(row_WOE), coefTbl.SE(row_WOE), ...
        coefTbl.tStat(row_WOE), coefTbl.pValue(row_WOE));
end

if any(row_WO)
    fprintf('WO vs BL:  estimate = %.6f, SE = %.6f, t = %.3f, p = %.4g\n', ...
        coefTbl.Estimate(row_WO), coefTbl.SE(row_WO), ...
        coefTbl.tStat(row_WO), coefTbl.pValue(row_WO));
end

%% ============================================================
% MODEL 2: ADD SESSION NUMBER (IMPORTANT FOR LEARNING)
% ============================================================

lme_all_sess = fitlme(LongTable, ...
    'Decel ~ Condition + SessionNum + (1|Mouse) + (1|SessionID)');

disp(' ')
disp('================ MODEL: CONDITION + SESSION ================')
disp(lme_all_sess)
disp(anova(lme_all_sess))

coefTbl2 = lme_all_sess.Coefficients;

row_AD2  = strcmp(coefTbl2.Name, 'Condition_AD');
row_WOE2 = strcmp(coefTbl2.Name, 'Condition_WOE');
row_WO2  = strcmp(coefTbl2.Name, 'Condition_WO');
row_Sess = strcmp(coefTbl2.Name, 'SessionNum');

fprintf('\n===== EFFECTS VS BASELINE (controlling for session) =====\n')

if any(row_AD2)
    fprintf('AD vs BL:  estimate = %.6f, p = %.4g\n', ...
        coefTbl2.Estimate(row_AD2), coefTbl2.pValue(row_AD2));
end

if any(row_WOE2)
    fprintf('WOE vs BL: estimate = %.6f, p = %.4g\n', ...
        coefTbl2.Estimate(row_WOE2), coefTbl2.pValue(row_WOE2));
end

if any(row_WO2)
    fprintf('WO vs BL:  estimate = %.6f, p = %.4g\n', ...
        coefTbl2.Estimate(row_WO2), coefTbl2.pValue(row_WO2));
end

if any(row_Sess)
    fprintf('Session effect: estimate = %.6f, p = %.4g\n', ...
        coefTbl2.Estimate(row_Sess), coefTbl2.pValue(row_Sess));
end

%% ============================================================
% MODEL 3: CONDITION x SESSION (LEARNING DIFFERENCES)
% ============================================================

lme_all_int = fitlme(LongTable, ...
    'Decel ~ Condition * SessionNum + (1|Mouse) + (1|SessionID)');

disp(' ')
disp('================ MODEL: CONDITION x SESSION ================')
disp(lme_all_int)
disp(anova(lme_all_int))

%% ============================================================
% OPTIONAL: RANDOM SLOPES (ADVANCED)
% ============================================================

try
    lme_rs = fitlme(LongTable, ...
        'Decel ~ Condition + SessionNum + (Condition|Mouse) + (1|SessionID)');

    disp(' ')
    disp('================ MODEL: RANDOM SLOPES ================')
    disp(lme_rs)
    disp(anova(lme_rs))
catch ME
    warning('Random slope model failed: %s', ME.message);
end

% predicted values for each condition
CondList = categorical({'BL','AD','WOE','WO'}, {'BL','AD','WOE','WO'});

PredTable = table;
PredTable.Mouse = categorical(repmat("dummy",4,1));
PredTable.SessionID = categorical(repmat("dummy__dummy",4,1));
PredTable.SessionNum = repmat(mean(LongTable.SessionNum,'omitnan'),4,1);
PredTable.Condition = CondList';

yhat = predict(lme_all_sess, PredTable);

figure; hold on
bar(1:4, yhat, 0.6, 'FaceColor',[0.8 0.8 0.8], 'EdgeColor','k')
plot(1:4, yhat, 'ko','MarkerFaceColor','k')

xticks(1:4)
xticklabels({'BL','AD','WOE','WO'})
ylabel('Model-estimated deceleration')
title('Mixed-effects model estimates')
yline(0,'k--')

%%
X = [1 2 3];
Y = [AnimalTable.AD_delta, AnimalTable.WOE_delta, AnimalTable.WO_delta];

figure; hold on
for i = 1:size(Y,1)
    plot(X, Y(i,:), '-o','Color',[0.7 0.7 0.7])
end

errorbar(X, mean(Y), std(Y)/sqrt(size(Y,1)), 'k-o','LineWidth',2)

xticks(X)
xticklabels({'AD','WOE','WO'})
yline(0,'k--')
ylabel('\Delta from baseline')
title('Adaptation trajectory')

%%
effects = [AD_est, WOE_est, WO_est];
errors  = [AD_se,  WOE_se,  WO_se];

figure; hold on
errorbar(1:3, effects, errors, 'ko','LineWidth',2)

xticks(1:3)
xticklabels({'AD','WOE','WO'})
yline(0,'k--')
ylabel('Effect vs BL')
title('Model coefficients')

%% 
% JUST AD AND WO

%% ============================================================
% MIXED EFFECTS MODEL: BL vs AD vs WO
% Main fixed effects of interest:
%   AD vs BL
%   WO vs BL
%
% Model:
%   Decel ~ Condition + SessionNum + (1|Mouse) + (1|SessionID)
% ============================================================

%---------------- BUILD LONG TABLE ----------------
LongRows = {};

for i = 1:height(SessionTable)

    mouse_i   = SessionTable.Mouse{i};
    sess_i    = SessionTable.Session{i};
    sessnum_i = SessionTable.SessionNum(i);

    % --- BL ---
    LongRows(end+1,:) = {mouse_i, sess_i, sessnum_i, 'BL', SessionTable.BL(i)}; %#ok<SAGROW>

    % --- AD ---
    LongRows(end+1,:) = {mouse_i, sess_i, sessnum_i, 'AD', SessionTable.AD(i)}; %#ok<SAGROW>

    % --- WO ---
    LongRows(end+1,:) = {mouse_i, sess_i, sessnum_i, 'WO', SessionTable.WO(i)}; %#ok<SAGROW>
end

LongTable = cell2table(LongRows, ...
    'VariableNames', {'Mouse','Session','SessionNum','Condition','Decel'});

LongTable = LongTable(~isnan(LongTable.Decel), :);

LongTable.Mouse     = categorical(LongTable.Mouse);
LongTable.Session   = categorical(LongTable.Session);
LongTable.Condition = categorical(LongTable.Condition);

% BL is reference
LongTable.Condition = reordercats(LongTable.Condition, {'BL','AD','WO'});

% nested session ID
LongTable.SessionID = categorical(strcat(string(LongTable.Mouse), "__", string(LongTable.Session)));

disp('================ LONG TABLE ================')
disp(LongTable)

%% ============================================================
% FIT MODEL
% ============================================================

lme_main = fitlme(LongTable, ...
    'Decel ~ Condition + SessionNum + (1|Mouse) + (1|SessionID)');

disp(' ')
disp('================ MODEL: BL vs AD vs WO ================')
disp(lme_main)

anovaTbl = anova(lme_main);
disp(anovaTbl)

coefTbl = lme_main.Coefficients;
disp(' ')
disp('================ FIXED EFFECTS ================')
disp(coefTbl)

%% ============================================================
% EXTRACT STATS TO REPORT
% ============================================================

row_AD   = strcmp(coefTbl.Name, 'Condition_AD');
row_WO   = strcmp(coefTbl.Name, 'Condition_WO');
row_Sess = strcmp(coefTbl.Name, 'SessionNum');
row_Int  = strcmp(coefTbl.Name, '(Intercept)');

BL_est = coefTbl.Estimate(row_Int);
BL_se  = coefTbl.SE(row_Int);

AD_est = coefTbl.Estimate(row_AD);
AD_se  = coefTbl.SE(row_AD);
AD_t   = coefTbl.tStat(row_AD);
AD_p   = coefTbl.pValue(row_AD);

WO_est = coefTbl.Estimate(row_WO);
WO_se  = coefTbl.SE(row_WO);
WO_t   = coefTbl.tStat(row_WO);
WO_p   = coefTbl.pValue(row_WO);

if any(row_Sess)
    Sess_est = coefTbl.Estimate(row_Sess);
    Sess_se  = coefTbl.SE(row_Sess);
    Sess_t   = coefTbl.tStat(row_Sess);
    Sess_p   = coefTbl.pValue(row_Sess);
else
    Sess_est = NaN; Sess_se = NaN; Sess_t = NaN; Sess_p = NaN;
end

fprintf('\n================ MAIN EFFECTS VS BASELINE ================\n');
fprintf('AD vs BL: estimate = %.6f, SE = %.6f, t = %.3f, p = %.4g\n', AD_est, AD_se, AD_t, AD_p);
fprintf('WO vs BL: estimate = %.6f, SE = %.6f, t = %.3f, p = %.4g\n', WO_est, WO_se, WO_t, WO_p);

if ~isnan(Sess_p)
    fprintf('SessionNum: estimate = %.6f, SE = %.6f, t = %.3f, p = %.4g\n', ...
        Sess_est, Sess_se, Sess_t, Sess_p);
end

%% ============================================================
% MODEL-ESTIMATED VALUES FOR BL / AD / WO
% ============================================================

CondList = categorical({'BL','AD','WO'}, {'BL','AD','WO'});

PredTable = table;
PredTable.Mouse     = categorical(repmat("dummy",3,1));
PredTable.SessionID = categorical(repmat("dummy__dummy",3,1));
PredTable.SessionNum = repmat(mean(LongTable.SessionNum,'omitnan'),3,1);
PredTable.Condition = CondList';

yhat = predict(lme_main, PredTable);

%% ============================================================
% PLOT 1: MODEL-ESTIMATED DECELERATION
% with AD vs BL and WO vs BL stats written below
% ============================================================

figure('Color','w','Position',[100 100 250 420]); 
hold on

bar(1:3, yhat, 0.62, 'FaceColor',[0.8 0.8 0.8], 'EdgeColor','k', 'LineWidth',1.2)
plot(1:3, yhat, 'ko', 'MarkerFaceColor','k', 'MarkerSize',8)

xticks(1:3)
xticklabels({'BL','AD','WO'})
ylabel('Model-estimated deceleration')
title('Mixed-effects model estimates')
yline(0,'k--','LineWidth',1)

box off
set(gca,'FontSize',12,'LineWidth',1.2)

yl = ylim;
yr = yl(2) - yl(1);

txt1 = sprintf('AD vs BL: \\beta = %.2f, SE = %.2f, t = %.2f, p = %.3g', AD_est, AD_se, AD_t, AD_p);
txt2 = sprintf('WO vs BL: \\beta = %.2f, SE = %.2f, t = %.2f, p = %.3g', WO_est, WO_se, WO_t, WO_p);

ylim([yl(1) - 0.22*yr, yl(2)])

text(0.70, yl(1) - 0.10*yr, txt1, ...
    'FontSize',10, 'Interpreter','tex', 'HorizontalAlignment','left');

text(0.70, yl(1) - 0.18*yr, txt2, ...
    'FontSize',10, 'Interpreter','tex', 'HorizontalAlignment','left');

%% ============================================================
% PLOT 2: MODEL COEFFICIENTS (EFFECT VS BL)
% with stats written below
% ============================================================

effects = [AD_est, WO_est];
errors  = [AD_se,  WO_se];

figure('Color','w','Position',[100 100 200 420]);
hold on

errorbar(1:2, effects, errors, 'ko', ...
    'LineWidth',2, 'MarkerFaceColor','w', 'MarkerSize',10)

yline(0,'k--','LineWidth',1)

xlim([0.5 2.5])
xticks(1:2)
xticklabels({'AD vs BL','WO vs BL'})
ylabel('Effect vs BL')
title('Model coefficients')

box off
set(gca,'FontSize',12,'LineWidth',1.2)

yl = ylim;
yr = yl(2) - yl(1);
ylim([yl(1) - 0.25*yr, yl(2)])

text(0.72, yl(1) - 0.11*yr, ...
    sprintf('AD: \\beta = %.2f, SE = %.2f, t = %.2f, p = %.3g', AD_est, AD_se, AD_t, AD_p), ...
    'FontSize',10, 'Interpreter','tex', 'HorizontalAlignment','left');

text(0.72, yl(1) - 0.20*yr, ...
    sprintf('WO: \\beta = %.2f, SE = %.2f, t = %.2f, p = %.3g', WO_est, WO_se, WO_t, WO_p), ...
    'FontSize',10, 'Interpreter','tex', 'HorizontalAlignment','left');

%% ============================================================
% PLOT 3: ANIMAL-LEVEL DELTA FROM BASELINE
% plots AD-BL and WO-BL
% ============================================================

Y = [AnimalTable.AD_delta, AnimalTable.WO_delta];
X = [1 2];

figure('Color','w','Position',[100 100 500 380]); 
hold on

for i = 1:size(Y,1)
    plot(X, Y(i,:), '-o', ...
        'Color',[0.72 0.72 0.72], ...
        'MarkerFaceColor',[0.72 0.72 0.72], ...
        'MarkerEdgeColor',[0.72 0.72 0.72], ...
        'LineWidth',1.1, ...
        'MarkerSize',6);
end

mY = mean(Y,1,'omitnan');
sY = std(Y,0,1,'omitnan') ./ sqrt(sum(~isnan(Y),1));

errorbar(X, mY, sY, 'k-o', ...
    'LineWidth',2.2, ...
    'MarkerFaceColor','k', ...
    'MarkerEdgeColor','k', ...
    'MarkerSize',7, ...
    'CapSize',0);

yline(0,'k--','LineWidth',1)

xlim([0.7 2.3])
xticks(X)
xticklabels({'AD-BL','WO-BL'})
ylabel('\Delta from baseline')
title('Animal-level delta from baseline')

box off
set(gca,'FontSize',12,'LineWidth',1.2)

% ============================================================
% OPTIONAL: MAIN EFFECT OF CONDITION FROM MODEL ANOVA
% ============================================================

disp(' ')
disp('================ ANOVA TABLE ================')
disp(anovaTbl)