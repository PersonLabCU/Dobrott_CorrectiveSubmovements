%% Replace all reaches aligned CSM with only the NoCSM template

an=fieldnames(ExcitRR_Sham);

thresh=0.16;
dirdir= 1;

ExcitRR_Ctrl= ExcitRR_Sham;

for a=1:length(an)
   
    CtrlReaches=ExcitRR_Ctrl.(an{a}).NoReReach;

    % NoCSM algined to T0
    CST0_CSM=AlignThreshCSM_T0(CtrlReaches,thresh,dirdir,200);
    CST0_CSM=stimvel_accel(CST0_CSM);

    % Ctrl to T1
    CST1_CSM=AlignThreshCSM_T1(CtrlReaches,thresh,dirdir,200);
    CST1_CSM=stimvel_accel(CST1_CSM);

    % Ctrl to T2
    CST2_CSM=AlignThreshCSM_T2(CtrlReaches,thresh,dirdir,200);
    CST2_CSM=stimvel_accel(CST2_CSM);

    % Ctrl to T3
    CST3_CSM=AlignThreshCSM_T3(CtrlReaches,thresh,dirdir,200);
    CST3_CSM=stimvel_accel(CST3_CSM);
    
    % Save to big structure 

     %Stimvel

     ExcitRR_Ctrl.(an{a}).CST0_CSM=CST0_CSM;
     ExcitRR_Ctrl.(an{a}).CST1_CSM=CST1_CSM;
     ExcitRR_Ctrl.(an{a}).CST2_CSM=CST2_CSM;
     ExcitRR_Ctrl.(an{a}).CST3_CSM=CST3_CSM;
        

end 

%% Now make big stimvel structures for CSM bins aligned to average CSM time for each bin 

T0_bigCSM=struct;
T1_bigCSM=struct;
T2_bigCSM= struct;
T3_bigCSM=struct;

ar=fieldnames(ExcitRR_Ctrl);
for a = 1:length(ar) 
   T0_bigCSM.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).ST0_CSM.Average;
   T1_bigCSM.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).ST1_CSM.Average;
   T2_bigCSM.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).ST2_CSM.Average;
   T3_bigCSM.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).ST3_CSM.Average;

    
end

T0_bigCSM.Average = nanmean(T0_bigCSM.stimvel,3);
T0_bigCSM.Std = nanstd(T0_bigCSM.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T0_bigCSM.Average),1); 
for nnum = 1:length(T0_bigCSM.Average)
    samp = squeeze(T0_bigCSM.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
T0_bigCSM.N = N;

T1_bigCSM.Average = nanmean(T1_bigCSM.stimvel,3);
T1_bigCSM.Std = nanstd(T1_bigCSM.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T1_bigCSM.Average),1); 
for nnum = 1:length(T1_bigCSM.Average)
    samp = squeeze(T1_bigCSM.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
T1_bigCSM.N = N;


T2_bigCSM.Average = nanmean(T2_bigCSM.stimvel,3);
T2_bigCSM.Std = nanstd(T2_bigCSM.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T2_bigCSM.Average),1); 
for nnum = 1:length(T2_bigCSM.Average)
    samp = squeeze(T2_bigCSM.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
T2_bigCSM.N = N;

%
T3_bigCSM.Average = nanmean(T3_bigCSM.stimvel,3);
T3_bigCSM.Std = nanstd(T3_bigCSM.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T3_bigCSM.Average),1); 
for nnum = 1:length(T3_bigCSM.Average)
    samp = squeeze(T3_bigCSM.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
T3_bigCSM.N = N;

%% Combine Stimvel and look at plots  - make big stimvel for control reaches that are aligned to the different CSM time bins 

CT0_bigCSM=struct;
CT1_bigCSM=struct;
CT2_bigCSM= struct;
CT3_bigCSM=struct;


ar=fieldnames(ExcitRR_Ctrl);
for a = 1:length(ar) 
   CT0_bigCSM.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).CST0_CSM.Average;
   CT1_bigCSM.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).CST1_CSM.Average;
   CT2_bigCSM.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).CST2_CSM.Average;
   CT3_bigCSM.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).CST3_CSM.Average;
    
end

CT0_bigCSM.Average = nanmean(CT0_bigCSM.stimvel,3);
CT0_bigCSM.Std = nanstd(CT0_bigCSM.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(CT0_bigCSM.Average),1); 
for nnum = 1:length(CT0_bigCSM.Average)
    samp = squeeze(CT0_bigCSM.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
CT0_bigCSM.N = N;

CT1_bigCSM.Average = nanmean(CT1_bigCSM.stimvel,3);
CT1_bigCSM.Std = nanstd(CT1_bigCSM.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(CT1_bigCSM.Average),1); 
for nnum = 1:length(CT1_bigCSM.Average)
    samp = squeeze(CT1_bigCSM.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
CT1_bigCSM.N = N;


CT2_bigCSM.Average = nanmean(CT2_bigCSM.stimvel,3);
CT2_bigCSM.Std = nanstd(CT2_bigCSM.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(CT2_bigCSM.Average),1); 
for nnum = 1:length(CT2_bigCSM.Average)
    samp = squeeze(CT2_bigCSM.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
CT2_bigCSM.N = N;

%
CT3_bigCSM.Average = nanmean(CT3_bigCSM.stimvel,3);
CT3_bigCSM.Std = nanstd(CT3_bigCSM.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(CT3_bigCSM.Average),1); 
for nnum = 1:length(CT3_bigCSM.Average)
    samp = squeeze(CT3_bigCSM.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
CT3_bigCSM.N = N;

%% %% ===== Make Final Master Structures for CSM and Control =====

Master_CSM = struct;
Master_CTRL = struct;

% ---- Combine the averaged bin traces into one matrix ----
% Dimensions: time × variable(velocity etc.) × binIndex (1–4)

Master_CSM.AllBins = cat(3, ...
    T0_bigCSM.Average, ...
    T1_bigCSM.Average, ...
    T2_bigCSM.Average, ...
    T3_bigCSM.Average);

Master_CTRL.AllBins = cat(3, ...
    CT0_bigCSM.Average, ...
    CT1_bigCSM.Average, ...
    CT2_bigCSM.Average, ...
    CT3_bigCSM.Average);

% ---- Compute the overall average across all bins ----
Master_CSM.Average = nanmean(Master_CSM.AllBins, 3);
Master_CSM.Std     = nanstd(Master_CSM.AllBins, 0, 3);

Master_CTRL.Average = nanmean(Master_CTRL.AllBins, 3);
Master_CTRL.Std     = nanstd(Master_CTRL.AllBins, 0, 3);

% ---- Compute N across all bins ----
N_csm = zeros(size(Master_CSM.Average,1),1);
for i = 1:length(N_csm)
    samp = squeeze(Master_CSM.AllBins(i,2,:));
    N_csm(i) = sum(~isnan(samp));
end
Master_CSM.N = N_csm;

N_ctrl = zeros(size(Master_CTRL.Average,1),1);
for i = 1:length(N_ctrl)
    samp = squeeze(Master_CTRL.AllBins(i,2,:));
    N_ctrl(i) = sum(~isnan(samp));
end
Master_CTRL.N = N_ctrl;


%% ===== Weighted Master Structures for CSM and Control =====

Master_CSM = struct;
Master_CTRL = struct;

% Collect bin averages
CSM_bins  = cat(3, T0_bigCSM.Average,  T1_bigCSM.Average,  T2_bigCSM.Average,  T3_bigCSM.Average);
CTRL_bins = cat(3, CT0_bigCSM.Average, CT1_bigCSM.Average, CT2_bigCSM.Average, CT3_bigCSM.Average);

% Collect counts per bin (number of reaches)
CSM_counts = [sum(~isnan(T0_bigCSM.Average(:,2))), ...
              sum(~isnan(T1_bigCSM.Average(:,2))), ...
              sum(~isnan(T2_bigCSM.Average(:,2))), ...
              sum(~isnan(T3_bigCSM.Average(:,2)))];

CTRL_counts = [sum(~isnan(CT0_bigCSM.Average(:,2))), ...
               sum(~isnan(CT1_bigCSM.Average(:,2))), ...
               sum(~isnan(CT2_bigCSM.Average(:,2))), ...
               sum(~isnan(CT3_bigCSM.Average(:,2)))];

% Save counts into Master structures
Master_CSM.N  = CSM_counts;
Master_CTRL.N = CTRL_counts;

% Normalize weights
CSM_w  = CSM_counts  ./ sum(CSM_counts);
CTRL_w = CTRL_counts ./ sum(CTRL_counts);

% Weighted averages
Master_CSM.Average  = sum(CSM_bins  .* permute(CSM_w,  [1 3 2]), 3);
Master_CTRL.Average = sum(CTRL_bins .* permute(CTRL_w, [1 3 2]), 3);

% Weighted STD
Master_CSM.Std = sqrt( sum((CSM_bins  - Master_CSM.Average).^2 .* permute(CSM_w, [1 3 2]), 3) );
Master_CTRL.Std = sqrt( sum((CTRL_bins - Master_CTRL.Average).^2 .* permute(CTRL_w, [1 3 2]), 3) );

%% ===== Plot Weighted CSM, Control, and the Difference =====
figure; clf; hold on;

timevec = Master_CSM.Average(:,1);   % if your structures contain time
dim = 8;                     % velocity/acc/pos index

CSM   = Master_CSM.Average(:,dim);
CTRL  = Master_CTRL.Average(:,dim);

CSM_sem  = Master_CSM.Std(:,dim) ./ sqrt(sum(Master_CSM.N));
CTRL_sem = Master_CTRL.Std(:,dim) ./ sqrt(sum(Master_CTRL.N));

% --- Plot Control with shading ---
fill([timevec; flipud(timevec)], ...
     [CTRL - CTRL_sem; flipud(CTRL + CTRL_sem)], ...
     [0.5 0.5 0.5], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
plot(timevec, CTRL, 'Color',[0.2 0.2 0.2], 'LineWidth', 2);

% --- Plot CSM with shading ---
fill([timevec; flipud(timevec)], ...
     [CSM - CSM_sem; flipud(CSM + CSM_sem)], ...
     CSMcolor, 'FaceAlpha', 0.25, 'EdgeColor','none');
plot(timevec, CSM, 'Color', CSMcolor, 'LineWidth', 2);

% --- Delta curve ---
Delta = CSM - CTRL;
plot(timevec, Delta, 'Color', [0.1 0.4 0.8], 'LineWidth', 2.2);

% --- 0 line ---
yline(0, '--k', 'LineWidth', 1);

% --- Formatting ---
xlabel('Time (s)');
ylabel('Outward Acceleration');
xlim([min(timevec) max(timevec)]);
set(gca, 'FontSize', 12, 'LineWidth', 1.4);

legend({'CTRL SEM','CTRL mean','CSM SEM','CSM mean','CSM - CTRL'}, ...
       'Location','Best','Box','off');

title('Weighted CSM vs Control and Delta');

%% ===== Z-score Difference Between CSM and CTRL =====
% dim = 2 already defined above
% Extract time series
dim=8;
CSM_curve  = Master_CSM.Average(:, dim);
CTRL_curve = Master_CTRL.Average(:, dim);

% Compute z-score of the difference curve relative to CTRL variability
Diff_curve = CSM_curve - CTRL_curve;

% Use CTRL std as baseline variability
baseline_std = Master_CTRL.Std(:, dim);

Zcurve = Diff_curve ./ baseline_std;

% Smooth (optional)
Zcurve_smooth = movmean(Zcurve, 5);

% ===== Add Heatmap Under PSTH =====
% Adjust main axis to make room
ax1 = gca;
ax1.Position(2) = ax1.Position(2) + 0.12;
ax1.Position(4) = ax1.Position(4) - 0.12;

% Create new axis for heatmap
ax2 = axes('Position', [ax1.Position(1), ax1.Position(2)-0.18, ax1.Position(3), 0.06]);

% Convert Z-curve to 2D for imagesc
Zmat = Zcurve_smooth';          % 1 × time
imagesc(timevec, 1, Zmat, 'Parent', ax2);

colormap(ax2, brewermap([], 'PrGn'));
caxis(ax2, [-2 2]);

% Remove ticks
set(ax2, 'YTick', [], 'XTick', []);

% Add colorbar
colorbar(ax2, 'Position', [0.92, ax2.Position(2), 0.02, ax2.Position(4)]);
xlabel(ax2, ''); ylabel(ax2, '');

title(ax1, 'Weighted CSM vs Control with Z-score Heatmap');
%% ===== Pointwise P-values and Z-scores Between Weighted CSM and CTRL =====
% new version

dim = 2;   % velocity / pos / acc dimension

CSM_curve  = Master_CSM.Average(:, dim);
CTRL_curve = Master_CTRL.Average(:, dim);

CSM_std  = Master_CSM.Std(:, dim);
CTRL_std = Master_CTRL.Std(:, dim);

ntime = length(CSM_curve);

% Match your original variable names and shapes
parray0  = NaN(ntime, 2);     % [pval, zval]
zarray0  = NaN(ntime, 1);     % mean z-score per timepoint
Zscores0 = NaN(ntime, 5);     % dummy structure (kept for compatibility)

for pnum = 1:ntime
    
    % --- Analog of a4b and b4b (experimental vs control) ---
    a4b = CSM_curve(pnum);    % experimental value
    b4b = CTRL_curve(pnum);   % control value
    
    % --- Approximate signed-rank using Gaussian z ---
    diff_val = a4b - b4b;
    pooled_std = CTRL_std(pnum);   % baseline variability
    
    if pooled_std == 0
        zval = 0;
        pval = 1;
    else
        zval = diff_val / pooled_std;
        pval = 2 * (1 - normcdf(abs(zval)));   % 2-sided
    end
    
    parray0(pnum, 1) = pval;
    parray0(pnum, 2) = zval;   % match your original "z-value" slot
    
    % Pointwise Z-score (same meaning as your previous zscores)
    zscores_point = zval;      % now a scalar per time point
    Zscores0(pnum, 1) = zscores_point;
    
    % Mean z-score = same as point z-score here
    zarray0(pnum, 1) = zscores_point;
end

% ===== Z-score Plot (same layout as your old version) =====
set(gcf, 'Position', [100, 100, 800, 600]);

% Main PSTH plot
StimvelPlot(Master_CSM, dim, CSMcolor, 1);
StimvelPlot(Master_CTRL, dim, 'k', 0);

yline(0, '--k');

xlim([-0.2 0.2]);
ylim('auto');
xlabel('Time from threshold crossing (s)');
ylabel('Outward Acceleration');

% Adjust axis for heatmap
ax1 = gca;
ax1.Position(2) = ax1.Position(2) + 0.1;
ax1.Position(4) = ax1.Position(4) - 0.1;

% ===== Heatmap of Z-scores Underneath =====
ax2 = axes('Position', [ax1.Position(1), ax1.Position(2) - 0.2, ax1.Position(3), 0.05]);

% 1×T heatmap (matches your zarray0 style)
imagesc(zarray0', 'Parent', ax2);
colormap(ax2, brewermap([], 'PrGn'));
caxis([-2, 2]);

% Colorbar
colorbar(ax2, 'Position', [0.92, ax2.Position(2), 0.02, ax2.Position(4)]);

set(ax2, 'XTick', [], 'YTick', [], 'YColor', 'none', 'XColor', 'none');


%% Plot stimvel for CSM

T0color=rgb('Thistle');
T1color=rgb('MediumOrchid');
T2color=rgb('SlateBlue');
T3color=rgb('DarkSlateBlue');
NRRcolor=rgb('OrangeRed');
stimcolor=rgb('DodgerBlue');

dim =5;

StimvelPlot(T0_bigCSM,dim,T0color,1)
StimvelPlot(T1_bigCSM,dim,T1color,0)
StimvelPlot(T2_bigCSM,dim,T2color,0)
StimvelPlot(T3_bigCSM,dim,T3color,0)
%StimvelPlot(NRR_bigCSM,dim,NRRcolor,0)

%xlim([-0.25 0.25])
xlim([-0.2 0.2])
%xticks(-0.2:0.01:0.2)
ylim('auto')
%ylim([-0.6 1.4])
%ylim([0 0.25])
%ylim([-20 20])
%pbaspect([1 2 1])
%ylim([0 0.3])
xline(0,'--')
yline(0,'--')

xlabel('Time from CSM(s)')
ylabel('Outward Position (cm)')
%ylabel('Outward Velocity (cm/s)')
%ylabel('Outward acceleration')
%% PLot all reaches together aligned to the different CSM times and then combined again 

CSMcolor=rgb('DarkViolet');

dim =5;

StimvelPlot(Master_CSM,dim,CSMcolor,1)
StimvelPlot(Master_CTRL,dim,'k',0)
ylim('auto')
xlim([-0.2 0.2])



%% Master CSM 
% Initialize arrays for p-values and Z-scores
parray0 = NaN(size(Master_CSM.AllBins, 1), 1);  % Store p-values
zarray0 =   NaN(size(Master_CSM.AllBins,1), 1)   ;    % Store Z-scores (same shape as stimvel)
Zscores0= NaN(size(Master_CSM.AllBins, 1), 5);

for pnum = 1:length(parray0)
    % Extract trial data for C_big (control) and T0_big (experimental)
    a4b = squeeze(Master_CSM.AllBins(pnum, 2, :));  % Experimental data (T0_big)
    b4b = squeeze(Master_CTRL.AllBins(pnum, 2, :));   % Control data (C_big), matching trial

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



%zscore plot 
% Adjust the figure to provide more vertical space
set(gcf, 'Position', [100, 100, 800, 600]); % [left, bottom, width, height]

dim = 2;
StimvelPlot(Master_CSM, dim,CSMcolor, 1);
StimvelPlot(Master_CTRL, dim, 'k', 0);

% Plot reference line
yline(0, '--k');

% Axis limits and labels
xlim([-0.2 0.2]);
%ylim([-0.5 1.5]);
%ylim([-0.05 0.27]);
%ylim([-20 15])
ylim('auto')
xlabel('Time from threshold crossing (s)');
%ylabel('Outward Velocity (cm/s)');
ylabel('Outward Acceleration');
%ylabel('Outward Position (cm)');

% Adjust the main plot's position to create space below
ax1 = gca;
ax1.Position(2) = ax1.Position(2) + 0.1;
ax1.Position(4) = ax1.Position(4) - 0.1;

% Add Z-score color bar below
ax2 = axes('Position', [ax1.Position(1), ax1.Position(2) - 0.2, ax1.Position(3), 0.05]);

% Plot Z-scores as a heatmap
imagesc(zarray0', 'Parent', ax2);  % Transpose to match time x trial layout
colormap(ax2, brewermap([], 'PrGn'));  % Diverging colormap
caxis([-3, 3]);  % Adjust the Z-score color limits
colorbar(ax2, 'Position', [0.92, ax2.Position(2), 0.02, ax2.Position(4)]);

% Hide ticks and labels for the Z-score axis
set(ax2, 'XTick', [], 'YTick', [], 'YColor', 'none', 'XColor', 'none');

%% ===== Master-level stats: comparing only Master_CSM.Average vs Master_CTRL.Average =====

% Extract the averaged traces (time × variable)
a4b = Master_CSM.Average(:,5);     % CSM master trace
b4b = Master_CTRL.Average(:,5);    % Control master trace

% --- Wilcoxon signed-rank test across time points ---
[p_value, ~, pstat] = signrank(a4b, b4b, 'method', 'approximate');

% Store p-value and Z-score
parray0 = p_value;
zarray0 = pstat.zval;

% --- Compute per-timepoint Z-score of CSM relative to control ---
mu_ctrl   = mean(b4b);        % single mean
sigma_ctrl = std(b4b);        % single std

Zscores0 = (a4b - mu_ctrl) ./ sigma_ctrl;    % vector, length = time

%% ===== Plotting =====

set(gcf, 'Position', [100, 100, 800, 600]);

dim = 5;
StimvelPlot(Master_CSM, dim, CSMcolor, 1);
StimvelPlot(Master_CTRL, dim, 'k', 0);

yline(0, '--k');

xlim([-0.2 0.2]);
ylim('auto');
xlabel('Time from threshold crossing (s)');
ylabel('Outward Acceleration');

% Move main axis up for z-score bar
ax1 = gca;
ax1.Position(2) = ax1.Position(2) + 0.1;
ax1.Position(4) = ax1.Position(4) - 0.1;

% ----- Z-score heatmap (single row) -----
ax2 = axes('Position', [ax1.Position(1), ax1.Position(2) - 0.2, ...
                        ax1.Position(3), 0.05]);

imagesc(Zscores0', 'Parent', ax2);   % 1×time row
colormap(ax2, brewermap([], 'PrGn'));
caxis([-3 3]);

colorbar(ax2, 'Position', [0.92, ax2.Position(2), 0.02, ax2.Position(4)]);

set(ax2, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');
%% T0
% Initialize arrays for p-values and Z-scores
parray0 = NaN(size(T0_bigCSM.stimvel, 1), 1);  % Store p-values
zarray0 =   NaN(size(T0_bigCSM.stimvel, 1), 1)   ;    % Store Z-scores (same shape as stimvel)
Zscores0= NaN(size(T0_bigCSM.stimvel, 1), 5);

for pnum = 1:length(parray0)
    % Extract trial data for C_big (control) and T0_big (experimental)
    a4b = squeeze(T0_bigCSM.stimvel(pnum, 2, :));  % Experimental data (T0_big)
    b4b = squeeze(CT0_bigCSM.stimvel(pnum, 2, :));   % Control data (C_big), matching trial

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

T0color=rgb('Thistle');
T1color=rgb('MediumOrchid');
T2color=rgb('SlateBlue');
T3color=rgb('DarkSlateBlue');

%zscore plot 
% Adjust the figure to provide more vertical space
set(gcf, 'Position', [100, 100, 800, 600]); % [left, bottom, width, height]

dim = 2;
% StimvelPlot(C_big, dim, 'k', 1);
% StimvelPlot(T0_big, dim, T0color, 0);
StimvelPlot(CT0_bigCSM, dim, 'k', 1);
StimvelPlot(T0_bigCSM, dim, T0color, 0);

% Plot reference line
yline(0, '--k');

% Axis limits and labels
xlim([-0.2 0.2]);
%ylim([-0.5 1.5]);
%ylim([-0.05 0.27]);
%ylim([-20 15])
ylim('auto')
xlabel('Time from threshold crossing (s)');
%ylabel('Outward Velocity (cm/s)');
ylabel('Outward Acceleration');
%ylabel('Outward Position (cm)');

% Adjust the main plot's position to create space below
ax1 = gca;
ax1.Position(2) = ax1.Position(2) + 0.1;
ax1.Position(4) = ax1.Position(4) - 0.1;

% Add Z-score color bar below
ax2 = axes('Position', [ax1.Position(1), ax1.Position(2) - 0.2, ax1.Position(3), 0.05]);

% Plot Z-scores as a heatmap
imagesc(zarray0', 'Parent', ax2);  % Transpose to match time x trial layout
colormap(ax2, brewermap([], 'PrGn'));  % Diverging colormap
caxis([-3, 3]);  % Adjust the Z-score color limits
colorbar(ax2, 'Position', [0.92, ax2.Position(2), 0.02, ax2.Position(4)]);

% Hide ticks and labels for the Z-score axis
set(ax2, 'XTick', [], 'YTick', [], 'YColor', 'none', 'XColor', 'none');



%% T1

% Initialize arrays for p-values and Z-scores
parray1 = NaN(size(T1_bigCSM.stimvel, 1), 1);  % Store p-values
zarray1 =   NaN(size(T1_bigCSM.stimvel, 1), 1)   ;    % Store Z-scores (same shape as stimvel)

for pnum = 1:length(parray1)
    % Extract trial data for C_big (control) and T0_big (experimental)
    a4b = squeeze(T1_bigCSM.stimvel(pnum, 2, :));  % Experimental data (T0_big)
    b4b = squeeze(CT1_bigCSM.stimvel(pnum, 2, :));   % Control data (C_big), matching trial

    % Calculate p-value using Wilcoxon signed-rank test (trial-wise)
    [parray1(pnum, 1), ~, pstat] = signrank(a4b, b4b, 'method', 'approximate');
    parray1(pnum, 2) = pstat.zval;  % Store Z-value from signrank

    % Calculate mean and standard deviation for the control trial (b4b)
    mu_b4b = mean(b4b);             % Mean of the current trial in C_big
    sigma_b4b = std(b4b);           % Standard deviation of the current trial in C_big

    % Calculate Z-scores for each time point in a4b relative to b4b stats
    zscores = (a4b - mu_b4b) ./ sigma_b4b;  % Point-wise Z-score
    mean_zscore=mean(zscores);
    zarray1(pnum,:)=mean_zscore;
end

%%

%zscore plot 
% Adjust the figure to provide more vertical space
set(gcf, 'Position', [100, 100, 800, 600]); % [left, bottom, width, height]

dim = 2;
StimvelPlot(CT1_bigCSM, dim, 'k', 1);
StimvelPlot(T1_bigCSM, dim, T1color, 0);

% Plot reference line
yline(0, '--k');

% Axis limits and labels
%xlim([-0.5 0.5]);
xlim([-0.2 0.2])
%ylim([-0.5 1.5])
ylim('auto')
%ylim([-20 20]);
%ylim([-0.05 0.27])
xlabel('Time from threshold crossing (s)');
ylabel('Outward Acceleration(cm/s)');
%ylabel('Outward Position(cm)');
%ylabel('Outward Velocity')

% Adjust the main plot's position to create space below
ax1 = gca;
ax1.Position(2) = ax1.Position(2) + 0.1;
ax1.Position(4) = ax1.Position(4) - 0.1;

% Add Z-score color bar below
ax2 = axes('Position', [ax1.Position(1), ax1.Position(2) - 0.2, ax1.Position(3), 0.05]);

% Plot Z-scores as a heatmap
imagesc(zarray1', 'Parent', ax2);  % Transpose to match time x trial layout
colormap(ax2, brewermap([], 'PrGn'));  % Diverging colormap
caxis([-3, 3]);  % Adjust the Z-score color limits
colorbar(ax2, 'Position', [0.92, ax2.Position(2), 0.02, ax2.Position(4)]);

% Hide ticks and labels for the Z-score axis
set(ax2, 'XTick', [], 'YTick', [], 'YColor', 'none', 'XColor', 'none');




%% T2

% Initialize arrays for p-values and Z-scores
parray2 = NaN(size(T2_bigCSM.stimvel, 1), 1);  % Store p-values
zarray2 =   NaN(size(T2_bigCSM.stimvel, 1), 1)   ;    % Store Z-scores (same shape as stimvel)

for pnum = 1:length(parray2)
    % Extract trial data for C_big (control) and T0_big (experimental)
    a4b = squeeze(T2_bigCSM.stimvel(pnum, 2, :));  % Experimental data (T0_big)
    b4b = squeeze(CT2_bigCSM.stimvel(pnum, 2, :));   % Control data (C_big), matching trial

    % Calculate p-value using Wilcoxon signed-rank test (trial-wise)
    [parray2(pnum, 1), ~, pstat] = signrank(a4b, b4b, 'method', 'approximate');
    parray2(pnum, 2) = pstat.zval;  % Store Z-value from signrank

    % Calculate mean and standard deviation for the control trial (b4b)
    mu_b4b = mean(b4b);             % Mean of the current trial in C_big
    sigma_b4b = std(b4b);           % Standard deviation of the current trial in C_big

    % Calculate Z-scores for each time point in a4b relative to b4b stats
    zscores = (a4b - mu_b4b) ./ sigma_b4b;  % Point-wise Z-score
    mean_zscore=mean(zscores);
    zarray2(pnum,:)=mean_zscore;
end

%%

%zscore plot 
% Adjust the figure to provide more vertical space
set(gcf, 'Position', [100, 100, 800, 600]); % [left, bottom, width, height]

dim = 2;
StimvelPlot(CT2_bigCSM, dim, 'k', 1);
StimvelPlot(T2_bigCSM, dim, T2color, 0);

% Plot reference line
yline(0, '--k');

% Axis limits and labels
xlim([-0.2 0.2]);
%ylim([-20 20]);
ylim('auto')
%ylim([-0.05 0.27]);
%ylim([-0.5 1.5])
ylabel('Outward Velocity')
%xlabel('Time from Pos. Thresh. (s)');
ylabel('Outward Acceleration(cm/s)');
%ylabel('Outward Position(cm)');
%
% Adjust the main plot's position to create space below
ax1 = gca;
ax1.Position(2) = ax1.Position(2) + 0.1;
ax1.Position(4) = ax1.Position(4) - 0.1;

% Add Z-score color bar below
ax2 = axes('Position', [ax1.Position(1), ax1.Position(2) - 0.2, ax1.Position(3), 0.05]);

% Plot Z-scores as a heatmap
imagesc(zarray2', 'Parent', ax2);  % Transpose to match time x trial layout
colormap(ax2, brewermap([], 'PrGn'));  % Diverging colormap
caxis([-3, 3]);  % Adjust the Z-score color limits
colorbar(ax2, 'Position', [0.92, ax2.Position(2), 0.02, ax2.Position(4)]);

% Hide ticks and labels for the Z-score axis
set(ax2, 'XTick', [], 'YTick', [], 'YColor', 'none', 'XColor', 'none');



%% T3

% Initialize arrays for p-values and Z-scores
parray3 = NaN(size(T3_bigCSM.stimvel, 1), 1);  % Store p-values
zarray3 =   NaN(size(T3_bigCSM.stimvel, 1), 1)   ;    % Store Z-scores (same shape as stimvel)

for pnum = 1:length(parray3)
    % Extract trial data for C_big (control) and T0_big (experimental)
    a4b = squeeze(T3_bigCSM.stimvel(pnum, 2, :));  % Experimental data (T0_big)
    b4b = squeeze(CT3_bigCSM.stimvel(pnum, 2, :));   % Control data (C_big), matching trial

    % Calculate p-value using Wilcoxon signed-rank test (trial-wise)
    [parray3(pnum, 1), ~, pstat] = signrank(a4b, b4b, 'method', 'approximate');
    parray3(pnum, 2) = pstat.zval;  % Store Z-value from signrank

    % Calculate mean and standard deviation for the control trial (b4b)
    mu_b4b = mean(b4b);             % Mean of the current trial in C_big
    sigma_b4b = std(b4b);           % Standard deviation of the current trial in C_big

    % Calculate Z-scores for each time point in a4b relative to b4b stats
    zscores = (a4b - mu_b4b) ./ sigma_b4b;  % Point-wise Z-score
    mean_zscore=mean(zscores);
    zarray3(pnum,:)=mean_zscore;
end

%%

%zscore plot 
% Adjust the figure to provide more vertical space
set(gcf, 'Position', [100, 100, 800, 600]); % [left, bottom, width, height]

dim = 2;
StimvelPlot(CT3_bigCSM, dim, 'k', 1);
StimvelPlot(T3_bigCSM, dim, T3color, 0);

% Plot reference line
yline(0, '--k');

% Axis limits and labels
xlim([-0.2 0.2]);
ylim('auto')
%ylim([-20 20]);
%ylim([-0.5 1.5]);
%ylim([0 0.25])
xlabel('Time from threshold crossing (s)');
ylabel('Outward Acceleration(cm/s)');
%ylabel('Outward Position(cm)');
%ylabel('Outward Velcity')

% Adjust the main plot's position to create space below
ax1 = gca;
ax1.Position(2) = ax1.Position(2) + 0.1;
ax1.Position(4) = ax1.Position(4) - 0.1;

% Add Z-score color bar below
ax2 = axes('Position', [ax1.Position(1), ax1.Position(2) - 0.2, ax1.Position(3), 0.05]);

% Plot Z-scores as a heatmap
imagesc(zarray3', 'Parent', ax2);  % Transpose to match time x trial layout
colormap(ax2, brewermap([], 'PrGn'));  % Diverging colormap
caxis([-3, 3]);  % Adjust the Z-score color limits
colorbar(ax2, 'Position', [0.92, ax2.Position(2), 0.02, ax2.Position(4)]);

% Hide ticks and labels for the Z-score axis
set(ax2, 'XTick', [], 'YTick', [], 'YColor', 'none', 'XColor', 'none');
%% 
% Create a new figure for Z-score heatmap
figure;
set(gcf, 'Position', [200, 200, 800, 300]); % Adjust figure size [left, bottom, width, height]

% Plot Z-scores as a heatmap
imagesc(zarray0');  % Transpose to match time x trial layout
colormap(brewermap([], 'PrGn'));  % Use a diverging colormap (red-blue)
caxis([-3, 3]);  % Set color scale for Z-scores (-3 to 3)

% Add color bar with increments of 0.5
c = colorbar;
c.Ticks = -3:0.5:3;  % Set tick marks at 0.5 increments
c.Label.String = 'Z-Score';  % Label for the color bar

% Hide default axis labels and ticks for clarity
set(gca, 'XTick', [], 'YTick', []);

% Optional: Add title to the heatmap
title('Z-Score Heatmap');
%% 
% Plot all z scores
time_vector = linspace(-0.5, 0.5, length(zarray0)); %

T0color=rgb('Thistle');
T1color=rgb('MediumOrchid');
T2color=rgb('SlateBlue');
T3color=rgb('DarkSlateBlue');

figure;
plot(time_vector, zarray0,'Color',T0color,'LineWidth',3);
hold on
plot(time_vector, zarray1,'Color',T1color,'LineWidth',3);
hold on
plot(time_vector, zarray2,'Color',T2color,'LineWidth',3);
hold on
plot(time_vector, zarray3,'Color',T3color,'LineWidth',3);
hold on
yline(2,'--')
yline(-2,'--')
xlim([-0.25 0.25])
ylim([-10 10])