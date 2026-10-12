%% Baseline versus Sham 
% % BSTIMOFF
an=fieldnames(ExcitRR_Ctrl);
thresh=0.1;
dirdir=1;
for a=1:length(an)
    % % BASELINE
    CtrlR= ExcitRR_Ctrl.(an{a}).CtrlReaches;
    StimR= ExcitRR.(an{a}).ReReaches;
    
     dim=7; %max outward velocity
     dir=1;

     %C= AlignPeakV(CtrlR,dim,1000);
     C=AlignThresh(CtrlR,thresh,dirdir,200);
     C=stimvel_accel(C);

     % S= AlignPeakV(StimR,dim,1000);  
    S=AlignStim(StimR,200)
     S=stimvel_accel(S);

    
    % ExcitRR_Ctrl.(an{a}).C_peakV=C;
    % ExcitRR.(an{a}).S_peakV=S;

    ExcitRR_Ctrl.(an{a}).C=C;
    ExcitRR.(an{a}).S=S;

end
%%
% % SHAM
an=fieldnames(ExcitRR_Sham);
thresh=0.1;
dirdir=1;
for a=1:length(an)
    CtrlR= ExcitRR_Sham.(an{a}).CtrlReaches;
   
     dim=7; %max outward velocity
     dir=1;
     %C= AlignPeakV(CtrlR,dim,1000);
     C=AlignThresh(CtrlR,thresh,dirdir,200)
     C=stimvel_accel(C);
   
    %ExcitRR_Sham.(an{a}).C_peakV=C;
    ExcitRR_Sham.(an{a}).C=C;

end

%% Align Stimvel Zscore 
AD_big = struct;
C_big = struct;

an=fieldnames(ExcitRR_Ctrl);

for a = 1:length(an)
   
    % AD_big.stimvel(:,:,a) = ExcitRR.(an{a}).S_peakV.Average;
    % C_big.stimvel(:,:,a) = ExcitRR_Ctrl.(an{a}).C_peakV.Average;

    AD_big.stimvel(:,:,a) = ExcitRR.(an{a}).S.Average;
    C_big.stimvel(:,:,a) = ExcitRR_Ctrl.(an{a}).C.Average;
end


C_big.Average = nanmean(C_big.stimvel,3);
C_big.Std = nanstd(C_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(C_big.Average),1); 
for nnum = 1:length(C_big.Average)
    samp = squeeze(C_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
C_big.N = N;

C_big=stimvel_accel(C_big);


AD_big.Average = nanmean(AD_big.stimvel,3);
AD_big.Std = nanstd(AD_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(AD_big.Average),1); 
for nnum = 1:length(AD_big.Average)
    samp = squeeze(AD_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
AD_big.N = N;

AD_big=stimvel_accel(AD_big);


%% Align Stimvel Zscore - SHAM CONTROLs 
Sham_big = struct;


an=fieldnames(ExcitRR_Sham);

for a = 1:length(an)
   
    %Sham_big.stimvel(:,:,a) = ExcitRR_Sham.(an{a}).C_peakV.Average;
    Sham_big.stimvel(:,:,a) = ExcitRR_Sham.(an{a}).C.Average;
  
end

%For structures aligned to CSM
Sham_big.Average = nanmean(Sham_big.stimvel,3);
Sham_big.Std = nanstd(Sham_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(Sham_big.Average),1); 
for nnum = 1:length(Sham_big.Average)
    samp = squeeze(Sham_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
Sham_big.N = N;

Sham_big=stimvel_accel(Sham_big);



%% Plot 
basecolor= rgb('Gray');
stimcolor=rgb('DodgerBlue');
shamcolor=rgb('Black');

dim = 2;

StimvelPlot(C_big,dim,basecolor,1);
StimvelPlot(AD_big,dim,stimcolor,0);
StimvelPlot(Sham_big,dim,shamcolor,0);
yline(0,'--')
ylim('auto')
xlim([-0.2 0.2])
%title('Max Valigned')
ylabel('Outward Velocity')
xlabel('Time from Mx V (s)')
xlabel('Time from Stim (s)')

%% Generate stimvel zscore
% Initialize arrays for p-values and Z-scores
parray0 = NaN(size(BL_CSM2_big.stimvel, 1), 1);  % Store p-values
zarray0 =   NaN(size(BL_CSM2_big.stimvel, 1), 1)   ;    % Store Z-scores (same shape as stimvel)
Zscores0= NaN(size(BL_CSM2_big.stimvel, 1), 3);

for pnum = 1:length(parray0)
    % Extract trial data for C_big (control) and T0_big (experimental)
    %a4b = squeeze(WOE_CSM2_big.stimvel(pnum, 8, :));  % Experimental data (T0_big)
    %b4b = squeeze(BL_CSM2_big.stimvel(pnum,8, :));   % Control data (C_big), matching trial
    % 
    %  a4b = squeeze(WOE_CSM3_big.stimvel(pnum, 5, :));  % Experimental data (T0_big)
    % b4b = squeeze(BL_CSM3_big.stimvel(pnum,5, :));   % Control data (C_big), matching trial


    a4b = squeeze(AD_CSM2_big.stimvel(pnum, 8, :));  % Experimental data (T0_big)
    b4b = squeeze(BL_CSM2_big.stimvel(pnum,8, :));   % Control data (C_big), matching trial

    % Calculate p-value using Wilcoxon signed-rank test (trial-wise)
    [parray0(pnum, 1), ~, pstat] = signrank(a4b, b4b, 'method', 'approximate');
    parray0(pnum, 2) = pstat.zval;  % Store Z-value from signrank

    % Calculate mean and standard deviation for the control trial (b4b)
    mu_b4b = mean(b4b);             % Mean of the current trial in C_big
    sigma_b4b = std(b4b);           % Standard deviation of the current trial in C_big

    % Calculate Z-scores for each time point in a4b relative to b4b stats
    zscores = (a4b - mu_b4b) ./ sigma_b4b;  % Point-wise Z-score
    Zscores0(pnum,:)=zscores;
    mean_zscore=mean(zscores);
    zarray0(pnum,:)=mean_zscore;
end

%%
BLcolor=rgb('Gray');
CSMcolor=rgb('Orchid');
%zscore plot 
% Adjust the figure to provide more vertical space
set(gcf, 'Position', [100, 100, 800, 600]); % [left, bottom, width, height]

dim = 8;
StimvelPlot(BL_CSM2_big, dim,BLcolor, 1);
% StimvelPlot(WOE_CSM2_big, dim,CSMcolor, 0);
StimvelPlot(AD_CSM2_big,dim,stimcolor,0);

% StimvelPlot(BL_CSM3_big, dim,BLcolor, 1);
% StimvelPlot(WOE_CSM3_big, dim,CSMcolor, 0);
% %StimvelPlot(AD_CSM3_big,dim,stimcolor,0);

% Plot reference line
yline(0, '--k');

% Axis limits and labels
ylim('auto')
%xlim([-0.2 0.2])
%ylim([-20 20])
xlabel('Time from CSM T2 (s)');

ylabel('Outward Acceleration');


% Adjust the main plot's position to create space below
ax1 = gca;
ax1.Position(2) = ax1.Position(2) + 0.1;
ax1.Position(4) = ax1.Position(4) - 0.1;

% Add Z-score color bar below
ax2 = axes('Position', [ax1.Position(1), ax1.Position(2) - 0.2, ax1.Position(3), 0.05]);

% Plot Z-scores as a heatmap
imagesc(zarray0', 'Parent', ax2);  % Transpose to match time x trial layout
colormap(ax2, brewermap([], 'PrGn'));  % Diverging colormap
caxis([-2, 2]);  % Adjust the Z-score color limits
colorbar(ax2, 'Position', [0.92, ax2.Position(2), 0.02, ax2.Position(4)]);

% Hide ticks and labels for the Z-score axis
set(ax2, 'XTick', [], 'YTick', [], 'YColor', 'none', 'XColor', 'none');


%%
CSMcolor = rgb('Orchid');

% zscore plot 
set(gcf, 'Position', [100, 100, 800, 600]); % [left, bottom, width, height]

dim = 8;

% Plot both datasets on the same axes
StimvelPlot(BL_CSM2_big, dim,BLcolor, 1);
hold on
StimvelPlot(WOE_CSM2_big, dim,CSMcolor, 0);
hold on
StimvelPlot(AD_CSM2_big,dim,stimcolor,0);
hold off;

% Plot reference line
yline(0, '--k');
xlim([-0.2 0.2])
ylim('auto')
% Axis limits and labels
ylim('auto')
xlabel('Time from MaxV (s)');
ylabel('Outward Acceleration');

% Adjust the main plot's position to make it shorter
ax1 = gca;
ax1.Position(2) = ax1.Position(2) + 0.15;  % Raise the plot up
ax1.Position(4) = ax1.Position(4) * 0.6;   % Shrink the height (e.g., 60% of original)

% Add Z-score color bar below
ax2 = axes('Position', [ax1.Position(1), ax1.Position(2) - 0.18, ax1.Position(3), 0.05]);

% Plot Z-scores as a heatmap
imagesc(zarray0', 'Parent', ax2);  % Transpose to match time x trial layout
colormap(ax2, brewermap([], 'PrGn'));  % Diverging colormap
caxis([-2, 2]);  % Adjust the Z-score color limits
colorbar(ax2, 'Position', [0.92, ax2.Position(2), 0.02, ax2.Position(4)]);

% Hide ticks and labels for the Z-score axis
set(ax2, 'XTick', [], 'YTick', [], 'YColor', 'none', 'XColor', 'none');


%% Generate stimvel zscore : BL vs SHAM


% Initialize arrays to store Z-scores
Zscores0 = NaN(size(BL_CSM2_big.stimvel, 1), size(BL_CSM2_big.stimvel, 3)); % [nUnits x nTrials]
zarray0 = NaN(size(BL_CSM2_big.stimvel, 1), 1);  % [nUnits x 1] mean z-score at T=8

for pnum = 1:size(BL_CSM2_big.stimvel, 1)
    % Extract timepoint 8 across trials
    sham_tp8 = squeeze(Sham_CSM2_big.stimvel(pnum, 5, :));  % SHAM reference
    bl_tp8 = squeeze(BL_CSM2_big.stimvel(pnum, 5, :));      % BL data to z-score

    % Compute SHAM reference stats
    mu_sham = mean(sham_tp8, 'omitnan');
    sigma_sham = std(sham_tp8, 'omitnan');

    % Compute z-scores of BL trials relative to SHAM distribution
    zscores = (bl_tp8 - mu_sham) ./ sigma_sham;

    % Store
    Zscores0(pnum, :) = zscores;
    zarray0(pnum) = mean(zscores, 'omitnan');
end

%% Plotting
BLcolor = rgb('Gray');
CSMcolor = rgb('Orchid');
shamcolor=rgb('Black');

% Adjust the figure to provide more vertical space
set(gcf, 'Position', [100, 100, 800, 600]); % [left, bottom, width, height]

dim = 5;
StimvelPlot(BL_CSM2_big, dim, BLcolor, 1);
StimvelPlot(Sham_CSM2_big, dim, shamcolor, 0);
StimvelPlot(AD_CSM2_big,dim,stimcolor,0);

% Plot reference line
yline(0, '--k');

% Axis limits and labels
ylim('auto')
xlabel('Time from MaxV (s)');
ylabel('Outward Acceleration');

% Adjust main plot position to create space below
ax1 = gca;
ax1.Position(2) = ax1.Position(2) + 0.1;
ax1.Position(4) = ax1.Position(4) - 0.1;

% Add Z-score heatmap below
ax2 = axes('Position', [ax1.Position(1), ax1.Position(2) - 0.2, ax1.Position(3), 0.05]);
imagesc(zarray0', 'Parent', ax2);  % Transpose to match expected orientation
colormap(ax2, brewermap([], 'PrGn'));  % Diverging colormap
caxis([-2, 2]);  % Adjust the Z-score color limits
colorbar(ax2, 'Position', [0.92, ax2.Position(2), 0.02, ax2.Position(4)]);

% Hide ticks and labels for the Z-score axis
set(ax2, 'XTick', [], 'YTick', [], 'YColor', 'none', 'XColor', 'none');
