%% Treat all ReReaches as T0 reaches - Stim off Reaches 
% Load in Structure with TimeBins already saved 
T0_divide=[];
thresh=0.16;
an=fieldnames(ExcitRR_Ctrl);
for a=1:length(an)
    
    t0r=ExcitRR_Ctrl.(an{a}).RRT0;
    t0r=accel_toreaches_jerk(t0r);
    [ReReach0, StimRR0, MaxVel0,ReachesWithoutSustainedJerk0] = CSMEarlyThresh_jerk_edit7(t0r,thresh);% only looking for max and no threshold 
    
    t1r=ExcitRR_Ctrl.(an{a}).RRT1;
    t1r=accel_toreaches_jerk(t1r);
    [ReReach1, StimRR1, MaxVel1,ReachesWithoutSustainedJerk1] = CSMEarlyThresh_jerk_edit7(t1r,thresh); % only looking for max and no threshold 
   
    t2r=ExcitRR_Ctrl.(an{a}).RRT2;
    t2r=accel_toreaches_jerk(t2r);
    [ReReach2, StimRR2, MaxVel2,ReachesWithoutSustainedJerk2] = CSMEarlyThresh_jerk_edit7(t2r,thresh);
   
    t3r=ExcitRR_Ctrl.(an{a}).RRT3;
    t3r=accel_toreaches_jerk(t3r);
    [ReReach3, StimRR3, MaxVel3,ReachesWithoutSustainedJerk3] = CSMEarlyThresh_jerk_edit7(t3r,thresh);
   
    
    T0_divide.(an{a}).ReReach0=ReReach0;
    T0_divide.(an{a}).StimRR0=StimRR0;
    T0_divide.(an{a}).NoJerk0=ReachesWithoutSustainedJerk0;

    T0_divide.(an{a}).ReReach1=ReReach1;
    T0_divide.(an{a}).StimRR1=StimRR1;
    T0_divide.(an{a}).NoJerk1=ReachesWithoutSustainedJerk1;

    T0_divide.(an{a}).ReReach2=ReReach2;
    T0_divide.(an{a}).StimRR2=StimRR2;
    T0_divide.(an{a}).NoJerk2=ReachesWithoutSustainedJerk2;
 
    T0_divide.(an{a}).ReReach3=ReReach3;
    T0_divide.(an{a}).StimRR3=StimRR3;
    T0_divide.(an{a}).NoJerk3=ReachesWithoutSustainedJerk3;


    % Save to ExcitRR_Sham
    % Run through all ReReaches - easier than doing all reaches separated
    AllRR=ExcitRR_Ctrl.(an{a}).ReReaches;
    AllRR=accel_toreaches_jerk(AllRR);
    [ReReachAll, StimRRAll, MaxVelAll,ReachesWithoutSustainedJerkAll] = CSMEarlyThresh_jerk_edit7(AllRR,thresh);
    ExcitRR_Ctrl.(an{a}).StimRRAll_thresh=StimRRAll;
    ExcitRR_Ctrl.(an{a}).ReReachAll_thresh=ReReachAll;

end

%% Treat all ReReaches as T0 reaches - Stim off Reaches 
% Load in Structure with TimeBins already saved 
T0_divide_Stim=[];
thresh=0.16;
an=fieldnames(ExcitRR);
for a=1:length(an)
    
    t0r=ExcitRR.(an{a}).RRT0;
    if isempty(t0r)
        ReReach0=[];
        StimRR0=[];
        MaxVel0=[];
        ReachesWithoutSustainedJerk0=[];
    else
        t0r=accel_toreaches_jerk(t0r);
        [ReReach0, StimRR0, MaxVel0,ReachesWithoutSustainedJerk0] = CSMEarlyThresh_jerk_edit7(t0r,thresh);
    end

    t1r=ExcitRR.(an{a}).RRT1;
    t1r=accel_toreaches_jerk(t1r);
    [ReReach1, StimRR1, MaxVel1,ReachesWithoutSustainedJerk1] = CSMEarlyThresh_jerk_edit7(t1r,thresh); % only looking for max and no threshold 
   
    t2r=ExcitRR.(an{a}).RRT2;
    t2r=accel_toreaches_jerk(t2r);
    [ReReach2, StimRR2, MaxVel2,ReachesWithoutSustainedJerk2] = CSMEarlyThresh_jerk_edit7(t2r,thresh);

    t3r=ExcitRR.(an{a}).RRT3;
    if isempty(t3r)
        ReReach3=[];
        StimRR3=[];
        MaxVel3=[];
        ReachesWithoutSustainedJerk3=[];
    else
        t3r=accel_toreaches_jerk(t3r);
        [ReReach3, StimRR3, MaxVel3,ReachesWithoutSustainedJerk3] = CSMEarlyThresh_jerk_edit7(t3r,thresh);
    end
    
    T0_divide_Stim.(an{a}).ReReach0=ReReach0;
    T0_divide_Stim.(an{a}).StimRR0=StimRR0;
    T0_divide_Stim.(an{a}).NoJerk0=ReachesWithoutSustainedJerk0;

    T0_divide_Stim.(an{a}).ReReach1=ReReach1;
    T0_divide_Stim.(an{a}).StimRR1=StimRR1;
    T0_divide_Stim.(an{a}).NoJerk1=ReachesWithoutSustainedJerk1;

    T0_divide_Stim.(an{a}).ReReach2=ReReach2;
    T0_divide_Stim.(an{a}).StimRR2=StimRR2;
    T0_divide_Stim.(an{a}).NoJerk2=ReachesWithoutSustainedJerk2;
 
    T0_divide_Stim.(an{a}).ReReach3=ReReach3;
    T0_divide_Stim.(an{a}).StimRR3=StimRR3;
    T0_divide_Stim.(an{a}).NoJerk3=ReachesWithoutSustainedJerk3;


    % Save to ExcitRR_Sham
    % Run through all ReReaches - easier than doing all reaches separated
    AllRR=ExcitRR.(an{a}).ReReaches;
    AllRR=accel_toreaches_jerk(AllRR);
    [ReReachAll, StimRRAll, MaxVelAll,ReachesWithoutSustainedJerkAll] = CSMEarlyThresh_jerk_edit7(AllRR,thresh);
    ExcitRR.(an{a}).StimRRAll_thresh=StimRRAll;
    ExcitRR.(an{a}).ReReachAll_thresh=ReReachAll;

end


%%  Stim off Reaches 
BigRR0=[]; BigRR1=[]; BigRR2=[]; BigRR3=[];

excitam=fieldnames(T0_divide);
for m = 1:length(excitam)
    t0 = T0_divide.(excitam{m}).StimRR0(:,1);
    BigRR0=[BigRR0;t0];
    
    t1 = T0_divide.(excitam{m}).StimRR1(:,1);
    BigRR1=[BigRR1;t1];
  
    t2 = T0_divide.(excitam{m}).StimRR2(:,1);
    BigRR2=[BigRR2;t2];
   
    t3 = T0_divide.(excitam{m}).StimRR3(:,1);
    BigRR3=[BigRR3;t3];

end
%%  Stim on Reaches 
BigRR0_stim = [];BigRR1_stim = [];BigRR2_stim = [];BigRR3_stim = [];

excitam = fieldnames(T0_divide_Stim);
for m = 1:length(excitam)
    currStruct = T0_divide_Stim.(excitam{m});
    
    % Check for existence and non-emptiness of StimRR0
    if isfield(currStruct, 'StimRR0') && ~isempty(currStruct.StimRR0)
        t0 = currStruct.StimRR0(:,1);
        BigRR0_stim = [BigRR0_stim; t0];
    end
    
    % Check for existence and non-emptiness of StimRR1
    if isfield(currStruct, 'StimRR1') && ~isempty(currStruct.StimRR1)
        t1 = currStruct.StimRR1(:,1);
        BigRR1_stim = [BigRR1_stim; t1];
    end
    
    % Check for existence and non-emptiness of StimRR2
    if isfield(currStruct, 'StimRR2') && ~isempty(currStruct.StimRR2)
        t2 = currStruct.StimRR2(:,1);
        BigRR2_stim = [BigRR2_stim; t2];
    end
    
    % Check for existence and non-emptiness of StimRR3
    if isfield(currStruct, 'StimRR3') && ~isempty(currStruct.StimRR3)
        t3 = currStruct.StimRR3(:,1);
        BigRR3_stim = [BigRR3_stim; t3];
    end
end

%% plot No Thresh vss thresh TO code 

nothreshcolor=rgb('black');
threshcolor=rgb('SeaGreen');
binwidth = 0.01;

subplot(2,2,1); hold on
histogram(BigRR0,'EdgeColor','none','FaceColor',nothreshcolor,'BinWidth', binwidth)

%xlim([-0.2 0])
%ylim([0 0.1])
xline(-0.02,'--')
xline(-0.08,'--')
yticks(0:50:ceil(max(ylim)))  % <-- Add this line
xlabel('Time from Positional Threshold (s)')
ylabel('Probability')
title('RRT0 reaches ')


subplot(2,2,2); hold on
histogram(BigRR1,'EdgeColor','none','FaceColor',nothreshcolor,'BinWidth', binwidth)

%xlim([-0.2 0])
%ylim([0 0.1])
xline(-0.02,'--')
xline(-0.08,'--')
yticks(0:50:ceil(max(ylim)))  % <-- Add this line
xlabel('Time from Positional Threshold (s)')
ylabel('Probability')
title('RRT1 reaches ')


subplot(2,2,3); hold on
histogram(BigRR2,'EdgeColor','none','FaceColor',nothreshcolor,'BinWidth', binwidth)
%xlim([-0.2 0])
%ylim([0 0.1])
xline(-0.02,'--')
xline(-0.08,'--')
yticks(0:50:ceil(max(ylim)))  % <-- Add this line
xlabel('Time from Positional Threshold (s)')
ylabel('Probability')
title('RRT2 reaches ')


subplot(2,2,4); hold on
histogram(BigRR3,'EdgeColor','none','FaceColor',nothreshcolor,'BinWidth', binwidth)
%xlim([-0.2 0])
%ylim([0 0.1])
xline(-0.02,'--')
xline(-0.08,'--')
yticks(0:50:ceil(max(ylim)))  % <-- Add this line
xlabel('Time from Positional Threshold (s)')
ylabel('Probability')
title('RRT3 reaches ')

%% plot No Thresh vss thresh TO code - Stim on Reaches 

nothreshcolor=rgb('black');
threshcolor=rgb('SeaGreen');
binwidth = 0.01;

subplot(2,2,1); hold on
histogram(BigRR0_stim,'EdgeColor','none','FaceColor',nothreshcolor,'BinWidth', binwidth)
%xlim([-0.2 0])
%ylim([0 0.1])
xline(-0.02,'--')
xline(-0.08,'--')
yticks(0:50:ceil(max(ylim)))  % <-- Add this line
xlabel('Time from Positional Threshold (s)')
ylabel('Probability')
title('RRT0 reaches ')


subplot(2,2,2); hold on
histogram(BigRR1_stim,'EdgeColor','none','FaceColor',nothreshcolor,'BinWidth', binwidth)

%xlim([-0.2 0])
%ylim([0 0.1])
xline(-0.02,'--')
xline(-0.08,'--')
yticks(0:50:ceil(max(ylim)))  % <-- Add this line
xlabel('Time from Positional Threshold (s)')
ylabel('Probability')
title('RRT1 reaches ')


subplot(2,2,3); hold on
histogram(BigRR2_stim,'EdgeColor','none','FaceColor',nothreshcolor,'BinWidth', binwidth)
%xlim([-0.2 0])
%ylim([0 0.1])
xline(-0.02,'--')
xline(-0.08,'--')
yticks(0:50:ceil(max(ylim)))  % <-- Add this line
xlabel('Time from Positional Threshold (s)')
ylabel('Probability')
title('RRT2 reaches ')


subplot(2,2,4); hold on
histogram(BigRR3_stim,'EdgeColor','none','FaceColor',nothreshcolor,'BinWidth', binwidth)
%xlim([-0.2 0])
%ylim([0 0.1])
xline(-0.02,'--')
xline(-0.08,'--')
yticks(0:50:ceil(max(ylim)))  % <-- Add this line
xlabel('Time from Positional Threshold (s)')
ylabel('Probability')
title('RRT3 reaches ')
%% Now filter through , get total number within pre thresh window - STIM OFF REACHES
for a = 1:length(an)
    RRtimes1= T0_divide.(an{a}).StimRR1;
    RRtimes2= T0_divide.(an{a}).StimRR2;
    RRtimes3= T0_divide.(an{a}).StimRR3;

    % Keep only values >= -0.08 and <= -0.02
    RRtimes1_filt = RRtimes1(RRtimes1 >= -0.08 & RRtimes1 <= -0.02);
    RRtimes2_filt = RRtimes2(RRtimes2 >= -0.08 & RRtimes2 <= -0.02);
    RRtimes3_filt = RRtimes3(RRtimes3 >= -0.08 & RRtimes3 <= -0.02);

    % (Optional) Save back the filtered values, if desired:
    ExcitRR_Ctrl.(an{a}).StimRR1_filt=  RRtimes1_filt;
    ExcitRR_Ctrl.(an{a}).StimRR2_filt=  RRtimes2_filt;
    ExcitRR_Ctrl.(an{a}).StimRR3_filt=  RRtimes3_filt;

    AllbinsRR=ExcitRR_Ctrl.(an{a}).AllBins;

    combinedData=[];
    combinedData = [AllbinsRR; RRtimes1_filt; RRtimes2_filt;RRtimes3_filt];
    ExcitRR_Ctrl.(an{a}).AllBins_prethresh=combinedData;

end

%% Now filter through , get total number within pre thresh window - STIM ON REACHES
for a = 1:length(an)
    RRtimes1= T0_divide_Stim.(an{a}).StimRR1;
    RRtimes2= T0_divide_Stim.(an{a}).StimRR2;
    RRtimes3= T0_divide_Stim.(an{a}).StimRR3;

    % Keep only values >= -0.08 and <= -0.02
    RRtimes1_filt = RRtimes1(RRtimes1 >= -0.08 & RRtimes1 <= -0.02);
    RRtimes2_filt = RRtimes2(RRtimes2 >= -0.08 & RRtimes2 <= -0.02);
    RRtimes3_filt = RRtimes3(RRtimes3 >= -0.08 & RRtimes3 <= -0.02);

    % (Optional) Save back the filtered values, if desired:
    ExcitRR.(an{a}).StimRR1_filt=  RRtimes1_filt;
    ExcitRR.(an{a}).StimRR2_filt=  RRtimes2_filt;
    ExcitRR.(an{a}).StimRR3_filt=  RRtimes3_filt;

    AllbinsRR=ExcitRR.(an{a}).AllBins;

    combinedData=[];
    combinedData = [AllbinsRR; RRtimes1_filt; RRtimes2_filt;RRtimes3_filt];
    ExcitRR.(an{a}).AllBins_prethresh=combinedData;

end
%% Plot All Groups Together with KS Test

% Define colors
stimcolor = rgb('CadetBlue'); % ExcitRR
ctrlcolor = rgb('LightSlateGray'); % ExcitRR_Ctrl
shamcolor = rgb('Black'); % ExcitRR_Sham (Dashed line)

% Initialize storage for KDEs
num_points = 100; % Number of points for KDE interpolation
xi_common = linspace(-0.1, 0.2, num_points); % Common x-axis

% -----------------------------
% Compute KDE for ExcitRR Group
% -----------------------------
all_f_excitrr = [];
excitam_excitrr = fieldnames(ExcitRR);
num_excitrr = length(excitam_excitrr);
all_data_excitrr = [];

for m = 1:num_excitrr
    allbins_excitrr = ExcitRR.(excitam_excitrr{m}).AllBins_prethresh;
    all_data_excitrr = [all_data_excitrr; allbins_excitrr];  % Collect for KS test
    [f_excitrr, xi_excitrr] = ksdensity(allbins_excitrr, 'Bandwidth', 0.008);
    f_interp_excitrr = interp1(xi_excitrr, f_excitrr, xi_common, 'linear', 0);
    all_f_excitrr = [all_f_excitrr; f_interp_excitrr];
end

% -----------------------------
% Compute KDE for ExcitRR_Ctrl Group
% -----------------------------
all_f_excitrr_ctrl = [];
excitam_excitrr_ctrl = fieldnames(ExcitRR_Ctrl);
num_excitrr_ctrl = length(excitam_excitrr_ctrl);
all_data_excitrr_ctrl = [];

for m = 1:num_excitrr_ctrl
    allbins_excitrr_ctrl = ExcitRR_Ctrl.(excitam_excitrr_ctrl{m}).AllBins_prethresh;
    all_data_excitrr_ctrl = [all_data_excitrr_ctrl; allbins_excitrr_ctrl];  % Collect for KS test
    [f_excitrr_ctrl, xi_excitrr_ctrl] = ksdensity(allbins_excitrr_ctrl, 'Bandwidth', 0.008);
    f_interp_excitrr_ctrl = interp1(xi_excitrr_ctrl, f_excitrr_ctrl, xi_common, 'linear', 0);
    all_f_excitrr_ctrl = [all_f_excitrr_ctrl; f_interp_excitrr_ctrl];
end

% -----------------------------
% Compute KDE for ExcitRR_Sham Group (Dashed Line Only, No Fill)
% -----------------------------
all_f_sham = [];
excitam_sham = fieldnames(ExcitRR_Sham);
num_sham = length(excitam_sham);
all_data_sham = [];

for m = 1:num_sham
    allbins_sham = ExcitRR_Sham.(excitam_sham{m}).AllBins_prethresh;
    all_data_sham = [all_data_sham; allbins_sham];  % Collect for KS test
    [f_sham, xi_sham] = ksdensity(allbins_sham, 'Bandwidth', 0.008);
    f_interp_sham = interp1(xi_sham, f_sham, xi_common, 'linear', 0);
    all_f_sham = [all_f_sham; f_interp_sham];
end

%% Run KS Test for comparisons
[h_excitrr_ctrl, p_excitrr_ctrl, ks_stat_excitrr_ctrl] = kstest2(all_data_excitrr, all_data_excitrr_ctrl);
[h_excitrr_sham, p_excitrr_sham, ks_stat_excitrr_sham] = kstest2(all_data_excitrr, all_data_sham);
[h_excitrr_ctrl_sham, p_excitrr_ctrl_sham, ks_stat_excitrr_ctrl_sham] = kstest2(all_data_excitrr_ctrl, all_data_sham);

% Display KS Test results in command window
fprintf('KS test between ExcitRR and ExcitRR_Ctrl:\n');
fprintf('  h = %d (1 = reject null hypothesis)\n', h_excitrr_ctrl);
fprintf('  p = %.4f\n', p_excitrr_ctrl);
fprintf('  KS statistic = %.4f\n\n', ks_stat_excitrr_ctrl);

fprintf('KS test between ExcitRR and ExcitRR_Sham:\n');
fprintf('  h = %d (1 = reject null hypothesis)\n', h_excitrr_sham);
fprintf('  p = %.4f\n', p_excitrr_sham);
fprintf('  KS statistic = %.4f\n\n', ks_stat_excitrr_sham);

fprintf('KS test between ExcitRR_Ctrl and ExcitRR_Sham:\n');
fprintf('  h = %d (1 = reject null hypothesis)\n', h_excitrr_ctrl_sham);
fprintf('  p = %.4f\n', p_excitrr_ctrl_sham);
fprintf('  KS statistic = %.4f\n\n', ks_stat_excitrr_ctrl_sham);

%% Plot All Groups Together
figure; hold on;

% Fill areas for each group
fill([xi_common, fliplr(xi_common)], [mean(all_f_excitrr, 1), fliplr(zeros(1, num_points))], ...
    stimcolor, 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % ExcitRR (Teal)

fill([xi_common, fliplr(xi_common)], [mean(all_f_excitrr_ctrl, 1), fliplr(zeros(1, num_points))], ...
    ctrlcolor, 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % ExcitRR_Ctrl (LightSlateGray)

% Plot mean KDE lines
plot(xi_common, mean(all_f_excitrr, 1), 'Color', stimcolor, 'LineWidth', 3); % ExcitRR (Teal)
plot(xi_common, mean(all_f_excitrr_ctrl, 1), 'Color', ctrlcolor, 'LineWidth', 3); % ExcitRR_Ctrl (LightSlateGray)

% Plot Sham group with dashed line, no fill
plot(xi_common, mean(all_f_sham, 1), '--', 'Color', shamcolor, 'LineWidth', 3); % Sham (Blue, Dashed)

% Set axis limits
xlim([-0.1 0.2]);

% Customize x-ticks
xticks([-0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2]);

% Labels and title
xlabel('Time from threshold crossing (s)');
ylabel('Density');
title('Average KDE across Animals (ExcitRR, ExcitRR_Ctrl, ExcitRR_Sham)');

% Annotate KS test results on plot
ymax = max([mean(all_f_excitrr, 1), mean(all_f_excitrr_ctrl, 1), mean(all_f_sham, 1)], [], 'all');
ks_text_1 = sprintf('KS (ExcitRR vs ExcitRR_Ctrl): p = %.4f, D = %.3f', p_excitrr_ctrl, ks_stat_excitrr_ctrl);
ks_text_2 = sprintf('KS (ExcitRR vs ExcitRR_Sham): p = %.4f, D = %.3f', p_excitrr_sham, ks_stat_excitrr_sham);
ks_text_3 = sprintf('KS (ExcitRR_Ctrl vs Sham): p = %.4f, D = %.3f', p_excitrr_ctrl_sham, ks_stat_excitrr_ctrl_sham);

% Display the KS test text on the plot
text(0.05, ymax * 0.9, ks_text_1, 'FontSize', 12, 'FontWeight', 'bold');
text(0.05, ymax * 0.85, ks_text_2, 'FontSize', 12, 'FontWeight', 'bold');
text(0.05, ymax * 0.8, ks_text_3, 'FontSize', 12, 'FontWeight', 'bold');

% Formatting
set(gca, 'FontSize', 12, 'LineWidth', 1);

hold off;
