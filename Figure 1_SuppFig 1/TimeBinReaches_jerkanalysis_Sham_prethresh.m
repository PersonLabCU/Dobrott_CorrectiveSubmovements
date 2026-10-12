%% Treat all ReReaches as T0 reaches 
% Load in Structure with TimeBins already saved 
T0_divide=[];
thresh=0.1; % change depending on the type of data 
ExcitRR_Sham=ExcitRR_Ctrl;
an=fieldnames(ExcitRR_Sham);
for a=1:length(an)
    
    t0r=ExcitRR_Sham.(an{a}).RRT0;
    t0r=accel_toreaches_jerk(t0r);
    [ReReach0, StimRR0, MaxVel0,ReachesWithoutSustainedJerk0] = CSMEarlyThresh_jerk_edit7(t0r,thresh);% only looking for max and no threshold 
    
    t1r=ExcitRR_Sham.(an{a}).RRT1;
    t1r=accel_toreaches_jerk(t1r);
    [ReReach1, StimRR1, MaxVel1,ReachesWithoutSustainedJerk1] = CSMEarlyThresh_jerk_edit7(t1r,thresh); % only looking for max and no threshold 
   
    t2r=ExcitRR_Sham.(an{a}).RRT2;
    t2r=accel_toreaches_jerk(t2r);
    [ReReach2, StimRR2, MaxVel2,ReachesWithoutSustainedJerk2] = CSMEarlyThresh_jerk_edit7(t2r,thresh);
   
    t3r=ExcitRR_Sham.(an{a}).RRT3;
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
    AllRR=ExcitRR_Sham.(an{a}).ReReaches;
    AllRR=accel_toreaches_jerk(AllRR);
    [ReReachAll, StimRRAll, MaxVelAll,ReachesWithoutSustainedJerkAll] = CSMEarlyThresh_jerk_edit7(AllRR,thresh);
    ExcitRR_Sham.(an{a}).StimRRAll_thresh=StimRRAll;
    ExcitRR_Sham.(an{a}).ReReachAll_thresh=ReReachAll;

end


%% 
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
%% plot No Thresh vss thresh TO code 

nothreshcolor=rgb('black');
threshcolor=rgb('SeaGreen');
binwidth = 0.01;

subplot(2,2,1); hold on
% histogram(BigRR0,'Normalization','probability','EdgeColor','none','FaceColor',nothreshcolor,'BinWidth', binwidth)
% hold on
% histogram(BigRR0_thresh,'Normalization','probability','EdgeColor','none','FaceColor',threshcolor,'BinWidth', binwidth)
histogram(BigRR0,'EdgeColor','none','FaceColor',nothreshcolor,'BinWidth', binwidth)

%xlim([-0.2 0])
%ylim([0 0.1])
xline(-0.01,'--')
xline(-0.05,'--')
yticks(0:50:ceil(max(ylim)))  % <-- Add this line
xlabel('Time from Positional Threshold (s)')
ylabel('Probability')
title('RRT0 reaches ')


subplot(2,2,2); hold on
% histogram(BigRR1,'Normalization','probability','EdgeColor','none','FaceColor',nothreshcolor,'BinWidth', binwidth)
% hold on
% histogram(BigRR1_thresh,'Normalization','probability','EdgeColor','none','FaceColor',threshcolor,'BinWidth', binwidth)
histogram(BigRR1,'EdgeColor','none','FaceColor',nothreshcolor,'BinWidth', binwidth)

%xlim([-0.2 0])
%ylim([0 0.1])
xline(-0.01,'--')
xline(-0.05,'--')
yticks(0:50:ceil(max(ylim)))  % <-- Add this line
xlabel('Time from Positional Threshold (s)')
ylabel('Probability')
title('RRT1 reaches ')


subplot(2,2,3); hold on
% histogram(BigRR2,'Normalization','probability','EdgeColor','none','FaceColor',nothreshcolor,'BinWidth', binwidth)
% hold on
% histogram(BigRR2_thresh,'Normalization','probability','EdgeColor','none','FaceColor',threshcolor,'BinWidth', binwidth)
histogram(BigRR2,'EdgeColor','none','FaceColor',nothreshcolor,'BinWidth', binwidth)
%xlim([-0.2 0])
%ylim([0 0.1])
xline(-0.01,'--')
xline(-0.05,'--')
yticks(0:50:ceil(max(ylim)))  % <-- Add this line
xlabel('Time from Positional Threshold (s)')
ylabel('Probability')
title('RRT2 reaches ')


subplot(2,2,4); hold on
% histogram(BigRR3,'Normalization','probability','EdgeColor','none','FaceColor',nothreshcolor,'BinWidth', binwidth)
% hold on
% histogram(BigRR3_thresh,'Normalization','probability','EdgeColor','none','FaceColor',threshcolor,'BinWidth', binwidth)
histogram(BigRR3,'EdgeColor','none','FaceColor',nothreshcolor,'BinWidth', binwidth)
%xlim([-0.2 0])
%ylim([0 0.1])
xline(-0.01,'--')
xline(-0.05,'--')
yticks(0:50:ceil(max(ylim)))  % <-- Add this line
xlabel('Time from Positional Threshold (s)')
ylabel('Probability')
title('RRT3 reaches ')
%% Now filter through , get total number within pre thresh window 
for a = 1:length(an)
    RRtimes1= T0_divide.(an{a}).StimRR1;
    RRtimes2= T0_divide.(an{a}).StimRR2;
    RRtimes3= T0_divide.(an{a}).StimRR3;

    % Keep only values >= -0.08 and <= -0.02
    RRtimes1_filt = RRtimes1(RRtimes1 >= -0.05 & RRtimes1 <= -0.01);
    RRtimes2_filt = RRtimes2(RRtimes2 >= -0.05 & RRtimes2 <= -0.01);
    RRtimes3_filt = RRtimes3(RRtimes3 >= -0.05 & RRtimes3 <= -0.01);

    % (Optional) Save back the filtered values, if desired:
    ExcitRR_Sham.(an{a}).StimRR1_filt=  RRtimes1_filt;
    ExcitRR_Sham.(an{a}).StimRR2_filt=  RRtimes2_filt;
    ExcitRR_Sham.(an{a}).StimRR3_filt=  RRtimes3_filt;

    AllbinsRR=ExcitRR_Sham.(an{a}).AllBins;

    combinedData=[];
    combinedData = [AllbinsRR; RRtimes1_filt; RRtimes2_filt;RRtimes3_filt];
    ExcitRR_Sham.(an{a}).AllBins_prethresh=combinedData;

end
%% plot KDE 
% -----------------------------
% Compute KDE for ExcitRR_Sham Group (Dashed Line Only, No Fill)
% -----------------------------
num_points = 100; % Number of points for KDE interpolation
xi_common = linspace(-0.1, 0.2, num_points); % Common x-axis for KDEs

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
%% Plot the KDEs with Full Area Under the Curve Filled
% -----------------------------
stimcolor = rgb('Gray');
figure; hold on;


% Fill area for ExcitRR_Sham group (blue)
fill([xi_common, fliplr(xi_common)], [mean(all_f_sham, 1), fliplr(zeros(1, num_points))], ...
    [0.5 0.5 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Plot mean KDE lines
plot(xi_common, mean(all_f_sham, 1), 'k', 'LineWidth', 3); % T0_divide

% Axis limits and labels
xlim([-0.1 0.2]);
xticks([-0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2]);
xlabel('Time from threshold crossing (s)');
ylabel('Density');
title('Average KDE across Animals (T0_divide vs ExcitRR_Sham)');

% Formatting
set(gca, 'FontSize', 12, 'LineWidth', 1);

% -----------------------------
% Add KS test results to the plot
% -----------------------------
ymax = max([mean(all_f_sham, 1), mean(all_f_sham, 1)], [], 'all');
%%  Find peaks and location of peaks 
% Compute mean KDE
mean_kde = mean(all_f_sham, 1);

% Find peaks in the KDE
[peakVals, peakLocs] = findpeaks(mean_kde, xi_common, ...
    'MinPeakProminence', 0.01); % Adjust this threshold as needed
% 
% % Plot detected peaks
plot(peakLocs, peakVals, 'ro', 'MarkerFaceColor', 'r');

% Compute distances between peaks (in seconds)
if length(peakLocs) >= 2
    peakDistances = diff(peakLocs);
    fprintf('Distances between peaks (s):\n');
    disp(peakDistances);

    % If you want an approximate frequency:
    mean_peak_distance = mean(peakDistances);
    fprintf('Mean distance between peaks: %.4f s\n', mean_peak_distance);
    fprintf('Estimated frequency: %.2f Hz\n', 1 / mean_peak_distance);
else
    fprintf('Not enough peaks detected to compute frequency.\n');
end


%% Find peaks and location of peaks 
% Compute mean KDE
mean_kde = mean(all_f_sham, 1);
%%
% Generate a 13 Hz sine wave
t = xi_common;  % x-axis: time (already used in your KDE)
amp = 0.05;     % adjust amplitude to match KDE y-axis scale
sine_wave = amp * sin(2 * pi * 13.85 * t);  % 13 Hz sine wave

% Offset the sine wave vertically to align with your KDE
% Example: offset it just above zero baseline, or near your curve
y_offset = 0.02;  % adjust based on your KDE height
sine_wave_shifted = sine_wave + y_offset;


% Find peaks in the KDE
[peakVals, peakLocs] = findpeaks(mean_kde, xi_common, ...
    'MinPeakProminence', 0.01); % Adjust this threshold as needed

% Plot detected peaks (on existing figure)
plot(peakLocs, peakVals, 'ro', 'MarkerFaceColor', 'r');
hold on

% Fill area for ExcitRR_Sham group (blue)
fill([xi_common, fliplr(xi_common)], [mean(all_f_sham, 1), fliplr(zeros(1, num_points))], ...
    [0.5 0.5 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Plot mean KDE lines
plot(xi_common, mean(all_f_sham, 1), 'k', 'LineWidth', 3); % T0_divide
hold on

% Plot the 13 Hz sine wave
plot(t, sine_wave_shifted, 'b--', 'LineWidth', 1.5);
hold on
% Label each peak with its time (in seconds)
for i = 1:length(peakLocs)
    text(peakLocs(i), peakVals(i) + 0.005, ...
        sprintf('%.3fs', peakLocs(i)), ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 10, 'Color', 'r');
end

% Compute distances between peaks (in seconds)
if length(peakLocs) >= 2
    peakDistances = diff(peakLocs);
    fprintf('Distances between peaks (s):\n');
    disp(peakDistances);

    % Approximate frequency
    mean_peak_distance = mean(peakDistances);
    fprintf('Mean distance between peaks: %.4f s\n', mean_peak_distance);
    fprintf('Estimated frequency: %.2f Hz\n', 1 / mean_peak_distance);
else
    fprintf('Not enough peaks detected to compute frequency.\n');
end
% Axis limits and labels
xlim([-0.1 0.2]);
xticks([-0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2]);
xline(-0.0636,'--')
xline(0.01212,'--')
xline(0.0818,'--')
xline(0.153,'--')
xlabel('Time from threshold crossing (s)');
ylabel('Density');
title('Average KDE across Animals (T0_divide vs ExcitRR_Sham)');

% Formatting
set(gca, 'FontSize', 12, 'LineWidth', 1);
%%
% Generate a 13 Hz sine wave
t = xi_common;  % x-axis: time (already used in your KDE)
amp = 0.05;     % adjust amplitude to match KDE y-axis scale
sine_wave = amp * sin(2 * pi * 13 * t);  % 13 Hz sine wave

% Offset the sine wave vertically to align with your KDE
% Example: offset it just above zero baseline, or near your curve
y_offset = 0.02;  % adjust based on your KDE height
sine_wave_shifted = sine_wave + y_offset;

% Plot the 13 Hz sine wave
plot(t, sine_wave_shifted, 'b--', 'LineWidth', 1.5);

% Optional: Add legend entry
legend('KDE', 'Detected Peaks', '13 Hz Sine Wave');


%% Different way to plot 
% -----------------------------
% Compute KDE for ExcitRR_Sham Group
% -----------------------------
num_points = 100; % Number of points for KDE interpolation
xi_common = linspace(-0.1, 0.2, num_points); % Common x-axis for KDEs

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

%% Plot KDE + Histogram using histogram()
csmcolor = rgb('DimGray');
figure; hold on;

% ----- Histogram -----
% histogram(all_data_sham, 40, 'Normalization', 'pdf', ...
%     'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

% ----- KDE Fill and Line -----
fill([xi_common, fliplr(xi_common)], ...
     [mean(all_f_sham, 1), fliplr(zeros(1, num_points))], ...
     csmcolor, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(xi_common, mean(all_f_sham, 1),'Color',csmcolor, 'LineWidth', 3);


% ----- Axis Labels & Formatting -----
xlim([-0.1 0.2]);
xticks([-0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2]);
xlabel('Time from threshold crossing (s)');
ylabel('Density');
title('Average KDE with Underlying Histogram (ExcitRR_Sham)');

set(gca, 'FontSize', 12, 'LineWidth', 1);
