
%% Need to have run the NoeffecRRtiming script to separate out reaches 
T0_divide=[];
thresh=0.16;
an=fieldnames(ExcitRR_Ctrl);
for a=1:length(an)
    t0r=ExcitRR_Ctrl.(an{a}).RRT0;
    t0r=accel_toreaches_jerk(t0r);
    [ReReach, StimRR, MaxVel,ReachesWithoutSustainedJerk] = CSMEarlyThresh_jerk_edit5(t0r,thresh);
    %[ReReach, StimRR, MaxVel] = ReReachFilter_accel_T0_edit(t0r);
    T0_divide.(an{a}).ReReach=ReReach;
    T0_divide.(an{a}).StimRR=StimRR;
    T0_divide.(an{a}).NoJerk=ReachesWithoutSustainedJerk;


end
%% 
BigRR=[];
ExcitMean=[];
ExcitStd=[];
excitam=fieldnames(T0_divide);
for m = 1:length(excitam)
    excitall = T0_divide.(excitam{m}).StimRR(:,1);
    BigRR=[BigRR;excitall];
    excitmean=mean(excitall); excitstd=std(excitall); 
    ExcitMean=[ExcitMean;excitmean];
    ExcitStd=[ExcitStd; excitstd]; 
  
end
%% plot 

num_bins=ceil(sqrt(length(BigRR)));
excit_mean=mean(BigRR) 
mean2=excit_mean;
excit_std=std(BigRR);
excitcolor=rgb('black');

num_bins=20;
figure; hold on
histogram(BigRR,num_bins,'Normalization','probability','EdgeColor','none','FaceColor',excitcolor)
hold on

xline(excit_mean, 'Color', excitcolor, 'LineWidth', 2);
%xlim([-0.2 0])
%ylim([0 0.1])
xlabel('Time from Positional Threshold (s)')
ylabel('Probability')
%% Plot bin edge values 
num_bins = ceil(sqrt(length(BigRR)));  % Calculate number of bins
excit_mean = mean(BigRR);
excit_std = std(BigRR);
excitcolor = rgb('black');

% Create histogram and capture bin edges
figure; hold on;
h = histogram(BigRR, num_bins, 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceColor', excitcolor);

% Get bin edges
bin_edges = h.BinEdges;  
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;  % Compute bin centers

% Plot vertical line at mean
xline(excit_mean, 'Color', excitcolor, 'LineWidth', 2);

% Set x-axis ticks at bin edges
xticks(bin_edges);  
xtickformat('%.3f');  % Format x-axis tick labels to 2 decimal places

% Labels
xlabel('Time from Positional Threshold (s)');
ylabel('Probability');
hold off;
%% Slightly different way to plot the histogram 
% Compute number of bins using Freedman-Diaconis rule (recommended)
q25 = prctile(BigRR, 25);  
q75 = prctile(BigRR, 75);
bin_width = 2 * (q75 - q25) / length(BigRR)^(1/3);
num_bins = round((max(BigRR) - min(BigRR)) / bin_width);

% Ensure num_bins is at least 1
num_bins = max(num_bins, 1);  

% Plot histogram
figure; hold on;
excitcolor = rgb('black');
h = histogram(BigRR, num_bins, 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceColor', excitcolor);

% Get bin edges
bin_edges = h.BinEdges;  
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

% Plot vertical line at mean
xline(mean(BigRR), 'Color', excitcolor, 'LineWidth', 2);

% Set x-axis ticks at bin edges with 3 decimal points
xticks(bin_edges);
xtickformat('%.3f');

% Labels
xlabel('Time from Positional Threshold (s)');
ylabel('Probability');
hold off;

%%
thresh=0.16;

dirdir= 1;
an=fieldnames(T0_divide);
for a=1:length(an)
    Rreaches=ExcitRR_Ctrl.(an{a}).RRT0;
    C=AlignThresh(Rreaches,thresh,dirdir,200); C=stimvel_accel(C); C=stimvel_accel_jerk(C);
    rrt=T0_divide.(an{a}).StimRR; rrt=rrt(:,1);
    

    % Timings taken from histo from all animals 
    TearlyReaches=[]; % CSM is occuring at the beginning of the reach 
    TthreshReaches=[]; % CSM is occurring closer to the threshold 
   
 
    
    % Bins for RR times 
    Tearly=find(rrt<0);
    Tthresh=find(rrt==0);
   

    rnum=fieldnames(Rreaches);
    
    % for all reaches with a CSM time that is earlier 
   for s=1:length(Tearly)
        idx0=Tearly(s);
        TearlyReaches.(rnum{idx0})= ExcitRR_Ctrl.(an{a}).RRT0.(rnum{idx0});  
    end 
    TearlyKin=BasicKinematics(TearlyReaches);
 
    STearly = AlignThresh(TearlyReaches,thresh,dirdir,200) ; % for heatmap
    STearly=stimvel_accel(STearly); STearly=stimvel_accel_jerk(STearly);
  

    for f=1:length(Tthresh)
        idx1=Tthresh(f);
        TthreshReaches.(rnum{idx1})= ExcitRR_Ctrl.(an{a}).RRT0.(rnum{idx1});  
    end 


     TthreshKin=BasicKinematics(TthreshReaches);
     STthresh = AlignThresh(TthreshReaches,thresh,dirdir,200);
     STthresh=stimvel_accel(STthresh);STthresh=stimvel_accel_jerk(STthresh);
    

    % Reaches 
     T0_divide.(an{a}).AllT0=Rreaches;
     T0_divide.(an{a}).EarlyReaches=TearlyReaches;
     T0_divide.(an{a}).ThreshReaches=TthreshReaches;
 
    tearlyfields=fieldnames(TearlyReaches);
    tthreshfields=fieldnames(TthreshReaches);

    rrnumearly=length(tearlyfields);
    rrnumthresh=length(tthreshfields);
 
    T0_divide.(an{a}).rrnumearly=rrnumearly;
    T0_divide.(an{a}).rrnumthresh=rrnumthresh;

     % Kinematics
     T0_divide.(an{a}).EarlyKin=TearlyKin;
     T0_divide.(an{a}).ThreshKin=TthreshKin;

     %Stimvel
     T0_divide.(an{a}).ST0=C;
     T0_divide.(an{a}).STearly=STearly;
     T0_divide.(an{a}).STthresh=STthresh;
 

     % success rate 
     [S_Tearly, F_Tearly] = SuccessFilter(TearlyReaches);
     [S_Tthresh, F_Tthresh] = SuccessFilter(TthreshReaches);

      T0_divide.(an{a}).successrate.early = length(fieldnames(S_Tearly)) / length(fieldnames(TearlyReaches));
      T0_divide.(an{a}).successrate.thresh= length(fieldnames(S_Tthresh)) / length(fieldnames(TthreshReaches));

end 

 %% Combine Stimvel and look at plots 
Early_big=struct;
Thresh_big=struct;
T0All_big=struct;

ar=fieldnames(T0_divide);
for a = 1:length(ar) 
  Early_big.stimvel(:,:,a) = T0_divide.(ar{a}).STearly.Average;
  Thresh_big.stimvel(:,:,a) = T0_divide.(ar{a}).STthresh.Average;
  T0All_big.stimvel(:,:,a) = T0_divide.(ar{a}).ST0.Average;
end

Early_big.Average = nanmean(Early_big.stimvel,3);
Early_big.Std = nanstd(Early_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(Early_big.Average),1); 
for nnum = 1:length(Early_big.Average)
    samp = squeeze(Early_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
Early_big.N = N;



T0All_big.Average = nanmean(T0All_big.stimvel,3);
T0All_big.Std = nanstd(T0All_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T0All_big.Average),1); 
for nnum = 1:length(T0All_big.Average)
    samp = squeeze(T0All_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
T0All_big.N = N;

Thresh_big.Average = nanmean(Thresh_big.stimvel,3);
Thresh_big.Std = nanstd(Thresh_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(Thresh_big.Average),1); 
for nnum = 1:length(Thresh_big.Average)
    samp = squeeze(Thresh_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
Thresh_big.N = N;



%% Plot stimvel
x1=0 ;x2=0.05;
%y_upper= 1.5; y_lower=-0.5;
%Add acceleration column to stimvel
T0color=rgb('Thistle');
Earlycolor=rgb('Red');
Threshcolor=rgb('DarkRed');


dim =5;

StimvelPlot(T0All_big,dim,T0color,1)
StimvelPlot(Early_big,dim,Earlycolor,0)
StimvelPlot(Thresh_big,dim,Threshcolor,0)

xlim([-0.2 0.2])
%xticks(-0.2:0.01:0.2)
ylim('auto')
%ylim([-0.8 1.5])
%ylim([0 0.27])
%ylim([-20 20])
%pbaspect([1 2 1])
%ylim([0 0.3])
xline(0,'--')
xline(0.05,'--')
yline(0,'--')

xlabel('Time from Threshold Crossing(s)')
%ylabel('Speed')
%ylabel('Outward Position (cm)')
ylabel('Outward Velocity (cm/s)')
%ylabel('Outward acceleration')
%% Separate out 0 CSM from StimRR
% For the CTRL reaches
excitam = fieldnames(ExcitRR_Ctrl);
for b=1:length(excitam)
    RR_Excit=ExcitRR_Ctrl.(excitam{b}).StimRR;
 
    nonZeroRows_ex = [];

    % Iterate over each row of the array
    for i = 1:size(RR_Excit, 1)
        % Check if the current row contains any zeros
        if any(RR_Excit(i,:) == 0)
            % If zeros are present, skip this row
            continue;
        end
        
        % Append the current row to nonZeroRows if no zeros are present
        nonZeroRows_ex = [nonZeroRows_ex; RR_Excit(i,:)];
    end
    
        % Display the resulting array
        disp(nonZeroRows_ex);
        ExcitRR_Ctrl.(excitam{b}).StimRREdit=nonZeroRows_ex;
end
%% Now combine TO divide with the StimRREdit - where all 0s have been removed

for a=1:length(an)
    combinedData=[];
    t0nums=T0_divide.(an{a}).StimRR(:,1);
    allnums=ExcitRR_Ctrl.(an{a}).StimRREdit(:,1);
    combinedData = [t0nums; allnums]; 
     ExcitRR_Ctrl.(an{a}).StimRRAll=combinedData;
end
%% Now plot combined histo 

BigRR_Ctrl=[];
Ctrl_mean=[];
Ctrl_std=[];
excitam=fieldnames(T0_divide);
for m = 1:length(excitam)
    excitall = T0_divide.(excitam{m}).StimRRAll(:,1);
    BigRR_Ctrl=[BigRR_Ctrl;excitall];
    excitmean=mean(excitall); excitstd=std(excitall); 
    Ctrl_mean=[Ctrl_mean;excitmean];
    Ctrl_std=[Ctrl_std; excitstd];   
end
%% Now for the stim RR

BigRR_Stim=[];
Stim_mean=[];
Stim_std=[];
an=fieldnames(ExcitRR);
for a = 1:length(an)
    excitall = ExcitRR.(an{a}).StimRR(:,1);
    BigRR_Stim=[BigRR_Stim;excitall];
    excitmean=mean(excitall); excitstd=std(excitall); 
    Stim_mean=[Stim_mean;excitmean];
    Stim_std=[Stim_std; excitstd];   
end
%% plot 

num_bins=ceil(sqrt(length(BigRR_Ctrl)));
excit_mean=mean(BigRR_Ctrl) 
mean2=excit_mean;
excit_std=std(BigRR_Ctrl);
ctrlcolor=rgb('black');

num_bins=25;
figure; hold on
histogram(BigRR_Ctrl,num_bins,'Normalization','pdf','EdgeColor','none','FaceColor',ctrlcolor)
hold on

%xline(excit_mean, 'Color', excitcolor, 'LineWidth', 2);
xline(0,'--')
%xlim([-0.2 0])
%ylim([0 0.1])
xlabel('Time from Positional Threshold (s)')
ylabel('Probability')


%%
%average ksdensity across all anims 
barcolor=rgb('DimGray');
[f2_comb,xi2_comb,bw]=ksdensity(BigRR_Ctrl,'Bandwidth',0.008);

ssmean=mean(BigRR_Ctrl);


% Plot the mean ksdensity on the same plot
figure; hold on
plot(xi2_comb, f2_comb,'LineWidth', 4,'Color',barcolor);

%xline(hf_mean, 'Color', RNstimcolor, 'LineWidth', 2);
xlim([-0.1 0.2])
xlabel('Time from Pos. Thresh')


%%
data = BigRR_Ctrl;  % Your data
num_binsall = 31;  % Number of bins for the histogram
num_binsstim=16;
stimcolor=rgb('DodgerBlue');
% Create a figure
figure; hold on;

% Plot the histogram
[counts, edges] = histcounts(data, num_binsall, 'Normalization', 'pdf');
[counts_stim, edges_stim] = histcounts(BigRR_Stim, num_binsstim, 'Normalization', 'pdf');

% Calculate bin centers
bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
bin_centers_stim = (edges_stim(1:end-1) + edges_stim(2:end)) / 2;

% Plot the histogram bars
hBar=bar(bin_centers, counts, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'BarWidth', 0.9);

% Hold on to overlay the ksdensity plot
hold on;

hBar1=bar(bin_centers_stim, counts_stim, 'FaceColor', [ 0.1172    0.5625    1.0000], 'EdgeColor', 'none', 'BarWidth', 0.9);

% Hold on to overlay the ksdensity plot
hold on;

% Plot the ksdensity on the same figure
[f2_comb, xi2_comb, bw] = ksdensity(data, 'Bandwidth', 0.008);
disp(['Default Bandwidth: ', num2str(bw)]);

hKsdensity = plot(xi2_comb, f2_comb, 'LineWidth', 3, 'Color', 'k');

hold on
% Plot the ksdensity on the same figure
[f2_comb_stim, xi2_comb_stim, bw] = ksdensity(BigRR_Stim, 'Bandwidth', 0.008);
disp(['Default Bandwidth: ', num2str(bw)]);

hKsdensity1 = plot(xi2_comb_stim, f2_comb_stim, 'LineWidth', 3, 'Color',stimcolor);

% Set axis limits
xlim([-0.1 0.2]);
%ylim([0, max(counts)*1.1]);

% Set custom x-ticks on the bottom axis
xticks([-0.1, -0.05,0, 0.05, 0.1, 0.15, 0.2]); % Define clear tick marks

set(gca, 'FontSize', 12);

% Add labels and title
xlabel('Time from threshold crossing (s)');
ylabel('Density');
title('Histogram and ksdensity plot');
legend([hBar, hKsdensity], {'Histogram', 'ksdensity'}, 'Location', 'best');
%% 
data = BigRR_Ctrl;  % Your data
num_bins = 31;  % Consistent number of bins
stimcolor = rgb('DodgerBlue');

% Define common bin edges
min_val = min([BigRR_Ctrl; BigRR_Stim]);  % Find the minimum value across both datasets
max_val = max([BigRR_Ctrl; BigRR_Stim]);  % Find the maximum value across both datasets
edges = linspace(min_val, max_val, num_bins + 1);  % Define common bin edges

% Compute histograms using the same bin edges
[counts, ~] = histcounts(data, edges, 'Normalization', 'pdf');
[counts_stim, ~] = histcounts(BigRR_Stim, edges, 'Normalization', 'pdf');

% Calculate bin centers
bin_centers = (edges(1:end-1) + edges(2:end)) / 2;

% Create figure
figure; hold on;

% Plot the histograms
hBar = bar(bin_centers, counts, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'BarWidth', 0.9);
hBar1 = bar(bin_centers, counts_stim, 'FaceColor', stimcolor, 'EdgeColor', 'none', 'BarWidth', 0.9,'FaceAlpha',0.5); % Adjust BarWidth if needed

% Compute and plot kernel density estimates
[f_ctrl, xi_ctrl] = ksdensity(data, 'Bandwidth', 0.008);
[f_stim, xi_stim] = ksdensity(BigRR_Stim, 'Bandwidth', 0.008);

hKsdensity = plot(xi_ctrl, f_ctrl, 'LineWidth', 3, 'Color', barcolor);
hKsdensity1 = plot(xi_stim, f_stim, 'LineWidth', 3, 'Color', stimcolor);

% Set axis limits
%xlim([min_val, max_val]);
xlim([-0.1 0.2])
% Customize x-ticks
xticks([-0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2]);

set(gca, 'FontSize', 12);

% Add labels and title
xlabel('Time from threshold crossing (s)');
ylabel('Density');
title('Histogram and ksdensity plot');

% Ensure x-ticks are set properly
ax = gca;
ax.XTick = [-0.1 -0.05 0 0.05 0.1 0.15 0.2];  % Manually set X-ticks if needed
ax.XTickLabelRotation = 0;  % Optional: Rotate labels if needed

% Increase tick mark size and make them thicker and black
ax.TickLength = [0.02, 0.02];  % Make the ticks longer
ax.XAxis.LineWidth = 1;  % Thicker line for the X-axis
ax.YAxis.LineWidth = 1;  % Thicker line for the Y-axis

% Set the color of tick marks to black
ax.XAxis.TickDirection = 'in';  % Make ticks point outwards
ax.YAxis.TickDirection = 'in';  % Make ticks point outwards

% Make tick labels larger
ax.XAxis.FontSize = 12;  % Increase font size for X-tick labels
ax.YAxis.FontSize = 12;  % Increase font size for Y-tick labels
%% KSdensity with area under the curve colored 
% Your data
data = BigRR_Ctrl;
num_bins = 31;  % Number of bins (no longer needed for histograms)
stimcolor = rgb('DodgerBlue');
ctrlcolor = rgb('Gray'); % Define control color

% Define common bin edges
min_val = min([BigRR_Ctrl; BigRR_Stim]);
max_val = max([BigRR_Ctrl; BigRR_Stim]);

% Compute kernel density estimates
[f_ctrl, xi_ctrl] = ksdensity(data, 'Bandwidth', 0.008);
[f_stim, xi_stim] = ksdensity(BigRR_Stim, 'Bandwidth', 0.008);

% [f_ctrl, xi_ctrl] = ksdensity(data);
% [f_stim, xi_stim] = ksdensity(BigRR_Stim);
% Create figure
figure; hold on;

% Fill under KDE curves
fill([xi_ctrl fliplr(xi_ctrl)], [f_ctrl zeros(size(f_ctrl))], ctrlcolor, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
fill([xi_stim fliplr(xi_stim)], [f_stim zeros(size(f_stim))], stimcolor, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

% Plot the density curves on top
plot(xi_ctrl, f_ctrl, 'LineWidth', 3, 'Color', ctrlcolor);
plot(xi_stim, f_stim, 'LineWidth', 3, 'Color', stimcolor);

% Set axis limits
xlim([-0.1 0.2]);

% Customize x-ticks
xticks([-0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2]);

% Axis properties
ax = gca;
ax.FontSize = 12;
ax.XTickLabelRotation = 0;
ax.TickLength = [0.02, 0.02];
ax.XAxis.LineWidth = 1;
ax.YAxis.LineWidth = 1;
ax.XAxis.TickDirection = 'in';
ax.YAxis.TickDirection = 'in';

% Labels and title
xlabel('Time from threshold crossing (s)');
ylabel('Density');
title('Kernel Density Estimate with Shaded Area');

hold off;
%% Average Ksdensity 
% Initialize storage for KDEs
num_points = 100; % Number of points for KDE interpolation
xi_common = linspace(-0.1, 0.2, num_points); % Common x-axis for KDEs

% -----------------------------
% Compute KDE for T0_divide.(excitam{m}).StimRRAll
% -----------------------------
all_f_T0 = []; % Store KDEs for T0_divide
excitam_T0 = fieldnames(T0_divide);
num_T0 = length(excitam_T0);

for m = 1:num_T0
    allbins_T0 = T0_divide.(excitam_T0{m}).StimRRAll;  % Extract data
    [f_T0, xi_T0] = ksdensity(allbins_T0, 'Bandwidth', 0.008);  % Compute KDE
    f_interp_T0 = interp1(xi_T0, f_T0, xi_common, 'linear', 0);  % Interpolate KDE
    all_f_T0 = [all_f_T0; f_interp_T0];  % Store interpolated KDE
end

% -----------------------------
% Compute KDE for ExcitRR.(an{a}).StimRR(:,1)
% -----------------------------
all_f_StimRR = []; % Store KDEs for ExcitRR
an = fieldnames(ExcitRR);
num_an = length(an);

for a = 1:num_an
    allbins_StimRR = ExcitRR.(an{a}).StimRR(:,1);  % Extract data
    [f_StimRR, xi_StimRR] = ksdensity(allbins_StimRR, 'Bandwidth', 0.008);  % Compute KDE
    f_interp_StimRR = interp1(xi_StimRR, f_StimRR, xi_common, 'linear', 0);  % Interpolate KDE
    all_f_StimRR = [all_f_StimRR; f_interp_StimRR];  % Store interpolated KDE
end

% Collect all data points across animals
all_data_T0 = [];
for m = 1:num_T0
    all_data_T0 = [all_data_T0; T0_divide.(excitam_T0{m}).StimRRAll];
end

all_data_StimRR = [];
for a = 1:num_an
    all_data_StimRR = [all_data_StimRR; ExcitRR.(an{a}).StimRR(:,1)];
end

% Run two-sample KS test
[h_ks, p_ks, ks_stat] = kstest2(all_data_T0, all_data_StimRR);

% Display results
fprintf('KS test result:\n');
fprintf('  h = %d (1 = reject null hypothesis)\n', h_ks);
fprintf('  p = %.4f\n', p_ks);
fprintf('  KS statistic = %.4f\n', ks_stat);
% -----------------------------
% Plot the KDEs with Full Area Under the Curve Filled
% -----------------------------
stimcolor=rgb('DodgerBlue');
figure; hold on;

% Fill area for T0_divide group (red)
fill([xi_common, fliplr(xi_common)], [mean(all_f_T0, 1), fliplr(zeros(1, num_points))], ...
    [0, 0, 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Red fill

% Fill area for ExcitRR group (blue)
fill([xi_common, fliplr(xi_common)], [mean(all_f_StimRR, 1), fliplr(zeros(1, num_points))], ...
    [0.1172,0.5625,1.0000], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Blue fill

% Plot mean KDE lines (only lines, no SEM)
plot(xi_common, mean(all_f_T0, 1), 'k', 'LineWidth', 3); % T0_divide (red)
plot(xi_common, mean(all_f_StimRR, 1),'Color',stimcolor, 'LineWidth', 3); % ExcitRR (blue)

% Set axis limits
xlim([-0.1 0.2]);
%ylim([0 35])
% Customize x-ticks
xticks([-0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2]);

% Labels and title
xlabel('Time from threshold crossing (s)');
ylabel('Density');
title('Average KDE across Animals (T0_divide vs ExcitRR)');

% Formatting
set(gca, 'FontSize', 12, 'LineWidth', 1);
% Format the result as a string
ks_text = sprintf('KS test: p = %.4f, D = %.3f', p_ks, ks_stat);

% Decide where to place it
% You can adjust the x and y position based on your y-limits
ymax = max([mean(all_f_T0, 1), mean(all_f_StimRR, 1)], [], 'all');
text_x = 0.05;  % Change this if it overlaps data
text_y = ymax * 0.9;

% Add to plot
text(text_x, text_y, ks_text, 'FontSize', 12, 'FontWeight', 'bold');

%% Sham vs. Ctrl
%% Average Ksdensity 
% Initialize storage for KDEs
num_points = 100; % Number of points for KDE interpolation
xi_common = linspace(-0.1, 0.2, num_points); % Common x-axis for KDEs

% -----------------------------
% Compute KDE for T0_divide.(excitam{m}).StimRRAll
% -----------------------------
all_f_T0 = []; % Store KDEs for T0_divide
excitam_T0 = fieldnames(T0_divide);
num_T0 = length(excitam_T0);

for m = 1:num_T0
    allbins_T0 = T0_divide.(excitam_T0{m}).StimRRAll;  % Extract data
    [f_T0, xi_T0] = ksdensity(allbins_T0, 'Bandwidth', 0.008);  % Compute KDE
    f_interp_T0 = interp1(xi_T0, f_T0, xi_common, 'linear', 0);  % Interpolate KDE
    all_f_T0 = [all_f_T0; f_interp_T0];  % Store interpolated KDE
end

% -----------------------------
% Compute KDE for ExcitRR_Sham.(an{a}).StimRR(:,1)
% -----------------------------
all_f_StimRR = []; % Store KDEs for ExcitRR_Sham
an = fieldnames(ExcitRR_Sham);
num_an = length(an);

for a = 1:num_an
    allbins_StimRR = ExcitRR_Sham.(an{a}).AllBins(:,1);  % Extract data
    [f_StimRR, xi_StimRR] = ksdensity(allbins_StimRR, 'Bandwidth', 0.008);  % Compute KDE
    f_interp_StimRR = interp1(xi_StimRR, f_StimRR, xi_common, 'linear', 0);  % Interpolate KDE
    all_f_StimRR = [all_f_StimRR; f_interp_StimRR];  % Store interpolated KDE
end

% -----------------------------
% Perform KS test on raw data
% -----------------------------
% Collect all data points across animals
all_data_T0 = [];
for m = 1:num_T0
    all_data_T0 = [all_data_T0; T0_divide.(excitam_T0{m}).StimRRAll];
end

all_data_StimRR = [];
for a = 1:num_an
    all_data_StimRR = [all_data_StimRR; ExcitRR_Sham.(an{a}).StimRR(:,1)];
end

% Run two-sample KS test
[h_ks, p_ks, ks_stat] = kstest2(all_data_T0, all_data_StimRR);
% Display results
fprintf('KS test result:\n');
fprintf('  h = %d (1 = reject null hypothesis)\n', h_ks);
fprintf('  p = %.4f\n', p_ks);
fprintf('  KS statistic = %.4f\n', ks_stat);

% -----------------------------
% Plot the KDEs with Full Area Under the Curve Filled
% -----------------------------
stimcolor = rgb('DodgerBlue');
figure; hold on;

% Fill area for T0_divide group (black)
fill([xi_common, fliplr(xi_common)], [mean(all_f_T0, 1), fliplr(zeros(1, num_points))], ...
    [0, 0, 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Fill area for ExcitRR_Sham group (blue)
fill([xi_common, fliplr(xi_common)], [mean(all_f_StimRR, 1), fliplr(zeros(1, num_points))], ...
    [0.1172, 0.5625, 1.0000], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Plot mean KDE lines
plot(xi_common, mean(all_f_T0, 1), 'k', 'LineWidth', 3); % T0_divide
plot(xi_common, mean(all_f_StimRR, 1), 'Color', stimcolor, 'LineWidth', 3); % ExcitRR_Sham

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
ymax = max([mean(all_f_T0, 1), mean(all_f_StimRR, 1)], [], 'all');
ks_text = sprintf('KS test: p = %.4f, D = %.3f', p_ks, ks_stat);
text(0.05, ymax * 0.9, ks_text, 'FontSize', 12, 'FontWeight', 'bold');
