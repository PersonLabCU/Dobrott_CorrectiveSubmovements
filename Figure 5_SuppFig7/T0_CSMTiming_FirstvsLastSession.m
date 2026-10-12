%% Need to have run the NoeffecRRtiming script to separate out reaches 
% First Session
T0_dividefirst=[];
thresh=0.16;
an=fieldnames(ExcitRR_Ctrl1);
for a=1:length(an)
    t0r=ExcitRR_Ctrl1.(an{a}).RRT0;
    t0r=accel_toreaches_jerk(t0r);
    [ReReach, StimRR, MaxVel,ReachesWithoutSustainedJerk] = CSMEarlyThresh_jerk_edit5(t0r,thresh);
    %[ReReach, StimRR, MaxVel] = ReReachFilter_accel_T0_edit(t0r);
    T0_dividefirst.(an{a}).ReReach=ReReach;
    T0_dividefirst.(an{a}).StimRR=StimRR;
    T0_dividefirst.(an{a}).NoJerk=ReachesWithoutSustainedJerk;

end
%%
%Last session
T0_dividelast=[];
thresh=0.16;
for a=1:length(an)
    t0r=ExcitRR_Ctrlend.(an{a}).RRT0;
    if isempty(t0r)
        t0r=[];
        ReReach=[];
        StimRR=[];
        NoJerk=[];
    else
    t0r=accel_toreaches_jerk(t0r);
    [ReReach, StimRR, MaxVel,ReachesWithoutSustainedJerk] = CSMEarlyThresh_jerk_edit5(t0r,thresh);
    end
    %[ReReach, StimRR, MaxVel] = ReReachFilter_accel_T0_edit(t0r);
    T0_dividelast.(an{a}).ReReach=ReReach;
    T0_dividelast.(an{a}).StimRR=StimRR;
    T0_dividelast.(an{a}).NoJerk=ReachesWithoutSustainedJerk;
   
end
%% 
BigRR_first=[];

excitam=fieldnames(T0_dividefirst);
for m = 1:length(excitam)
    excitall = T0_dividefirst.(excitam{m}).StimRR(:,1);
    BigRR_first=[BigRR_first;excitall];
    
end
%% 
BigRR_last=[];

excitam=fieldnames(T0_dividelast);
for m = 1:4
    excitall = T0_dividelast.(excitam{m}).StimRR(:,1);
    BigRR_last=[BigRR_last;excitall];
    
end

%% For the CTRL first session reaches

excitam = fieldnames(ExcitRR_Ctrl1);
for b=1:length(excitam)
    RR_Excit=ExcitRR_Ctrl1.(excitam{b}).StimRR;
 
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
        ExcitRR_Ctrl1.(excitam{b}).StimRREdit=nonZeroRows_ex;
end

%% For the CTRL first session reaches
excitam = fieldnames(ExcitRR_Ctrlend);
for b=1:length(excitam)
    RR_Excit=ExcitRR_Ctrlend.(excitam{b}).StimRR;
 
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
        ExcitRR_Ctrlend.(excitam{b}).StimRREdit=nonZeroRows_ex;
end

%% 
BigRR_first1=[];

excitam=fieldnames(ExcitRR_Ctrl1);
for m =1:length(excitam)
    allbins=[];
    excitall = ExcitRR_Ctrl1.(excitam{m}).StimRREdit(:,1);
    t0bins=T0_dividefirst.(excitam{m}).StimRR(:,1);
    allbins=[excitall;t0bins];
    ExcitRR_Ctrl1.(excitam{m}).AllBins=allbins;
    BigRR_first1=[BigRR_first1;excitall];
    
end
%% 
BigRR_last1 = [];

excitam = fieldnames(ExcitRR_Ctrlend);
for m = 1:4
    allbins = [];
    
    % Extract data
    excitall = ExcitRR_Ctrlend.(excitam{m}).StimRREdit(:,1);
    t0bins = T0_dividelast.(excitam{m}).StimRR(:,1);
    

    % Concatenate filtered data
    allbins = [excitall; t0bins];
    
    % Store results
    ExcitRR_Ctrlend.(excitam{m}).AllBins = allbins;
    BigRR_last1 = [BigRR_last1; excitall];
end
%% just change the last animal who has no T0
ExcitRR_Ctrlend.EX028_strongstim.AllBins=ExcitRR_Ctrlend.EX028_strongstim.StimRREdit(:,1);
%% Combine T0s and all the rest 
BigRR_firstAll=[]; BigRR_lastAll=[];
BigRR_firstAll=[BigRR_first;BigRR_first1];
BigRR_lastAll=[BigRR_last;BigRR_last1];

%% Now plot 
%num_bins=ceil(sqrt(length(BigRR_firstAll)));
% First vs. last 
num_bins = 17;  % Consistent number of bins
firstcolor = rgb('Black');
lastcolor=rgb('Red');
% Define common bin edges
min_val = min([BigRR_firstAll; BigRR_lastAll]);  % Find the minimum value across both datasets
max_val = max([BigRR_firstAll; BigRR_lastAll]);  % Find the maximum value across both datasets
edges = linspace(min_val, max_val, num_bins + 1);  % Define common bin edges

% Compute histograms using the same bin edges
[counts, ~] = histcounts(BigRR_firstAll, edges, 'Normalization', 'pdf');
[counts_stim, ~] = histcounts(BigRR_lastAll, edges, 'Normalization', 'pdf');

% Calculate bin centers
bin_centers = (edges(1:end-1) + edges(2:end)) / 2;

% Create figure
figure; hold on;

% Plot the histograms
hBar = bar(bin_centers, counts, 'FaceColor',firstcolor, 'EdgeColor', 'none', 'BarWidth', 0.9,'FaceAlpha',0.50);
hBar1 = bar(bin_centers, counts_stim, 'FaceColor', lastcolor, 'EdgeColor', 'none', 'BarWidth', 0.9,'FaceAlpha',0.50); % Adjust BarWidth if needed

% Compute and plot kernel density estimates
[f_ctrl, xi_ctrl] = ksdensity(BigRR_firstAll, 'Bandwidth', 0.015);
[f_stim, xi_stim] = ksdensity(BigRR_lastAll, 'Bandwidth', 0.01);

hKsdensity = plot(xi_ctrl, f_ctrl, 'LineWidth', 3, 'Color', firstcolor);
hKsdensity1 = plot(xi_stim, f_stim, 'LineWidth', 3, 'Color', lastcolor);

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

%% All 3 on the same plot 
% Now plot 
% First vs. last 
num_bins = 20;  % Consistent number of bins
firstcolor = rgb('Black');
lastcolor=rgb('OrangeRed');
shamcolor=rgb('DimGray');
% Define common bin edges
min_val = min([BigRR_firstAll; BigRR_lastAll;BigRR_Sham]);  % Find the minimum value across both datasets
max_val = max([BigRR_firstAll; BigRR_lastAll;BigRR_Sham]);  % Find the maximum value across both datasets
edges = linspace(min_val, max_val, num_bins + 1);  % Define common bin edges

% Compute histograms using the same bin edges
[counts, ~] = histcounts(BigRR_firstAll, edges, 'Normalization', 'pdf');
[counts_stim, ~] = histcounts(BigRR_lastAll, edges, 'Normalization', 'pdf');
[counts_sham, ~] = histcounts(BigRR_Sham, edges, 'Normalization', 'pdf');

% Calculate bin centers
bin_centers = (edges(1:end-1) + edges(2:end)) / 2;

% Create figure
figure; hold on;

% Plot the histograms
hBar = bar(bin_centers, counts, 'FaceColor',firstcolor, 'EdgeColor', 'none', 'BarWidth', 0.9,'FaceAlpha',0.75);
hBar1 = bar(bin_centers, counts_stim, 'FaceColor', lastcolor, 'EdgeColor', 'none', 'BarWidth', 0.9,'FaceAlpha',0.50); % Adjust BarWidth if needed
hBar2 = bar(bin_centers, counts_sham, 'FaceColor', shamcolor, 'EdgeColor', 'none', 'BarWidth', 0.9,'FaceAlpha',0.50);
% Compute and plot kernel density estimates
[f_ctrl, xi_ctrl] = ksdensity(BigRR_firstAll, 'Bandwidth', 0.008);
[f_stim, xi_stim] = ksdensity(BigRR_lastAll, 'Bandwidth', 0.008);
[f_sham, xi_sham] = ksdensity(BigRR_Sham, 'Bandwidth', 0.008);

hKsdensity = plot(xi_ctrl, f_ctrl, 'LineWidth', 3, 'Color', firstcolor);
hKsdensity1 = plot(xi_stim, f_stim, 'LineWidth', 3, 'Color', lastcolor);
hKsdensity2 = plot(xi_sham, f_sham, 'LineWidth', 3, 'Color',shamcolor);

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

%% KSdensity plots
%% Now plot filled ksdensity
num_bins = 17;  % Number of bins (not needed for KDE)
firstcolor = rgb('Black');
lastcolor = rgb('Red');

% Compute kernel density estimates
[f_first, xi_first] = ksdensity(BigRR_firstAll, 'Bandwidth', 0.015);
[f_last, xi_last] = ksdensity(BigRR_lastAll, 'Bandwidth', 0.01);

% Create figure
figure; hold on;

% Fill under KDE curves
fill([xi_first fliplr(xi_first)], [f_first zeros(size(f_first))], firstcolor, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
fill([xi_last fliplr(xi_last)], [f_last zeros(size(f_last))], lastcolor, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

% Plot the density curves on top
plot(xi_first, f_first, 'LineWidth', 3, 'Color', firstcolor);
plot(xi_last, f_last, 'LineWidth', 3, 'Color', lastcolor);

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
%% Sham and last 
%% Sham vs. CtrlEnd - Filled ksdensity plot

shamcolor = rgb('DimGray');
ctrlEndcolor = rgb('OrangeRed'); % Assuming CtrlEnd is represented by lastcolor

% Compute kernel density estimates
[f_sham, xi_sham] = ksdensity(BigRR_Sham, 'Bandwidth', 0.008);
[f_ctrlEnd, xi_ctrlEnd] = ksdensity(BigRR_lastAll, 'Bandwidth', 0.008);

% Create figure
figure; hold on;

% Fill under KDE curves
fill([xi_sham fliplr(xi_sham)], [f_sham zeros(size(f_sham))], shamcolor, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
fill([xi_ctrlEnd fliplr(xi_ctrlEnd)], [f_ctrlEnd zeros(size(f_ctrlEnd))], ctrlEndcolor, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

% Plot the density curves on top
plot(xi_sham, f_sham, 'LineWidth', 3, 'Color', shamcolor);
plot(xi_ctrlEnd, f_ctrlEnd, 'LineWidth', 3, 'Color', ctrlEndcolor);

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
title('Sham vs. CtrlEnd Kernel Density Estimate');

hold off;
%% average ks density 
% Initialize storage for KDEs
num_points = 100; % Number of points for KDE interpolation
all_f = []; % Store KDEs
xi_common = linspace(-0.1, 0.2, num_points); % Common x-axis for all KDEs

excitam = fieldnames(ExcitRR_Ctrlend);
num_animals = length(excitam);

for m = 1:num_animals
    % Extract data
    allbins = ExcitRR_Ctrlend.(excitam{m}).AllBins;
    
    % Compute KDE
    [f, xi] = ksdensity(allbins, 'Bandwidth', 0.008);
    
    % Interpolate KDE onto common x-axis
    f_interp = interp1(xi, f, xi_common, 'linear', 0);
    
    % Store interpolated KDE
    all_f = [all_f; f_interp];
end

% Compute mean and standard error across animals
mean_f = mean(all_f, 1);
sem_f = std(all_f, 0, 1) / sqrt(num_animals);

% Plot mean KDE
figure; hold on;
fill([xi_common, fliplr(xi_common)], [mean_f - sem_f, fliplr(mean_f + sem_f)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(xi_common, mean_f, 'k', 'LineWidth', 3);

% Set axis limits
xlim([-0.1 0.2]);

% Customize x-ticks
xticks([-0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2]);

% Labels and title
xlabel('Time from threshold crossing (s)');
ylabel('Density');
title('Average KDE across Animals');

% Formatting
set(gca, 'FontSize', 12, 'LineWidth', 1);
%%
% Initialize storage for KDEs
num_points = 100; % Number of points for KDE interpolation
xi_common = linspace(-0.1, 0.2, num_points); % Common x-axis

% -----------------------------
% Compute KDE for Ctrlend Group
% -----------------------------
all_f_ctrl = []; % Store KDEs for Ctrlend
excitam_ctrl = fieldnames(ExcitRR_Ctrlend);
num_ctrl = length(excitam_ctrl);

for m = 1:num_ctrl
    allbins_ctrl = ExcitRR_Ctrlend.(excitam_ctrl{m}).AllBins;
    [f_ctrl, xi_ctrl] = ksdensity(allbins_ctrl, 'Bandwidth', 0.008);
    f_interp_ctrl = interp1(xi_ctrl, f_ctrl, xi_common, 'linear', 0);
    all_f_ctrl = [all_f_ctrl; f_interp_ctrl];
end

% Compute mean and SEM for Ctrlend
mean_f_ctrl = mean(all_f_ctrl, 1);
sem_f_ctrl = std(all_f_ctrl, 0, 1) / sqrt(num_ctrl);

% -----------------------------
% Compute KDE for Sham Group
% -----------------------------
all_f_sham = []; % Store KDEs for Sham
excitam_sham = fieldnames(ExcitRR_Sham);
num_sham = length(excitam_sham);

for m = 1:num_sham
    allbins_sham = ExcitRR_Sham.(excitam_sham{m}).AllBins;
    [f_sham, xi_sham] = ksdensity(allbins_sham, 'Bandwidth', 0.008);
    f_interp_sham = interp1(xi_sham, f_sham, xi_common, 'linear', 0);
    all_f_sham = [all_f_sham; f_interp_sham];
end

% Compute mean and SEM for Sham
mean_f_sham = mean(all_f_sham, 1);
sem_f_sham = std(all_f_sham, 0, 1) / sqrt(num_sham);

% -----------------------------
% Plot the KDEs with Fill
% -----------------------------
figure; hold on;

% Fill area for Sham group (gray)
fill([xi_common, fliplr(xi_common)], [mean_f_sham - sem_f_sham, fliplr(mean_f_sham + sem_f_sham)], ...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Light gray fill

% Fill area for Ctrlend group (black)
fill([xi_common, fliplr(xi_common)], [mean_f_ctrl - sem_f_ctrl, fliplr(mean_f_ctrl + sem_f_ctrl)], ...
    [0, 0, 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Black fill

% Plot mean KDE lines
plot(xi_common, mean_f_sham, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 3); % Sham (gray)
plot(xi_common, mean_f_ctrl, 'k', 'LineWidth', 3); % Ctrlend (black)

% Set axis limits
xlim([-0.1 0.2]);

% Customize x-ticks
xticks([-0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2]);

% Labels and title
xlabel('Time from threshold crossing (s)');
ylabel('Density');
title('Average KDE across Animals (Ctrlend vs Sham)');

% Formatting
set(gca, 'FontSize', 12, 'LineWidth', 1);
%%
% Initialize storage for KDEs
num_points = 100; % Number of points for KDE interpolation
xi_common = linspace(-0.1, 0.2, num_points); % Common x-axis

% -----------------------------
% Compute KDE for Ctrlend Group
% -----------------------------
all_f_ctrl = []; % Store KDEs for Ctrlend
excitam_ctrl = fieldnames(ExcitRR_Ctrlend);
num_ctrl = length(excitam_ctrl);

for m = 1:num_ctrl
    allbins_ctrl = ExcitRR_Ctrlend.(excitam_ctrl{m}).AllBins;
    [f_ctrl, xi_ctrl] = ksdensity(allbins_ctrl, 'Bandwidth', 0.008);
    f_interp_ctrl = interp1(xi_ctrl, f_ctrl, xi_common, 'linear', 0);
    all_f_ctrl = [all_f_ctrl; f_interp_ctrl];
end

% -----------------------------
% Compute KDE for Sham Group
% -----------------------------
all_f_sham = []; % Store KDEs for Sham
excitam_sham = fieldnames(ExcitRR_Sham);
num_sham = length(excitam_sham);

for m = 1:num_sham
    allbins_sham = ExcitRR_Sham.(excitam_sham{m}).AllBins;
    [f_sham, xi_sham] = ksdensity(allbins_sham, 'Bandwidth', 0.008);
    f_interp_sham = interp1(xi_sham, f_sham, xi_common, 'linear', 0);
    all_f_sham = [all_f_sham; f_interp_sham];
end

% -----------------------------
% Plot the KDEs with Fill
% -----------------------------
figure; hold on;

% Fill area for Sham group (gray)
fill(xi_common, mean(all_f_sham, 1), ...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Light gray fill

% Fill area for Ctrlend group (black)
fill(xi_common, mean(all_f_ctrl, 1), ...
    [0, 0, 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Black fill

% Plot mean KDE lines (only lines, no SEM)
plot(xi_common, mean(all_f_sham, 1), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 3); % Sham (gray)
plot(xi_common, mean(all_f_ctrl, 1), 'k', 'LineWidth', 3); % Ctrlend (black)

% Set axis limits
xlim([-0.1 0.2]);

% Customize x-ticks
xticks([-0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2]);

% Labels and title
xlabel('Time from threshold crossing (s)');
ylabel('Density');
title('Average KDE across Animals (Ctrlend vs Sham)');

% Formatting
set(gca, 'FontSize', 12, 'LineWidth', 1);
%%
% Initialize storage for KDEs
num_points = 100; % Number of points for KDE interpolation
xi_common = linspace(-0.1, 0.2, num_points); % Common x-axis

% -----------------------------
% Compute KDE for Ctrlend Group
% -----------------------------
all_f_ctrl = []; % Store KDEs for Ctrlend
excitam_ctrl = fieldnames(ExcitRR_Ctrlend);
num_ctrl = length(excitam_ctrl);

for m = 1:num_ctrl
    allbins_ctrl = ExcitRR_Ctrlend.(excitam_ctrl{m}).AllBins;
    [f_ctrl, xi_ctrl] = ksdensity(allbins_ctrl, 'Bandwidth', 0.008);
    f_interp_ctrl = interp1(xi_ctrl, f_ctrl, xi_common, 'linear', 0);
    all_f_ctrl = [all_f_ctrl; f_interp_ctrl];
end

% -----------------------------
% Compute KDE for Sham Group
% -----------------------------
all_f_sham = []; % Store KDEs for Sham
excitam_sham = fieldnames(ExcitRR_Sham);
num_sham = length(excitam_sham);

for m = 1:num_sham
    allbins_sham = ExcitRR_Sham.(excitam_sham{m}).AllBins;
    [f_sham, xi_sham] = ksdensity(allbins_sham, 'Bandwidth', 0.008);
    f_interp_sham = interp1(xi_sham, f_sham, xi_common, 'linear', 0);
    all_f_sham = [all_f_sham; f_interp_sham];
end

% -----------------------------
% Compute KDE for Ctrlfirst Group (ExcitRR_ctrlfirst)
% -----------------------------
all_f_ctrlfirst = []; % Store KDEs for Ctrlfirst
excitam_ctrlfirst = fieldnames(ExcitRR_Ctrl1);
num_ctrlfirst = length(excitam_ctrlfirst);

for m = 1:num_ctrlfirst
    allbins_ctrlfirst = ExcitRR_Ctrl1.(excitam_ctrlfirst{m}).AllBins;
    [f_ctrlfirst, xi_ctrlfirst] = ksdensity(allbins_ctrlfirst, 'Bandwidth', 0.008);
    f_interp_ctrlfirst = interp1(xi_ctrlfirst, f_ctrlfirst, xi_common, 'linear', 0);
    all_f_ctrlfirst = [all_f_ctrlfirst; f_interp_ctrlfirst];
end

% -----------------------------
% Plot the KDEs with Full Area Under the Curve Filled
% -----------------------------
figure; hold on;

% Fill area for Sham group (gray)
fill([xi_common, fliplr(xi_common)], [mean(all_f_sham, 1), fliplr(zeros(1, num_points))], ...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Light gray fill

% Fill area for Ctrlend group (black)
fill([xi_common, fliplr(xi_common)], [mean(all_f_ctrl, 1), fliplr(zeros(1, num_points))], ...
    [0, 0, 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Black fill

% Fill area for Ctrlfirst group (blue)
fill([xi_common, fliplr(xi_common)], [mean(all_f_ctrlfirst, 1), fliplr(zeros(1, num_points))], ...
    [0, 0, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Blue fill

% Plot mean KDE lines (only lines, no SEM)
plot(xi_common, mean(all_f_sham, 1), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 3); % Sham (gray)
plot(xi_common, mean(all_f_ctrl, 1), 'k', 'LineWidth', 3); % Ctrlend (black)
plot(xi_common, mean(all_f_ctrlfirst, 1), 'b', 'LineWidth', 3); % Ctrlfirst (blue)

% Set axis limits
xlim([-0.1 0.2]);

% Customize x-ticks
xticks([-0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2]);

% Labels and title
xlabel('Time from threshold crossing (s)');
ylabel('Density');
title('Average KDE across Animals (Ctrlend, Sham, Ctrlfirst)');

% Formatting
set(gca, 'FontSize', 12, 'LineWidth', 1);
%% added with the KS desnity test 
%% Code for Comparing Ctrlend, Sham, and Ctrlfirst with KS Test
% Initialize storage for KDEs
num_points = 100; % Number of points for KDE interpolation
xi_common = linspace(-0.1, 0.2, num_points); % Common x-axis

% -----------------------------
% Compute KDE for Ctrlend Group
% -----------------------------
all_f_ctrl = []; % Store KDEs for Ctrlend
excitam_ctrl = fieldnames(ExcitRR_Ctrlend);
num_ctrl = length(excitam_ctrl);
all_data_ctrl = [];

for m = 1:num_ctrl
    allbins_ctrl = ExcitRR_Ctrlend.(excitam_ctrl{m}).AllBins;
    all_data_ctrl = [all_data_ctrl; allbins_ctrl];  % Collect for KS test
    [f_ctrl, xi_ctrl] = ksdensity(allbins_ctrl, 'Bandwidth', 0.008);
    f_interp_ctrl = interp1(xi_ctrl, f_ctrl, xi_common, 'linear', 0);
    all_f_ctrl = [all_f_ctrl; f_interp_ctrl];
end

% -----------------------------
% Compute KDE for Sham Group
% -----------------------------
all_f_sham = []; % Store KDEs for Sham
excitam_sham = fieldnames(ExcitRR_Sham);
num_sham = length(excitam_sham);
all_data_sham = [];

for m = 1:num_sham
    allbins_sham = ExcitRR_Sham.(excitam_sham{m}).AllBins;
    all_data_sham = [all_data_sham; allbins_sham];  % Collect for KS test
    [f_sham, xi_sham] = ksdensity(allbins_sham, 'Bandwidth', 0.008);
    f_interp_sham = interp1(xi_sham, f_sham, xi_common, 'linear', 0);
    all_f_sham = [all_f_sham; f_interp_sham];
end

% -----------------------------
% Compute KDE for Ctrlfirst Group (ExcitRR_Ctrl1)
% -----------------------------
all_f_ctrlfirst = []; % Store KDEs for Ctrlfirst
excitam_ctrlfirst = fieldnames(ExcitRR_Ctrl1);
num_ctrlfirst = length(excitam_ctrlfirst);
all_data_ctrlfirst = [];

for m = 1:num_ctrlfirst
    allbins_ctrlfirst = ExcitRR_Ctrl1.(excitam_ctrlfirst{m}).AllBins;
    all_data_ctrlfirst = [all_data_ctrlfirst; allbins_ctrlfirst];  % Collect for KS test
    [f_ctrlfirst, xi_ctrlfirst] = ksdensity(allbins_ctrlfirst, 'Bandwidth', 0.008);
    f_interp_ctrlfirst = interp1(xi_ctrlfirst, f_ctrlfirst, xi_common, 'linear', 0);
    all_f_ctrlfirst = [all_f_ctrlfirst; f_interp_ctrlfirst];
end

% -----------------------------
% Run KS Test for comparisons
% -----------------------------
[h_ctrl_sham, p_ctrl_sham, ks_stat_ctrl_sham] = kstest2(all_data_ctrl, all_data_sham);
[h_ctrl_ctrlfirst, p_ctrl_ctrlfirst, ks_stat_ctrl_ctrlfirst] = kstest2(all_data_ctrl, all_data_ctrlfirst);
[h_sham_ctrlfirst, p_sham_ctrlfirst, ks_stat_sham_ctrlfirst] = kstest2(all_data_sham, all_data_ctrlfirst);

% Display KS Test results in command window
fprintf('KS test between Ctrlend and Sham:\n');
fprintf('  h = %d (1 = reject null hypothesis)\n', h_ctrl_sham);
fprintf('  p = %.4f\n', p_ctrl_sham);
fprintf('  KS statistic = %.4f\n\n', ks_stat_ctrl_sham);

fprintf('KS test between Ctrlend and Ctrlfirst:\n');
fprintf('  h = %d (1 = reject null hypothesis)\n', h_ctrl_ctrlfirst);
fprintf('  p = %.4f\n', p_ctrl_ctrlfirst);
fprintf('  KS statistic = %.4f\n\n', ks_stat_ctrl_ctrlfirst);

fprintf('KS test between Sham and Ctrlfirst:\n');
fprintf('  h = %d (1 = reject null hypothesis)\n', h_sham_ctrlfirst);
fprintf('  p = %.4f\n', p_sham_ctrlfirst);
fprintf('  KS statistic = %.4f\n', ks_stat_sham_ctrlfirst);

% -----------------------------
% Plot the KDEs with Full Area Under the Curve Filled
% -----------------------------
figure; hold on;

% Fill areas for the groups
fill([xi_common, fliplr(xi_common)], [mean(all_f_sham, 1), fliplr(zeros(1, num_points))], ...
    [0.5, 0.5, 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Sham (gray)
fill([xi_common, fliplr(xi_common)], [mean(all_f_ctrl, 1), fliplr(zeros(1, num_points))], ...
    [0, 0, 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Ctrlend (black)
fill([xi_common, fliplr(xi_common)], [mean(all_f_ctrlfirst, 1), fliplr(zeros(1, num_points))], ...
    [0, 0, 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Ctrlfirst (blue)

% Plot mean KDE lines
plot(xi_common, mean(all_f_sham, 1), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 3); % Sham (gray)
plot(xi_common, mean(all_f_ctrl, 1), 'k', 'LineWidth', 3); % Ctrlend (black)
plot(xi_common, mean(all_f_ctrlfirst, 1), 'b', 'LineWidth', 3); % Ctrlfirst (blue)

% Set axis limits
xlim([-0.1 0.2]);

% Customize x-ticks
xticks([-0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2]);

% Labels and title
xlabel('Time from threshold crossing (s)');
ylabel('Density');
title('Average KDE across Animals (Ctrlend, Sham, Ctrlfirst)');

% Annotate KS test results on plot
ymax = max([mean(all_f_sham, 1), mean(all_f_ctrl, 1), mean(all_f_ctrlfirst, 1)], [], 'all');
ks_text_1 = sprintf('KS (Ctrlend vs Sham): p = %.4f, D = %.3f', p_ctrl_sham, ks_stat_ctrl_sham);
ks_text_2 = sprintf('KS (Ctrlend vs Ctrlfirst): p = %.4f, D = %.3f', p_ctrl_ctrlfirst, ks_stat_ctrl_ctrlfirst);
ks_text_3 = sprintf('KS (Sham vs Ctrlfirst): p = %.4f, D = %.3f', p_sham_ctrlfirst, ks_stat_sham_ctrlfirst);

% Display the KS test text on the plot
text(0.05, ymax * 0.9, ks_text_1, 'FontSize', 12, 'FontWeight', 'bold');
text(0.05, ymax * 0.85, ks_text_2, 'FontSize', 12, 'FontWeight', 'bold');
text(0.05, ymax * 0.8, ks_text_3, 'FontSize', 12, 'FontWeight', 'bold');

% Formatting
set(gca, 'FontSize', 12, 'LineWidth', 1);

%%
% Now plot ksdensity for all three datasets
num_bins = 20;  % Consistent number of bins
firstcolor = rgb('Black');
lastcolor = rgb('OrangeRed');
shamcolor = rgb('DimGray');

% Define common bin edges
min_val = min([BigRR_firstAll; BigRR_lastAll; BigRR_Sham]);  % Find the minimum value across all datasets
max_val = max([BigRR_firstAll; BigRR_lastAll; BigRR_Sham]);  % Find the maximum value across all datasets
edges = linspace(min_val, max_val, num_bins + 1);  % Define common bin edges

% Compute kernel density estimates
[f_ctrl, xi_ctrl] = ksdensity(BigRR_firstAll, 'Bandwidth', 0.008);
[f_stim, xi_stim] = ksdensity(BigRR_lastAll, 'Bandwidth', 0.008);
[f_sham, xi_sham] = ksdensity(BigRR_Sham, 'Bandwidth', 0.008);

% Create figure
figure; hold on;

% Plot the filled ksdensity estimates
fill(xi_ctrl, f_ctrl, firstcolor, 'FaceAlpha', 0.75, 'EdgeColor', 'none');  % Plot first data (Ctrl)
fill(xi_stim, f_stim, lastcolor, 'FaceAlpha', 0.50, 'EdgeColor', 'none');  % Plot last data (Stim)
fill(xi_sham, f_sham, shamcolor, 'FaceAlpha', 0.50, 'EdgeColor', 'none');  % Plot Sham data

% Set axis limits
xlim([-0.1 0.2]);

% Customize x-ticks
xticks([-0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2]);

set(gca, 'FontSize', 12);

% Add labels and title
xlabel('Time from threshold crossing (s)');
ylabel('Density');
title('Filled ksdensity plot');

% Ensure x-ticks are set properly
ax = gca;
ax.XTick = [-0.1 -0.05 0 0.05 0.1 0.15 0.2];  % Manually set X-ticks if needed
ax.XTickLabelRotation = 0;  % Optional: Rotate labels if needed

% Increase tick mark size and make them thicker and black
ax.TickLength = [0.02, 0.02];  % Make the ticks longer
ax.XAxis.LineWidth = 1;  % Thicker line for the X-axis
ax.YAxis.LineWidth = 1;  % Thicker line for the Y-axis

% Set the color of tick marks to black
ax.XAxis.TickDirection = 'in';  % Make ticks point inwards
ax.YAxis.TickDirection = 'in';  % Make ticks point inwards

% Make tick labels larger
ax.XAxis.FontSize = 12;  % Increase font size for X-tick labels
ax.YAxis.FontSize = 12;  % Increase font size for Y-tick labels
%% Sham vs Ctrl
% Initialize storage for KDEs
num_points = 100; % Number of points for KDE interpolation
xi_common = linspace(-0.1, 0.2, num_points); % Common x-axis for KDEs

% -----------------------------
% Compute KDE for T0_divide StimRRAll
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
% Compute KDE for ExcitRR_Sham AllBins
% -----------------------------
all_f_Sham = []; % Store KDEs for ExcitRR_Sham
excitam_Sham = fieldnames(ExcitRR_Sham);
num_Sham = length(excitam_Sham);

for m = 1:num_Sham
    allbins_Sham = ExcitRR_Sham.(excitam_Sham{m}).AllBins;  % Extract data
    [f_Sham, xi_Sham] = ksdensity(allbins_Sham, 'Bandwidth', 0.008);  % Compute KDE
    f_interp_Sham = interp1(xi_Sham, f_Sham, xi_common, 'linear', 0);  % Interpolate KDE
    all_f_Sham = [all_f_Sham; f_interp_Sham];  % Store interpolated KDE
end

% -----------------------------
% Plot the KDEs with Full Area Under the Curve Filled
% -----------------------------
shamcolor=rgb('DimGray');
ctrlcolor=rgb('Black');
figure; hold on;

% Fill area for T0_divide group (red)
fill([xi_common, fliplr(xi_common)], [mean(all_f_T0, 1), fliplr(zeros(1, num_points))], ...
    [0, 0, 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Red fill

% Fill area for ExcitRR_Sham group (blue)
fill([xi_common, fliplr(xi_common)], [mean(all_f_Sham, 1), fliplr(zeros(1, num_points))], ...
    [0.4102, 0.4102, 0.4102], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Blue fill

% Plot mean KDE lines (only lines, no SEM)
plot(xi_common, mean(all_f_T0, 1), 'k', 'LineWidth', 3); % T0_divide (red)
plot(xi_common, mean(all_f_Sham, 1),'Color',shamcolor, 'LineWidth', 3); % ExcitRR_Sham (blue)

% Set axis limits
xlim([-0.1 0.2]);

% Customize x-ticks
xticks([-0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2]);

% Labels and title
xlabel('Time from threshold crossing (s)');
ylabel('Density');
title('Average KDE across Animals (T0_divide vs ExcitRR_Sham)');

% Formatting
set(gca, 'FontSize', 12, 'LineWidth', 1);
%% Just sham
% Initialize storage for KDEs
num_points = 100; % Number of points for KDE interpolation
xi_common = linspace(-0.1, 0.2, num_points); % Common x-axis for KDEs

% -----------------------------
% Compute KDE for ExcitRR_Sham AllBins
% -----------------------------
all_f_Sham = []; % Store KDEs for ExcitRR_Sham
excitam_Sham = fieldnames(ExcitRR_Sham);
num_Sham = length(excitam_Sham);

for m = 1:num_Sham
    allbins_Sham = ExcitRR_Sham.(excitam_Sham{m}).AllBins;  % Extract data
    [f_Sham, xi_Sham] = ksdensity(allbins_Sham, 'Bandwidth', 0.008);  % Compute KDE
    f_interp_Sham = interp1(xi_Sham, f_Sham, xi_common, 'linear', 0);  % Interpolate KDE
    all_f_Sham = [all_f_Sham; f_interp_Sham];  % Store interpolated KDE
end

% -----------------------------
% Plot the KDE with Full Area Under the Curve Filled for Sham
% -----------------------------
shamcolor = rgb('DimGray');
figure; hold on;

% Fill area for ExcitRR_Sham group (gray)
fill([xi_common, fliplr(xi_common)], [mean(all_f_Sham, 1), fliplr(zeros(1, num_points))], ...
    shamcolor, 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Gray fill for the area under the curve

% Plot mean KDE line for Sham group
plot(xi_common, mean(all_f_Sham, 1), 'Color', shamcolor, 'LineWidth', 3); % ExcitRR_Sham (gray)

% Set axis limits
xlim([-0.1 0.2]);

% Customize x-ticks
xticks([-0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2]);

% Labels and title
xlabel('Time from threshold crossing (s)');
ylabel('Density');
title('Average KDE for ExcitRR_Sham');

% Formatting
set(gca, 'FontSize', 12, 'LineWidth', 1);

%% Code for First vs. Last Stim Reaches (with KS Test)
% Initialize storage for KDEs
num_points = 100; % Number of points for KDE interpolation
xi_common = linspace(-0.1, 0.2, num_points); % Common x-axis for KDEs

% -----------------------------
% Compute KDE for ExcitRR.(...) StimRR(:,1)
% -----------------------------
all_f_T0 = []; % Store KDEs for ExcitRR
excitam_T0 = fieldnames(ExcitRR);
num_T0 = length(excitam_T0);
all_data_T0 = [];

for m = 1:num_T0
    allbins_T0 = ExcitRR.(excitam_T0{m}).StimRR(:,1);  % Extract data
    all_data_T0 = [all_data_T0; allbins_T0];           % Collect for KS test
    [f_T0, xi_T0] = ksdensity(allbins_T0, 'Bandwidth', 0.008);
    f_interp_T0 = interp1(xi_T0, f_T0, xi_common, 'linear', 0);
    all_f_T0 = [all_f_T0; f_interp_T0];
end

% -----------------------------
% Compute KDE for ExcitRR_end.(...) StimRR(:,1)
% -----------------------------
all_f_Sham = []; % Store KDEs for ExcitRR_end
excitam_Sham = fieldnames(ExcitRR_end);
num_Sham = length(excitam_Sham);
all_data_Sham = [];

for m = 1:num_Sham
    allbins_Sham = ExcitRR_end.(excitam_Sham{m}).StimRR(:,1);  % Extract data
    all_data_Sham = [all_data_Sham; allbins_Sham];             % Collect for KS test
    [f_Sham, xi_Sham] = ksdensity(allbins_Sham, 'Bandwidth', 0.008);
    f_interp_Sham = interp1(xi_Sham, f_Sham, xi_common, 'linear', 0);
    all_f_Sham = [all_f_Sham; f_interp_Sham];
end

% -----------------------------
% Run two-sample KS test
% -----------------------------
[h_ks, p_ks, ks_stat] = kstest2(all_data_T0, all_data_Sham);

% Display results in command window
fprintf('KS test result:\n');
fprintf('  h = %d (1 = reject null hypothesis)\n', h_ks);
fprintf('  p = %.4f\n', p_ks);
fprintf('  KS statistic = %.4f\n', ks_stat);

% -----------------------------
% Plot the KDEs with Full Area Under the Curve Filled
% -----------------------------
stim1color = rgb('DeepSkyBlue');
stim2color = rgb('DarkBlue');
figure; hold on;

% Fill areas
fill([xi_common, fliplr(xi_common)], [mean(all_f_T0, 1), fliplr(zeros(1, num_points))], ...
    stim1color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([xi_common, fliplr(xi_common)], [mean(all_f_Sham, 1), fliplr(zeros(1, num_points))], ...
    stim2color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Plot mean KDE lines
plot(xi_common, mean(all_f_T0, 1), 'Color', stim1color, 'LineWidth', 3);
plot(xi_common, mean(all_f_Sham, 1), 'Color', stim2color, 'LineWidth', 3);

% Axis settings
xlim([-0.1 0.2]);
xticks([-0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2]);
xlabel('Time from threshold crossing (s)');
ylabel('Density');
title('Average KDE across Animals (First vs Last Stim Reaches)');

% KS test annotation on plot
ymax = max([mean(all_f_T0, 1), mean(all_f_Sham, 1)], [], 'all');
ks_text = sprintf('KS Test: p = %.4f, D = %.3f', p_ks, ks_stat);
text(0.05, ymax * 0.9, ks_text, 'FontSize', 12, 'FontWeight', 'bold');

% Formatting
set(gca, 'FontSize', 12, 'LineWidth', 1);
