
%% Need to have run the NoeffecRRtiming script to separate out reaches 
T0_divide=[];
thresh=0.1;
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
%Last session
T0_dividestim=[];
thresh=0.1;
for a=1:length(an)
    t0r=ExcitRR.(an{a}).RRT0;
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
    T0_dividestim.(an{a}).ReReach=ReReach;
    T0_dividestim.(an{a}).StimRR=StimRR;
    T0_dividestim.(an{a}).NoJerk=ReachesWithoutSustainedJerk;
   
end
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
%% For the Stim reaches 
excitam = fieldnames(ExcitRR);
for b=1:length(excitam)
    RR_Excit=ExcitRR.(excitam{b}).StimRR;
 
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
        ExcitRR.(excitam{b}).StimRREdit=nonZeroRows_ex;
end

%% 
BigRR_ctrl=[];

excitam=fieldnames(ExcitRR_Ctrl);
for m =1:length(excitam)
    allbins=[];
    excitall = ExcitRR_Ctrl.(excitam{m}).StimRREdit(:,1);
    t0bins=T0_divide.(excitam{m}).StimRR(:,1);
    allbins=[excitall;t0bins];
    ExcitRR_Ctrl.(excitam{m}).AllBins=allbins;
    BigRR_ctrl=[BigRR_ctrl;excitall];
    
end
%% 
BigRR_stim=[];

excitam = fieldnames(ExcitRR);
for m = 1:length(excitam)
    allbins = [];
    
    % Extract data
    excitall = ExcitRR.(excitam{m}).StimRREdit(:,1);
    
    % Check if t0bins exists and is non-empty
    if isfield(T0_dividestim, excitam{m}) && ~isempty(T0_dividestim.(excitam{m}).StimRR)
        t0bins = T0_dividestim.(excitam{m}).StimRR(:,1);
    else
        t0bins = []; % Empty case
    end

    % Concatenate filtered data
    allbins = [excitall; t0bins];

    % Store results
    ExcitRR.(excitam{m}).AllBins = allbins;
    BigRR_stim = [BigRR_stim; excitall];
end
%% 
% Define colors
stimcolor = rgb('Teal'); % ExcitRR
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

for m = 1:num_excitrr
    allbins_excitrr = ExcitRR.(excitam_excitrr{m}).AllBins;
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

for m = 1:num_excitrr_ctrl
    allbins_excitrr_ctrl = ExcitRR_Ctrl.(excitam_excitrr_ctrl{m}).AllBins;
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

for m = 1:num_sham
    allbins_sham = ExcitRR_Sham.(excitam_sham{m}).AllBins;
    [f_sham, xi_sham] = ksdensity(allbins_sham, 'Bandwidth', 0.008);
    f_interp_sham = interp1(xi_sham, f_sham, xi_common, 'linear', 0);
    all_f_sham = [all_f_sham; f_interp_sham];
end

%% Plot All Groups Together

% Define colors
stimcolor = rgb('Teal'); % ExcitRR
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

for m = 1:num_excitrr
    allbins_excitrr = ExcitRR.(excitam_excitrr{m}).AllBins;
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

for m = 1:num_excitrr_ctrl
    allbins_excitrr_ctrl = ExcitRR_Ctrl.(excitam_excitrr_ctrl{m}).AllBins;
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

for m = 1:num_sham
    allbins_sham = ExcitRR_Sham.(excitam_sham{m}).AllBins;
    [f_sham, xi_sham] = ksdensity(allbins_sham, 'Bandwidth', 0.008);
    f_interp_sham = interp1(xi_sham, f_sham, xi_common, 'linear', 0);
    all_f_sham = [all_f_sham; f_interp_sham];
end

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

% Formatting
set(gca, 'FontSize', 12, 'LineWidth', 1);

%legend({'ExcitRR (Teal)', 'ExcitRR_Ctrl (LightSlateGray)', 'ExcitRR_Sham (Dashed)'}, 'Location', 'NorthEast');

hold off;

%% edited to add KS tests 
%% Plot All Groups Together with KS Test

% Define colors
stimcolor = rgb('Teal'); % ExcitRR
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
    allbins_excitrr = ExcitRR.(excitam_excitrr{m}).AllBins;
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
    allbins_excitrr_ctrl = ExcitRR_Ctrl.(excitam_excitrr_ctrl{m}).AllBins;
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
    allbins_sham = ExcitRR_Sham.(excitam_sham{m}).AllBins;
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
