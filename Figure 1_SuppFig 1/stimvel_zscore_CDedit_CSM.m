%% T0
% Initialize arrays for p-values and Z-scores
parray0 = NaN(size(T0_big.stimvel, 1), 1);  % Store p-values
zarray0 =   NaN(size(T0_big.stimvel, 1), 1)   ;    % Store Z-scores (same shape as stimvel)
Zscores0= NaN(size(T0_big.stimvel, 1), 5);

for pnum = 1:length(parray0)
    % Extract trial data for C_big (control) and T0_big (experimental)
    a4b = squeeze(T0_bigCSM.stimvel(pnum, 8, :));  % Experimental data (T0_big)
    b4b = squeeze(CT0_bigCSM.stimvel(pnum, 8, :));   % Control data (C_big), matching trial

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

dim = 8;
% StimvelPlot(C_big, dim, 'k', 1);
% StimvelPlot(T0_big, dim, T0color, 0);
StimvelPlot(CT0_bigCSM, dim, 'k', 1);
StimvelPlot(T0_bigCSM, dim, T0color, 0);

% Plot reference line
yline(0, '--k');

% Axis limits and labels
xlim([-0.5 0.5]);
%ylim([-0.5 1.5]);
%ylim([-0.05 0.27]);
ylim([-20 20])
xlabel('Time from Pos. Thresh. (s)');
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
parray1 = NaN(size(T1_big.stimvel, 1), 1);  % Store p-values
zarray1 =   NaN(size(T1_big.stimvel, 1), 1)   ;    % Store Z-scores (same shape as stimvel)

for pnum = 1:length(parray1)
    % Extract trial data for C_big (control) and T0_big (experimental)
    a4b = squeeze(T1_bigCSM.stimvel(pnum, 8, :));  % Experimental data (T0_big)
    b4b = squeeze(CT1_bigCSM.stimvel(pnum, 8, :));   % Control data (C_big), matching trial

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

dim = 8;
StimvelPlot(CT1_bigCSM, dim, 'k', 1);
StimvelPlot(T1_bigCSM, dim, T1color, 0);

% Plot reference line
yline(0, '--k');

% Axis limits and labels
xlim([-0.5 0.5]);
%ylim([-0.5 1.5])
%ylim('auto')
ylim([-20 20]);
%ylim([-0.05 0.27])
xlabel('Time from Pos. Thresh. (s)');
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
parray2 = NaN(size(T2_big.stimvel, 1), 1);  % Store p-values
zarray2 =   NaN(size(T2_big.stimvel, 1), 1)   ;    % Store Z-scores (same shape as stimvel)

for pnum = 1:length(parray2)
    % Extract trial data for C_big (control) and T0_big (experimental)
    a4b = squeeze(T2_bigCSM.stimvel(pnum, 8, :));  % Experimental data (T0_big)
    b4b = squeeze(CT2_bigCSM.stimvel(pnum, 8, :));   % Control data (C_big), matching trial

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

dim = 8;
StimvelPlot(CT2_bigCSM, dim, 'k', 1);
StimvelPlot(T2_bigCSM, dim, T2color, 0);

% Plot reference line
yline(0, '--k');

% Axis limits and labels
xlim([-0.5 0.5]);
ylim([-20 20]);
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
parray3 = NaN(size(T3_big.stimvel, 1), 1);  % Store p-values
zarray3 =   NaN(size(T3_big.stimvel, 1), 1)   ;    % Store Z-scores (same shape as stimvel)

for pnum = 1:length(parray3)
    % Extract trial data for C_big (control) and T0_big (experimental)
    a4b = squeeze(T3_bigCSM.stimvel(pnum, 8, :));  % Experimental data (T0_big)
    b4b = squeeze(CT3_bigCSM.stimvel(pnum, 8, :));   % Control data (C_big), matching trial

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

dim = 8;
StimvelPlot(CT3_bigCSM, dim, 'k', 1);
StimvelPlot(T3_bigCSM, dim, T3color, 0);

% Plot reference line
yline(0, '--k');

% Axis limits and labels
xlim([-0.5 0.5]);
ylim([-20 20]);
%ylim([-0.5 1.5]);
%ylim([0 0.25])
xlabel('Time from Pos. Thresh. (s)');
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