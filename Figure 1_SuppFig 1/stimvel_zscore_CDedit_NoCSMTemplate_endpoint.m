%% Replace all reaches aligned CSM with only the NoCSM template

an=fieldnames(ExcitRR_Sham);

thresh=0.16;
dirdir= 1;

ExcitRR_Ctrl= ExcitRR_Sham;

for a=1:length(an)
   
    CtrlReaches=ExcitRR_Ctrl.(an{a}).NoReReach;
    T0Reaches=ExcitRR_Ctrl.(an{a}).RRT0;
    T1Reaches=ExcitRR_Ctrl.(an{a}).RRT1;
    T2Reaches=ExcitRR_Ctrl.(an{a}).RRT2;
    T3Reaches=ExcitRR_Ctrl.(an{a}).RRT3;

    % NoCSM aligned to endpoint 
    %NoCSM_end=AlignEndpoint_max1(CtrlReaches,500);
    NoCSM_end=AlignEndpoint(CtrlReaches,500);
    NoCSM_end=stimvel_accel(NoCSM_end);

    %T0 aligned to endpoint
    %T0_end=AlignEndpoint_max1(T0Reaches,500);
    T0_end=AlignEndpoint(T0Reaches,500);
    T0_end=stimvel_accel(T0_end);

    %T1 aligned to endpoint
    %T1_end=AlignEndpoint_max1(T1Reaches,500);
    T1_end=AlignEndpoint(T1Reaches,500);
    T1_end=stimvel_accel(T1_end);

    %T2 aligned to endpoint
    %T2_end=AlignEndpoint_max1(T2Reaches,500);
    T2_end=AlignEndpoint(T2Reaches,500);
    T2_end=stimvel_accel(T2_end);

    %T3 aligned to endpoint
    %T3_end=AlignEndpoint_max1(T3Reaches,500);
    T3_end=AlignEndpoint(T3Reaches,500);
    T3_end=stimvel_accel(T3_end);

    % Save to big structure

    %Stimvel

    ExcitRR_Ctrl.(an{a}).NoCSM_end=NoCSM_end;
    ExcitRR_Ctrl.(an{a}).T0_end=T0_end;
    ExcitRR_Ctrl.(an{a}).T1_end=T1_end;
    ExcitRR_Ctrl.(an{a}).T2_end=T2_end;
    ExcitRR_Ctrl.(an{a}).T3_end=T3_end;


end 

%% Now make big stimvel structures for CSM bins aligned to average CSM time for each bin 
NoCSM_bigCSM=struct;
T0_bigCSM=struct;
T1_bigCSM=struct;
T2_bigCSM= struct;
T3_bigCSM=struct;

ar=fieldnames(ExcitRR_Ctrl);
for a = 1:length(ar) 
   NoCSM_bigCSM.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).NoCSM_end.Average;
   T0_bigCSM.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).T0_end.Average;
   T1_bigCSM.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).T1_end.Average;
   T2_bigCSM.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).T2_end.Average;
   T3_bigCSM.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).T3_end.Average;

    
end

NoCSM_bigCSM.Average = nanmean(NoCSM_bigCSM.stimvel,3);
NoCSM_bigCSM.Std = nanstd(NoCSM_bigCSM.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(NoCSM_bigCSM.Average),1); 
for nnum = 1:length(NoCSM_bigCSM.Average)
    samp = squeeze(NoCSM_bigCSM.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
NoCSM_bigCSM.N = N;


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




%% Plot stimvel for CSM

T0color=rgb('Thistle');
T1color=rgb('MediumOrchid');
T2color=rgb('SlateBlue');
T3color=rgb('DarkSlateBlue');
NRRcolor=rgb('Black');
stimcolor=rgb('DodgerBlue');

dim =2;

StimvelPlot(T0_bigCSM,dim,T0color,1)
StimvelPlot(T1_bigCSM,dim,T1color,0)
StimvelPlot(T2_bigCSM,dim,T2color,0)
StimvelPlot(T3_bigCSM,dim,T3color,0)
StimvelPlot(NoCSM_bigCSM,dim,NRRcolor,0)

%xlim([-0.25 0.25])
yline(0.267,'--') % Average NoCSM template endpint 
ylim([0.1 0.3])
xlim([-0.2 0.1])
%xticks(-0.2:0.01:0.2)
ylim('auto')
%ylim([-0.6 1.4])
%ylim([0 0.25])
%ylim([-20 20])
pbaspect([1 2 1])
%ylim([0 0.3])
xline(0,'--')
yline(0,'--')


xlabel('Time from CSM(s)')
ylabel('Outward Position (cm)')
%ylabel('Outward Velocity (cm/s)')
%ylabel('Outward acceleration')





%% Narrower plot 
% Plot stimvel for endpoint

figure('Color','w','Position',[300 300 350 300]); % <<< narrower width

T0color=rgb('Thistle');
T1color=rgb('MediumOrchid');
T2color=rgb('SlateBlue');
T3color=rgb('DarkSlateBlue');
NRRcolor=rgb('Black');
stimcolor=rgb('DodgerBlue');

dim = 2;

hold on
StimvelPlot(T0_bigCSM,dim,T0color,0)
StimvelPlot(T1_bigCSM,dim,T1color,0)
StimvelPlot(T2_bigCSM,dim,T2color,0)
StimvelPlot(T3_bigCSM,dim,T3color,0)
StimvelPlot(NoCSM_bigCSM,dim,NRRcolor,0)

% --- axis limits ---
xlim([-0.2 0.1])
ylim([0.1 0.3])

% --- reference lines ---
yline(0.267,'--') % Average NoCSM template endpoint
xline(0,'--')
yline(0,'--')

% --- labels ---
xlabel('Time from CSM(s)')
ylabel('Outward Position (cm)')

% optional: tighter margins
set(gca,'LooseInset',max(get(gca,'TightInset'),0.02))
box off
%% T1

% Initialize arrays for p-values and Z-scores
parray1 = NaN(size(T1_bigCSM.stimvel, 1), 1);  % Store p-values
zarray1 =   NaN(size(T1_bigCSM.stimvel, 1), 1)   ;    % Store Z-scores (same shape as stimvel)

for pnum = 1:length(parray1)
    % Extract trial data for C_big (control) and T0_big (experimental)
    a4b = squeeze(T1_bigCSM.stimvel(pnum, 2, :));  % Experimental data (T0_big)
    b4b = squeeze(NoCSM_bigCSM.stimvel(pnum, 2, :));   % Control data (C_big), matching trial

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
StimvelPlot(NoCSM_bigCSM, dim, 'k', 1);
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
    b4b = squeeze(NoCSM_bigCSM.stimvel(pnum, 2, :));   % Control data (C_big), matching trial

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
StimvelPlot(NoCSM_bigCSM, dim, 'k', 1);
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