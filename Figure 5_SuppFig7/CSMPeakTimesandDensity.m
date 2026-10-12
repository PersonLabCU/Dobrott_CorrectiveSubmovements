%% PLot peak times and density 
%% -----------------------------
%% Control only: Scatter plot for all four sessions

segments = {'first5','middle5','last5'};
sessions = {'1','2','3','End'}; % session names

% Combine all Control PeakInfo structs
PeakInfo_Control = {PeakInfo_Ctrl1, PeakInfo_Ctrl2, PeakInfo_Ctrl3, PeakInfo_CtrlEnd};

figure; hold on;

% Colors for sessions
session_colors = lines(4);  % MATLAB default colormap for 4 sessions

% Marker shapes for segments
segment_markers = containers.Map({'first5','middle5','last5'}, {'o','s','^'});

% Legend tracking
legend_entries = containers.Map;

for sess_idx = 1:4
    PeakInfo = PeakInfo_Control{sess_idx};
    color = session_colors(sess_idx,:);
    
    for s = 1:length(segments)
        x = PeakInfo.(segments{s}).peak_loc;
        y = PeakInfo.(segments{s}).peak_val;
        m = segment_markers(segments{s});
        
        combo = ['Session' sessions{sess_idx} '_' segments{s}];
        if ~isKey(legend_entries, combo)
            h = plot(x, y, m, 'MarkerSize', 10, 'MarkerFaceColor', color, 'MarkerEdgeColor', 'k', 'DisplayName', combo);
            legend_entries(combo) = h;
        else
            plot(x, y, m, 'MarkerSize', 10, 'MarkerFaceColor', color, 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
        end
    end
end

xlabel('Peak Time (s)');
ylabel('Peak Density');
title('Control: KDE Peaks for All Sessions');
legend('Location','northwest');
set(gca,'FontSize',12,'LineWidth',1);
grid on;
hold off;

%% -----------------------------
%% Control + Stim: Scatter plot for all four sessions
%% -----------------------------
segments = {'first5','middle5','last5'};
sessions = {'1','2','3','End'};

% Control and Stim structs
PeakInfo_Control = {PeakInfo_Ctrl1, PeakInfo_Ctrl2, PeakInfo_Ctrl3, PeakInfo_CtrlEnd};
PeakInfo_Stim    = {PeakInfo_Stim1, PeakInfo_Stim2, PeakInfo_Stim3, PeakInfo_StimEnd};

figure; hold on;

% Colors: Control vs Stim
group_colors = containers.Map({'Control','Stim'}, {[0 0.4470 0.7410],[0.8500 0.3250 0.0980]}); % blue vs orange
% Marker shapes for segments
segment_markers = containers.Map({'first5','middle5','last5'}, {'o','s','^'});

% Legend tracking
legend_entries = containers.Map;

% Plot Control
for sess_idx = 1:4
    PeakInfo = PeakInfo_Control{sess_idx};
    for s = 1:length(segments)
        x = PeakInfo.(segments{s}).peak_loc;
        y = PeakInfo.(segments{s}).peak_val;
        m = segment_markers(segments{s});
        c = group_colors('Control');
        
        combo = ['Ctrl_S' sessions{sess_idx} '_' segments{s}];
        if ~isKey(legend_entries, combo)
            h = plot(x, y, m, 'MarkerSize', 10, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'DisplayName', combo);
            legend_entries(combo) = h;
        else
            plot(x, y, m, 'MarkerSize', 10, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
        end
    end
end

% Plot Stim
for sess_idx = 1:4
    PeakInfo = PeakInfo_Stim{sess_idx};
    for s = 1:length(segments)
        x = PeakInfo.(segments{s}).peak_loc;
        y = PeakInfo.(segments{s}).peak_val;
        m = segment_markers(segments{s});
        c = group_colors('Stim');
        
        combo = ['Stim_S' sessions{sess_idx} '_' segments{s}];
        if ~isKey(legend_entries, combo)
            h = plot(x, y, m, 'MarkerSize', 10, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'DisplayName', combo);
            legend_entries(combo) = h;
        else
            plot(x, y, m, 'MarkerSize', 10, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
        end
    end
end

xlabel('Peak Time (s)');
ylabel('Peak Density');
title('KDE Peaks: Control and Stim Across All Sessions');
legend('Location','northwest');
set(gca,'FontSize',12,'LineWidth',1);
grid on;
hold off;

segments = {'first5','middle5','last5'};
sessions = {'1','2','3','End'};
%%
% Control and Stim structs
PeakInfo_Control = {PeakInfo_Ctrl1, PeakInfo_Ctrl2, PeakInfo_Ctrl3, PeakInfo_CtrlEnd};
PeakInfo_Stim    = {PeakInfo_Stim1, PeakInfo_Stim2, PeakInfo_Stim3, PeakInfo_StimEnd};

figure; hold on;

% Colors for sessions (same for Control and Stim, can adjust shades if desired)
session_colors = lines(4); % 4 different colors for sessions

% Marker shapes for segments
segment_markers = containers.Map({'first5','middle5','last5'}, {'o','s','^'});

% Legend tracking
legend_entries = containers.Map;

% --- Plot Control ---
for sess_idx = 1:4
    PeakInfo = PeakInfo_Control{sess_idx};
    c = session_colors(sess_idx,:); % color by session
    
    for s = 1:length(segments)
        x = PeakInfo.(segments{s}).peak_loc;
        y = PeakInfo.(segments{s}).peak_val;
        m = segment_markers(segments{s});
        
        combo = ['Ctrl_S' sessions{sess_idx} '_' segments{s}];
        if ~isKey(legend_entries, combo)
            h = plot(x, y, m, 'MarkerSize', 10, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'DisplayName', combo);
            legend_entries(combo) = h;
        else
            plot(x, y, m, 'MarkerSize', 10, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
        end
    end
end

% --- Plot Stim ---
for sess_idx = 1:4
    PeakInfo = PeakInfo_Stim{sess_idx};
    c = session_colors(sess_idx,:); % same session color for Stim
    
    for s = 1:length(segments)
        x = PeakInfo.(segments{s}).peak_loc;
        y = PeakInfo.(segments{s}).peak_val;
        m = segment_markers(segments{s});
        
        combo = ['Stim_S' sessions{sess_idx} '_' segments{s}];
        if ~isKey(legend_entries, combo)
            h = plot(x, y, m, 'MarkerSize', 10, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'DisplayName', combo);
            legend_entries(combo) = h;
        else
            plot(x, y, m, 'MarkerSize', 10, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
        end
    end
end

xlabel('Peak Time (s)');
ylabel('Peak Density');
xlim([0.04 0.1])
ylim([0 40])
title('KDE Peaks: Control and Stim Across Sessions');
legend('Location','northwest');
set(gca,'FontSize',12,'LineWidth',1);
grid on;
hold off;

%%
% Control and Stim structs
PeakInfo_Control = {PeakInfo_Ctrl1, PeakInfo_Ctrl2, PeakInfo_Ctrl3, PeakInfo_CtrlEnd};
PeakInfo_Stim    = {PeakInfo_Stim1, PeakInfo_Stim2, PeakInfo_Stim3, PeakInfo_StimEnd};

figure; hold on;

% Colors for sessions (same for Control and Stim)
session_colors = lines(4); % 4 different colors for sessions

% Marker shapes for segments
segment_markers = containers.Map({'first5','middle5','last5'}, {'o','s','^'});

% Legend tracking
legend_entries = containers.Map;

% --- Plot Control ---
for sess_idx = 1:4
    PeakInfo = PeakInfo_Control{sess_idx};
    c = session_colors(sess_idx,:); % color by session
    
    for s = 1:length(segments)
        x = PeakInfo.(segments{s}).peak_loc;
        y = PeakInfo.(segments{s}).peak_val;
        m = segment_markers(segments{s});
        
        combo = ['Ctrl_S' sessions{sess_idx} '_' segments{s}];
        if ~isKey(legend_entries, combo)
            h = plot(x, y, m, 'MarkerSize', 10, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'DisplayName', combo);
            legend_entries(combo) = h;
        else
            plot(x, y, m, 'MarkerSize', 10, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
        end
    end
end

% % --- Plot Stim ---
% for sess_idx = 1:4
%     PeakInfo = PeakInfo_Stim{sess_idx};
%     c = session_colors(sess_idx,:); % same session color for Stim
% 
%     for s = 1:length(segments)
%         x = PeakInfo.(segments{s}).peak_loc;
%         y = PeakInfo.(segments{s}).peak_val;
%         m = segment_markers(segments{s});
% 
%         combo = ['Stim_S' sessions{sess_idx} '_' segments{s}];
%         if ~isKey(legend_entries, combo)
%             h = plot(x, y, m, 'MarkerSize', 10, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'DisplayName', combo);
%             legend_entries(combo) = h;
%         else
%             plot(x, y, m, 'MarkerSize', 10, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
%         end
%     end
% end

xlabel('Peak Time (s)');
ylabel('Peak Density');
xlim([0.04 0.1])
ylim([0 20])
title('KDE Peaks: Control and Stim Across Sessions');

% Move legend outside the plot
legend('Location','eastoutside');  % puts legend to the right
set(gca,'FontSize',12,'LineWidth',1);
grid on;
hold off;

%%
%%
segments = {'first5','last5'}; % only first and last
sessions = {'1','2','3','End'};

% Control and Stim structs
PeakInfo_Control = {PeakInfo_Ctrl1, PeakInfo_Ctrl2, PeakInfo_Ctrl3, PeakInfo_CtrlEnd};
PeakInfo_Stim    = {PeakInfo_Stim1, PeakInfo_Stim2, PeakInfo_Stim3, PeakInfo_StimEnd};

figure; hold on;

% Colors for sessions (same for Control and Stim)
session_colors = lines(4); % 4 colors for 4 sessions

% Marker shapes for segments
segment_markers = containers.Map({'first5','last5'}, {'o','^'});

% --- Plot Control ---
for sess_idx = 1:4
    PeakInfo = PeakInfo_Control{sess_idx};
    c = session_colors(sess_idx,:); % color by session
    
    x_vals = zeros(1,2);
    y_vals = zeros(1,2);

    for s = 1:length(segments)
        x = PeakInfo.(segments{s}).peak_loc;
        y = PeakInfo.(segments{s}).peak_val;
        m = segment_markers(segments{s});
        
        % Store for connecting line
        x_vals(s) = x;
        y_vals(s) = y;
        
        % Plot point
        plot(x, y, m, 'MarkerSize', 10, ...
            'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', ...
            'DisplayName', ['Ctrl_S' sessions{sess_idx} '_' segments{s}]);
    end

    % Connect first5 → last5
    plot(x_vals, y_vals, '-', 'Color', c, 'LineWidth', 1.5, 'HandleVisibility','off');
end

xlabel('Peak Time (s)');
ylabel('Peak Density');
xlim([0.04 0.1]);
ylim([0 20]);
title('KDE Peaks: Control Across Sessions (First5 vs Last5)');
legend('Location','eastoutside');
set(gca, 'FontSize',12,'LineWidth',1);
grid on;
hold off;
