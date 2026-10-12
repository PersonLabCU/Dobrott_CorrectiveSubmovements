
% --- Helper function to format subplots ---
function formatPlot()
    set(gca, 'YDir', 'normal');
    hold on;
    ylims = ylim;
    line([1 1], ylims, 'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 1.2); % vertical line at 1 ms
    rectangle('Position', [50, ylims(1), 50, ylims(2)-ylims(1)], 'EdgeColor', 'r', 'LineWidth', 1.5, 'LineStyle', '--'); % box from 50–100
    line([0 0], ylims, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2); % vertical line at 0 ms
end