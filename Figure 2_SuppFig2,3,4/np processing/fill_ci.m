%% Function to Plot Confidence Interval Shading
function fill_ci(time_axis, ci, color)
    fill([time_axis, fliplr(time_axis)], [ci(1, :), fliplr(ci(2, :))], color, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
end