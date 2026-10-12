%% Plot example reachstart and show where stim threshold is crossed 
ExcitRR_Ctrl=ExcitRR_Sham;
an=fieldnames(ExcitRR_Ctrl);


kinc=ExcitRR_Ctrl.EX010.RRT1;
fieldnames(kinc);
%%
kind=kinc.R447.kin;
specific_x=0.16;

endpt=max(kind(:,2));

col = kind(:,1) - kind(1,1);

%%

endptx=max(kind(:,2));
endpty=max(kind(:,3));
col = kind(:,1) - kind(1,1); % Define col based on your data

cmap=colormap(flipud(slanCM('gem',length(col))));
% Plot each segment with the corresponding color
hold on
for i = 1:length(col)-1
    plot(kind(i:i+1,2), kind(i:i+1,3), 'LineWidth', 4, 'Color', cmap(i,:));
end

% Highlight points based on specific_y
h = findobj(gca, 'Type', 'Line');
xData = kind(:,2);
yData = kind(:,3);
tolerance = 1e-3;
idx = abs(xData - specific_x) < tolerance;
y_coords = yData(idx);

p = plot(specific_x * ones(size(y_coords)),y_coords, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'red');
xline(0.16, '--');
xline(endptx)
yline(endpty)

% Adjust x-axis ticks
% xticks = get(gca, 'XTick');
% %new_xticks = xticks - kind(1,1);
% set(gca, 'XTickLabel', new_xticks);

% Set aspect ratio and labels
ylim([0 0.25])
pbaspect([2 1.2 1]);
xlabel('Outward Position');
ylabel('Upward Position');

colorbar; % Add a colorbar to show the mapping of col
%% Different method for plotting 
% Preprocess data: Only include points from x >= 0
start_idx = find(kind(:,2) >= 0, 1, 'first');
kind = kind(start_idx:end, :);

% Extract position and time (or any variable to color by)
x = kind(:,2);               % Outward position
y = kind(:,3);               % Upward position
col = kind(:,1) - kind(1,1); % Time relative to starting point

% Normalize col to [0, 1]
col_norm = (col - min(col)) / (max(col) - min(col));

% Prepare data for surface (line with color gradient)
xLine = [x(1:end-1)'; x(2:end)'];
yLine = [y(1:end-1)'; y(2:end)'];
zLine = zeros(size(xLine)); % Keep z constant for 2D plotting
cLine = [col_norm(1:end-1)'; col_norm(2:end)']; % Colors for each segment

% Plot
figure;
hold on;
%colormap(flipud(slanCM('dense', 256))); % Or use your preferred colormap
%colormap(flipud(slanCM('tempo', 256))); % Or use your preferred colormap
colormap(flipud(slanCM('deep', 256))); % Or use your preferred colormap
surface(xLine, yLine, zLine, cLine, ...
        'EdgeColor', 'interp', ...
        'FaceColor', 'none', ...
        'LineWidth', 4);

colorbar;

% Highlight specific x location
tolerance = 1e-3;
idx = abs(x - specific_x) < tolerance;
y_coords = y(idx);
plot(specific_x * ones(size(y_coords)), y_coords, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'red');

% Reference lines
xline(0.16, '--');
xline(max(x));
yline(max(y));

% Axis formatting
ylim([0 0.25]);
xlim([0 0.3])
%pbaspect([2 1.2 1]);
xlabel('Outward Position');
ylabel('Upward Position');
title('Time-Colored Reach Trajectory');



%%

% endptx=max(kind(:,2));
% endpty=max(kind(:,4));
col = kind(:,1) - kind(1,1); % Define col based on your data
%cmap = colormap(cool(length(col))); % Create a colormap with the same length as col

cmap=colormap(flipud(slanCM('dense',length(col))));
% Plot each segment with the corresponding color
hold on
for i = 1:length(col)-1
    plot(kind(i:i+1,2), kind(i:i+1,4), 'LineWidth', 4, 'Color', cmap(i,:));
end

% Highlight points based on specific_y
h = findobj(gca, 'Type', 'Line');
xData = kind(:,2);
yData = kind(:,4);
tolerance = 1e-3;
idx = abs(xData - specific_x) < tolerance;
y_coords = yData(idx);

p = plot(specific_x * ones(size(y_coords)),y_coords, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'red');
xline(0.16, '--');
xline(endptx)
yline(endpty)

% Adjust x-axis ticks
% xticks = get(gca, 'XTick');
% %new_xticks = xticks - kind(1,1);
% set(gca, 'XTickLabel', new_xticks);

% Set aspect ratio and labels
%ylim([0 0.25])
ylim([0 0.1])
pbaspect([2 0.75 1]);
xlabel('Outward Position');
ylabel('Upward Position');

colorbar; % Add a colorbar to show the mapping of col
%%

x = kind(:,2);
y = kind(:,4);
col_dense = linspace(min(col), max(col), 1000); % Denser col values
x_dense = interp1(col, x, col_dense);
y_dense = interp1(col, y, col_dense);

cmap = colormap(flipud(slanCM('dense', length(col_dense)))); % Adjust for new length
hold on
for i = 1:length(col_dense)-1
    plot(x_dense(i:i+1), y_dense(i:i+1), 'LineWidth', 4, 'Color', cmap(i,:));
end


%%
plot(kind(:,2),kind(:,3),'LineWidth',4)

