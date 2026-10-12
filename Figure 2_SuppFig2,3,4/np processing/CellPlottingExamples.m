% Example PCs, MFs, and Nuclear cells 

% First plot to find example cells

%% 
%plotting

epo = 500;

evnt = "end_point";


%folder_path="C:\Users\Courtney\Documents\NeuropixelsResults\Cell Catalog\050622_SAM012\endpoint_raster";
%folder_path= "C:\Users\Courtney\Documents\NeuropixelsResults\Cell Catalog\051022_SAM012\Vmax_raster";
%for c = 1:length(data.units)
%for c = filtered_modulated_indices
for c=15

    figure; hold on
    
    subplot(2,1,1); hold on
    %h = title(['cluster' ' ' num2str(data.units(c).cluster_ID) ' ' 'depth' ' ' num2str(data.units(c).depth)]) ;
    h = title(['cluster' ' ' num2str(data.units(c).cluster_ID) ' ' 'depth' ' ' num2str(data.units(c).neurontype)]) ;
    set(h,'interpreter','none')
    
    % Plots Raster from binary matrix
    tick_length = 1; % length of tick marks for spikes
    spike_tick_width = 0.5; % width of lines marking spikes
    event_tick_width = 2; % width of lines marking secondary events
    
    
    rastr = data.units(c).(evnt).Raster(:,(1000-epo):(2000-(1000-epo)));
    [y_co,x_co] = find(rastr);
    x_co = x_co - epo;
    line([x_co';x_co'],[y_co'-(tick_length/2);y_co'+(tick_length/2)],'Color', [0 0 0], 'LineWidth', spike_tick_width);
    line([zeros(1,size(rastr,1));zeros(1,size(rastr,1))],[(1:1:size(rastr,1))-.5;(1:1:size(rastr,1))+.5],'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);
    axis ij;
    ylim([0,size(rastr,1)])
    xlim([-epo,epo])
    ylabel('Trial')
    
    
    subplot(2,1,2); hold on
    title('peth')
    
    psth_smooth = data.units(c).(evnt).PSTH_smooth;
    psth_raw = data.units(c).(evnt).PSTH_raw;
    td1 = length(psth_raw)-length(psth_smooth)+1;
    td2 = td1/2;
    plot((-epo:epo),psth_smooth(((1000-epo)-td2):((2000-(1000-epo))-td2)),'k','LineWidth',3)
    xlim([-epo,epo])
%     ylim([0,100])
  
    xlabel('Time from endpoint')
    ylabel('FR (Hz)')
%     filename=sprintf('endpt_raster%d.svg',c);
%     %filename=sprintf('vmax_raster%d.svg',c);
%     full_path=fullfile(folder_path,filename);
%     saveas(gcf,full_path);
    
    %close;
    
if mod(c,10) == 0
    drawnow
    input("press enter when done looking at figures")
    close all
end

end

%% editable figure producer
epo = 500;
evnt = "end_point";

 folder_path = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\Figures\ForKat"

for c = 1

    % Use painters renderer for editable vector output
    figure('Renderer','painters'); 

    % ---------------- Raster ----------------
    subplot(2,1,1); hold on
    title(['cluster ' num2str(data.units(c).cluster_ID) ...
           '  type ' num2str(data.units(c).neurontype)], ...
           'Interpreter','none')

    tick_length = 1;      
    spike_tick_width = 0.5;
    event_tick_width = 2;

    rastr = data.units(c).(evnt).Raster(:,(1000-epo):(2000-(1000-epo)));
    [y_co,x_co] = find(rastr);
    x_co = x_co - epo;

    % Plot spikes as individual short lines (separate objects in Illustrator)
    for i = 1:numel(x_co)
        plot([x_co(i) x_co(i)], ...
             [y_co(i)-tick_length/2 y_co(i)+tick_length/2], ...
             'k-', 'LineWidth', spike_tick_width);
    end

    % Plot event alignment ticks
    for trial = 1:size(rastr,1)
        plot([0 0], [trial-0.5 trial+0.5], ...
             'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);
    end

    axis ij;
    ylim([0, size(rastr,1)])
    xlim([-epo, epo])
    ylabel('Trial')

    % ---------------- PSTH ----------------
    subplot(2,1,2); hold on
    title('peth')

    psth_smooth = data.units(c).(evnt).PSTH_smooth;
    psth_raw = data.units(c).(evnt).PSTH_raw;
    td1 = length(psth_raw) - length(psth_smooth) + 1;
    td2 = td1/2;

    plot((-epo:epo), ...
        psth_smooth(((1000-epo)-td2):((2000-(1000-epo))-td2)), ...
        'k','LineWidth',3)

    xlim([-epo, epo])
    xlabel('Time from endpoint')
    ylabel('FR (Hz)')

    % ---------------- Save as SVG ----------------
    filename = sprintf('endpt_rasterexample%d.svg', c);
    full_path = fullfile(folder_path, filename);

    % Best for Illustrator: print with -painters
    print(gcf, full_path, '-dsvg', '-painters');

    % If on MATLAB R2020a+, this is also good:
    % exportgraphics(gcf, full_path, 'ContentType','vector');

    if mod(c,10) == 0
        drawnow
        input("press enter when done looking at figures")
        close all
    end
end

%% Smoothed version
epo = 500;
evnt = "end_point";

 folder_path = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\Figures\ForKat\050522_SAM009";

% Gaussian filter parameters
smoother = 5; % standard deviation
gauss_filter = fspecial('gaussian', [1, 6 * smoother], smoother);

for c = 8

    figure('Renderer','painters'); 

    % ---------------- Raster ----------------
    subplot(2,1,1); hold on
    title(['cluster ' num2str(data.units(c).cluster_ID) ...
           '  type ' num2str(data.units(c).neurontype)], ...
           'Interpreter','none')

    tick_length = 1;      
    spike_tick_width = 0.5;
    event_tick_width = 2;

    rastr = data.units(c).(evnt).Raster(:,(1000-epo):(2000-(1000-epo)));
    [y_co,x_co] = find(rastr);
    x_co = x_co - epo;

    for i = 1:numel(x_co)
        plot([x_co(i) x_co(i)], ...
             [y_co(i)-tick_length/2 y_co(i)+tick_length/2], ...
             'k-', 'LineWidth', spike_tick_width);
    end

    % for trial = 1:size(rastr,1)
    %     plot([0 0], [trial-0.5 trial+0.5], ...
    %          'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);
    % end

    axis ij;
    ylim([0, size(rastr,1)])
    xlim([-epo, epo])
    ylabel('Trial')

    % ---------------- PSTH ----------------
    subplot(2,1,2); hold on
    title('peth')

    psth_raw = data.units(c).(evnt).PSTH_raw;
    td1 = length(psth_raw) - length(data.units(c).(evnt).PSTH_smooth) + 1;
    td2 = td1/2;

    % Apply Gaussian smoothing to the raw PSTH
    psth_gaussian = conv(psth_raw, gauss_filter, 'same');

    % Crop the smoothed PSTH for plotting
    time_window = ((1000-epo)-td2):((2000-(1000-epo))-td2);
    plot((-epo:epo), psth_gaussian(time_window), ...
         'k', 'LineWidth', 3)
    
    ylim([0 70])
    xlim([-epo, epo])
    xlabel('Time from endpoint')
    ylabel('FR (Hz)')

    % ---------------- Save as SVG ----------------
    filename = sprintf('endpt_rasterNC%d.svg', c);
    full_path = fullfile(folder_path, filename);
    print(gcf, full_path, '-dsvg', '-painters');

    if mod(c,10) == 0
        drawnow
        input("press enter when done looking at figures")
        close all
    end
end

%% plot example kinematics 
%Reaches_align=AlignEndpoint(Reaches,500);

C= AlignEndpoint_seconds(Reaches,1000);
%%
dim=5;

% Make a wide, short figure
figW = 600;   % width (px)
figH = 300;   % height (px) -> very short
%figure('Renderer','painters','Position',[100 100 figW figH]);

StimvelPlot(C,dim,'k',1)
ylim([-0.5 2.5])
xlim([-0.5 0.5])
yline(0,'--k')
ylabel('Outward Velocity(cm/s)')
xlabel('Time from endpoint (s)')

%% Behavior plotting 

CSMcolor=rgb('Orchid');
%zscore plot 
% Adjust the figure to provide more vertical space
set(gcf, 'Position', [100, 100, 600, 300]); % [left, bottom, width, height]
hold on
dim = 5;
StimvelPlot(C,dim,'k',0);
ylim([-0.5 2.5])
xlim([-0.5 0.5])
yline(0,'--k')
ylabel('Outward Velocity(cm/s)')
xlabel('Time from endpoint (s)')


% Adjust the main plot's position to create space below
ax1 = gca;
ax1.Position(2) = ax1.Position(2) + 0.1;
ax1.Position(4) = ax1.Position(4) - 0.1;

