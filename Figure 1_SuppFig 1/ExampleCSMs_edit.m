%% Get CSM reaches 
names=fieldnames(ExcitRR_Sham);
T0Reaches=ExcitRR_Sham.(names{1}).RRT0;
T1Reaches=ExcitRR_Sham.(names{1}).RRT1;
T2Reaches=ExcitRR_Sham.(names{1}).RRT2;
T3Reaches=ExcitRR_Sham.(names{4}).RRT3;

NoCSMReaches=ExcitRR_Sham.(names{2}).NoReReach;

% Now do stimvel
thresh = .16; dirdir = 1;
T0=AlignThresh(T0Reaches,thresh,dirdir,500);
%T0=AlignThresh_speed(T0Reaches,thresh,dirdir,500);
T0=stimvel_accel(T0);
T0=stimvel_accel_jerk(T0);

T1=AlignThresh(T1Reaches,thresh,dirdir,500);
%T1=AlignThresh_speed(T1Reaches,thresh,dirdir,500);
T1=stimvel_accel(T1);
T1=stimvel_accel_jerk(T1);

T2=AlignThresh(T2Reaches,thresh,dirdir,500);
%T2=AlignThresh_speed(T2Reaches,thresh,dirdir,500);
T2=stimvel_accel(T2);
T2=stimvel_accel_jerk(T2);

T3=AlignThresh(T3Reaches,thresh,dirdir,500);
%T3=AlignThresh_speed(T3Reaches,thresh,dirdir,500);
T3=stimvel_accel(T3);
T3=stimvel_accel_jerk(T3);

NRR=AlignThresh(NoCSMReaches,thresh,dirdir,500);
%NRR=AlignThresh_speed(NoCSMReaches,thresh,dirdir,500);
NRR=stimvel_accel(NRR);
NRR=stimvel_accel_jerk(NRR);

T0fields=fieldnames(T0Reaches);
T1fields=fieldnames(T1Reaches);
T2fields=fieldnames(T2Reaches);
T3fields=fieldnames(T3Reaches);
NoCSMfields=fieldnames(NoCSMReaches);

%% Now want to go thrugh and plot position traces?? Need to find good examples of ReReaches - preferably in the outward position 



folder_path="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\Example CSMs\T0 CSM Examples";
%figure; hold on
for k=1:length(T0fields)
%for k=1:10
    RR=T0.stimvel(:,8,k);
    RR_smooth = smoothdata(RR, 'gaussian', 6);  % 7-point Gaussian window

    times=T0.stimvel(:,1,2);
    
    sz=20;
    figure;
    plot(times,RR_smooth,'Color',[0.5977    0.1953    0.7969],'LineWidth',2)
    % hold on
    % plot(times,RR,'Color',[0 0 0],'LineWidth',3)
    yline(0,'--')
    xlabel('Time from threshold crossing (s)')
    ylabel('Outward Acceleration (cm)')
    %title(Sucfields(k)) ;
    title(sprintf('Reach #%d', k));
    %ylabel('Upward Position (cm)')
    %ylabel('Lateral Position (cm)')
    %ylim([0 0.1])
    %ylim([0 0.25])
    ylim('auto')
    xlim([-0.1 0.2])
    
  
    %filename=sprintf('SucStimReach%d.tif',k);
    % filename=sprintf('FailStimReach%d.tif',k);
    % full_path=fullfile(folder_path,filename);
    % saveas(gcf,full_path);
    % 
    % close;
    if mod(k,10) == 0
    drawnow
    input("press enter when done looking at figures")
    close all
    end


end

%%
folder_path = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Neuron Reviews\Example CSMs\T0 CSM Examples";

%for k = 1:length(T0fields)
for k=17
    % Extract each signal
    pos = T0.stimvel(:,2,k);       % position
    vel = T0.stimvel(:,5,k);       % velocity
    acc = T0.stimvel(:,8,k);       % acceleration
    jerk = T0.stimvel(:,11,k);     % jerk
    
    % Smooth each signal
    pos_smooth  = smoothdata(pos, 'gaussian', 6);
    vel_smooth  = smoothdata(vel, 'gaussian', 6);
    acc_smooth  = smoothdata(acc, 'gaussian', 6);
    jerk_smooth = smoothdata(jerk, 'gaussian', 6);

    % Time vector
    times = T0.stimvel(:,1,2);

    % Create figure
    figure('Position',[100 100 400 800]); % taller figure for column layout

    % Plot position
    subplot(4,1,1)
    plot(times, pos_smooth, 'b', 'LineWidth', 2)
    ylabel('Position')
    yline(0,'--')
    xlim([-0.2 0.2])
    title(sprintf('Reach #%d', k))

    % Plot velocity
    subplot(4,1,2)
    plot(times, vel_smooth, 'r', 'LineWidth', 2)
    ylabel('Velocity')
    yline(0,'--')
    xlim([-0.2 0.2])

    % Plot acceleration
    subplot(4,1,3)
    plot(times, acc_smooth, 'm', 'LineWidth', 2)
    ylabel('Acceleration')
    yline(0,'--')
    xlim([-0.2 0.2])

    % Plot jerk
    subplot(4,1,4)
    plot(times, jerk_smooth, 'k', 'LineWidth', 2)
    ylabel('Jerk')
    yline(0,'--')
    xlabel('Time from threshold crossing (s)')
    xlim([-0.2 0.2])

    % Optional: save figure
    % filename = sprintf('Reach_%d.tif', k);
    % full_path = fullfile(folder_path, filename);
    % saveas(gcf, full_path);

    % Pause every 10 reaches for review
    if mod(k,10) == 0
        drawnow
        input("Press enter when done looking at figures")
        close all
    end
end

%%


%for k = 1:length(T1fields)
for k=18
    % Extract each signal
    pos = T1.stimvel(:,2,k);       % position
    vel = T1.stimvel(:,5,k);       % velocity
    acc = T1.stimvel(:,8,k);       % acceleration
    jerk = T1.stimvel(:,11,k);     % jerk
    
    % Smooth each signal
    pos_smooth  = smoothdata(pos, 'gaussian', 6);
    vel_smooth  = smoothdata(vel, 'gaussian', 6);
    acc_smooth  = smoothdata(acc, 'gaussian', 6);
    jerk_smooth = smoothdata(jerk, 'gaussian', 6);

    % Time vector
    times = T1.stimvel(:,1,2);

    % Create figure
    figure('Position',[100 100 400 800]); % taller figure for column layout

    % Plot position
    subplot(4,1,1)
    plot(times, pos_smooth, 'b', 'LineWidth', 2)
    ylabel('Position')
    yline(0,'--')
    xlim([-0.2 0.2])
    title(sprintf('Reach #%d', k))

    % Plot velocity
    subplot(4,1,2)
    plot(times, vel_smooth, 'r', 'LineWidth', 2)
    ylabel('Velocity')
    yline(0,'--')
    xlim([-0.2 0.2])

    % Plot acceleration
    subplot(4,1,3)
    plot(times, acc_smooth, 'm', 'LineWidth', 2)
    ylabel('Acceleration')
    yline(0,'--')
    xlim([-0.2 0.2])

    % Plot jerk
    subplot(4,1,4)
    plot(times, jerk_smooth, 'k', 'LineWidth', 2)
    ylabel('Jerk')
    yline(0,'--')
    xlabel('Time from threshold crossing (s)')
    xlim([-0.2 0.2])

    % Optional: save figure
    % filename = sprintf('Reach_%d.tif', k);
    % full_path = fullfile(folder_path, filename);
    % saveas(gcf, full_path);

    % Pause every 10 reaches for review
    if mod(k,10) == 0
        drawnow
        input("Press enter when done looking at figures")
        close all
    end
end
%% T2 reaches 

%for k = 1:length(T2fields)
for k=20
    % Extract each signal
    pos = T2.stimvel(:,2,k);       % position
    vel = T2.stimvel(:,5,k);       % velocity
    acc = T2.stimvel(:,8,k);       % acceleration
    jerk = T2.stimvel(:,11,k);     % jerk
    
    % Smooth each signal
    pos_smooth  = smoothdata(pos, 'gaussian', 6);
    vel_smooth  = smoothdata(vel, 'gaussian', 6);
    acc_smooth  = smoothdata(acc, 'gaussian', 6);
    jerk_smooth = smoothdata(jerk, 'gaussian', 6);

    % Time vector
    times = T2.stimvel(:,1,2);

    % Create figure
    figure('Position',[100 100 400 800]); % taller figure for column layout

    % Plot position
    subplot(4,1,1)
    plot(times, pos_smooth, 'b', 'LineWidth', 2)
    ylabel('Position')
    yline(0,'--')
    xlim([-0.2 0.2])
    title(sprintf('Reach #%d', k))

    % Plot velocity
    subplot(4,1,2)
    plot(times, vel_smooth, 'r', 'LineWidth', 2)
    ylabel('Velocity')
    yline(0,'--')
    xlim([-0.2 0.2])

    % Plot acceleration
    subplot(4,1,3)
    plot(times, acc_smooth, 'm', 'LineWidth', 2)
    ylabel('Acceleration')
    yline(0,'--')
    xlim([-0.2 0.2])

    % Plot jerk
    subplot(4,1,4)
    plot(times, jerk_smooth, 'k', 'LineWidth', 2)
    ylabel('Jerk')
    yline(0,'--')
    xlabel('Time from threshold crossing (s)')
    xlim([-0.2 0.2])

    % Optional: save figure
    % filename = sprintf('Reach_%d.tif', k);
    % full_path = fullfile(folder_path, filename);
    % saveas(gcf, full_path);

    % Pause every 10 reaches for review
    if mod(k,10) == 0
        drawnow
        input("Press enter when done looking at figures")
        close all
    end
end


%% T3 reaches 

%for k = 1:length(T3fields)
for k=62

    % Extract each signal
    pos = T3.stimvel(:,2,k);       % position
    vel = T3.stimvel(:,5,k);       % velocity
    acc = T3.stimvel(:,8,k);       % acceleration
    jerk = T3.stimvel(:,11,k);     % jerk
    
    % Smooth each signal
    pos_smooth  = smoothdata(pos, 'gaussian', 6);
    vel_smooth  = smoothdata(vel, 'gaussian', 6);
    acc_smooth  = smoothdata(acc, 'gaussian', 6);
    jerk_smooth = smoothdata(jerk, 'gaussian', 6);

    % Time vector
    times = T3.stimvel(:,1,2);

    % Create figure
    figure('Position',[100 100 400 800]); % taller figure for column layout

    % Plot position
    subplot(4,1,1)
    plot(times, pos_smooth, 'b', 'LineWidth', 2)
    ylabel('Position')
    yline(0,'--')
    xlim([-0.2 0.3])
    title(sprintf('Reach #%d', k))

    % Plot velocity
    subplot(4,1,2)
    plot(times, vel_smooth, 'r', 'LineWidth', 2)
    ylabel('Velocity')
    yline(0,'--')
    xlim([-0.2 0.3])

    % Plot acceleration
    subplot(4,1,3)
    plot(times, acc_smooth, 'm', 'LineWidth', 2)
    ylabel('Acceleration')
    yline(0,'--')
    xlim([-0.2 0.3])

    % Plot jerk
    subplot(4,1,4)
    plot(times, jerk_smooth, 'k', 'LineWidth', 2)
    ylabel('Jerk')
    yline(0,'--')
    xlabel('Time from threshold crossing (s)')
    xlim([-0.2 0.3])

    % Optional: save figure
    % filename = sprintf('Reach_%d.tif', k);
    % full_path = fullfile(folder_path, filename);
    % saveas(gcf, full_path);

    % Pause every 10 reaches for review
    if mod(k,10) == 0
        drawnow
        input("Press enter when done looking at figures")
        close all
    end
end


%% NoCSM Reaches 

%for k = 1:length(NoCSMfields)
for k=33
    % Extract each signal
    pos = NRR.stimvel(:,2,k);       % position
    vel = NRR.stimvel(:,5,k);       % velocity
    acc = NRR.stimvel(:,8,k);       % acceleration
    jerk = NRR.stimvel(:,11,k);     % jerk
    
    % Smooth each signal
    pos_smooth  = smoothdata(pos, 'gaussian', 6);
    vel_smooth  = smoothdata(vel, 'gaussian', 6);
    acc_smooth  = smoothdata(acc, 'gaussian', 6);
    jerk_smooth = smoothdata(jerk, 'gaussian', 6);

    % Time vector
    times = NRR.stimvel(:,1,2);

    % Create figure
    figure('Position',[100 100 400 800]); % taller figure for column layout

    % Plot position
    subplot(4,1,1)
    plot(times, pos_smooth, 'b', 'LineWidth', 2)
    ylabel('Position')
    yline(0,'--')
    xlim([-0.2 0.2])
    title(sprintf('Reach #%d', k))

    % Plot velocity
    subplot(4,1,2)
    plot(times, vel_smooth, 'r', 'LineWidth', 2)
    ylabel('Velocity')
    yline(0,'--')
    xlim([-0.2 0.2])

    % Plot acceleration
    subplot(4,1,3)
    plot(times, acc_smooth, 'm', 'LineWidth', 2)
    ylabel('Acceleration')
    yline(0,'--')
    xlim([-0.2 0.2])

    % Plot jerk
    subplot(4,1,4)
    plot(times, jerk_smooth, 'k', 'LineWidth', 2)
    ylabel('Jerk')
    yline(0,'--')
    xlabel('Time from threshold crossing (s)')
    xlim([-0.2 0.2])

    % Optional: save figure
    % filename = sprintf('Reach_%d.tif', k);
    % full_path = fullfile(folder_path, filename);
    % saveas(gcf, full_path);

    % Pause every 10 reaches for review
    if mod(k,10) == 0
        drawnow
        input("Press enter when done looking at figures")
        close all
    end
end
%% Now want to go thrugh and plot position traces?? Need to find good examples of ReReaches - preferably in the outward position 



folder_path="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\CSM for ControlAnimals\CSMIndividualReaches_animals\EX002\T2CSMs" ;
%figure; hold on
%for k=1:length(T2fields)
for k=1:10
    RR=T2.stimvel(:,8,k);
    RR_smooth = smoothdata(RR, 'gaussian', 6);  % 7-point Gaussian window

    times=T2.stimvel(:,1,2);
    
    sz=20;
    figure;
    plot(times,RR_smooth,'Color',[0.5977    0.1953    0.7969],'LineWidth',2)
    % hold on
    % plot(times,RR,'Color',[0 0 0],'LineWidth',3)
    yline(0,'--')
    xlabel('Time from threshold crossing (s)')
    ylabel('Outward Acceleration (cm)')
    %title(Sucfields(k)) ;
    title(sprintf('Reach #%d', k));
    %ylabel('Upward Position (cm)')
    %ylabel('Lateral Position (cm)')
    %ylim([0 0.1])
    %ylim([0 0.25])
    ylim('auto')
    xlim([-0.1 0.2])
    
  
    %filename=sprintf('SucStimReach%d.tif',k);
    % filename=sprintf('FailStimReach%d.tif',k);
    % full_path=fullfile(folder_path,filename);
    % saveas(gcf,full_path);
    % 
    % close;
    if mod(k,10) == 0
    drawnow
    input("press enter when done looking at figures")
    close all
    end


end
 %%


folder_path="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\CSM for ControlAnimals\CSMIndividualReaches_animals\EX002\T2CSMs" ;
%figure; hold on
for k=1:length(T1fields)
%for k=1:10
    RR=T1.stimvel(:,8,k);
    RR_smooth = smoothdata(RR, 'gaussian', 6);  % 7-point Gaussian window

    times=T2.stimvel(:,1,2);
    
    sz=20;
    figure;
    plot(times,RR_smooth,'Color',[0.5977    0.1953    0.7969],'LineWidth',2)
    % hold on
    % plot(times,RR,'Color',[0 0 0],'LineWidth',3)
    yline(0,'--')
    xlabel('Time from threshold crossing (s)')
    ylabel('Outward Acceleration (cm)')
    %title(Sucfields(k)) ;
    title(sprintf('Reach #%d', k));
    %ylabel('Upward Position (cm)')
    %ylabel('Lateral Position (cm)')
    %ylim([0 0.1])
    %ylim([0 0.25])
    ylim('auto')
    xlim([-0.1 0.2])
    
  
    %filename=sprintf('SucStimReach%d.tif',k);
    % filename=sprintf('FailStimReach%d.tif',k);
    % full_path=fullfile(folder_path,filename);
    % saveas(gcf,full_path);
    % 
    % close;
    if mod(k,10) == 0
    drawnow
    input("press enter when done looking at figures")
    close all
    end


end
%% Now create an example plot
% T1idx=[5,6,11,18]';
% T2idx=[7,8,10,16,19]';

T0idx=[3,17,24,36]';
T1idx=[5,6,11,18]';
T2idx=[10,11,13,20]';
T3idx=[10,18,19,22]';

%T0 reaches 
figure; hold on
% Generate blue shades for T0
nT0 = length(T0idx);
greens = [linspace(0.2, 0.4, nT1)', linspace(0.4, 0.7, nT1)', linspace(0.8, 1.0, nT1)'];


for h = 1:length(T0idx)
    j = T0idx(h);
    RR=T0.stimvel(:,8,j);
    RR_smooth = smoothdata(RR, 'gaussian', 8);  % 7-point Gaussian window

    times=T0.stimvel(:,1,2);

    %plot(times,RR_smooth,'Color',[0 0 0],'LineWidth',2)
    plot(times, RR_smooth, 'Color', greens(h,:), 'LineWidth', 2)
    yline(0,'--')
    ylim('auto')
    xlim([-0.15 0.2])
    title('T0 reaches')

end


figure; hold on
% Generate blue shades for T1
nT1 = length(T1idx);
blues = [linspace(0.2, 0.4, nT1)', linspace(0.4, 0.7, nT1)', linspace(0.8, 1.0, nT1)'];


for h = 1:length(T1idx)
    j = T1idx(h);
    RR=T1.stimvel(:,8,j);
    RR_smooth = smoothdata(RR, 'gaussian', 8);  % 7-point Gaussian window

    times=T1.stimvel(:,1,2);

    %plot(times,RR_smooth,'Color',[0 0 0],'LineWidth',2)
    plot(times, RR_smooth, 'Color', blues(h,:), 'LineWidth', 2)
    yline(0,'--')
    ylim('auto')
    xlim([-0.15 0.2])
    title('T1 reaches')

end


% Generate distinct purple shades for T2idx
nT2 = length(T2idx);
purples = [linspace(0.4, 0.8, nT2)', linspace(0, 0.3, nT2)', linspace(0.6, 0.9, nT2)'];

figure; hold on
for i = 1:length(T2idx)
    k = T2idx(i);
    RR=T2.stimvel(:,8,k);
    RR_smooth = smoothdata(RR, 'gaussian', 8);  % 7-point Gaussian window
    times=T2.stimvel(:,1,2);

 
    %plot(times,RR_smooth,'Color',[0 0 0],'LineWidth',2)
    plot(times, RR_smooth, 'Color', purples(i,:), 'LineWidth', 2)

    yline(0,'--')
    xlabel('Time from threshold crossing (s)')
    ylabel('Outward Acceleration (cm)')
    ylim('auto')
    xlim([-0.15 0.2])
    title('T2 reaches ')


end



% Generate distinct purple shades for T3idx
nT3 = length(T3idx);
purples = [linspace(0.4, 0.8, nT2)', linspace(0, 0.3, nT2)', linspace(0.6, 0.9, nT2)'];

figure; hold on
for i = 1:length(T3idx)
    k = T3idx(i);
    RR=T3.stimvel(:,8,k);
    RR_smooth = smoothdata(RR, 'gaussian', 8);  % 7-point Gaussian window
    times=T3.stimvel(:,1,2);

 
    %plot(times,RR_smooth,'Color',[0 0 0],'LineWidth',2)
    plot(times, RR_smooth, 'Color', purples(i,:), 'LineWidth', 2)

    yline(0,'--')
    xlabel('Time from threshold crossing (s)')
    ylabel('Outward Acceleration (cm)')
    ylim('auto')
    xlim([-0.15 0.2])
    title('T3 reaches ')


end

%%
T1idx = [5,6,11,18]';
T2idx = [7,8,10,16,19]';

allidx = [T2idx; T1idx];
nTraces = length(allidx);

% Create purple color gradient: from light lavender to deep purple
purples = [linspace(0.85, 0.4, nTraces)', linspace(0.7, 0.1, nTraces)', linspace(1.0, 0.6, nTraces)'];

figure; hold on

% Get shared time vector
times = T2.stimvel(:,1,2);
t0_idx = find(times >= 0, 1, 'first');

for i = 1:nTraces
    idx = allidx(i);

    % Determine which struct to use (T1 or T2)
    if ismember(idx, T2idx)
        RR = T2.stimvel(:,8,idx);
    else
        RR = T1.stimvel(:,8,idx);
    end

    RR_smooth = smoothdata(RR, 'gaussian', 8);
    plot(times, RR_smooth, 'Color', purples(i,:), 'LineWidth', 2)

    % Zero crossing detection
    s_post = RR_smooth(t0_idx:end);
    t_post = times(t0_idx:end);
    sign_change_idx = find(diff(sign(s_post)) ~= 0, 1, 'first');

    if ~isempty(sign_change_idx) && sign_change_idx < length(t_post)
        x1 = t_post(sign_change_idx);
        x2 = t_post(sign_change_idx + 1);
        y1 = s_post(sign_change_idx);
        y2 = s_post(sign_change_idx + 1);
        x_cross = x1 - y1 * (x2 - x1) / (y2 - y1);

        plot(x_cross, 0, 'mo', 'MarkerSize', 6, 'LineWidth', 1.5)
    else
        fprintf('No zero crossing found for trace %d\n', idx);
    end
end

% Plot formatting
yline(0, '--')
xline(0,'--')
xlabel('Time from threshold crossing (s)')
ylabel('Outward Acceleration (cm)')
xlim([-0.15 0.2])
ylim('auto')

%%
% T1idx = [5,6,11,18]';
% T2idx = [7,8,10,16,19]';


T0idx=[3,17,24,36]';
T1idx=[5,6,11,18]';
T2idx=[10,11,13,20]';
T3idx=[10,18,19,22]';

allidx = [T2idx; T1idx];
nTraces = length(allidx);

% Create purple color gradient: from light lavender to deep purple
purples = [linspace(0.85, 0.4, nTraces)', linspace(0.7, 0.1, nTraces)', linspace(1.0, 0.6, nTraces)'];

% Create figure with subplots
%figure;
figure('Position', [100, 100, 600, 800]);
% Get shared time vector
times = T2.stimvel(:,1,2);
t0_idx = find(times >= 0, 1, 'first');

for i = 1:nTraces
    idx = allidx(i);

    % Determine which struct to use (T1 or T2)
    if ismember(idx, T2idx)
        RR_pos = T2.stimvel(:,2,idx);  % Outward position (2nd column)
        RR_accel = T2.stimvel(:,8,idx); % Outward acceleration (8th column)
    else
        RR_pos = T1.stimvel(:,2,idx);  % Outward position (2nd column)
        RR_accel = T1.stimvel(:,8,idx); % Outward acceleration (8th column)
    end

    % Smooth data
    RR_pos_smooth = smoothdata(RR_pos, 'gaussian', 8);
    RR_accel_smooth = smoothdata(RR_accel, 'gaussian', 8);

    % Plot position (subplot 1)
    subplot(2, 1, 1); % First subplot for position
    plot(times, RR_pos_smooth, 'Color', purples(i,:), 'LineWidth', 2);
    hold on;

    % Plot acceleration (subplot 2)
    subplot(2, 1, 2); % Second subplot for acceleration
    plot(times, RR_accel_smooth, 'Color', purples(i,:), 'LineWidth', 2);
    hold on;

    % Zero crossing detection for acceleration (after threshold)
    s_post = RR_accel_smooth(t0_idx:end);  % Acceleration after threshold
    t_post = times(t0_idx:end);            % Time after threshold
    sign_change_idx = find(diff(sign(s_post)) ~= 0, 1, 'first');  % First zero crossing

    if ~isempty(sign_change_idx) && sign_change_idx < length(t_post)
        x1 = t_post(sign_change_idx);
        x2 = t_post(sign_change_idx + 1);
        y1 = s_post(sign_change_idx);
        y2 = s_post(sign_change_idx + 1);
        x_cross = x1 - y1 * (x2 - x1) / (y2 - y1);  % Interpolation for zero crossing

        % Find the index closest to the zero-crossing time in the position data
        [~, idx_cross] = min(abs(times - x_cross));

        % Plot the zero crossing on both the position plot and the acceleration plot
        subplot(2, 1, 1); % Position subplot
        plot(x_cross, RR_pos_smooth(idx_cross), 'mo', 'MarkerSize', 6, 'LineWidth', 1.5);

        subplot(2, 1, 2); % Acceleration subplot
        plot(x_cross, 0, 'mo', 'MarkerSize', 6, 'LineWidth', 1.5); % Zero crossing on accel plot
    else
        fprintf('No zero crossing found for trace %d\n', idx);
    end

    % Mark the threshold crossing on the acceleration plot
    threshold_time = 0;
    [~, idx_threshold] = min(abs(times - threshold_time));
    subplot(2, 1, 2); % Acceleration subplot
    %plot(times(idx_threshold), RR_accel_smooth(idx_threshold), 'go', 'MarkerSize', 6, 'LineWidth', 1.5); % Threshold crossing (green)
end

% Format position plot
subplot(2, 1, 1);
yline(0, '--');
xlabel('Time from threshold crossing (s)');
ylabel('Outward Position (cm)');
xlim([-0.15 0.2]);
ylim('auto');
title('Outward Position');

% Format acceleration plot
subplot(2, 1, 2);
yline(0, '--');
xlabel('Time from threshold crossing (s)');
ylabel('Outward Acceleration (cm/s^2)');
xlim([-0.15 0.2]);
ylim('auto');
title('Outward Acceleration');

%% Now plotting all subtypes: 
% Define indices for each time point
T0idx = [3,17,24,36]';
T1idx = [5,6,11,18]';
T2idx = [10,11,13,20]';
T3idx = [10,18,19,22]';

all_structs = {'T0','T1','T2','T3'};
all_idx = {T0idx, T1idx, T2idx, T3idx};

% Colors
nTracesMax = max(cellfun(@length, all_idx));
purples = [linspace(0.85, 0.4, nTracesMax)', linspace(0.7, 0.1, nTracesMax)', linspace(1.0, 0.6, nTracesMax)'];

% Loop over each timepoint figure
for t = 1:length(all_structs)
    figName = all_structs{t};
    idxList = all_idx{t};
    nTraces = length(idxList);

    % Create figure
    figure('Name', figName, 'NumberTitle', 'off', 'Position', [100, 100, 600, 800]);
    
    % Time vector (assuming all T structs have the same time vector)
    times = eval([figName '.stimvel(:,1,2)']); 
    t0_idx = find(times >= 0, 1, 'first');

    % Loop over traces
    for i = 1:nTraces
        idx = idxList(i);

        % Extract signals
        RR_pos  = eval([figName '.stimvel(:,2,idx)']);
        RR_vel  = eval([figName '.stimvel(:,5,idx)']);
        RR_acc  = eval([figName '.stimvel(:,8,idx)']);
        RR_jerk = eval([figName '.stimvel(:,11,idx)']);

        % Smooth signals
        pos_smooth  = smoothdata(RR_pos,'gaussian',6);
        vel_smooth  = smoothdata(RR_vel,'gaussian',6);
        acc_smooth  = smoothdata(RR_acc,'gaussian',6);
        jerk_smooth = smoothdata(RR_jerk,'gaussian',6);

        % Zero-crossing detection in acceleration
        s_post = acc_smooth(t0_idx:end);
        t_post = times(t0_idx:end);
        sign_change_idx = find(diff(sign(s_post)) ~= 0, 1, 'first');

        if ~isempty(sign_change_idx) && sign_change_idx < length(t_post)
            % Interpolated zero crossing
            x1 = t_post(sign_change_idx);
            x2 = t_post(sign_change_idx+1);
            y1 = s_post(sign_change_idx);
            y2 = s_post(sign_change_idx+1);
            x_cross = x1 - y1*(x2-x1)/(y2-y1);

            % Find index closest to zero-crossing
            [~, idx_cross] = min(abs(times - x_cross));
        else
            x_cross = NaN;
            idx_cross = NaN;
        end

        % Plot subplots
        % Position
        subplot(4,1,1); hold on;
        plot(times, pos_smooth, 'Color', purples(i,:), 'LineWidth', 2);
        if ~isnan(idx_cross)
            plot(x_cross, pos_smooth(idx_cross), 'mo', 'MarkerSize',6,'LineWidth',1.5);
        end

        % Velocity
        subplot(4,1,2); hold on;
        plot(times, vel_smooth, 'Color', purples(i,:), 'LineWidth', 2);
        if ~isnan(idx_cross)
            plot(x_cross, vel_smooth(idx_cross), 'mo', 'MarkerSize',6,'LineWidth',1.5);
        end

        % Acceleration
        subplot(4,1,3); hold on;
        plot(times, acc_smooth, 'Color', purples(i,:), 'LineWidth', 2);
        if ~isnan(idx_cross)
            plot(x_cross, 0, 'mo', 'MarkerSize',6,'LineWidth',1.5);
        end
        yline(0,'--');

        % Jerk
        subplot(4,1,4); hold on;
        plot(times, jerk_smooth, 'Color', purples(i,:), 'LineWidth', 2);
        if ~isnan(idx_cross)
            plot(x_cross, jerk_smooth(idx_cross), 'mo', 'MarkerSize',6,'LineWidth',1.5);
        end
    end

    % Format subplots
    subplot(4,1,1);
    ylabel('Position (cm)'); xlim([-0.15 0.2]); title([figName ' - Position']);

    subplot(4,1,2);
    ylabel('Velocity (cm/s)'); xlim([-0.15 0.2]); title([figName ' - Velocity']);

    subplot(4,1,3);
    ylabel('Acceleration (cm/s^2)'); xlim([-0.15 0.2]); title([figName ' - Acceleration']);

    subplot(4,1,4);
    ylabel('Jerk (cm/s^3)'); xlabel('Time from threshold crossing (s)'); xlim([-0.15 0.2]); title([figName ' - Jerk']);
end

%% Different version without points 

T0idx = [3,17,24,36]';
T1idx = [5,6,11,18]';
T2idx = [10,11,13,20]';
%T3idx = [10,18,19,22]'; % from animal 1 or 2 I think ?? 
T3idx=[41,48,59,62]'; % from animal 4
NoCSMidx=[2,4,9,15,40];

all_structs = {'T0','T1','T2','T3','NRR'};
all_idx = {T0idx, T1idx, T2idx, T3idx, NoCSMidx};

% Colors
nTracesMax = max(cellfun(@length, all_idx));
purples = [linspace(0.85, 0.4, nTracesMax)', linspace(0.7, 0.1, nTracesMax)', linspace(1.0, 0.6, nTracesMax)'];

% Loop over each timepoint figure
for t = 1:length(all_structs)
    figName = all_structs{t};
    idxList = all_idx{t};
    nTraces = length(idxList);

    % Create figure (narrower width)
    figure('Name', figName, 'NumberTitle', 'off', 'Position', [100, 100, 400, 800]);
    
    % Time vector (assuming all T structs have the same time vector)
    times = eval([figName '.stimvel(:,1,2)']); 
    t0_idx = find(times >= 0, 1, 'first');

    % Loop over traces
    for i = 1:nTraces
        idx = idxList(i);

        % Extract signals
        RR_pos  = eval([figName '.stimvel(:,2,idx)']);
        RR_vel  = eval([figName '.stimvel(:,5,idx)']);
        RR_acc  = eval([figName '.stimvel(:,8,idx)']);
        RR_jerk = eval([figName '.stimvel(:,11,idx)']);

        % Smooth signals
        % pos_smooth  = smoothdata(RR_pos,'gaussian',6);
        % vel_smooth  = smoothdata(RR_vel,'gaussian',6);
        % acc_smooth  = smoothdata(RR_acc,'gaussian',6);
        % jerk_smooth = smoothdata(RR_jerk,'gaussian',6);
        
        % Look at plots without smoothing 
        pos_smooth  = RR_pos;
        vel_smooth  = RR_vel;
        acc_smooth  = RR_acc;
        jerk_smooth = RR_jerk;

        % Plot subplots
        % Position
        subplot(4,1,1); hold on;
        plot(times, pos_smooth, 'Color', purples(i,:), 'LineWidth', 2);
        yline(0,'--');
        ylim([0 0.3])

        % Velocity
        subplot(4,1,2); hold on;
        plot(times, vel_smooth, 'Color', purples(i,:), 'LineWidth', 2);
        yline(0,'--');
        %ylim([-1 2])
        %ylim([0 3])
            
        % Acceleration
        subplot(4,1,3); hold on;
        plot(times, acc_smooth, 'Color', purples(i,:), 'LineWidth', 2);
        yline(0,'--');
        %ylim([-40 40])

        % Jerk
        subplot(4,1,4); hold on;
        plot(times, jerk_smooth, 'Color', purples(i,:), 'LineWidth', 2);
        yline(0,'--');
        %ylim ([-2200 2000])

    end

    % % Format subplots
    subplot(4,1,1);
    ylabel('Position (cm)'); xlim([-0.15 0.3]); title([figName ' - Position']);

    subplot(4,1,2);
    ylabel('Velocity (cm/s)'); xlim([-0.15 0.3]); title([figName ' - Velocity']);

    subplot(4,1,3);
    ylabel('Acceleration (cm/s^2)'); xlim([-0.15 0.3]); title([figName ' - Acceleration']);

    subplot(4,1,4);
    ylabel('Jerk (cm/s^3)'); xlabel('Time from threshold crossing (s)'); xlim([-0.15 0.3]); title([figName ' - Jerk']);


    % Format subplots
    % subplot(4,1,1);
    % ylabel('Position (cm)'); xlim([-0.15 0.2]); title([figName ' - Position']);
    % 
    % subplot(4,1,2);
    % ylabel('Velocity (cm/s)'); xlim([-0.15 0.2]); title([figName ' - Velocity']);
    % % 
    % subplot(4,1,3);
    % ylabel('Acceleration (cm/s^2)'); xlim([-0.15 0.3]); title([figName ' - Acceleration']);
    % 
    % subplot(4,1,4);
    % ylabel('Jerk (cm/s^3)'); xlabel('Time from threshold crossing (s)'); xlim([-0.15 0.3]); title([figName ' - Jerk']);
end


%% Saving animation
% T1idx = [5,6,11,18]'; 
% T2idx = [7,8,10,16,19]';

T1idx = [5,18]'; 
T2idx = [7,8,19]';
% 
% T1idx = [6]';
% T2idx = [16]';

allidx = [T2idx; T1idx];
nTraces = length(allidx);

% Color before and after zero crossing
beforeColor = [0.6 0.6 1];  % Light purple
afterColor = [0.4 0 0.6];   % Deep purple

% Get shared time vector
times = T2.stimvel(:,1,2);
t0_idx = find(times >= 0, 1, 'first');

% Create figure
%figure('Position', [100, 100, 600, 800]);
fig = figure('Position', [100, 100, 600, 800], 'Color', 'w');  % White background everywhere


% Setup subplots
subplot(2, 1, 1);
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
hold on;
xlabel('Time from threshold crossing (s)');
ylabel('Outward Position (cm)');
title('Outward Position');
xlim([-0.15 0.2]);
ylim([0 3.0]);
yline(0, '--');

subplot(2, 1, 2);
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
hold on;
xlabel('Time from threshold crossing (s)');
ylabel('Outward Acceleration (cm/s^2)');
title('Outward Acceleration');
xlim([-0.15 0.2]);
ylim([-500 400]);
yline(0, '--');
%
% Create and open video file

saveFolder = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\CSM for ControlAnimals\CSMIndividualReaches_animals";  % <-- change to your path
filename = fullfile(saveFolder, 'ReachAnimation.mp4');
v = VideoWriter(filename, 'MPEG-4');

v.FrameRate = 120;  % Can adjust for speed
open(v);

for i = 1:nTraces
    idx = allidx(i);

    % Choose T1 or T2
    if ismember(idx, T2idx)
        RR_pos = T2.stimvel(:,2,idx);
        RR_accel = T2.stimvel(:,8,idx);
    else
        RR_pos = T1.stimvel(:,2,idx);
        RR_accel = T1.stimvel(:,8,idx);
    end

    % Smooth
    RR_pos_smooth = smoothdata(RR_pos, 'gaussian', 8);
    RR_accel_smooth = smoothdata(RR_accel, 'gaussian', 8);

    RR_pos_smooth = RR_pos_smooth * 10;  % Convert position to centimeters
    RR_accel_smooth= RR_accel_smooth * 10;  % Convert acceleration to centimeters (cm/s^2)


    % Find zero-crossing after t = 0
    s_post = RR_accel_smooth(t0_idx:end);
    t_post = times(t0_idx:end);
    sign_change_idx = find(diff(sign(s_post)) ~= 0, 1, 'first');

    if ~isempty(sign_change_idx) && sign_change_idx < length(t_post)
        x1 = t_post(sign_change_idx);
        x2 = t_post(sign_change_idx + 1);
        y1 = s_post(sign_change_idx);
        y2 = s_post(sign_change_idx + 1);
        x_cross = x1 - y1 * (x2 - x1) / (y2 - y1);
        [~, idx_cross] = min(abs(times - x_cross));
    else
        x_cross = inf;
        idx_cross = length(times);
        fprintf('No zero crossing found for trace %d\n', idx);
    end

    % Animate line segment by segment
    for t = 2:length(times)
        current_time = times(t);
        if current_time > 0 && current_time > x_cross
            c = afterColor;
        else
            c = beforeColor;
        end

        subplot(2,1,1);
        plot(times(t-1:t), RR_pos_smooth(t-1:t), '-', 'Color', c, 'LineWidth', 2);

        subplot(2,1,2);
        plot(times(t-1:t), RR_accel_smooth(t-1:t), '-', 'Color', c, 'LineWidth', 2);

        % Write this frame to video
        frame = getframe(gcf);
        writeVideo(v, frame);

        pause(0.005);  % Optional visual pause
    end

    % Mark zero crossing (optional)
    if isfinite(x_cross)
        subplot(2,1,1);
        plot(x_cross, RR_pos_smooth(idx_cross), 'mo', 'MarkerSize', 6, 'LineWidth', 1.5, 'MarkerFaceColor','m');

        subplot(2,1,2);
        plot(x_cross, 0, 'mo', 'MarkerSize', 6, 'LineWidth', 1.5,'MarkerFaceColor','m');

        % Save this marking as a frame
        frame = getframe(gcf);
        writeVideo(v, frame);
    end

    pause(0.005);  % Delay between traces (can be removed)
   % pause(0.1);  % Delay between traces (can be removed)
end

% Close video
close(v);
disp('Animation saved to ReachAnimation.mp4');
%% Faster version

%% Saving animation
T1idx = [5,6,11,18]'; 
T2idx = [7,8,10,16,19]';
% 
% T1idx = [6]';
% T2idx = [16]';

allidx = [T2idx; T1idx];
nTraces = length(allidx);

% Color before and after zero crossing
beforeColor = [0.6 0.6 1];  % Light purple
afterColor = [0.4 0 0.6];   % Deep purple

% Get shared time vector
times = T2.stimvel(:,1,2);
t0_idx = find(times >= 0, 1, 'first');

% Create figure
fig = figure('Position', [100, 100, 600, 800], 'Color', 'w');  % White background everywhere

% Setup subplots
subplot(2, 1, 1);
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
hold on;
xlabel('Time from threshold crossing (s)');
ylabel('Outward Position (cm)');
title('Outward Position');
xlim([-0.15 0.2]);
ylim([0 3.0]);
yline(0, '--');

subplot(2, 1, 2);
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
hold on;
xlabel('Time from threshold crossing (s)');
ylabel('Outward Acceleration (cm/s^2)');
title('Outward Acceleration');
xlim([-0.15 0.2]);
ylim([-500 400]);
yline(0, '--');

% Create and open video file
saveFolder = "C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\CSM for ControlAnimals\CSMIndividualReaches_animals";  
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end
filename = fullfile(saveFolder, 'ReachAnimation.mp4');
v = VideoWriter(filename, 'MPEG-4');

v.FrameRate = 30;   % faster playback (can set to 60 if you want)
open(v);

for i = 1:nTraces
    idx = allidx(i);

    % Choose T1 or T2
    if ismember(idx, T2idx)
        RR_pos = T2.stimvel(:,2,idx);
        RR_accel = T2.stimvel(:,8,idx);
    else
        RR_pos = T1.stimvel(:,2,idx);
        RR_accel = T1.stimvel(:,8,idx);
    end

    % Smooth
    RR_pos_smooth = smoothdata(RR_pos, 'gaussian', 8);
    RR_accel_smooth = smoothdata(RR_accel, 'gaussian', 8);

    RR_pos_smooth = RR_pos_smooth * 10;   % Convert position to cm
    RR_accel_smooth= RR_accel_smooth * 10; % Convert accel to cm/s^2

    % Find zero-crossing after t = 0
    s_post = RR_accel_smooth(t0_idx:end);
    t_post = times(t0_idx:end);
    sign_change_idx = find(diff(sign(s_post)) ~= 0, 1, 'first');

    if ~isempty(sign_change_idx) && sign_change_idx < length(t_post)
        x1 = t_post(sign_change_idx);
        x2 = t_post(sign_change_idx + 1);
        y1 = s_post(sign_change_idx);
        y2 = s_post(sign_change_idx + 1);
        x_cross = x1 - y1 * (x2 - x1) / (y2 - y1);
        [~, idx_cross] = min(abs(times - x_cross));
    else
        x_cross = inf;
        idx_cross = length(times);
        fprintf('No zero crossing found for trace %d\n', idx);
    end

    % -----------------------------
    % Auto-compute skip factor so each trace lasts ~2 seconds
    traceDuration = 2;                        % desired seconds per trace
    nFramesPerTrace = traceDuration * v.FrameRate;
    skip = ceil(length(times) / nFramesPerTrace);
    % -----------------------------

    % Animate line segment by segment (with skipping)
    for t = 2:skip:length(times)
        current_time = times(t);
        if current_time > 0 && current_time > x_cross
            c = afterColor;
        else
            c = beforeColor;
        end

        subplot(2,1,1);
        plot(times(t-1:t), RR_pos_smooth(t-1:t), '-', 'Color', c, 'LineWidth', 2);

        subplot(2,1,2);
        plot(times(t-1:t), RR_accel_smooth(t-1:t), '-', 'Color', c, 'LineWidth', 2);

        % Write this frame to video
        frame = getframe(gcf);
        writeVideo(v, frame);
    end

    % Mark zero crossing (optional)
    if isfinite(x_cross)
        subplot(2,1,1);
        plot(x_cross, RR_pos_smooth(idx_cross), 'mo', 'MarkerSize', 6, ...
             'LineWidth', 1.5, 'MarkerFaceColor', 'm');

        subplot(2,1,2);
        plot(x_cross, 0, 'mo', 'MarkerSize', 6, ...
             'LineWidth', 1.5, 'MarkerFaceColor', 'm');

        % Save this marking as a frame
        frame = getframe(gcf);
        writeVideo(v, frame);
    end
end

% Close video
close(v);
disp('Animation saved to ReachAnimation.mp4');



%% Now just get a single image
% Create a figure with white background
fig = figure('Position', [100, 100, 600, 800], 'Color', 'w');

% Define the indices
T1idx = [5,6,11,18]';
T2idx = [7,8,10,16,19]';
allidx = [T2idx; T1idx];
nTraces = length(allidx);

% Create a color gradient for the traces (light lavender to deep purple)
purples = [linspace(0.85, 0.4, nTraces)', linspace(0.7, 0.1, nTraces)', linspace(1.0, 0.6, nTraces)'];

% Get the shared time vector
times = T2.stimvel(:,1,2);
t0_idx = find(times >= 0, 1, 'first');

% Create the position and acceleration plots
for i = 1:nTraces
    idx = allidx(i);

    % Determine which struct to use (T1 or T2)
    if ismember(idx, T2idx)
        RR_pos = T2.stimvel(:,2,idx);  % Outward position (2nd column)
        RR_accel = T2.stimvel(:,8,idx); % Outward acceleration (8th column)
    else
        RR_pos = T1.stimvel(:,2,idx);  % Outward position (2nd column)
        RR_accel = T1.stimvel(:,8,idx); % Outward acceleration (8th column)
    end

    % Smooth data
    RR_pos_smooth = smoothdata(RR_pos, 'gaussian', 8);
    RR_accel_smooth = smoothdata(RR_accel, 'gaussian', 8);

    % Convert to centimeters
    RR_pos_smooth_cm = RR_pos_smooth * 10;
    RR_accel_smooth_cm = RR_accel_smooth * 10;

    % Find zero-crossing after t = 0
    s_post = RR_accel_smooth_cm(t0_idx:end);  % Acceleration after threshold
    t_post = times(t0_idx:end);                % Time after threshold
    sign_change_idx = find(diff(sign(s_post)) ~= 0, 1, 'first');  % First zero crossing

    if ~isempty(sign_change_idx) && sign_change_idx < length(t_post)
        x1 = t_post(sign_change_idx);
        x2 = t_post(sign_change_idx + 1);
        y1 = s_post(sign_change_idx);
        y2 = s_post(sign_change_idx + 1);
        x_cross = x1 - y1 * (x2 - x1) / (y2 - y1);  % Interpolation for zero crossing
        [~, idx_cross] = min(abs(times - x_cross));
    else
        x_cross = inf;
        idx_cross = length(times);
        fprintf('No zero crossing found for trace %d\n', idx);
    end

    % Plot before and after zero-crossing with same color gradient
    subplot(2, 1, 1); % Position plot
    before_cross_idx = times < x_cross;
    after_cross_idx = times >= x_cross;

    % Plot the position before the crossing (gradient color)
    plot(times(before_cross_idx), RR_pos_smooth_cm(before_cross_idx), 'Color', purples(i,:), 'LineWidth', 2);
    hold on;
    % Plot the position after the crossing (same gradient color)
    plot(times(after_cross_idx), RR_pos_smooth_cm(after_cross_idx), 'Color', purples(i,:), 'LineWidth', 2);

    subplot(2, 1, 2); % Acceleration plot
    % Plot the acceleration before the crossing (gradient color)
    plot(times(before_cross_idx), RR_accel_smooth_cm(before_cross_idx), 'Color', purples(i,:), 'LineWidth', 2);
    hold on;
    % Plot the acceleration after the crossing (same gradient color)
    plot(times(after_cross_idx), RR_accel_smooth_cm(after_cross_idx), 'Color', purples(i,:), 'LineWidth', 2);

    % Plot the zero-crossing point on both plots
    if isfinite(x_cross)
        subplot(2, 1, 1);
        plot(x_cross, RR_pos_smooth_cm(idx_cross), 'mo', 'MarkerSize', 6, 'LineWidth', 1.5);
        subplot(2, 1, 2);
        plot(x_cross, 0, 'mo', 'MarkerSize', 6, 'LineWidth', 1.5); % Zero crossing on accel plot
    end
end

% Format position plot
subplot(2, 1, 1);
yline(0, '--');
xlabel('Time from threshold crossing (s)');
ylabel('Outward Position (cm)');
xlim([-0.15 0.2]);
ylim([0 3.0]);  % Set the y-limit for position (0 to 3.0 cm)
title('Outward Position');

% Format acceleration plot
subplot(2, 1, 2);
yline(0, '--');
xlabel('Time from threshold crossing (s)');
ylabel('Outward Acceleration (cm/s^2)');
xlim([-0.15 0.2]);
ylim([-600 400]);  % Set the y-limit for acceleration (-600 to 400 cm/s^2)
title('Outward Acceleration');

%% Now plot single reach over time 

T2idx=[7,8,10,16,19]';

time=T2.stimvel(:,1,1);
outpos=T2.stimvel(:,2,16);
uppos=T2.stimvel(:,3,16);

figure; hold on
plot(outpos, uppos)
%%
T2idx = [7, 8, 10, 16, 19]';

time = T2.stimvel(:,1,1);
outpos = T2.stimvel(:,2,16);
uppos = T2.stimvel(:,3,16);

% Create colormap
cmap = jet(length(time));

% Normalize time for coloring
time_norm = (time - min(time)) / (max(time) - min(time));
color_idx = round(time_norm * (length(cmap)-1)) + 1;

figure; hold on
for i = 1:length(time)-1
    plot(outpos(i:i+1), uppos(i:i+1), 'Color', cmap(color_idx(i), :), 'LineWidth', 2);
end
colormap(jet);
colorbar;
title('Trajectory Colored by Time');
xlabel('Outward Position');
ylabel('Upward Position');
