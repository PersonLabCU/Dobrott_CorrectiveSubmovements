%% CODE FOR LOOKING AT CSM TIMES AND PLOTTING KDES FOR FIRST, MIDDLE AND LAST TRIALS 

%% Get RR times in order - First control session
% First have to get the times for all the T0
T0_divide=[];
thresh=0.16;
an=fieldnames(ExcitRR_Ctrl1);
for a=1:length(an)
    t0r=ExcitRR_Ctrl1.(an{a}).RRT0;

        % Skip animals that have no T0 reaches
    if isempty(t0r)
        T0_divide.(an{a}).ReReach = [];
        T0_divide.(an{a}).StimRR = [];
        T0_divide.(an{a}).NoJerk = [];
        continue
    end

    t0r=accel_toreaches_jerk(t0r);
    [ReReach, StimRR, MaxVel,ReachesWithoutSustainedJerk] = CSMEarlyThresh_jerk_edit5(t0r,thresh);
    %[ReReach, StimRR, MaxVel] = ReReachFilter_accel_T0_edit(t0r);
    T0_divide.(an{a}).ReReach=ReReach;
    T0_divide.(an{a}).StimRR=StimRR;
    T0_divide.(an{a}).NoJerk=ReachesWithoutSustainedJerk;
end
%% Replace 0s in StimRR with the corresponding times from T0_divide.StimRR - FIRST SESSION CONTROL 
for b = 1:length(an)
    RR_Excit = ExcitRR_Ctrl1.(an{b}).StimRR(:,1);

   % Skip animals with no T0 times
    if ~isfield(T0_divide, an{b}) || isempty(T0_divide.(an{b}).StimRR)
        ExcitRR_Ctrl1.(an{b}).StimRR_replace = RR_Excit; % keep as is
        continue
    end

    T0_vals  = T0_divide.(an{b}).StimRR;  % contains replacement times
    
    % Find linear indices of all zero entries in RR_Excit
    zero_idx = find(RR_Excit == 0);
    
    % Make sure the number of replacement values matches the number of zeros
    n_replace = min(length(zero_idx), numel(T0_vals));
    
    % Replace 0's in order with the times from T0_divide
    RR_Excit(zero_idx(1:n_replace)) = T0_vals(1:n_replace);
    
    % Save the edited version
    ExcitRR_Ctrl1.(an{b}).StimRR_replace = RR_Excit;
end


%% Get RR times in order - second session
% First have to get the times for all the T0
T0_divide_2=[];
thresh=0.16;
an=fieldnames(ExcitRR_Ctrl2);
for a=1:length(an)
    t0r=ExcitRR_Ctrl2.(an{a}).RRT0;

        % Skip animals that have no T0 reaches
    if isempty(t0r)
        T0_divide_2.(an{a}).ReReach = [];
        T0_divide_2.(an{a}).StimRR = [];
        T0_divide_2.(an{a}).NoJerk = [];
        continue
    end

    t0r=accel_toreaches_jerk(t0r);
    [ReReach, StimRR, MaxVel,ReachesWithoutSustainedJerk] = CSMEarlyThresh_jerk_edit5(t0r,thresh);
    %[ReReach, StimRR, MaxVel] = ReReachFilter_accel_T0_edit(t0r);
    T0_divide_2.(an{a}).ReReach=ReReach;
    T0_divide_2.(an{a}).StimRR=StimRR;
    T0_divide_2.(an{a}).NoJerk=ReachesWithoutSustainedJerk;
end
%% Replace 0s in StimRR with the corresponding times from T0_divide.StimRR - -Second session
for b = 1:length(an)
    RR_Excit = ExcitRR_Ctrl2.(an{b}).StimRR(:,1);

   % Skip animals with no T0 times
    if ~isfield(T0_divide_2, an{b}) || isempty(T0_divide_2.(an{b}).StimRR)
        ExcitRR_Ctrl2.(an{b}).StimRR_replace = RR_Excit; % keep as is
        continue
    end

    T0_vals  = T0_divide_2.(an{b}).StimRR;  % contains replacement times
    
    % Find linear indices of all zero entries in RR_Excit
    zero_idx = find(RR_Excit == 0);
    
    % Make sure the number of replacement values matches the number of zeros
    n_replace = min(length(zero_idx), numel(T0_vals));
    
    % Replace 0's in order with the times from T0_divide
    RR_Excit(zero_idx(1:n_replace)) = T0_vals(1:n_replace);
    
    % Save the edited version
    ExcitRR_Ctrl2.(an{b}).StimRR_replace = RR_Excit;
end

%% Get RR times in order -third session
% First have to get the times for all the T0
T0_divide_3=[];
thresh=0.16;
an=fieldnames(ExcitRR_Ctrl3);
for a=1:length(an)
    t0r=ExcitRR_Ctrl3.(an{a}).RRT0;

        % Skip animals that have no T0 reaches
    if isempty(t0r)
        T0_divide_3.(an{a}).ReReach = [];
        T0_divide_3.(an{a}).StimRR = [];
        T0_divide_3.(an{a}).NoJerk = [];
        continue
    end

    t0r=accel_toreaches_jerk(t0r);
    [ReReach, StimRR, MaxVel,ReachesWithoutSustainedJerk] = CSMEarlyThresh_jerk_edit5(t0r,thresh);
    %[ReReach, StimRR, MaxVel] = ReReachFilter_accel_T0_edit(t0r);
    T0_divide_3.(an{a}).ReReach=ReReach;
    T0_divide_3.(an{a}).StimRR=StimRR;
    T0_divide_3.(an{a}).NoJerk=ReachesWithoutSustainedJerk;
end
%% Replace 0s in StimRR with the corresponding times from T0_divide.StimRR - FIRST SESSION CONTROL -Second session
for b = 1:length(an)
    RR_Excit = ExcitRR_Ctrl3.(an{b}).StimRR(:,1);

   % Skip animals with no T0 times
    if ~isfield(T0_divide_3, an{b}) || isempty(T0_divide_3.(an{b}).StimRR)
        ExcitRR_Ctrl3.(an{b}).StimRR_replace = RR_Excit; % keep as is
        continue
    end

    T0_vals  = T0_divide_3.(an{b}).StimRR;  % contains replacement times
    
    % Find linear indices of all zero entries in RR_Excit
    zero_idx = find(RR_Excit == 0);
    
    % Make sure the number of replacement values matches the number of zeros
    n_replace = min(length(zero_idx), numel(T0_vals));
    
    % Replace 0's in order with the times from T0_divide
    RR_Excit(zero_idx(1:n_replace)) = T0_vals(1:n_replace);
    
    % Save the edited version
    ExcitRR_Ctrl3.(an{b}).StimRR_replace = RR_Excit;
end
%% Get RR times in order - LAST CONTROL SESSION
T0_divide_ctrlend = [];

an = fieldnames(ExcitRR_Ctrlend);
for a = 1:length(an)
    t0r = ExcitRR_Ctrlend.(an{a}).RRT0;

    % Skip animals that have no T0 reaches
    if isempty(t0r)
        T0_divide_ctrlend.(an{a}).ReReach = [];
        T0_divide_ctrlend.(an{a}).StimRR = [];
        T0_divide_ctrlend.(an{a}).NoJerk = [];
        continue
    end

    % Process animals with T0 reaches
    t0r = accel_toreaches_jerk(t0r);
    [ReReach, StimRR, MaxVel, ReachesWithoutSustainedJerk] = ...
        CSMEarlyThresh_jerk_edit5(t0r, thresh);

    T0_divide_ctrlend.(an{a}).ReReach = ReReach;
    T0_divide_ctrlend.(an{a}).StimRR = StimRR;
    T0_divide_ctrlend.(an{a}).NoJerk = ReachesWithoutSustainedJerk;
end

%% Replace 0s in StimRR with the corresponding times from T0_divide_ctrlend.StimRR - LAST SESSION CONTROL
for b = 1:length(an)
    RR_Excit = ExcitRR_Ctrlend.(an{b}).StimRR(:,1);

    % Skip animals with no T0 times
    if ~isfield(T0_divide_ctrlend, an{b}) || isempty(T0_divide_ctrlend.(an{b}).StimRR)
        ExcitRR_Ctrlend.(an{b}).StimRR_replace = RR_Excit; % keep as is
        continue
    end

    T0_vals = T0_divide_ctrlend.(an{b}).StimRR;

    % Find linear indices of all zero entries in RR_Excit
    zero_idx = find(RR_Excit == 0);

    % Make sure the number of replacement values matches the number of zeros
    n_replace = min(length(zero_idx), numel(T0_vals));

    % Replace 0's in order with the times from T0_divide
    RR_Excit(zero_idx(1:n_replace)) = T0_vals(1:n_replace);

    % Save the edited version
    ExcitRR_Ctrlend.(an{b}).StimRR_replace = RR_Excit;
end




%% Now plot CSM times 

groups_to_split = {'ExcitRR_Ctrl1', 'ExcitRR_Ctrlend','ExcitRR_Ctrl2','ExcitRR_Ctrl3'}; % leave Sham alone

for g = 1:length(groups_to_split)
    group = groups_to_split{g};
    animals = fieldnames(eval(group));
    
    for a = 1:length(animals)
        dat = eval([group '.(animals{a}).StimRR_replace']);
        
        % Skip if not enough reaches
        if numel(dat) < 5
            warning('%s: not enough reaches for animal %s', group, animals{a});
            continue
        end
        
        % First 5
        first5 = dat(1:min(5,numel(dat)));
        
        % Middle 5
        mid_idx = floor(numel(dat)/2)-2 : floor(numel(dat)/2)+2;
        mid_idx = mid_idx(mid_idx >= 1 & mid_idx <= numel(dat)); % bounds check
        middle5 = dat(mid_idx);
        
        % Last 5
        last5 = dat(max(1,end-4):end);
        
        % Save back into structure
        eval([group '.(animals{a}).StimRR_first5 = first5;']);
        eval([group '.(animals{a}).StimRR_middle5 = middle5;']);
        eval([group '.(animals{a}).StimRR_last5 = last5;']);
    end
end



%%

num_points = 100;
xi_common = linspace(-0.1, 0.2, num_points);

% Colors
stimcolor = rgb('OrangeRed'); % ExcitRR
ctrlcolor = rgb('Orange');    % ExcitRR_Ctrl
shamcolor = rgb('Black');     % Sham (dashed)


% Step 3: KDE for all reaches combined

figure; hold on;

% ExcitRR_Ctrl1 (first session)
all_data_ctrl1 = [];
animals_ctrl1 = fieldnames(ExcitRR_Ctrl1);
for a = 1:length(animals_ctrl1)
    dat = ExcitRR_Ctrl1.(animals_ctrl1{a}).StimRR_replace;
    all_data_ctrl1 = [all_data_ctrl1; dat];
end
[f_ctrl1, xi_ctrl1] = ksdensity(all_data_ctrl1, 'Bandwidth', 0.008);
f_interp_ctrl1 = interp1(xi_ctrl1, f_ctrl1, xi_common, 'linear', 0);
plot(xi_common, f_interp_ctrl1, 'Color', ctrlcolor, 'LineWidth', 2);

% ExcitRR_Ctrlend (last session)
all_data_ctrlend = [];
animals_ctrlend = fieldnames(ExcitRR_Ctrlend);
for a = 1:length(animals_ctrlend)
    dat = ExcitRR_Ctrlend.(animals_ctrlend{a}).StimRR_replace;
    all_data_ctrlend = [all_data_ctrlend; dat];
end
[f_ctrlend, xi_ctrlend] = ksdensity(all_data_ctrlend, 'Bandwidth', 0.008);
f_interp_ctrlend = interp1(xi_ctrlend, f_ctrlend, xi_common, 'linear', 0);
plot(xi_common, f_interp_ctrlend, 'Color', stimcolor, 'LineWidth', 2);

% ExcitRR_Sham
all_data_sham = [];
animals_sham = fieldnames(ExcitRR_Sham);
for a = 1:length(animals_sham)
    dat = ExcitRR_Sham.(animals_sham{a}).AllBins_prethresh;
    all_data_sham = [all_data_sham; dat];
end
[f_sham, xi_sham] = ksdensity(all_data_sham, 'Bandwidth', 0.008);
f_interp_sham = interp1(xi_sham, f_sham, xi_common, 'linear', 0);
plot(xi_common, f_interp_sham, '--', 'Color', shamcolor, 'LineWidth', 2);

title('All reaches combined');
xlabel('Time from threshold crossing (s)');
ylabel('Density');
xlim([-0.1 0.2]);
ylim([ 0 14])
set(gca, 'FontSize', 12, 'LineWidth', 1);
hold off;


%% ============================================================
% BUILD 60-POINT MEAN CURVE ACROSS ANIMALS
%
% Sessions:
%   ExcitRR_Ctrl1   = first session
%   ExcitRR_Ctrl2   = second session
%   ExcitRR_Ctrl3   = third session
%   ExcitRR_Ctrlend = last session
%
% Within each session:
%   1:5   = first 5 reaches
%   6:10  = middle 5 reaches
%   11:15 = last 5 reaches
%
% Total x-axis:
%   1:15   = Ctrl1
%   16:30  = Ctrl2
%   31:45  = Ctrl3
%   46:60  = Ctrlend
%
% One point = mean across animals for that reach
% Error bar = SEM across animals
%% ============================================================

group_names = {'ExcitRR_Ctrl1', 'ExcitRR_Ctrl2', 'ExcitRR_Ctrl3', 'ExcitRR_Ctrlend'};
nSessions = numel(group_names);
nPerChunk = 5;
nPerSession = 15;
nTotal = 60;

% Load groups
G = cell(1, nSessions);
for s = 1:nSessions
    G{s} = eval(group_names{s});
end

% Only keep animals present in all four groups
animals = fieldnames(G{1});
for s = 2:nSessions
    animals = intersect(animals, fieldnames(G{s}));
end

nAnimals = numel(animals);

% Preallocate: rows = animals, cols = 60 reaches
all_data = nan(nAnimals, nTotal);

% ---------------- Reconstruct each 15-reach session block ----------------
for a = 1:nAnimals
    animal_name = animals{a};

    for s = 1:nSessions
        this_group = G{s};

        block_start = (s-1)*nPerSession + 1;
        block_end   = s*nPerSession;

        % Initialize this 15-reach session block
        session_vec = nan(1, nPerSession);

        % ---- First 5 ----
        if isfield(this_group.(animal_name), 'StimRR_first5')
            first5 = this_group.(animal_name).StimRR_first5;
            first5 = first5(:)'; % force row
            n_use = min(numel(first5), nPerChunk);
            session_vec(1:n_use) = first5(1:n_use);
        end

        % ---- Middle 5 ----
        if isfield(this_group.(animal_name), 'StimRR_middle5')
            middle5 = this_group.(animal_name).StimRR_middle5;
            middle5 = middle5(:)'; % force row
            n_use = min(numel(middle5), nPerChunk);
            session_vec(6:5+n_use) = middle5(1:n_use);
        end

        % ---- Last 5 ----
        if isfield(this_group.(animal_name), 'StimRR_last5')
            last5 = this_group.(animal_name).StimRR_last5;
            last5 = last5(:)'; % force row
            n_use = min(numel(last5), nPerChunk);
            session_vec(11:10+n_use) = last5(1:n_use);
        end

        % Save this session block into the full 60-reach matrix
        all_data(a, block_start:block_end) = session_vec;
    end
end

% ---------------- Mean and SEM across animals ----------------
mean_all = mean(all_data, 1, 'omitnan');
sem_all  = std(all_data, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(all_data), 1));

x_all = 1:nTotal;

% ---------------- Check where data are missing ----------------
n_per_reach = sum(~isnan(all_data), 1);

disp('Number of animals contributing to each reach:');
disp(n_per_reach);

disp('Indices with no data at all:');
disp(find(n_per_reach == 0));

% ---------------- Exponential fit ----------------
% Single exponential:
% y = A * exp(B*x) + C

exp_model = @(p, x) p(1) .* exp(p(2) .* x) + p(3);

valid_idx = ~isnan(mean_all) & ~isnan(sem_all);
x_fit = x_all(valid_idx)';
y_fit = mean_all(valid_idx)';

p_fit = [];
t_dense = linspace(min(x_fit), max(x_fit), 500)';
y_exp_dense = nan(size(t_dense));

try
    A_guess = max(y_fit) - min(y_fit);
    B_guess = -0.02;
    C_guess = min(y_fit);
    p0 = [A_guess, B_guess, C_guess];

    lb = [-Inf, -Inf, -Inf];
    ub = [ Inf,  Inf,  Inf];

    opts = optimoptions('lsqcurvefit', 'Display', 'off');
    p_fit = lsqcurvefit(exp_model, p0, x_fit, y_fit, lb, ub, opts);

    y_exp_dense = exp_model(p_fit, t_dense);
catch
    warning('Exponential fit failed.');
end

%% ---------------- Plot ----------------
figure('Color','w','Position',[120 120 1300 450]); hold on

% Plot only valid mean points with SEM
errorbar(x_all(valid_idx), mean_all(valid_idx), sem_all(valid_idx), 'ko', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 5, ...
    'MarkerFaceColor', 'w', ...
    'CapSize', 0);

% Connect only valid points
plot(x_all(valid_idx), mean_all(valid_idx), 'k-', 'LineWidth', 1.2);

% Exponential fit
if ~isempty(p_fit)
    plot(t_dense, y_exp_dense, 'r-', 'LineWidth', 2.5);
end

% Session dividers
xline(15.5, ':k', 'LineWidth', 1.2);
xline(30.5, ':k', 'LineWidth', 1.2);
xline(45.5, ':k', 'LineWidth', 1.2);

xlabel('Reach number');
ylabel('Mean CSM time (s)');
title('Mean CSM time across animals across four sessions');
xlim([1 60]);
xticks(1:5:60);
box off;

legend({'Mean \pm SEM', 'Mean curve', 'Exponential fit'}, 'Location', 'best');

% ---------------- Session labels ----------------
yl = ylim;
text(7.5,  yl(2), 'Ctrl1',   'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(22.5, yl(2), 'Ctrl2',   'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(37.5, yl(2), 'Ctrl3',   'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(52.5, yl(2), 'Ctrlend', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% ---------------- Print fit results ----------------
fprintf('\n================ SINGLE EXPONENTIAL FIT =================\n');
if ~isempty(p_fit)
    fprintf('Model: y = A * exp(B*x) + C\n');
    fprintf('A = %.6f\n', p_fit(1));
    fprintf('B = %.6f\n', p_fit(2));
    fprintf('C = %.6f\n', p_fit(3));
else
    fprintf('Exponential fit failed.\n');
end


%% ============================================================
% PLOT FIRST, SECOND, AND FINAL SESSION (45 TOTAL REACHES)
%% ============================================================

group_names = {'ExcitRR_Ctrl1', 'ExcitRR_Ctrl2', 'ExcitRR_Ctrlend'};
nSessions = numel(group_names);

nPerChunk = 5;
nPerSession = 15;
nTotal = nSessions * nPerSession;

% Load groups
G = cell(1, nSessions);
for s = 1:nSessions
    G{s} = eval(group_names{s});
end

% Only animals present in all sessions
animals = fieldnames(G{1});
for s = 2:nSessions
    animals = intersect(animals, fieldnames(G{s}));
end

nAnimals = numel(animals);

% Preallocate
all_data = nan(nAnimals, nTotal);

% ---------------- Build data matrix ----------------
for a = 1:nAnimals
    animal = animals{a};

    for s = 1:nSessions
        this_group = G{s};

        block_start = (s-1)*nPerSession + 1;

        % Initialize 15-length session
        session_vec = nan(1, nPerSession);

        % ---- First 5 ----
        if isfield(this_group.(animal), 'StimRR_first5')
            v = this_group.(animal).StimRR_first5;
            v = v(:)';
            session_vec(1:min(5,numel(v))) = v(1:min(5,numel(v)));
        end

        % ---- Middle 5 ----
        if isfield(this_group.(animal), 'StimRR_middle5')
            v = this_group.(animal).StimRR_middle5;
            v = v(:)';
            session_vec(6:5+min(5,numel(v))) = v(1:min(5,numel(v)));
        end

        % ---- Last 5 ----
        if isfield(this_group.(animal), 'StimRR_last5')
            v = this_group.(animal).StimRR_last5;
            v = v(:)';
            session_vec(11:10+min(5,numel(v))) = v(1:min(5,numel(v)));
        end

        % Store
        all_data(a, block_start:block_start+14) = session_vec;
    end
end

% ---------------- Mean + SEM ----------------
mean_all = mean(all_data, 1, 'omitnan');
sem_all  = std(all_data, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(all_data),1));

%x = 1:nTotal;
x = 0:(nTotal-1);

% ---------------- Exponential fit ----------------
exp_model = @(p,x) p(1).*exp(p(2).*x) + p(3);

valid = ~isnan(mean_all);
x_fit = x(valid)';
y_fit = mean_all(valid)';

p_fit = [];
t_dense = linspace(min(x_fit), max(x_fit), 500)';
y_dense = nan(size(t_dense));

try
    p0 = [range(y_fit), -0.02, min(y_fit)];
    p_fit = lsqcurvefit(exp_model, p0, x_fit, y_fit, [], [], ...
        optimoptions('lsqcurvefit','Display','off'));

    y_dense = exp_model(p_fit, t_dense);
catch
    warning('Fit failed');
end

% ---------------- Plot ----------------
figure('Color','w','Position',[200 200 1100 450]); hold on
% 
% % Mean ± SEM
errorbar(x(valid), mean_all(valid), sem_all(valid), 'ko', ...
    'LineWidth', 1.5, 'MarkerFaceColor','w','CapSize',0);

% Connect points
plot(x(valid), mean_all(valid), 'k-', 'LineWidth', 1.2);

% plot(x(valid), mean_all(valid), 'ko', ...
%     'MarkerFaceColor','w', ...
%     'MarkerSize',6);

% Exponential fit
if ~isempty(p_fit)
    plot(t_dense, y_dense, 'r-', 'LineWidth', 2.5);
end

% Session dividers
xline(15.5, '--k');
xline(30.5, '--k');

xlabel('Reach number');
ylabel('Mean CSM time (s)');
title('Mean CSM time across animals (Ctrl1, Ctrl2, Ctrlend)');
xlim([1 45]);
%xlim([0 45])
%xticks(1:5:45);
ylim([-0.1 0.2])
box off;

legend({'Mean \pm SEM','Mean','Exp fit'}, 'Location','best');

% ---------------- Labels ----------------
yl = ylim;
text(7.5,  yl(2), 'Ctrl1',  'HorizontalAlignment','center','VerticalAlignment','bottom');
text(22.5, yl(2), 'Ctrl2',  'HorizontalAlignment','center','VerticalAlignment','bottom');
text(37.5, yl(2), 'Ctrlend','HorizontalAlignment','center','VerticalAlignment','bottom');

% ---------------- Print fit ----------------
fprintf('\n==== EXPONENTIAL FIT ====\n');
if ~isempty(p_fit)
    fprintf('y = %.5f * exp(%.5f*x) + %.5f\n', p_fit(1), p_fit(2), p_fit(3));
end

%% Slightly different version
%% ============================================================
% PLOT FIRST, SECOND, AND FINAL SESSION (45 TOTAL REACHES)
% USING BI-EXPONENTIAL FIT VIA fit_biexponential3
%
% Sessions:
%   Ctrl1   = first session
%   Ctrl2   = second session
%   Ctrlend = final session
%
% Within each session:
%   1:5   = first 5 reaches
%   6:10  = middle 5 reaches
%   11:15 = last 5 reaches
%
% Total x-axis:
%   1:15   = Ctrl1
%   16:30  = Ctrl2
%   31:45  = Ctrlend
%
% One point = mean across animals for that reach
%% ============================================================

group_names = {'ExcitRR_Ctrl1', 'ExcitRR_Ctrl2', 'ExcitRR_Ctrlend'};
nSessions   = numel(group_names);

nPerChunk   = 5;
nPerSession = 15;
nTotal      = nSessions * nPerSession;

% Load groups
G = cell(1, nSessions);
for s = 1:nSessions
    G{s} = eval(group_names{s});
end

% Keep only animals present in all three groups
animals = fieldnames(G{1});
for s = 2:nSessions
    animals = intersect(animals, fieldnames(G{s}));
end

nAnimals = numel(animals);

% Preallocate full data matrix
all_data = nan(nAnimals, nTotal);

% ---------------- Build animal x reach matrix ----------------
for a = 1:nAnimals
    animal_name = animals{a};

    for s = 1:nSessions
        this_group = G{s};
        block_start = (s-1)*nPerSession + 1;

        % Initialize this 15-reach session block
        session_vec = nan(1, nPerSession);

        % ---- First 5 ----
        if isfield(this_group.(animal_name), 'StimRR_first5')
            v = this_group.(animal_name).StimRR_first5;
            v = v(:)';   % force row
            n_use = min(numel(v), nPerChunk);
            session_vec(1:n_use) = v(1:n_use);
        end

        % ---- Middle 5 ----
        if isfield(this_group.(animal_name), 'StimRR_middle5')
            v = this_group.(animal_name).StimRR_middle5;
            v = v(:)';   % force row
            n_use = min(numel(v), nPerChunk);
            session_vec(6:5+n_use) = v(1:n_use);
        end

        % ---- Last 5 ----
        if isfield(this_group.(animal_name), 'StimRR_last5')
            v = this_group.(animal_name).StimRR_last5;
            v = v(:)';   % force row
            n_use = min(numel(v), nPerChunk);
            session_vec(11:10+n_use) = v(1:n_use);
        end

        % Store this session block
        all_data(a, block_start:block_start+nPerSession-1) = session_vec;
    end
end

% ---------------- Mean across animals per reach ----------------
mean_all = mean(all_data, 1, 'omitnan');
x = 1:nTotal;

% ---------------- Bi-exponential fit ----------------
valid = ~isnan(mean_all);
t_fit = x(valid)';
y_fit_data = mean_all(valid)';

[tau_fast, tau_slow, fit_params, y_fit, t_dense, y_fit_dense] = fit_biexponential3(t_fit, y_fit_data);

% ---------------- Plot ----------------
figure('Color','w','Position',[200 200 1100 450]); 
hold on

% % Mean points only
% plot(x(valid), mean_all(valid), 'ko', ...
%     'MarkerFaceColor', 'w', ...
%     'MarkerSize', 6, ...
%     'DisplayName', 'Mean across animals');
plot(x(valid), mean_all(valid), 'ko', ...
    'MarkerFaceColor','k', ...
    'MarkerSize',6);

% Smooth bi-exponential fit
plot(t_dense, y_fit_dense, 'r-', 'LineWidth', 2.5, ...
    'DisplayName', sprintf('Bi-exp fit: \\tau_{fast}=%.3f, \\tau_{slow}=%.3f', tau_fast, tau_slow));

% Session dividers
xline(15.5, '--k', 'LineWidth', 1.2);
xline(30.5, '--k', 'LineWidth', 1.2);
yline(0.0581,'k')

xlabel('Reach number');
ylabel('Mean CSM time (s)');
title('Mean CSM time across animals (Ctrl1, Ctrl2, Ctrlend)');

xlim([1 45]);
ylim([-0.1 0.2])
%xticks(1:5:45);
xticks(5:5:45)
xticklabels(string(5:5:45))
box off;
legend('Location', 'best');

% Session labels
yl = ylim;
text(7.5,  yl(2), 'Ctrl1',   'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(22.5, yl(2), 'Ctrl2',   'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(37.5, yl(2), 'Ctrlend', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% ---------------- Print fit results ----------------
fprintf('\n================ BI-EXPONENTIAL FIT =================\n');
fprintf('A_fast   = %.6f\n', fit_params.A_fast);
fprintf('tau_fast = %.6f\n', fit_params.tau_fast);
fprintf('A_slow   = %.6f\n', fit_params.A_slow);
fprintf('tau_slow = %.6f\n', fit_params.tau_slow);
fprintf('C        = %.6f\n', fit_params.C);

if isfield(fit_params, 'R2')
    fprintf('R2       = %.6f\n', fit_params.R2);
end
if isfield(fit_params, 'RMSE')
    fprintf('RMSE     = %.6f\n', fit_params.RMSE);
end

%% Single exponential 

%% Slightly different version
%% ============================================================
% PLOT FIRST, SECOND, AND FINAL SESSION (45 TOTAL REACHES)
% USING SINGLE EXPONENTIAL FIT
%
% Sessions:
%   Ctrl1   = first session
%   Ctrl2   = second session
%   Ctrlend = final session
%
% Within each session:
%   1:5   = first 5 reaches
%   6:10  = middle 5 reaches
%   11:15 = last 5 reaches
%
% Total x-axis:
%   1:15   = Ctrl1
%   16:30  = Ctrl2
%   31:45  = Ctrlend
%
% One point = mean across animals for that reach
%% ============================================================

group_names = {'ExcitRR_Ctrl1', 'ExcitRR_Ctrl2', 'ExcitRR_Ctrlend'};
nSessions   = numel(group_names);

nPerChunk   = 5;
nPerSession = 15;
nTotal      = nSessions * nPerSession;

% Load groups
G = cell(1, nSessions);
for s = 1:nSessions
    G{s} = eval(group_names{s});
end

% Keep only animals present in all three groups
animals = fieldnames(G{1});
for s = 2:nSessions
    animals = intersect(animals, fieldnames(G{s}));
end

nAnimals = numel(animals);

% Preallocate full data matrix
all_data = nan(nAnimals, nTotal);

% ---------------- Build animal x reach matrix ----------------
for a = 1:nAnimals
    animal_name = animals{a};

    for s = 1:nSessions
        this_group = G{s};
        block_start = (s-1)*nPerSession + 1;

        % Initialize this 15-reach session block
        session_vec = nan(1, nPerSession);

        % ---- First 5 ----
        if isfield(this_group.(animal_name), 'StimRR_first5')
            v = this_group.(animal_name).StimRR_first5;
            v = v(:)';   % force row
            n_use = min(numel(v), nPerChunk);
            session_vec(1:n_use) = v(1:n_use);
        end

        % ---- Middle 5 ----
        if isfield(this_group.(animal_name), 'StimRR_middle5')
            v = this_group.(animal_name).StimRR_middle5;
            v = v(:)';   % force row
            n_use = min(numel(v), nPerChunk);
            session_vec(6:5+n_use) = v(1:n_use);
        end

        % ---- Last 5 ----
        if isfield(this_group.(animal_name), 'StimRR_last5')
            v = this_group.(animal_name).StimRR_last5;
            v = v(:)';   % force row
            n_use = min(numel(v), nPerChunk);
            session_vec(11:10+n_use) = v(1:n_use);
        end

        % Store this session block
        all_data(a, block_start:block_start+nPerSession-1) = session_vec;
    end
end

% ---------------- Mean across animals per reach ----------------
mean_all = mean(all_data, 1, 'omitnan');
x = 1:nTotal;

% ---------------- Single exponential fit ----------------
valid = ~isnan(mean_all);
t_fit = x(valid)';
y_fit_data = mean_all(valid)';

% Model: y = A * exp(B*x) + C
exp_model = @(p, t) p(1) .* exp(p(2) .* t) + p(3);

% Initial guesses
A_guess = y_fit_data(1) - y_fit_data(end);
B_guess = -0.02;
C_guess = y_fit_data(end);

p0 = [A_guess, B_guess, C_guess];
lb = [-Inf, -Inf, -Inf];
ub = [ Inf,  Inf,  Inf];

opts = optimoptions('lsqcurvefit', 'Display', 'off');
p_fit = lsqcurvefit(exp_model, p0, t_fit, y_fit_data, lb, ub, opts);

A = p_fit(1);
B = p_fit(2);
C = p_fit(3);

% Time constant
tau = inf;
if B ~= 0
    tau = -1 / B;
end

% Fitted values
y_fit = exp_model(p_fit, t_fit);
t_dense = linspace(min(t_fit), max(t_fit), 300)';
y_fit_dense = exp_model(p_fit, t_dense);

% Goodness of fit
residuals = y_fit_data - y_fit;
SS_res = sum(residuals.^2);
SS_tot = sum((y_fit_data - mean(y_fit_data)).^2);
R2 = 1 - SS_res / SS_tot;
RMSE = sqrt(mean(residuals.^2));

% ---------------- Plot ----------------
figure('Color','w','Position',[200 200 1100 450]); 
hold on

% Mean points only
plot(x(valid), mean_all(valid), 'ko', ...
    'MarkerFaceColor','k', ...
    'MarkerSize',6);

% Smooth single exponential fit
plot(t_dense, y_fit_dense, 'r-', 'LineWidth', 2.5, ...
    'DisplayName', sprintf('Single exp fit: \\tau=%.3f', tau));

% Session dividers
xline(15.5, '--k', 'LineWidth', 1.2);
xline(30.5, '--k', 'LineWidth', 1.2);
yline(0.0581,'k')

xlabel('Reach number');
ylabel('Mean CSM time (s)');
title('Mean CSM time across animals (Ctrl1, Ctrl2, Ctrlend)');

xlim([1 45]);
ylim([-0.1 0.2])
xticks(5:5:45)
xticklabels(string(5:5:45))
box off;
legend('Location', 'best');

% Session labels
yl = ylim;
text(7.5,  yl(2), 'Ctrl1',   'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(22.5, yl(2), 'Ctrl2',   'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(37.5, yl(2), 'Ctrlend', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% ---------------- Print fit results ----------------
fprintf('\n================ SINGLE EXPONENTIAL FIT =================\n');
fprintf('A   = %.6f\n', A);
fprintf('B   = %.6f\n', B);
fprintf('C   = %.6f\n', C);
fprintf('tau = %.6f\n', tau);
fprintf('Equation: y = %.6f * exp(%.6f*x) + %.6f\n', A, B, C);
fprintf('R2  = %.6f\n', R2);
fprintf('RMSE = %.6f\n', RMSE);


%% Step 4: KDE for first, middle, last 5 reaches - first session (Ctrl1)

figure; hold on;

segments = {'first5', 'middle5', 'last5'};
colors = lines(3); % 3 distinct colors

for s = 1:length(segments)
    all_segment = [];
    for a = 1:length(animals_ctrl1)
        dat = ExcitRR_Ctrl1.(animals_ctrl1{a}).(['StimRR_' segments{s}]);
        if ~isempty(dat)
            all_segment = [all_segment; dat];
        end
    end
    [f_seg, xi_seg] = ksdensity(all_segment, 'Bandwidth', 0.008);
    f_interp_seg = interp1(xi_seg, f_seg, xi_common, 'linear', 0);
    plot(xi_common, f_interp_seg, 'Color', colors(s,:), 'LineWidth', 2, 'DisplayName', segments{s});
end

title('First Session: First, Middle, Last 5 Reaches');
xlabel('Time from threshold crossing (s)');
ylabel('Density');
xlim([-0.1 0.2]);
%ylim([ 0 14])
legend('show');
set(gca, 'FontSize', 12, 'LineWidth', 1);
hold off;


%% Step 4: KDE for first, middle, last 5 reaches - first session (Ctrl1)
% saving and plotting KDE times 
figure; hold on;

segments = {'first5', 'middle5', 'last5'};
colors = lines(3); % 3 distinct colors
xi_common = linspace(-0.1, 0.2, 500); % define common support for interpolation

PeakInfo_Ctrl1 = struct(); % initialize structure

for s = 1:length(segments)
    all_segment = [];
    for a = 1:length(animals_ctrl1)
        dat = ExcitRR_Ctrl1.(animals_ctrl1{a}).(['StimRR_' segments{s}]);
        if ~isempty(dat)
            all_segment = [all_segment; dat];
        end
    end

    % --- Compute KDE ---
    [f_seg, xi_seg] = ksdensity(all_segment, 'Bandwidth', 0.008);

    % --- Interpolate to common x-axis ---
    f_interp_seg = interp1(xi_seg, f_seg, xi_common, 'linear', 0);

    % --- Find peak location and density ---
    [peak_val, idx_peak] = max(f_interp_seg);
    peak_loc = xi_common(idx_peak);

    % --- Save results in structure ---
    PeakInfo_Ctrl1.(segments{s}).peak_loc = peak_loc;
    PeakInfo_Ctrl1.(segments{s}).peak_val = peak_val;
    PeakInfo_Ctrl1.(segments{s}).xi = xi_common;
    PeakInfo_Ctrl1.(segments{s}).f = f_interp_seg;

    % --- Plot KDE ---
    plot(xi_common, f_interp_seg, 'Color', colors(s,:), 'LineWidth', 2, ...
        'DisplayName', segments{s});

    % --- Mark and label the peak ---
    plot(peak_loc, peak_val, 'o', 'MarkerEdgeColor', colors(s,:), ...
        'MarkerFaceColor', colors(s,:), 'HandleVisibility', 'off');
    text(peak_loc, peak_val + 0.05*peak_val, sprintf('%.3f', peak_loc), ...
        'Color', colors(s,:), 'FontSize', 10, 'HorizontalAlignment', 'center');
end

title('First Session: First, Middle, Last 5 Reaches');
xlabel('Time from threshold crossing (s)');
ylabel('Density');
xlim([-0.1 0.2]);
legend('show');
set(gca, 'FontSize', 12, 'LineWidth', 1);
hold off;

% --- Display summary table in console ---
disp(struct2table(PeakInfo_Ctrl1));

%% Step 5: KDE for first, middle, last 5 reaches - last session (Ctrlend)
%% -----------------------------
figure; hold on;

for s = 1:length(segments)
    all_segment = [];
    for a = 1:length(animals_ctrlend)
        dat = ExcitRR_Ctrlend.(animals_ctrlend{a}).(['StimRR_' segments{s}]);
        if ~isempty(dat)
            all_segment = [all_segment; dat];
        end
    end
    [f_seg, xi_seg] = ksdensity(all_segment, 'Bandwidth', 0.008);
    f_interp_seg = interp1(xi_seg, f_seg, xi_common, 'linear', 0);
    plot(xi_common, f_interp_seg, 'Color', colors(s,:), 'LineWidth', 2, 'DisplayName', segments{s});
end

title('Last Session: First, Middle, Last 5 Reaches');
xlabel('Time from threshold crossing (s)');
ylabel('Density');
xlim([-0.1 0.2])
ylim([ 0 14]);
legend('show');
set(gca, 'FontSize', 12, 'LineWidth', 1);
hold off;
%%
%% Step 4: KDE for first, middle, last 5 reaches - last session (CtrlEnd)

figure; hold on;

segments = {'first5', 'middle5', 'last5'};
colors = lines(3); % 3 distinct colors
xi_common = linspace(-0.1, 0.2, 500); % define common support for interpolation

PeakInfo_CtrlEnd = struct(); % initialize structure

for s = 1:length(segments)
    all_segment = [];
    for a = 1:length(animals_ctrlend)
        dat = ExcitRR_Ctrlend.(animals_ctrlend{a}).(['StimRR_' segments{s}]);
        if ~isempty(dat)
            all_segment = [all_segment; dat];
        end
    end

    % --- Compute KDE ---
    [f_seg, xi_seg] = ksdensity(all_segment, 'Bandwidth', 0.008);

    % --- Interpolate to common x-axis ---
    f_interp_seg = interp1(xi_seg, f_seg, xi_common, 'linear', 0);

    % --- Find peak location and density ---
    [peak_val, idx_peak] = max(f_interp_seg);
    peak_loc = xi_common(idx_peak);

    % --- Save results in structure ---
    PeakInfo_CtrlEnd.(segments{s}).peak_loc = peak_loc;
    PeakInfo_CtrlEnd.(segments{s}).peak_val = peak_val;
    PeakInfo_CtrlEnd.(segments{s}).xi = xi_common;
    PeakInfo_CtrlEnd.(segments{s}).f = f_interp_seg;

    % --- Plot KDE ---
    plot(xi_common, f_interp_seg, 'Color', colors(s,:), 'LineWidth', 2, ...
        'DisplayName', segments{s});

    % --- Mark and label the peak ---
    plot(peak_loc, peak_val, 'o', 'MarkerEdgeColor', colors(s,:), ...
        'MarkerFaceColor', colors(s,:), 'HandleVisibility', 'off');
    text(peak_loc, peak_val + 0.05*peak_val, sprintf('%.3f', peak_loc), ...
        'Color', colors(s,:), 'FontSize', 10, 'HorizontalAlignment', 'center');
end

title('Last Session: First, Middle, Last 5 Reaches');
xlabel('Time from threshold crossing (s)');
ylabel('Density');
xlim([-0.1 0.2])
ylim([0 14]);
legend('show');
set(gca, 'FontSize', 12, 'LineWidth', 1);
hold off;

% --- Display summary table in console ---
disp(struct2table(PeakInfo_CtrlEnd));

%% NOW adding sham 

figure; hold on;

segments = {'first5', 'middle5', 'last5'};
colors = lines(3); % 3 distinct colors for first/middle/last

% Ctrl1 KDEs
for s = 1:length(segments)
    all_segment = [];
    for a = 1:length(animals_ctrl1)
        dat = ExcitRR_Ctrl1.(animals_ctrl1{a}).(['StimRR_' segments{s}]);
        if ~isempty(dat)
            all_segment = [all_segment; dat];
        end
    end
    [f_seg, xi_seg] = ksdensity(all_segment, 'Bandwidth', 0.008);
    f_interp_seg = interp1(xi_seg, f_seg, xi_common, 'linear', 0);
    plot(xi_common, f_interp_seg, 'Color', colors(s,:), 'LineWidth', 2, 'DisplayName', segments{s});
end

% Sham KDE (all reaches combined)
all_data_sham = [];
for a = 1:length(animals_sham)
    dat = ExcitRR_Sham.(animals_sham{a}).AllBins_prethresh;
    all_data_sham = [all_data_sham; dat];
end
[f_sham, xi_sham] = ksdensity(all_data_sham, 'Bandwidth', 0.008);
f_interp_sham = interp1(xi_sham, f_sham, xi_common, 'linear', 0);
plot(xi_common, f_interp_sham, '--', 'Color', shamcolor, 'LineWidth', 2, 'DisplayName', 'Sham');

title('First Session: First, Middle, Last 5 Reaches + Sham');
xlabel('Time from threshold crossing (s)');
ylabel('Density');
xlim([-0.1 0.2]);
legend('show');
set(gca, 'FontSize', 12, 'LineWidth', 1);
hold off;

%% 
figure; hold on;

% Ctrlend KDEs
for s = 1:length(segments)
    all_segment = [];
    for a = 1:length(animals_ctrlend)
        dat = ExcitRR_Ctrlend.(animals_ctrlend{a}).(['StimRR_' segments{s}]);
        if ~isempty(dat)
            all_segment = [all_segment; dat];
        end
    end
    [f_seg, xi_seg] = ksdensity(all_segment, 'Bandwidth', 0.008);
    f_interp_seg = interp1(xi_seg, f_seg, xi_common, 'linear', 0);
    plot(xi_common, f_interp_seg, 'Color', colors(s,:), 'LineWidth', 2, 'DisplayName', segments{s});
end

% Sham KDE (all reaches combined)
plot(xi_common, f_interp_sham, '--', 'Color', shamcolor, 'LineWidth', 2, 'DisplayName', 'Sham');

title('Last Session: First, Middle, Last 5 Reaches + Sham');
xlabel('Time from threshold crossing (s)');
ylabel('Density');
xlim([-0.1 0.2]);
ylim([0 14])
legend('show');
set(gca, 'FontSize', 12, 'LineWidth', 1);
hold off;


%% STIM REACHES KDE PLOT

num_points = 100;
xi_common = linspace(-0.1, 0.2, num_points);

% Colors
stimcolor = rgb('OrangeRed'); % ExcitRR
ctrlcolor = rgb('Orange');    % ExcitRR_Ctrlend (last session)
shamcolor = rgb('Black');     % Sham (dashed)

figure; hold on;

% -----------------------------
% ExcitRR (first session)
% -----------------------------
all_data_first = [];
animals_first = fieldnames(ExcitRR);
for a = 1:length(animals_first)
    dat = ExcitRR.(animals_first{a}).StimRR(:,1);
    all_data_first = [all_data_first; dat];
end
[f_first, xi_first] = ksdensity(all_data_first, 'Bandwidth', 0.008);
f_interp_first = interp1(xi_first, f_first, xi_common, 'linear', 0);
plot(xi_common, f_interp_first, 'Color', stimcolor, 'LineWidth', 2, 'DisplayName', 'ExcitRR (First Session)');

% -----------------------------
% ExcitRR_Ctrlend (last session)
% -----------------------------
all_data_last = [];
animals_last = fieldnames(ExcitRR_end);
for a = 1:length(animals_last)
    dat = ExcitRR_end.(animals_last{a}).StimRR(:,1);
    all_data_last = [all_data_last; dat];
end
[f_last, xi_last] = ksdensity(all_data_last, 'Bandwidth', 0.008);
f_interp_last = interp1(xi_last, f_last, xi_common, 'linear', 0);
plot(xi_common, f_interp_last, 'Color', ctrlcolor, 'LineWidth', 2, 'DisplayName', 'ExcitRR_Ctrlend (Last Session)');

% -----------------------------
% Sham group
% -----------------------------
all_data_sham = [];
animals_sham = fieldnames(ExcitRR_Sham);
for a = 1:length(animals_sham)
    dat = ExcitRR_Sham.(animals_sham{a}).AllBins_prethresh;
    all_data_sham = [all_data_sham; dat];
end
[f_sham, xi_sham] = ksdensity(all_data_sham, 'Bandwidth', 0.008);
f_interp_sham = interp1(xi_sham, f_sham, xi_common, 'linear', 0);
plot(xi_common, f_interp_sham, '--', 'Color', shamcolor, 'LineWidth', 2, 'DisplayName', 'Sham');

% -----------------------------
% Figure formatting
% -----------------------------
title('Stim Reaches: First vs Last Session + Sham');
xlabel('Time from threshold crossing (s)');
ylabel('Density');
xlim([-0.1 0.2]);

legend('show');
set(gca, 'FontSize', 12, 'LineWidth', 1);
hold off;

%% -----------------------------
%% STIMULATED REACHES: First Session + Sham
%% -----------------------------
figure; hold on;

segments = {'first5', 'middle5', 'last5'};
colors = lines(3); % Colors for first/middle/last

% Ctrl1 (first session) KDEs
for s = 1:length(segments)
    all_segment = [];
    for a = 1:length(animals_ctrl1)
        dat = ExcitRR.(animals_ctrl1{a}).StimRR(:,1);  % Take only first column (CSM times)
        % Extract first/middle/last 5 for this animal
        switch segments{s}
            case 'first5'
                dat_seg = dat(1:min(5,numel(dat)));
            case 'middle5'
                mid_idx = floor(numel(dat)/2)-2 : floor(numel(dat)/2)+2;
                mid_idx = mid_idx(mid_idx >= 1 & mid_idx <= numel(dat));
                dat_seg = dat(mid_idx);
            case 'last5'
                dat_seg = dat(max(1,end-4):end);
        end
        if ~isempty(dat_seg)
            all_segment = [all_segment; dat_seg];
        end
    end
    [f_seg, xi_seg] = ksdensity(all_segment, 'Bandwidth', 0.008);
    f_interp_seg = interp1(xi_seg, f_seg, xi_common, 'linear', 0);
    plot(xi_common, f_interp_seg, 'Color', colors(s,:), 'LineWidth', 2, 'DisplayName', segments{s});
end

% Sham KDE (all reaches combined)
all_data_sham = [];
for a = 1:length(animals_sham)
    dat = ExcitRR_Sham.(animals_sham{a}).AllBins_prethresh;
    all_data_sham = [all_data_sham; dat];
end
[f_sham, xi_sham] = ksdensity(all_data_sham, 'Bandwidth', 0.008);
f_interp_sham = interp1(xi_sham, f_sham, xi_common, 'linear', 0);
plot(xi_common, f_interp_sham, '--', 'Color', shamcolor, 'LineWidth', 2, 'DisplayName', 'Sham');

title('First Session Stim Reaches: First, Middle, Last 5 + Sham');
xlabel('Time from threshold crossing (s)');
ylabel('Density');
xlim([-0.1 0.2]);
legend('show');
set(gca, 'FontSize', 12, 'LineWidth', 1);
hold off;

%% -----------------------------
%% STIMULATED REACHES: Last Session + Sham
%% -----------------------------
figure; hold on;

% Ctrlend (last session) KDEs
for s = 1:length(segments)
    all_segment = [];
    for a = 1:length(animals_ctrlend)
        dat = ExcitRR_end.(animals_ctrlend{a}).StimRR(:,1);  % Take only first column (CSM times)
        % Extract first/middle/last 5
        switch segments{s}
            case 'first5'
                dat_seg = dat(1:min(5,numel(dat)));
            case 'middle5'
                mid_idx = floor(numel(dat)/2)-2 : floor(numel(dat)/2)+2;
                mid_idx = mid_idx(mid_idx >= 1 & mid_idx <= numel(dat));
                dat_seg = dat(mid_idx);
            case 'last5'
                dat_seg = dat(max(1,end-4):end);
        end
        if ~isempty(dat_seg)
            all_segment = [all_segment; dat_seg];
        end
    end
    [f_seg, xi_seg] = ksdensity(all_segment, 'Bandwidth', 0.008);
    f_interp_seg = interp1(xi_seg, f_seg, xi_common, 'linear', 0);
    plot(xi_common, f_interp_seg, 'Color', colors(s,:), 'LineWidth', 2, 'DisplayName', segments{s});
end

% Sham KDE (all reaches combined)
plot(xi_common, f_interp_sham, '--', 'Color', shamcolor, 'LineWidth', 2, 'DisplayName', 'Sham');

title('Last Session Stim Reaches: First, Middle, Last 5 + Sham');
xlabel('Time from threshold crossing (s)');
ylabel('Density');
xlim([-0.1 0.2]);

legend('show');
set(gca, 'FontSize', 12, 'LineWidth', 1);
hold off;
%%
%% -----------------------------
%% First Session Stim Reaches: First, Middle, Last 5 + Sham
%% -----------------------------
figure; hold on;

segments = {'first5', 'middle5', 'last5'};
colors = lines(3); % Colors for first/middle/last
xi_common = linspace(-0.1, 0.2, 500); % common support for interpolation
PeakInfo_Stim1 = struct(); % store peaks for first session

% Ctrl1 (first session) KDEs
for s = 1:length(segments)
    all_segment = [];
    for a = 1:length(animals_ctrl1)
        dat = ExcitRR.(animals_ctrl1{a}).StimRR(:,1);  % Take first column (CSM times)
        % Extract first/middle/last 5 for this animal
        switch segments{s}
            case 'first5'
                dat_seg = dat(1:min(5,numel(dat)));
            case 'middle5'
                mid_idx = floor(numel(dat)/2)-2 : floor(numel(dat)/2)+2;
                mid_idx = mid_idx(mid_idx >= 1 & mid_idx <= numel(dat));
                dat_seg = dat(mid_idx);
            case 'last5'
                dat_seg = dat(max(1,end-4):end);
        end
        if ~isempty(dat_seg)
            all_segment = [all_segment; dat_seg];
        end
    end
    
    % --- Compute KDE ---
    [f_seg, xi_seg] = ksdensity(all_segment, 'Bandwidth', 0.008);
    f_interp_seg = interp1(xi_seg, f_seg, xi_common, 'linear', 0);
    
    % --- Find peak location and density ---
    [peak_val, idx_peak] = max(f_interp_seg);
    peak_loc = xi_common(idx_peak);
    
    % --- Save results ---
    PeakInfo_Stim1.(segments{s}).peak_loc = peak_loc;
    PeakInfo_Stim1.(segments{s}).peak_val = peak_val;
    PeakInfo_Stim1.(segments{s}).xi = xi_common;
    PeakInfo_Stim1.(segments{s}).f = f_interp_seg;
    
    % --- Plot KDE ---
    plot(xi_common, f_interp_seg, 'Color', colors(s,:), 'LineWidth', 2, 'DisplayName', segments{s});
    
    % --- Mark peak ---
    plot(peak_loc, peak_val, 'o', 'MarkerEdgeColor', colors(s,:), 'MarkerFaceColor', colors(s,:), 'HandleVisibility','off');
    text(peak_loc, peak_val + 0.05*peak_val, sprintf('%.3f', peak_loc), 'Color', colors(s,:), 'FontSize', 10, 'HorizontalAlignment','center');
end

% Sham KDE (all reaches combined)
all_data_sham = [];
for a = 1:length(animals_sham)
    dat = ExcitRR_Sham.(animals_sham{a}).AllBins_prethresh;
    all_data_sham = [all_data_sham; dat];
end
[f_sham, xi_sham] = ksdensity(all_data_sham, 'Bandwidth', 0.008);
f_interp_sham = interp1(xi_sham, f_sham, xi_common, 'linear', 0);
plot(xi_common, f_interp_sham, '--', 'Color', shamcolor, 'LineWidth', 2, 'DisplayName', 'Sham');

title('First Session Stim Reaches: First, Middle, Last 5 + Sham');
xlabel('Time from threshold crossing (s)');
ylabel('Density');
xlim([-0.1 0.2]);
legend('show');
set(gca, 'FontSize', 12, 'LineWidth', 1);
hold off;

disp(struct2table(PeakInfo_Stim1));


%% -----------------------------
%% Last Session Stim Reaches: First, Middle, Last 5 + Sham
%% -----------------------------
figure; hold on;

PeakInfo_StimEnd = struct(); % store peaks for last session

% Ctrlend (last session) KDEs
for s = 1:length(segments)
    all_segment = [];
    for a = 1:length(animals_ctrlend)
        dat = ExcitRR_end.(animals_ctrlend{a}).StimRR(:,1);  % Take first column (CSM times)
        % Extract first/middle/last 5
        switch segments{s}
            case 'first5'
                dat_seg = dat(1:min(5,numel(dat)));
            case 'middle5'
                mid_idx = floor(numel(dat)/2)-2 : floor(numel(dat)/2)+2;
                mid_idx = mid_idx(mid_idx >= 1 & mid_idx <= numel(dat));
                dat_seg = dat(mid_idx);
            case 'last5'
                dat_seg = dat(max(1,end-4):end);
        end
        if ~isempty(dat_seg)
            all_segment = [all_segment; dat_seg];
        end
    end
    
    % --- Compute KDE ---
    [f_seg, xi_seg] = ksdensity(all_segment, 'Bandwidth', 0.008);
    f_interp_seg = interp1(xi_seg, f_seg, xi_common, 'linear', 0);
    
    % --- Find peak location and density ---
    [peak_val, idx_peak] = max(f_interp_seg);
    peak_loc = xi_common(idx_peak);
    
    % --- Save results ---
    PeakInfo_StimEnd.(segments{s}).peak_loc = peak_loc;
    PeakInfo_StimEnd.(segments{s}).peak_val = peak_val;
    PeakInfo_StimEnd.(segments{s}).xi = xi_common;
    PeakInfo_StimEnd.(segments{s}).f = f_interp_seg;
    
    % --- Plot KDE ---
    plot(xi_common, f_interp_seg, 'Color', colors(s,:), 'LineWidth', 2, 'DisplayName', segments{s});
    
    % --- Mark peak ---
    plot(peak_loc, peak_val, 'o', 'MarkerEdgeColor', colors(s,:), 'MarkerFaceColor', colors(s,:), 'HandleVisibility','off');
    text(peak_loc, peak_val + 0.05*peak_val, sprintf('%.3f', peak_loc), 'Color', colors(s,:), 'FontSize', 10, 'HorizontalAlignment','center');
end

% Sham KDE (reuse previously computed f_interp_sham)
plot(xi_common, f_interp_sham, '--', 'Color', shamcolor, 'LineWidth', 2, 'DisplayName', 'Sham');

title('Last Session Stim Reaches: First, Middle, Last 5 + Sham');
xlabel('Time from threshold crossing (s)');
ylabel('Density');
xlim([-0.1 0.2]);
legend('show');
set(gca, 'FontSize', 12, 'LineWidth', 1);
hold off;

disp(struct2table(PeakInfo_StimEnd));

%% Now plotting those peak times and densities 


segments = {'first5','middle5','last5'};
session_labels = {'First Session','Last Session'};

% Initialize arrays
peak_times = [];
peak_densities = [];
group_labels = {};  % control vs stim
session_group = {}; % first vs last session
segment_group = {}; % first/middle/last 5

% --- Control first session ---
for s = 1:length(segments)
    peak_times(end+1) = PeakInfo_Ctrl1.(segments{s}).peak_loc;
    peak_densities(end+1) = PeakInfo_Ctrl1.(segments{s}).peak_val;
    group_labels{end+1} = 'Control';
    session_group{end+1} = 'First';
    segment_group{end+1} = segments{s};
end

% --- Control last session ---
for s = 1:length(segments)
    peak_times(end+1) = PeakInfo_CtrlEnd.(segments{s}).peak_loc;
    peak_densities(end+1) = PeakInfo_CtrlEnd.(segments{s}).peak_val;
    group_labels{end+1} = 'Control';
    session_group{end+1} = 'Last';
    segment_group{end+1} = segments{s};
end

% --- Stim first session ---
for s = 1:length(segments)
    peak_times(end+1) = PeakInfo_Stim1.(segments{s}).peak_loc;
    peak_densities(end+1) = PeakInfo_Stim1.(segments{s}).peak_val;
    group_labels{end+1} = 'Stim';
    session_group{end+1} = 'First';
    segment_group{end+1} = segments{s};
end

% --- Stim last session ---
for s = 1:length(segments)
    peak_times(end+1) = PeakInfo_StimEnd.(segments{s}).peak_loc;
    peak_densities(end+1) = PeakInfo_StimEnd.(segments{s}).peak_val;
    group_labels{end+1} = 'Stim';
    session_group{end+1} = 'Last';
    segment_group{end+1} = segments{s};
end

% Create table for plotting
T_peaks = table(peak_times', peak_densities', group_labels', session_group', segment_group', ...
    'VariableNames', {'PeakTime','PeakDensity','Group','Session','Segment'});


%% Plot peak times and densities

figure; hold on;

% Use different colors for Control vs Stim
group_colors = containers.Map({'Control','Stim'},{[0 0.4470 0.7410],[0.8500 0.3250 0.0980]});
marker_shapes = containers.Map({'First','Last'},{'o','s'}); % first vs last session

for i = 1:height(T_peaks)
    x = T_peaks.PeakTime(i);
    y = T_peaks.PeakDensity(i);
    c = group_colors(T_peaks.Group{i});
    m = marker_shapes(T_peaks.Session{i});
    plot(x, y, m, 'MarkerSize',8, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k');
end

xlabel('Peak Time (s)');
ylabel('Peak Density');
title('KDE Peaks: Time vs Density');
legend({'Control First','Control Last','Stim First','Stim Last'}, 'Location','northwest');
set(gca, 'FontSize',12, 'LineWidth',1);
hold off;
%%
%% -----------------------------
%% Scatter plot: Peak Time vs Peak Density
%% -----------------------------
figure; hold on;

% Colors for Control vs Stim
group_colors = containers.Map({'Control','Stim'}, {[0 0.4470 0.7410], [0.8500 0.3250 0.0980]});
% Marker shapes for First vs Last session
marker_shapes = containers.Map({'First','Last'}, {'o','s'});

% Keep track of plotted combinations for legend
plotted_combos = containers.Map;

for i = 1:height(T_peaks)
    x = T_peaks.PeakTime(i);
    y = T_peaks.PeakDensity(i);
    c = group_colors(T_peaks.Group{i});
    m = marker_shapes(T_peaks.Session{i});
    
    % For legend: only plot first occurrence
    combo = [T_peaks.Group{i} '_' T_peaks.Session{i}];
    if ~isKey(plotted_combos, combo)
        h = plot(x, y, m, 'MarkerSize', 10, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'DisplayName', combo);
        plotted_combos(combo) = h;
    else
        plot(x, y, m, 'MarkerSize', 10, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
    end
end

xlabel('Peak Time (s)');
ylabel('Peak Density');
title('KDE Peaks: Time vs Density for First, Middle, Last 5 Reaches');
legend('Location','northwest');
set(gca, 'FontSize',12, 'LineWidth',1);
grid on;
hold off;

%%
%% -----------------------------
%% Scatter plot: Control only with session and segment encoding
%% -----------------------------
T_control = T_peaks(strcmp(T_peaks.Group,'Control'), :);  % filter for Control

figure; hold on;

% Colors for session
session_colors = containers.Map({'First','Last'}, {[0 0.4470 0.7410],[0 0.745 0.227]}); % blue vs green

% Marker shapes for segment
segment_markers = containers.Map({'first5','middle5','last5'}, {'o','s','^'}); % circle, square, triangle

% Keep track of legend entries
legend_entries = containers.Map;

for i = 1:height(T_control)
    x = T_control.PeakTime(i);
    y = T_control.PeakDensity(i);
    c = session_colors(T_control.Session{i});
    m = segment_markers(T_control.Segment{i});
    
    % Create unique label for legend
    combo = [T_control.Session{i} '_' T_control.Segment{i}];
    if ~isKey(legend_entries, combo)
        h = plot(x, y, m, 'MarkerSize', 10, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'DisplayName', combo);
        legend_entries(combo) = h;
    else
        plot(x, y, m, 'MarkerSize', 10, 'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'HandleVisibility','off');
    end
end

xlabel('Peak Time (s)');
ylabel('Peak Density');
title('Control KDE Peaks: Time vs Density for First, Middle, Last 5 Reaches');
legend('Location','northwest');
set(gca, 'FontSize',12, 'LineWidth',1);
grid on;
hold off;


%% T2- first session

% Initialize arrays for p-values and Z-scores
parray2 = NaN(size(T2_big.stimvel, 1), 1);  % Store p-values
zarray2 =   NaN(size(T2_big.stimvel, 1), 1)   ;    % Store Z-scores (same shape as stimvel)

for pnum = 1:length(parray2)
    % Extract trial data for C_big (control) and T0_big (experimental)
    a4b = squeeze(CT2_big.stimvel(pnum, 5, :));  % Experimental data (T0_big)
    b4b = squeeze(ShamT2_big.stimvel(pnum, 5, :));   % Control data (C_big), matching trial

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
T2color=rgb('SlateBlue');
%zscore plot 
% Adjust the figure to provide more vertical space
set(gcf, 'Position', [100, 100, 800, 600]); % [left, bottom, width, height]

dim = 5;
StimvelPlot(CT2_big, dim,T2color, 1);
StimvelPlot(ShamT2_big, dim,'k', 0);

% Plot reference line
yline(0, '--k');

% Axis limits and labels
xlim([-0.2 0.2]);
%ylim([-20 20]);
ylim('auto')
%ylim([-0.05 0.27]);
%ylim([-0.5 1.5])
ylabel('Outward Velocity (cm/s)')
%xlabel('Time from Pos. Thresh. (s)');
%ylabel('Outward Acceleration(cm/s)');
%ylabel('Outward Position(cm)');
xlabel('Time from threshold crossing (s)')
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

%% T2 last session

% Initialize arrays for p-values and Z-scores
parray2 = NaN(size(T2_big.stimvel, 1), 1);  % Store p-values
zarray2 =   NaN(size(T2_big.stimvel, 1), 1)   ;    % Store Z-scores (same shape as stimvel)

for pnum = 1:length(parray2)
    % Extract trial data for C_big (control) and T0_big (experimental)
    a4b = squeeze(CT2end_big.stimvel(pnum, 5, :));  % Experimental data (T0_big)
    b4b = squeeze(ShamT2_big.stimvel(pnum, 5, :));   % Control data (C_big), matching trial

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
T2color=rgb('SlateBlue');
%zscore plot 
% Adjust the figure to provide more vertical space
set(gcf, 'Position', [100, 100, 800, 600]); % [left, bottom, width, height]

dim = 5;
StimvelPlot(CT2end_big, dim,T2color, 1);
StimvelPlot(ShamT2_big, dim,'k', 0);

% Plot reference line
yline(0, '--k');

% Axis limits and labels
xlim([-0.2 0.2]);
%ylim([-20 20]);
ylim('auto')
%ylim([-0.05 0.27]);
%ylim([-0.5 1.5])
ylabel('Outward Velocity (cm/s)')
%xlabel('Time from Pos. Thresh. (s)');
%ylabel('Outward Acceleration(cm/s)');
%ylabel('Outward Position(cm)');
xlabel('Time from threshold crossing (s)')
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

