%% Find min and max of accel during CSM time window ( 100ms after stim ) 

%% Replace all reaches aligned CSM with only the NoCSM template

an=fieldnames(ExcitRR_Sham);

thresh=0.16;
dirdir= 1;

ExcitRR_Ctrl= ExcitRR_Sham;

for a=1:length(an)
   
    CtrlReaches=ExcitRR_Ctrl.(an{a}).NoReReach;

    % NoCSM algined to T0
    CST0_CSM=AlignThreshCSM_T0(CtrlReaches,thresh,dirdir,200);
    CST0_CSM=stimvel_accel(CST0_CSM);

    % Ctrl to T1
    CST1_CSM=AlignThreshCSM_T1(CtrlReaches,thresh,dirdir,200);
    CST1_CSM=stimvel_accel(CST1_CSM);

    % Ctrl to T2
    CST2_CSM=AlignThreshCSM_T2(CtrlReaches,thresh,dirdir,200);
    CST2_CSM=stimvel_accel(CST2_CSM);

    % Ctrl to T3
    CST3_CSM=AlignThreshCSM_T3(CtrlReaches,thresh,dirdir,200);
    CST3_CSM=stimvel_accel(CST3_CSM);
    
    % Save to big structure 

     %Stimvel

     ExcitRR_Ctrl.(an{a}).CST0_CSM=CST0_CSM;
     ExcitRR_Ctrl.(an{a}).CST1_CSM=CST1_CSM;
     ExcitRR_Ctrl.(an{a}).CST2_CSM=CST2_CSM;
     ExcitRR_Ctrl.(an{a}).CST3_CSM=CST3_CSM;
        

end
%%
%T0



an=fieldnames(ExcitRR_Ctrl);

for a=1:length(an)
    stimvel=ExcitRR_Ctrl.(an{a}).ST0_CSM.stimvel;
    Accelmin=[];
    Accelmax=[];
   stimvelsize=size(stimvel);
   num=(stimvelsize(3));

    for b=1:num
        reach=stimvel(:,:,b);
        decel_window=reach(11:21,:); % 100 to 0
        accel_window=reach(21:31,:); % 0 to 100

        % decel_window=reach(11:21,:); % 100 to 0
        % accel_window=reach(21:36,:); % 0 to 100
       
        accelmin=min(decel_window(:,8));
        accelmax=max(accel_window(:,8));
        Accelmin=[Accelmin;accelmin];
        Accelmax=[Accelmax;accelmax];
    end

    ExcitRR_Ctrl.(an{a}).AccelMinT0=Accelmin;
    ExcitRR_Ctrl.(an{a}).AccelMaxT0=Accelmax;
    clear Accelmin Accelmax 
end

%%
%T1
for a=1:length(an)
    stimvel=ExcitRR_Ctrl.(an{a}).ST1_CSM.stimvel;
    Accelmin=[];
    Accelmax=[];
   stimvelsize=size(stimvel);
   num=(stimvelsize(3));

    for b=1:num
        reach=stimvel(:,:,b);
        decel_window=reach(11:21,:);
        accel_window=reach(21:31,:);


        % decel_window=reach(11:21,:); % 100 to 0
        % accel_window=reach(21:36,:); % 0 to 100
       
        accelmin=min(decel_window(:,8));
        accelmax=max(accel_window(:,8));


        Accelmin=[Accelmin;accelmin];
        Accelmax=[Accelmax;accelmax];
    end

    ExcitRR_Ctrl.(an{a}).AccelMinT1=Accelmin;
    ExcitRR_Ctrl.(an{a}).AccelMaxT1=Accelmax;
    clear Accelmin Accelmax 
end

%% 
%T2
for a=1:length(an)
    stimvel=ExcitRR_Ctrl.(an{a}).ST2_CSM.stimvel;
    Accelmin=[];
    Accelmax=[];
   stimvelsize=size(stimvel);
   num=(stimvelsize(3));

    for b=1:num
        reach=stimvel(:,:,b);
        decel_window=reach(11:21,:);
        accel_window=reach(21:31,:);


        % decel_window=reach(11:21,:); % 100 to 0
        % accel_window=reach(21:36,:); % 0 to 100
       
        accelmin=min(decel_window(:,8));
        accelmax=max(accel_window(:,8));
        Accelmin=[Accelmin;accelmin];
        Accelmax=[Accelmax;accelmax];
    end

    ExcitRR_Ctrl.(an{a}).AccelMinT2=Accelmin;
    ExcitRR_Ctrl.(an{a}).AccelMaxT2=Accelmax;
    clear Accelmin Accelmax 
end

%% 
%T3
for a=1:length(an)
    stimvel=ExcitRR_Ctrl.(an{a}).ST3_CSM.stimvel;
    Accelmin=[];
    Accelmax=[];
   stimvelsize=size(stimvel);
   num=(stimvelsize(3));

    for b=1:num
        reach=stimvel(:,:,b);
        decel_window=reach(11:21,:);
        accel_window=reach(21:31,:);


        % decel_window=reach(11:21,:); % 100 to 0
        % accel_window=reach(21:36,:); % 0 to 100
       
        accelmin=min(decel_window(:,8));
        accelmax=max(accel_window(:,8));
        Accelmin=[Accelmin;accelmin];
        Accelmax=[Accelmax;accelmax];
    end

    ExcitRR_Ctrl.(an{a}).AccelMinT3=Accelmin;
    ExcitRR_Ctrl.(an{a}).AccelMaxT3=Accelmax;
    clear Accelmin Accelmax 
end




%% Different analysis approach 

% Estimate CSM gain across subtypes
% Gain = slope of AccelMax predicted by |AccelMin|

an = fieldnames(ExcitRR_Ctrl);

subtypes = {'T0','T1','T2','T3'};

T0color = rgb('Thistle');
T1color = rgb('MediumOrchid');
T2color = rgb('SlateBlue');
T3color = rgb('DarkSlateBlue');

subtype_colors = {T0color, T1color, T2color, T3color};

minReaches = 10;   % exclude animal/subtype combos with too few reaches

GainResults = struct;

for s = 1:length(subtypes)

    subtype = subtypes{s};

    gains = NaN(length(an),1);
    intercepts = NaN(length(an),1);
    p_gain = NaN(length(an),1);
    r2 = NaN(length(an),1);
    n_reaches = NaN(length(an),1);

    for a = 1:length(an)

        min_field = ['AccelMin' subtype];
        max_field = ['AccelMax' subtype];

        Accelmin = abs(ExcitRR_Ctrl.(an{a}).(min_field) * 10);
        Accelmax = ExcitRR_Ctrl.(an{a}).(max_field) * 10;

        valid_idx = ~isnan(Accelmin) & ~isnan(Accelmax);
        Accelmin = Accelmin(valid_idx);
        Accelmax = Accelmax(valid_idx);

        n_reaches(a) = length(Accelmin);

        if n_reaches(a) < minReaches
            continue
        end

        mdl = fitlm(Accelmin, Accelmax);

        intercepts(a) = mdl.Coefficients.Estimate(1);
        gains(a) = mdl.Coefficients.Estimate(2);
        p_gain(a) = mdl.Coefficients.pValue(2);
        r2(a) = mdl.Rsquared.Ordinary;

    end

    GainResults.(subtype).animal = an;
    GainResults.(subtype).gain = gains;
    GainResults.(subtype).intercept = intercepts;
    GainResults.(subtype).p_gain = p_gain;
    GainResults.(subtype).r2 = r2;
    GainResults.(subtype).n_reaches = n_reaches;

    fprintf('\n%s\n', subtype);
    fprintf('Animals included: %d / %d\n', sum(~isnan(gains)), length(an));
    fprintf('Mean gain: %.3f ± %.3f SEM\n', ...
        mean(gains, 'omitnan'), std(gains, 'omitnan') ./ sqrt(sum(~isnan(gains))));
end



%% Plot regression fits for each subtype

figure; hold on

xs = linspace(-100, 500, 100);

for s = 1:length(subtypes)

    subtype = subtypes{s};
    col = subtype_colors{s};

    gains = GainResults.(subtype).gain;
    intercepts = GainResults.(subtype).intercept;

    valid = ~isnan(gains) & ~isnan(intercepts);

    fit_values = NaN(sum(valid), length(xs));

    idx = find(valid);

    for i = 1:length(idx)

        a = idx(i);

        yfit = intercepts(a) + gains(a) .* xs;
        fit_values(i,:) = yfit;

        plot(xs, yfit, 'Color', [col 0.25], 'LineWidth', 2)

    end

    mean_fit = mean(fit_values, 1, 'omitnan');
    sem_fit = std(fit_values, 0, 1, 'omitnan') ./ sqrt(size(fit_values,1));

    fill([xs fliplr(xs)], ...
         [mean_fit + sem_fit, fliplr(mean_fit - sem_fit)], ...
         col, 'FaceAlpha', 0.15, 'EdgeColor', 'none')

    plot(xs, mean_fit, 'Color', col, 'LineWidth', 4)

end

plot([0 500], [0 500], '--k', 'LineWidth', 1.5)

xlim([0 500])
ylim([0 500])

xlabel('Peak Decel magnitude')
ylabel('CSM Peak Accel')
title('Regression-derived CSM gain by subtype')
legend(subtypes, 'Location', 'best')
set(gcf, 'Position', [100 100 500 450])



%% Plot gain across CSM subtypes with r and p values

figure; hold on

all_gain_mat = NaN(length(an), length(subtypes));

for s = 1:length(subtypes)
    subtype = subtypes{s};
    all_gain_mat(:,s) = GainResults.(subtype).gain;
end

% connect same animal across subtypes
for a = 1:length(an)
    plot(1:length(subtypes), all_gain_mat(a,:), '-', ...
        'Color', [0.75 0.75 0.75], 'LineWidth', 1)
end

for s = 1:length(subtypes)

    subtype = subtypes{s};
    col = subtype_colors{s};
    gains = all_gain_mat(:,s);

    scatter(ones(length(gains),1)*s, gains, 35, col, 'filled', ...
        'jitter', 'on', 'jitterAmount', 0.12, ...
        'MarkerFaceAlpha', 0.65)

    mean_gain = mean(gains, 'omitnan');
    sem_gain = std(gains, 'omitnan') ./ sqrt(sum(~isnan(gains)));

    errorbar(s, mean_gain, sem_gain, 'o', ...
        'Color', col, ...
        'MarkerFaceColor', col, ...
        'MarkerEdgeColor', col, ...
        'LineWidth', 2, ...
        'MarkerSize', 10, ...
        'CapSize', 10)

    % compute animal-level r values for this subtype
    r_vals = NaN(length(an),1);

    for a = 1:length(an)

        min_field = ['AccelMin' subtype];
        max_field = ['AccelMax' subtype];

        Accelmin = abs(ExcitRR_Ctrl.(an{a}).(min_field) * 10);
        Accelmax = ExcitRR_Ctrl.(an{a}).(max_field) * 10;

        valid_idx = ~isnan(Accelmin) & ~isnan(Accelmax);
        Accelmin = Accelmin(valid_idx);
        Accelmax = Accelmax(valid_idx);

        if length(Accelmin) < minReaches
            continue
        end

        r_vals(a) = corr(Accelmin, Accelmax);

    end

    mean_r = mean(r_vals, 'omitnan');

    % Fisher z-transform and t-test against zero
    z_vals = 0.5 * log((1 + r_vals) ./ (1 - r_vals));
    z_vals = z_vals(~isnan(z_vals) & ~isinf(z_vals));

    if length(z_vals) >= 2
        [~, p_r] = ttest(z_vals);
    else
        p_r = NaN;
    end

    % print r and p above each subtype
    text(s, mean_gain + sem_gain + 0.05, ...
        sprintf('r = %.2f\np = %.3f', mean_r, p_r), ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'Color', col)

end

xlim([0.5 length(subtypes)+0.5])
xticks(1:length(subtypes))
xticklabels(subtypes)

ylabel('Gain: slope of CSM accel vs preceding decel')
xlabel('CSM subtype')
title('CSM gain across subtypes')

set(gcf, 'Position', [100 100 450 400])


%% Gain across ALL CSMs (equal weighting of animals)

figure; hold on

RRcolor = rgb('Black');

xs = linspace(-100,500,100);

fit_values = [];
slopes = [];
intercepts = [];
r_vals = [];

for a = 1:length(an)

    Accelmin = [];
    Accelmax = [];

    % concatenate T0-T3
    for t = 0:3

        min_field = sprintf('AccelMinT%d',t);
        max_field = sprintf('AccelMaxT%d',t);

        if isfield(ExcitRR_Ctrl.(an{a}),min_field)

            Accelmin = [Accelmin;
                abs(ExcitRR_Ctrl.(an{a}).(min_field))*10];

            Accelmax = [Accelmax;
                ExcitRR_Ctrl.(an{a}).(max_field)*10];

        end

    end

    valid = ~isnan(Accelmin) & ~isnan(Accelmax);

    Accelmin = Accelmin(valid);
    Accelmax = Accelmax(valid);

    if length(Accelmin) < minReaches
        continue
    end

    mdl = fitlm(Accelmin,Accelmax);

    slope = mdl.Coefficients.Estimate(2);
    intercept = mdl.Coefficients.Estimate(1);

    slopes = [slopes; slope];
    intercepts = [intercepts; intercept];

    [r,p] = corr(Accelmin,Accelmax);

    r_vals = [r_vals; r];

    yfit = intercept + slope*xs;

    fit_values = [fit_values; yfit];

    % individual animal
    plot(xs,yfit,'Color',[0.7 0.7 0.7],'LineWidth',1.5)

end

% Mean regression across animals

mean_fit = mean(fit_values,1);
sem_fit = std(fit_values,0,1)./sqrt(size(fit_values,1));

fill([xs fliplr(xs)],...
    [mean_fit+sem_fit fliplr(mean_fit-sem_fit)],...
    RRcolor,...
    'FaceAlpha',0.2,...
    'EdgeColor','none')

plot(xs,mean_fit,...
    'Color',RRcolor,...
    'LineWidth',5)

plot([0 500],[0 500],'--k')

% Statistics

mean_slope = mean(slopes);
sem_slope = std(slopes)/sqrt(length(slopes));

mean_r = mean(r_vals);

z = atanh(r_vals);
[~,p_group] = ttest(z);

text(250,425,...
    sprintf(['Slope = %.2f \\pm %.2f\n' ...
             'Mean r = %.2f\n' ...
             'p = %.3f'],...
             mean_slope,sem_slope,...
             mean_r,p_group),...
    'FontWeight','bold',...
    'FontSize',11,...
    'BackgroundColor','w')

xlim([0 500])
ylim([0 500])

xlabel('Peak Deceleration (|cm/s^2|)')
ylabel('Peak Corrective Acceleration (cm/s^2)')
title('Gain across all CSMs')

set(gca,'Box','off')
set(gca,'TickDir','in')
set(gcf, 'Position', [100 100 500 450])

%% Subtracting from template reaches 

%% Residual acceleration analysis: CSM reach - subtype-matched NoCSM template

% Assumes you already created:
% ExcitRR_Ctrl.(animal).ST0_CSM.stimvel
% ExcitRR_Ctrl.(animal).CST0_CSM.stimvel
% ...
% through T3

an = fieldnames(ExcitRR_Ctrl);

subtypes = {'T0','T1','T2','T3'};

CSM_fields  = {'ST0_CSM','ST1_CSM','ST2_CSM','ST3_CSM'};
CTRL_fields = {'CST0_CSM','CST1_CSM','CST2_CSM','CST3_CSM'};

for s = 1:length(subtypes)

    subtype = subtypes{s};
    CSM_field = CSM_fields{s};
    CTRL_field = CTRL_fields{s};

    for a = 1:length(an)

        % CSM reaches
        csm_stimvel = ExcitRR_Ctrl.(an{a}).(CSM_field).stimvel;

        % Subtype-matched NoCSM reaches
        ctrl_stimvel = ExcitRR_Ctrl.(an{a}).(CTRL_field).stimvel;

        % Average NoCSM acceleration template
        % channel 8 = acceleration
        ctrl_template_accel = mean(ctrl_stimvel(:,8,:),3,'omitnan');

        ResidAccelMin = [];
        ResidAccelMax = [];

        stimvelsize = size(csm_stimvel);
        num = stimvelsize(3);

        for b = 1:num

            reach = csm_stimvel(:,:,b);

            % Raw CSM acceleration
            csm_accel = reach(:,8);

            % Residual acceleration = CSM acceleration - NoCSM template
            resid_accel = csm_accel - ctrl_template_accel;

            % Optional smoothing
            resid_accel = smoothdata(resid_accel,'sgolay',5);

            % Same windows as before
            decel_window = resid_accel(11:21); % -100 to 0 ms
            accel_window = resid_accel(21:31); % 0 to 100 ms

            % Residual peaks
            resid_decel_min = min(decel_window);
            resid_accel_max = max(accel_window);

            ResidAccelMin = [ResidAccelMin; resid_decel_min];
            ResidAccelMax = [ResidAccelMax; resid_accel_max];

        end

        % Save residual metrics
        ExcitRR_Ctrl.(an{a}).(['ResidAccelMin' subtype]) = ResidAccelMin;
        ExcitRR_Ctrl.(an{a}).(['ResidAccelMax' subtype]) = ResidAccelMax;

    end
end


%% Plot residual decel/accel correlations by subtype and animal

an = fieldnames(ExcitRR_Ctrl);

subtypes = {'T0','T1','T2','T3'};

T0color = rgb('Thistle');
T1color = rgb('MediumOrchid');
T2color = rgb('SlateBlue');
T3color = rgb('DarkSlateBlue');

subtype_colors = {T0color,T1color,T2color,T3color};

minReaches = 10;

GainResults_Resid = struct;

figure; hold on

for s = 1:length(subtypes)

    subtype = subtypes{s};
    col = subtype_colors{s};

    gains = NaN(length(an),1);
    intercepts = NaN(length(an),1);
    r_vals = NaN(length(an),1);
    p_vals = NaN(length(an),1);
    n_reaches = NaN(length(an),1);

    for a = 1:length(an)

        xfield = ['ResidAccelMin' subtype];
        yfield = ['ResidAccelMax' subtype];

        % Use absolute residual decel magnitude
        x = abs(ExcitRR_Ctrl.(an{a}).(xfield) * 10);
        y = ExcitRR_Ctrl.(an{a}).(yfield) * 10;

        valid = ~isnan(x) & ~isnan(y);
        x = x(valid);
        y = y(valid);

        n_reaches(a) = length(x);

        if length(x) < minReaches
            continue
        end

        % Correlation
        [r_vals(a),p_vals(a)] = corr(x,y,'rows','complete');

        % Linear fit
        pf = polyfit(x,y,1);

        gains(a) = pf(1);
        intercepts(a) = pf(2);

        xs = linspace(0,500,100);
        yfit = polyval(pf,xs);

        plot(xs,yfit,...
            'Color',[col 0.25],...
            'LineWidth',1.5)

    end

    % Save
    GainResults_Resid.(subtype).gain = gains;
    GainResults_Resid.(subtype).intercept = intercepts;
    GainResults_Resid.(subtype).r = r_vals;
    GainResults_Resid.(subtype).p = p_vals;
    GainResults_Resid.(subtype).n_reaches = n_reaches;

    % Mean fit across animals
    valid_animals = ~isnan(gains) & ~isnan(intercepts);

    xs = linspace(0,500,100);
    fit_values = NaN(sum(valid_animals),length(xs));

    idx = find(valid_animals);

    for i = 1:length(idx)
        aa = idx(i);
        fit_values(i,:) = intercepts(aa) + gains(aa).*xs;
    end

    mean_fit = mean(fit_values,1,'omitnan');
    sem_fit = std(fit_values,0,1,'omitnan') ./ sqrt(size(fit_values,1));

    fill([xs fliplr(xs)],...
        [mean_fit + sem_fit fliplr(mean_fit - sem_fit)],...
        col,...
        'FaceAlpha',0.15,...
        'EdgeColor','none')

    plot(xs,mean_fit,...
        'Color',col,...
        'LineWidth',4)

    fprintf('\n%s residual accel correlation\n',subtype);
    fprintf('Animals included: %d / %d\n',sum(valid_animals),length(an));
    fprintf('Mean slope: %.3f ± %.3f SEM\n', ...
        mean(gains,'omitnan'), ...
        std(gains,'omitnan') ./ sqrt(sum(valid_animals)));
    fprintf('Mean r: %.3f\n',mean(r_vals,'omitnan'));

end

plot([0 500],[0 500],'--k','LineWidth',1.5)

xlabel('Residual peak decel magnitude')
ylabel('Residual CSM peak accel')
title('Residual acceleration coupling by CSM subtype')

xlim([0 500])
ylim([0 500])

legend(subtypes,'Location','best')
box off
set(gca,'FontSize',12)
set(gcf,'Position',[100 100 550 450])

%% Residual acceleration analysis: CSM reach - subtype-matched NoCSM template
% Peak methods:
% 'raw'        = min/max after smoothing
% 'top3'       = mean of 3 strongest decel/accel samples
% 'percentile' = 5th/95th percentile
% 'peakarea'   = mean around local peak +/- 2 samples

peak_method = 'top3';   % recommended starting point

an = fieldnames(ExcitRR_Ctrl);

subtypes = {'T0','T1','T2','T3'};

CSM_fields  = {'ST0_CSM','ST1_CSM','ST2_CSM','ST3_CSM'};
CTRL_fields = {'CST0_CSM','CST1_CSM','CST2_CSM','CST3_CSM'};

for s = 1:length(subtypes)

    subtype = subtypes{s};
    CSM_field = CSM_fields{s};
    CTRL_field = CTRL_fields{s};

    for a = 1:length(an)

        csm_stimvel = ExcitRR_Ctrl.(an{a}).(CSM_field).stimvel;
        ctrl_stimvel = ExcitRR_Ctrl.(an{a}).(CTRL_field).Average;

        % subtype-matched NoCSM acceleration template
        ctrl_template_accel =ctrl_stimvel(:,8),3,'omitnan';

        ResidAccelMin = [];
        ResidAccelMax = [];

        num = size(csm_stimvel,3);

        for b = 1:num

            reach = csm_stimvel(:,:,b);

            % residual acceleration = CSM accel - subtype-matched NoCSM template
            resid_accel = reach(:,8) - ctrl_template_accel;

            % smooth residual acceleration
            resid_accel = smoothdata(resid_accel,'sgolay',5);

            % windows
            decel_window = resid_accel(11:21);  % before CSM
            accel_window = resid_accel(21:31);  % after CSM

            switch peak_method

                case 'raw'
                    resid_decel_min = min(decel_window);
                    resid_accel_max = max(accel_window);

                case 'top3'
                    resid_decel_min = mean(mink(decel_window,3));
                    resid_accel_max = mean(maxk(accel_window,3));

                case 'percentile'
                    resid_decel_min = prctile(decel_window,5);
                    resid_accel_max = prctile(accel_window,95);

                case 'peakarea'
                    [~,iMin] = min(decel_window);
                    idxMin = max(1,iMin-2):min(length(decel_window),iMin+2);
                    resid_decel_min = mean(decel_window(idxMin));

                    [~,iMax] = max(accel_window);
                    idxMax = max(1,iMax-2):min(length(accel_window),iMax+2);
                    resid_accel_max = mean(accel_window(idxMax));

            end

            ResidAccelMin = [ResidAccelMin; resid_decel_min];
            ResidAccelMax = [ResidAccelMax; resid_accel_max];

        end

        ExcitRR_Ctrl.(an{a}).(['ResidAccelMin' subtype]) = ResidAccelMin;
        ExcitRR_Ctrl.(an{a}).(['ResidAccelMax' subtype]) = ResidAccelMax;

    end
end

%% Plot residual decel/accel correlations by subtype and animal
% Uses:
% ResidAccelMinT0-T3
% ResidAccelMaxT0-T3

an = fieldnames(ExcitRR_Ctrl);

subtypes = {'T0','T1','T2','T3'};

T0color = rgb('Thistle');
T1color = rgb('MediumOrchid');
T2color = rgb('SlateBlue');
T3color = rgb('DarkSlateBlue');

subtype_colors = {T0color,T1color,T2color,T3color};

minReaches = 10;

GainResults_Resid = struct;

figure; hold on

for s = 1:length(subtypes)

    subtype = subtypes{s};
    col = subtype_colors{s};

    gains = NaN(length(an),1);
    intercepts = NaN(length(an),1);
    r_vals = NaN(length(an),1);
    p_vals = NaN(length(an),1);
    n_reaches = NaN(length(an),1);

    for a = 1:length(an)
        % 
        % xfield = ['ResidAccelMin' subtype];
        % yfield = ['ResidAccelMax' subtype];
         xfield = ['AccelMin' subtype];
        yfield = ['AccelMax' subtype];

        x = abs(ExcitRR_Ctrl.(an{a}).(xfield) * 10);
        y = ExcitRR_Ctrl.(an{a}).(yfield) * 10;

        valid = ~isnan(x) & ~isnan(y) & ...
                ~isinf(x) & ~isinf(y);

        x = x(valid);
        y = y(valid);

        % Optional: remove extreme outliers
        valid_range = x >= 0 & x <= 500 & y >= -100 & y <= 500;
        x = x(valid_range);
        y = y(valid_range);

        n_reaches(a) = length(x);

        if n_reaches(a) < minReaches
            continue
        end

        [r_vals(a),p_vals(a)] = corr(x,y,'rows','complete');

        pf = polyfit(x,y,1);

        gains(a) = pf(1);
        intercepts(a) = pf(2);

        xs = linspace(0,500,100);
        yfit = polyval(pf,xs);

        % Individual animal fit
        plot(xs,yfit,...
            'Color',[col 0.25],...
            'LineWidth',1.5,...
            'HandleVisibility','off')

    end

    GainResults_Resid.(subtype).animal = an;
    GainResults_Resid.(subtype).gain = gains;
    GainResults_Resid.(subtype).intercept = intercepts;
    GainResults_Resid.(subtype).r = r_vals;
    GainResults_Resid.(subtype).p = p_vals;
    GainResults_Resid.(subtype).n_reaches = n_reaches;

    % Mean animal fit for subtype

    valid_animals = ~isnan(gains) & ~isnan(intercepts);

    xs = linspace(0,500,100);
    fit_values = NaN(sum(valid_animals),length(xs));

    idx = find(valid_animals);

    for i = 1:length(idx)
        aa = idx(i);
        fit_values(i,:) = intercepts(aa) + gains(aa).*xs;
    end

    mean_fit = mean(fit_values,1,'omitnan');
    sem_fit = std(fit_values,0,1,'omitnan') ./ sqrt(size(fit_values,1));

    fill([xs fliplr(xs)],...
        [mean_fit + sem_fit fliplr(mean_fit - sem_fit)],...
        col,...
        'FaceAlpha',0.15,...
        'EdgeColor','none',...
        'HandleVisibility','off')

    plot(xs,mean_fit,...
        'Color',col,...
        'LineWidth',4,...
        'DisplayName',subtype)

    % Print subtype summary

    z_vals = atanh(r_vals);
    z_vals = z_vals(~isnan(z_vals) & ~isinf(z_vals));

    if length(z_vals) >= 2
        [~,p_group] = ttest(z_vals);
        mean_r = tanh(mean(z_vals));
    else
        p_group = NaN;
        mean_r = NaN;
    end

    fprintf('\n%s residual acceleration coupling\n',subtype);
    fprintf('Animals included: %d / %d\n',sum(valid_animals),length(an));
    fprintf('Mean slope: %.3f ± %.3f SEM\n', ...
        mean(gains,'omitnan'), ...
        std(gains,'omitnan') ./ sqrt(sum(valid_animals)));
    fprintf('Mean r Fisher-transformed: %.3f\n',mean_r);
    fprintf('Group-level p from Fisher z t-test: %.4f\n',p_group);

end

%% Overall pooled regression across all residual CSMs

Xall = [];
Yall = [];

for s = 1:length(subtypes)

    subtype = subtypes{s};

    for a = 1:length(an)

        xfield = ['ResidAccelMin' subtype];
        yfield = ['ResidAccelMax' subtype];

        x = abs(ExcitRR_Ctrl.(an{a}).(xfield) * 10);
        y = ExcitRR_Ctrl.(an{a}).(yfield) * 10;

        valid = ~isnan(x) & ~isnan(y) & ...
                ~isinf(x) & ~isinf(y);

        x = x(valid);
        y = y(valid);

        valid_range = x >= 0 & x <= 500 & y >= -100 & y <= 500;
        x = x(valid_range);
        y = y(valid_range);

        Xall = [Xall; x];
        Yall = [Yall; y];

    end
end

[r_all,p_all] = corr(Xall,Yall,'rows','complete');

mdl_all = fitlm(Xall,Yall);

xs_all = linspace(0,500,100)';
[y_all,CI_all] = predict(mdl_all,xs_all);

fill([xs_all; flipud(xs_all)],...
     [CI_all(:,1); flipud(CI_all(:,2))],...
     [0.5 0.5 0.5],...
     'FaceAlpha',0.18,...
     'EdgeColor','none',...
     'HandleVisibility','off')

plot(xs_all,y_all,...
    'k',...
    'LineWidth',5,...
    'DisplayName','All CSMs')

R2_all = mdl_all.Rsquared.Ordinary;

text(260,430,...
    sprintf('All CSMs\nr = %.2f\np = %.4f\nR^2 = %.2f', ...
    r_all,p_all,R2_all),...
    'FontWeight','bold',...
    'FontSize',11,...
    'BackgroundColor','w')

plot([0 500],[0 500],'--k','LineWidth',1.2,'HandleVisibility','off')

xlabel('Residual peak decel magnitude')
ylabel('Residual peak accel')
title('Template-subtracted residual acceleration coupling')

xlim([0 500])
ylim([-100 500])

legend('Location','best')

box off
set(gca,'FontSize',12)
set(gca,'TickDir','in')
set(gcf,'Position',[100 100 550 450])


%% Trying something new

%% Template-subtracted peak decel / peak accel analysis
% For each animal and CSM subtype:
% 1. Make subtype-matched NoCSM acceleration template
% 2. Subtract template from each CSM reach
% 3. Find residual decel and residual accel peaks
% 4. Plot correlations

ExcitRR_Ctrl = ExcitRR_Sham;
an = fieldnames(ExcitRR_Ctrl);

subtypes = {'T0','T1','T2','T3'};
CSM_fields  = {'ST0_CSM','ST1_CSM','ST2_CSM','ST3_CSM'};
CTRL_fields = {'CST0_CSM','CST1_CSM','CST2_CSM','CST3_CSM'};

T0color = rgb('Thistle');
T1color = rgb('MediumOrchid');
T2color = rgb('SlateBlue');
T3color = rgb('DarkSlateBlue');
subtype_colors = {T0color,T1color,T2color,T3color};

minReaches = 10;

% Extract residual peaks

for s = 1:length(subtypes)

    subtype = subtypes{s};

    for a = 1:length(an)

        csm_stimvel  = ExcitRR_Ctrl.(an{a}).(CSM_fields{s}).stimvel;
        ctrl_stimvel = ExcitRR_Ctrl.(an{a}).(CTRL_fields{s}).stimvel;

        % Average NoCSM template for this animal/subtype
        template_accel = mean(ctrl_stimvel(:,8,:),3,'omitnan');

        ResidDecelPeak = [];
        ResidAccelPeak = [];

        for b = 1:size(csm_stimvel,3)

            reach_accel = csm_stimvel(:,8,b);

            % Subtract template
            resid_accel = reach_accel - template_accel;

            % Smooth after subtraction
            resid_accel = smoothdata(resid_accel,'sgolay',5);

            % Windows
            decel_window = resid_accel(11:21);
            accel_window = resid_accel(21:31);

            % Find peaks from residual trace
            decel_peak = min(decel_window);  % most negative residual decel
            accel_peak = max(accel_window);  % most positive residual accel

            ResidDecelPeak = [ResidDecelPeak; decel_peak];
            ResidAccelPeak = [ResidAccelPeak; accel_peak];

        end

        ExcitRR_Ctrl.(an{a}).(['ResidDecelPeak' subtype]) = ResidDecelPeak;
        ExcitRR_Ctrl.(an{a}).(['ResidAccelPeak' subtype]) = ResidAccelPeak;

    end
end

% Plot residual peak correlations

figure; hold on

all_x = [];
all_y = [];

for s = 1:length(subtypes)

    subtype = subtypes{s};
    col = subtype_colors{s};

    for a = 1:length(an)

        x = abs(ExcitRR_Ctrl.(an{a}).(['ResidDecelPeak' subtype])) * 10;
        y = ExcitRR_Ctrl.(an{a}).(['ResidAccelPeak' subtype]) * 10;

        valid = ~isnan(x) & ~isnan(y) & ~isinf(x) & ~isinf(y);
        x = x(valid);
        y = y(valid);

        if length(x) < minReaches
            continue
        end

        all_x = [all_x; x];
        all_y = [all_y; y];

        % Per-animal regression line
        pf = polyfit(x,y,1);
        xs = linspace(min(x),max(x),50);
        ys = polyval(pf,xs);

        plot(xs,ys,...
            'Color',[col 0.35],...
            'LineWidth',1.5,...
            'HandleVisibility','off')

    end
end

%% Overall regression across all reaches

[r_all,p_all] = corr(all_x,all_y,'rows','complete');

mdl_all = fitlm(all_x,all_y);

xs = linspace(min(all_x),max(all_x),200)';
[yhat,CI] = predict(mdl_all,xs);

fill([xs; flipud(xs)],...
     [CI(:,1); flipud(CI(:,2))],...
     [0.5 0.5 0.5],...
     'FaceAlpha',0.2,...
     'EdgeColor','none',...
     'HandleVisibility','off')

plot(xs,yhat,'k','LineWidth',4,'DisplayName','All CSMs')

R2 = mdl_all.Rsquared.Ordinary;

text(max(all_x)*0.55,max(all_y)*0.9,...
    sprintf('All reaches\nr = %.2f\np = %.4f\nR^2 = %.2f',r_all,p_all,R2),...
    'FontSize',11,...
    'FontWeight','bold',...
    'BackgroundColor','w')

xlabel('Template-subtracted decel peak magnitude')
ylabel('Template-subtracted accel peak')
title('Residual acceleration peaks after subtracting NoCSM template')

legend('Location','best')
box off
set(gca,'FontSize',12)
set(gca,'TickDir','in')
set(gcf,'Position',[100 100 550 450])

%% taking out outliers 

%% Residual peak gain analysis with outlier removal

an = fieldnames(ExcitRR_Ctrl);

subtypes = {'T0','T1','T2','T3'};

T0color = rgb('Thistle');
T1color = rgb('MediumOrchid');
T2color = rgb('SlateBlue');
T3color = rgb('DarkSlateBlue');

subtype_colors = {T0color,T1color,T2color,T3color};

minReaches = 10;
maxVal = 600;

GainResults = struct;

% Calculate gain/correlation for each animal and subtype

for s = 1:length(subtypes)

    subtype = subtypes{s};

    gains = NaN(length(an),1);
    intercepts = NaN(length(an),1);
    p_gain = NaN(length(an),1);
    r2 = NaN(length(an),1);
    r_vals = NaN(length(an),1);
    p_vals = NaN(length(an),1);
    n_reaches = NaN(length(an),1);

    for a = 1:length(an)

        min_field = ['ResidDecelPeak' subtype];
        max_field = ['ResidAccelPeak' subtype];

        Accelmin = abs(ExcitRR_Ctrl.(an{a}).(min_field)) * 10;
        Accelmax = ExcitRR_Ctrl.(an{a}).(max_field) * 10;

        valid_idx = ~isnan(Accelmin) & ~isnan(Accelmax) & ...
                    ~isinf(Accelmin) & ~isinf(Accelmax) & ...
                    Accelmin >= 0 & Accelmin <= maxVal & ...
                    Accelmax >= 0 & Accelmax <= maxVal;

        Accelmin = Accelmin(valid_idx);
        Accelmax = Accelmax(valid_idx);

        n_reaches(a) = length(Accelmin);

        if n_reaches(a) < minReaches
            continue
        end

        mdl = fitlm(Accelmin,Accelmax);

        intercepts(a) = mdl.Coefficients.Estimate(1);
        gains(a) = mdl.Coefficients.Estimate(2);
        p_gain(a) = mdl.Coefficients.pValue(2);
        r2(a) = mdl.Rsquared.Ordinary;

        [r_vals(a),p_vals(a)] = corr(Accelmin,Accelmax,'rows','complete');

    end

    GainResults.(subtype).animal = an;
    GainResults.(subtype).gain = gains;
    GainResults.(subtype).intercept = intercepts;
    GainResults.(subtype).p_gain = p_gain;
    GainResults.(subtype).r2 = r2;
    GainResults.(subtype).r = r_vals;
    GainResults.(subtype).p = p_vals;
    GainResults.(subtype).n_reaches = n_reaches;

end


%% Plot regression fits for each subtype

figure; hold on

xs = linspace(0, maxVal, 100);

for s = 1:length(subtypes)

    subtype = subtypes{s};
    col = subtype_colors{s};

    gains = GainResults.(subtype).gain;
    intercepts = GainResults.(subtype).intercept;

    valid = ~isnan(gains) & ~isnan(intercepts);

    fit_values = NaN(sum(valid), length(xs));
    idx = find(valid);

    for i = 1:length(idx)

        a = idx(i);

        yfit = intercepts(a) + gains(a) .* xs;
        fit_values(i,:) = yfit;

        plot(xs, yfit, ...
            'Color', [col 0.25], ...
            'LineWidth', 2, ...
            'HandleVisibility','off')

    end

    mean_fit = mean(fit_values, 1, 'omitnan');
    sem_fit = std(fit_values, 0, 1, 'omitnan') ./ sqrt(size(fit_values,1));

    fill([xs fliplr(xs)], ...
         [mean_fit + sem_fit, fliplr(mean_fit - sem_fit)], ...
         col, ...
         'FaceAlpha', 0.15, ...
         'EdgeColor', 'none', ...
         'HandleVisibility','off')

    plot(xs, mean_fit, ...
        'Color', col, ...
        'LineWidth', 4, ...
        'DisplayName', subtype)

end

plot([0 maxVal], [0 maxVal], '--k', 'LineWidth', 1.5, 'HandleVisibility','off')

xlim([0 maxVal])
ylim([0 maxVal])

xlabel('Template-subtracted peak decel magnitude')
ylabel('Template-subtracted peak accel')
title('Residual CSM gain by subtype')
legend(subtypes, 'Location', 'best')
box off
set(gca,'TickDir','in')
set(gcf, 'Position', [100 100 500 450])


%% Plot gain across CSM subtypes with r and p values

figure; hold on

all_gain_mat = NaN(length(an), length(subtypes));

for s = 1:length(subtypes)
    subtype = subtypes{s};
    all_gain_mat(:,s) = GainResults.(subtype).gain;
end

% connect same animal across subtypes
for a = 1:length(an)
    plot(1:length(subtypes), all_gain_mat(a,:), '-', ...
        'Color', [0.75 0.75 0.75], ...
        'LineWidth', 1)
end

for s = 1:length(subtypes)

    subtype = subtypes{s};
    col = subtype_colors{s};

    gains = all_gain_mat(:,s);

    scatter(ones(length(gains),1)*s, gains, 35, col, 'filled', ...
        'jitter', 'on', ...
        'jitterAmount', 0.12, ...
        'MarkerFaceAlpha', 0.65)

    mean_gain = mean(gains, 'omitnan');
    sem_gain = std(gains, 'omitnan') ./ sqrt(sum(~isnan(gains)));

    errorbar(s, mean_gain, sem_gain, 'o', ...
        'Color', col, ...
        'MarkerFaceColor', col, ...
        'MarkerEdgeColor', col, ...
        'LineWidth', 2, ...
        'MarkerSize', 10, ...
        'CapSize', 10)

    r_vals = GainResults.(subtype).r;

    mean_r = mean(r_vals, 'omitnan');

    z_vals = atanh(r_vals);
    z_vals = z_vals(~isnan(z_vals) & ~isinf(z_vals));

    if length(z_vals) >= 2
        [~, p_r] = ttest(z_vals);
    else
        p_r = NaN;
    end

    text(s, mean_gain + sem_gain + 0.05, ...
        sprintf('r = %.2f\np = %.3f', mean_r, p_r), ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'Color', col)

end

xlim([0.5 length(subtypes)+0.5])
xticks(1:length(subtypes))
xticklabels(subtypes)

ylabel('Gain: slope of residual accel vs residual decel')
xlabel('CSM subtype')
title('Residual CSM gain across subtypes')

box off
set(gca,'TickDir','in')
set(gcf, 'Position', [100 100 450 400])


%% Gain across ALL CSMs, equal weighting of animals

figure; hold on

RRcolor = rgb('Black');

xs = linspace(0,maxVal,100);

fit_values = [];
slopes = [];
intercepts_all = [];
r_vals_all = [];

for a = 1:length(an)

    Accelmin = [];
    Accelmax = [];

    for t = 0:3

        min_field = sprintf('ResidDecelPeakT%d',t);
        max_field = sprintf('ResidAccelPeakT%d',t);

        if isfield(ExcitRR_Ctrl.(an{a}),min_field)

            Accelmin = [Accelmin; ...
                abs(ExcitRR_Ctrl.(an{a}).(min_field)) * 10];

            Accelmax = [Accelmax; ...
                ExcitRR_Ctrl.(an{a}).(max_field) * 10];

        end

    end

    valid = ~isnan(Accelmin) & ~isnan(Accelmax) & ...
            ~isinf(Accelmin) & ~isinf(Accelmax) & ...
            Accelmin >= 0 & Accelmin <= maxVal & ...
            Accelmax >= 0 & Accelmax <= maxVal;

    Accelmin = Accelmin(valid);
    Accelmax = Accelmax(valid);

    if length(Accelmin) < minReaches
        continue
    end

    mdl = fitlm(Accelmin,Accelmax);

    slope = mdl.Coefficients.Estimate(2);
    intercept = mdl.Coefficients.Estimate(1);

    slopes = [slopes; slope];
    intercepts_all = [intercepts_all; intercept];

    r = corr(Accelmin,Accelmax,'rows','complete');
    r_vals_all = [r_vals_all; r];

    yfit = intercept + slope*xs;

    fit_values = [fit_values; yfit];

    plot(xs,yfit,...
        'Color',[0.7 0.7 0.7],...
        'LineWidth',1.5)

end

mean_fit = mean(fit_values,1,'omitnan');
sem_fit = std(fit_values,0,1,'omitnan') ./ sqrt(size(fit_values,1));

fill([xs fliplr(xs)],...
    [mean_fit+sem_fit fliplr(mean_fit-sem_fit)],...
    RRcolor,...
    'FaceAlpha',0.2,...
    'EdgeColor','none')

plot(xs,mean_fit,...
    'Color',RRcolor,...
    'LineWidth',5)

plot([0 maxVal],[0 maxVal],'--k','LineWidth',1.2)

mean_slope = mean(slopes,'omitnan');
sem_slope = std(slopes,'omitnan') ./ sqrt(sum(~isnan(slopes)));

mean_r = mean(r_vals_all,'omitnan');

z = atanh(r_vals_all);
z = z(~isnan(z) & ~isinf(z));

if length(z) >= 2
    [~,p_group] = ttest(z);
else
    p_group = NaN;
end

text(maxVal*0.5,maxVal*0.85,...
    sprintf(['Slope = %.2f \\pm %.2f\n' ...
             'Mean r = %.2f\n' ...
             'p = %.3f'],...
             mean_slope,sem_slope,...
             mean_r,p_group),...
    'FontWeight','bold',...
    'FontSize',11,...
    'BackgroundColor','w')

xlim([0 500])
ylim([0 500])

xlabel('Template-subtracted peak decel magnitude')
ylabel('Template-subtracted peak accel')
title('Residual gain across all CSMs')

box off
set(gca,'TickDir','in')
set(gcf, 'Position', [100 100 500 450])
