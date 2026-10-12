%%
%% Use ReReach Filter - Control Reaches 
ExcitRR_Ctrl=[]; % Savethis structure to use for graphing rereach timing in sep. script
for a=1:length(an)

    Creaches=ChR2.(an{a}).Control;

 
    % ReReach
    %[ReReaches,StimRR,MaxVel]=ReReachFilter_accel_ctrlreaches(Creaches); % currently looking within 200ms 
    [ReReaches,StimRR,MaxVel]=ReReachFilter_accel_ctrlreaches_edit(Creaches);% new code 
    RRKin=BasicKinematics(ReReaches);

    C = AlignThresh(ReReaches,thresh,dirdir,200);
    C=stimvel_accel(C);

    % find difference in fieldnames to get stim re-reaches without RR
    fields1 = fieldnames(Creaches);
    fields2 = fieldnames(ReReaches);

    % Find field names present in structure1 but not in structure2
    fields_only_in_structure1 = setdiff(fields1, fields2);
    NoReach=struct;
    for g=1:length(fields_only_in_structure1)
        if isempty(fields_only_in_structure1)
            NoReach=[];
        else
        NoReach.(fields_only_in_structure1{g})=Creaches.(fields_only_in_structure1{g});
        end
    end
    
    if isempty(NoReach)
        NRkin=[];
        CN=[];
    else
        NRKin=BasicKinematics(NoReach);
        CN = AlignThresh(NoReach,thresh,dirdir,500);
        CN=stimvel_accel(CN); %causing issues in the code 
    end 

    [S_C, F_C] = SuccessFilter(ReReaches);
    [S_S, F_S] = SuccessFilter(NoReach);
    
    ExcitRR_Ctrl.(an{a}).successrate.ReReach = length(fieldnames(S_C)) / length(fieldnames(Creaches));
    ExcitRR_Ctrl.(an{a}).successrate.NoReReach = length(fieldnames(S_S)) / length(fieldnames(Creaches));
    ExcitRR_Ctrl.(an{a}).CtrlReaches=Creaches;
    ExcitRR_Ctrl.(an{a}).C=C;
    ExcitRR_Ctrl.(an{a}).ReReaches=ReReaches;
    ExcitRR_Ctrl.(an{a}).RRKin=RRKin;
    ExcitRR_Ctrl.(an{a}).NoReReach=NoReach;
    ExcitRR_Ctrl.(an{a}).NRKin=NRKin;
    ExcitRR_Ctrl.(an{a}).CN=CN;
    ExcitRR_Ctrl.(an{a}).StimRR=StimRR;
    ExcitRR_Ctrl.(an{a}).MaxVel=MaxVel;
    % get % of reaches with a rereach
    cnum=length(fields1);
    rrnum=length(fields2);
    percentrr=rrnum/cnum;
    ExcitRR_Ctrl.(an{a}).Cnum=cnum;
    ExcitRR_Ctrl.(an{a}).RRnum=rrnum;
    ExcitRR_Ctrl.(an{a}).PercentRR=percentrr;
  
end 

%% Use ReReach Filter- STIM reaches 
ExcitRR=[]; % Savethis structure to use for graphing rereach timing in sep. script
for a=1:length(an)
    Sreaches=ChR2.(an{a}).Stim;
  

    % Stim Re
    [ReReaches,StimRR,MaxVel]=ReReachFilter_accel(Sreaches);
    RRKin=BasicKinematics(ReReaches);
    
    S = AlignStim(ReReaches,200);
    S=stimvel_accel(S);
    % find difference in fieldnames to get stim re-reaches without RR
    fields1 = fieldnames(Sreaches);
    fields2 = fieldnames(ReReaches);

    % Find field names present in structure1 but not in structure2
    fields_only_in_structure1 = setdiff(fields1, fields2);
    NoReach=struct;
    for g=1:length(fields_only_in_structure1)
        if isempty(fields_only_in_structure1)
            NoReach=[];
        else
        NoReach.(fields_only_in_structure1{g})=Sreaches.(fields_only_in_structure1{g});
        end
    end
    
    if isempty(NoReach)
        NRkin=[];
        SN=[];
    else
        NRKin=BasicKinematics(NoReach);
        SN = AlignStim(NoReach,500);
        %SN=stimvel_accel(SN); %causing issues in the code 
    end 

    [S_C, F_C] = SuccessFilter(ReReaches);
    [S_S, F_S] = SuccessFilter(NoReach);
    
    ExcitRR.(an{a}).successrate.ReReach = length(fieldnames(S_C)) / length(fieldnames(Sreaches));
    ExcitRR.(an{a}).successrate.NoReReach = length(fieldnames(S_S)) / length(fieldnames(Sreaches));
    ExcitRR.(an{a}).ReReaches=ReReaches;
    ExcitRR.(an{a}).RRKin=RRKin;
    ExcitRR.(an{a}).S=S;
    ExcitRR.(an{a}).NoReReach=NoReach;
    ExcitRR.(an{a}).NRKin=NRKin;
    ExcitRR.(an{a}).SN=SN;
    ExcitRR.(an{a}).StimRR=StimRR;
    ExcitRR.(an{a}).MaxVel=MaxVel;
    snum=length(fields1);
    rrnum=length(fields2);
    percentrr=rrnum/snum;
    ExcitRR.(an{a}).Snum=snum;
    ExcitRR.(an{a}).RRnum=rrnum;
    ExcitRR.(an{a}).PercentRR=percentrr;
end 

%% Find min and max of accel during CSM time window ( 100ms after stim ) 
an=fieldnames(ChR2);
for a=1:length(an)
    stimvel=ExcitRR.(an{a}).S.stimvel;
    Accelmin=[];
    Accelmax=[];
   stimvelsize=size(stimvel);
   num=(stimvelsize(3));

    for b=1:num
        reach=stimvel(:,:,b);
        window=reach(21:31,:);% row numbers for time of stim on to 100ms after - only for reaches with 200ms 
        %window=reach(51:61,:);% row numbers for time of stim on to 100ms after - only for reaches with 200ms 
        accelmin=min(window(:,8));
        accelmax=max(window(:,8));
        Accelmin=[Accelmin;accelmin];
        Accelmax=[Accelmax;accelmax];
    end

    ExcitRR.(an{a}).AccelMin=Accelmin;
    ExcitRR.(an{a}).AccelMax=Accelmax;
    clear Accelmin Accelmax 
end

%%
% CSM Gain error 
% Have to have an ExcitRR structure with Min and Max 

% For IntARN Data 
%  Plotting
an=fieldnames(ExcitRR);
vd_corr = []; 
vd_corr_p = []; 

vd_corr_rand = [];
vd_corr_p_rand = [];

fits = [];

RRcolor = rgb('DodgerBlue');
%RRcolor = rgb('Teal');

% For pooling all data
all_accelmin = [];
all_accelmax = [];

% for plotting all linear correlations
figure; hold on

for a = 1:length(an)
    
    % Take absolute value of AccelMin
    Accelmin = abs(ExcitRR.(an{a}).AccelMin * 10);
    Accelmax = ExcitRR.(an{a}).AccelMax * 10;

    % accumulate for global correlation
    all_accelmin = [all_accelmin; Accelmin];
    all_accelmax = [all_accelmax; Accelmax];

    % compute correlation coefficient and p-values
    [vd_corr(a), vd_corr_p(a)] = corr(Accelmin, Accelmax);
    
    % shuffle control for ReReach
    corr_temp = NaN(1,1000); 
    p_temp = NaN(1,1000);
    for j = 1:1000
        rand_inx = randsample(length(Accelmin), length(Accelmin));
        [corr_temp(j), p_temp(j)] = corr(Accelmin(rand_inx), Accelmax);
    end
    vd_corr_rand(a) = mean(corr_temp);
    vd_corr_p_rand(a) = mean(p_temp);

    % linear fit
    pf = polyfit(Accelmin, Accelmax, 1);
    fits(a,:) = pf;

    % for plotting all linear correlations
    %xs = [0 400];  % now using absolute values
    %xs = [0 500];  % now using absolute values
    xs = [-100 500];  % now using absolute values
    pv = polyval(pf, xs);

    if vd_corr_p(a) <= 0.05
        %col = [0.1172 0.5625 1.0000];
        col = [0 0.5 0.5];
    else
        col = [0.5 0.5 0.5];
    end
    plot(xs, pv, 'Color', col, 'LineWidth', 1.5)

    % Add individual correlation coefficient as text
    midx = mean(xs);
    midy = mean(pv);
    text(midx, midy, sprintf('r=%.2f', vd_corr(a)), ...
        'FontSize', 8, 'Color', col(1:3), 'HorizontalAlignment', 'center')
end


% --- Summary scatter with Fisher z-transform analysis ---

% Fisher z-transform correlation stats
z_corr = 0.5 * log((1 + vd_corr) ./ (1 - vd_corr));  % Fisher z-transform
mean_z = mean(z_corr);
sem_z = std(z_corr) / sqrt(length(z_corr));
mean_r = (exp(2*mean_z) - 1) / (exp(2*mean_z) + 1);  % Back-transform
[h, p_val] = ttest(z_corr);

% Print summary
fprintf('Number of significant animals (p<=0.05): %d out of %d\n', sum(vd_corr_p <= 0.05), length(vd_corr));
fprintf('Mean correlation (untransformed): %.3f ± %.3f (SEM)\n', mean(vd_corr), std(vd_corr)/sqrt(length(vd_corr)));
fprintf('Mean correlation (Fisher z-transformed): %.3f ± %.3f (SEM)\n', mean_r, sem_z);
fprintf('Significance of average correlation (Fisher z t-test): p = %.4f\n', p_val);
fprintf('Mean shuffle correlation: %.3f\n', mean(vd_corr_rand));


% --- Plot mean linear fit with error shading ---
%xs = linspace(0, 400, 100);  % x-range for positive absolute values
xs = linspace(-100, 500, 100);  % x-range for positive absolute values
fit_values = NaN(length(an), length(xs));

for i = 1:length(an)
    fit_values(i, :) = polyval(fits(i,:), xs);
end

mean_fit = mean(fit_values, 1);
sem_fit = std(fit_values, 0, 1) / sqrt(size(fit_values, 1));

% fill([xs fliplr(xs)], ...
%      [mean_fit + sem_fit, fliplr(mean_fit - sem_fit)], ...
%      [0.1172 0.5625 1.0000], 'FaceAlpha', 0.2, 'EdgeColor', 'none');

fill([xs fliplr(xs)], ...
     [mean_fit + sem_fit, fliplr(mean_fit - sem_fit)], ...
     [0 0.5 0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'none');

%plot(xs, mean_fit, 'Color', [0.1172 0.5625 1.0000], 'LineWidth', 4)
plot(xs, mean_fit, 'Color', [0 0.5 0.5], 'LineWidth', 4)

% --- Plot unity line ---
%plot([0 400], [0 400], '--k', 'LineWidth', 1.5)
%plot([0 500], [0 500], '--k', 'LineWidth', 1.5)
plot([-100 500], [-100 500], '--k', 'LineWidth', 1.5)

% --- Compute and plot overall correlation coefficient ---
[r_all, p_all] = corr(all_accelmin, all_accelmax);

text_x = 300;  % updated x for label positioning
text_y = polyval(mean(fits), text_x);
text(text_x, text_y, sprintf('mean r = %.2f', r_all), ...
    'FontSize', 12, 'Color', [0.1172 0.5625 1.0000], ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')



% --- Final plot formatting ---
% xlim([0 400])
% ylim([0 400])
% xlim([0 500])
% ylim([0 500])
xlim([-100 500])
ylim([-100 500])
%axis square
xlabel('Peak Decel (|value|)')
ylabel('CSM Peak Accel')

% --- Summary scatter ---


% --- Summary scatter plot ---
figure; hold on
scatter(ones(length(an),1), vd_corr, 15, RRcolor, 'filled', ...
    'jitter', 'on', 'jitterAmount', .15, 'MarkerFaceAlpha', .5)
scatter(ones(length(an),1)*1.5, vd_corr_rand, 30, [.8 .8 .8], 'filled', ...
    'jitter', 'on', 'jitterAmount', .15)
scatter(1, mean(vd_corr), 170, RRcolor, 'filled')
scatter(1.5, mean(vd_corr_rand), 170, [.6 .6 .6], 'filled')
%ylim([-0.7 .1])
xlim([0.75 1.75])

set(gca, 'XTick', [1 1.5], 'XTickLabel', {'Real', 'Shuffle'})
ylabel('Correlation coefficient (r)')
set(gcf, 'Position', [100 100 300 400])

% --- Add text for mean Fisher transformed correlation and p-value ---
text(1, 0.05, sprintf('mean z(r) = %.2f\np = %.3f', mean_r, p_val), ...
    'HorizontalAlignment', 'center', 'FontSize', 10, ...
    'FontWeight', 'bold', 'Color', RRcolor)
%% Slightly different method 
% Correlation analysis per animal and across animals

vd_corr = []; 
vd_corr_p = []; 
vd_corr_rand = [];
vd_corr_p_rand = [];
fits = [];

%RRcolor = rgb('DodgerBlue');
RRcolor = rgb('Teal');

% For pooling all data
all_accelmin = [];
all_accelmax = [];

% for plotting all linear correlations
figure; hold on

for a = 1:length(an)
    
    Accelmin = ExcitRR.(an{a}).AccelMin * 10;
    Accelmax = ExcitRR.(an{a}).AccelMax * 10;

    % accumulate for global correlation
    all_accelmin = [all_accelmin; Accelmin];
    all_accelmax = [all_accelmax; Accelmax];

    % compute correlation coefficient and p-values
    [vd_corr(a), vd_corr_p(a)] = corr(Accelmin, Accelmax);
    
    % shuffle control
    corr_temp = NaN(1,1000); 
    p_temp = NaN(1,1000);
    for j = 1:1000
        rand_inx = randsample(length(Accelmin), length(Accelmin));
        [corr_temp(j), p_temp(j)] = corr(Accelmin(rand_inx), Accelmax);
    end
    vd_corr_rand(a) = mean(corr_temp);
    vd_corr_p_rand(a) = mean(p_temp);

    % linear fit
    pf = polyfit(Accelmin, Accelmax, 1);
    fits(a,:) = pf;

    % for plotting all linear correlations
   % xs = [0 -400];  % flipped x-axis range
    xs = [-100 500];  % now using absolute values
    pv = polyval(pf, xs);

    if vd_corr_p(a) <= 0.05
        %col = [0.1172 0.5625 1.0000];
        col = [0 0.5 0.5];
    else
        col = [0.5 0.5 0.5];
    end
    plot(xs, pv, 'Color', col, 'LineWidth', 1.5)

    % Add individual correlation coefficient as text
    midx = mean(xs);
    midy = mean(pv);
    text(midx, midy, sprintf('r=%.2f', vd_corr(a)), ...
        'FontSize', 8, 'Color', col(1:3), 'HorizontalAlignment', 'center')
end

% --- Plot mean linear fit with error shading ---
%xs = linspace(-400, 0, 100);  % smooth x-range
xs = linspace(-100, 500, 100);  % x-range for positive absolute values
fit_values = NaN(length(an), length(xs));
for i = 1:length(an)
    fit_values(i, :) = polyval(fits(i,:), xs);
end
mean_fit = mean(fit_values, 1);
sem_fit = std(fit_values, 0, 1) / sqrt(size(fit_values, 1));

% Plot shaded error
% fill([xs fliplr(xs)], ...
%      [mean_fit + sem_fit, fliplr(mean_fit - sem_fit)], ...
%      [0.1172 0.5625 1.0000], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% plot(xs, mean_fit, 'Color', [0.1172 0.5625 1.0000], 'LineWidth', 4)

fill([xs fliplr(xs)], ...
     [mean_fit + sem_fit, fliplr(mean_fit - sem_fit)], ...
     [0 0.5 0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(xs, mean_fit, 'Color', [0 0.5 0.5], 'LineWidth', 4)

% plot unity line
plot([-100 500], [-100 500], '--k', 'LineWidth', 1.5)
% --- Overall correlation across all animals ---
[r_all, p_all] = corr(all_accelmin, all_accelmax);

% Add overall r to plot
text_x = -200;
text_y = polyval(mean(fits), text_x);
% % --- Add text for mean Fisher transformed correlation and p-value ---
text(1, 0.05, sprintf('mean z(r) = %.2f\np = %.3f', mean_r, p_val), ...
    'HorizontalAlignment', 'center', 'FontSize', 10, ...
    'FontWeight', 'bold', 'Color', RRcolor)

xlim([-400 0])
ylim([-100 500])
set(gca, 'XDir', 'reverse')
ylabel('CSM peak accel')
xlabel('Peak Decel')

% --- Summary scatter with Fisher z-transform analysis ---

% Fisher z-transform correlation stats
z_corr = 0.5 * log((1 + vd_corr) ./ (1 - vd_corr));  % Fisher z-transform
mean_z = mean(z_corr);
sem_z = std(z_corr) / sqrt(length(z_corr));
mean_r = (exp(2*mean_z) - 1) / (exp(2*mean_z) + 1);  % Back-transform
[h, p_val] = ttest(z_corr);

% Print summary
fprintf('Number of significant animals (p<=0.05): %d out of %d\n', sum(vd_corr_p <= 0.05), length(vd_corr));
fprintf('Mean correlation (untransformed): %.3f ± %.3f (SEM)\n', mean(vd_corr), std(vd_corr)/sqrt(length(vd_corr)));
fprintf('Mean correlation (Fisher z-transformed): %.3f ± %.3f (SEM)\n', mean_r, sem_z);
fprintf('Significance of average correlation (Fisher z t-test): p = %.4f\n', p_val);
fprintf('Mean shuffle correlation: %.3f\n', mean(vd_corr_rand));

% --- Summary scatter plot ---
figure; hold on
scatter(ones(length(an),1), vd_corr, 15, RRcolor, 'filled', ...
    'jitter', 'on', 'jitterAmount', .15, 'MarkerFaceAlpha', .5)
scatter(ones(length(an),1)*1.5, vd_corr_rand, 30, [.8 .8 .8], 'filled', ...
    'jitter', 'on', 'jitterAmount', .15)
scatter(1, mean(vd_corr), 170, RRcolor, 'filled')
scatter(1.5, mean(vd_corr_rand), 170, [.6 .6 .6], 'filled')
ylim([-0.7 .1])
xlim([0.75 1.75])

set(gca, 'XTick', [1 1.5], 'XTickLabel', {'Real', 'Shuffle'})
ylabel('Correlation coefficient (r)')
set(gcf, 'Position', [100 100 300 400])

% --- Add text for mean Fisher transformed correlation and p-value ---
text(1, 0.05, sprintf('mean z(r) = %.2f\np = %.3f', mean_r, p_val), ...
    'HorizontalAlignment', 'center', 'FontSize', 10, ...
    'FontWeight', 'bold', 'Color', RRcolor)

%% Correlation gain comparison

% Your correlation coefficients from two conditions
r_condition1 = vd_corr1;  % e.g., 50 values from condition 1
r_condition2 = vd_corr;  % e.g., 42 values from condition 2

% Fisher r-to-z transform
z_condition1 = atanh(r_condition1);
z_condition2 = atanh(r_condition2);

% Perform two-sample t-test (unequal sample sizes, unpaired)
[~, p, ~, stats] = ttest2(z_condition1, z_condition2);

% Report results
fprintf('Two-sample t-test on Fisher z-transformed r-values:\n');
fprintf('t(%d) = %.2f, p = %.3g\n', stats.df, stats.tstat, p);
