% Success Rate Comparison : CSM versus NoCSM
an=fieldnames(ExcitRR_Sham);
SuccessRate=[];

for a=1:length(an)
    CSMreaches=ExcitRR_Sham.(an{a}).ReReaches;
    NoCSMreaches=ExcitRR_Sham.(an{a}).NoReReach;

    [S_CSM, F_CSM] = SuccessFilter(CSMreaches);
    [S_NoCSM, F_NoCSM] = SuccessFilter(NoCSMreaches);

    SuccessRate.(an{a}).successrate.csm = length(fieldnames(S_CSM)) / length(fieldnames(CSMreaches));
    SuccessRate.(an{a}).successrate.nocsm = length(fieldnames(S_NoCSM)) / length(fieldnames(NoCSMreaches));
end
%%
% Initialize
an = fieldnames(ExcitRR_Sham);
SuccessRate = struct();
All_CSM_Success = [];
All_CSM_Total = [];
All_NoCSM_Success = [];
All_NoCSM_Total = [];

% Calculate success rates per animal
for a = 1:length(an)
    CSMreaches = ExcitRR_Sham.(an{a}).ReReaches;
    NoCSMreaches = ExcitRR_Sham.(an{a}).NoReReach;

    [S_CSM, ~] = SuccessFilter(CSMreaches);
    [S_NoCSM, ~] = SuccessFilter(NoCSMreaches);

    % Success rates
    n_CSM_success = length(fieldnames(S_CSM));
    n_CSM_total = length(fieldnames(CSMreaches));
    n_NoCSM_success = length(fieldnames(S_NoCSM));
    n_NoCSM_total = length(fieldnames(NoCSMreaches));

    SuccessRate.(an{a}).successrate.csm = n_CSM_success / n_CSM_total;
    SuccessRate.(an{a}).successrate.nocsm = n_NoCSM_success / n_NoCSM_total;

    % Store for group-level analysis
    All_CSM_Success(end+1) = n_CSM_success;
    All_CSM_Total(end+1) = n_CSM_total;
    All_NoCSM_Success(end+1) = n_NoCSM_success;
    All_NoCSM_Total(end+1) = n_NoCSM_total;
end

%% Group-level success rates and confidence intervals (Bootstrap)
nBoot = 10000;

% Bootstrap for CSM
CSM_rates = All_CSM_Success ./ All_CSM_Total;
boot_csm = bootstrp(nBoot, @mean, CSM_rates);
ci_csm = prctile(boot_csm, [2.5 97.5]);

% Bootstrap for NoCSM
NoCSM_rates = All_NoCSM_Success ./ All_NoCSM_Total;
boot_nocsm = bootstrp(nBoot, @mean, NoCSM_rates);
ci_nocsm = prctile(boot_nocsm, [2.5 97.5]);

% Group means
mean_csm = mean(CSM_rates);
mean_nocsm = mean(NoCSM_rates);

%% Plot
figure;
hold on;
bar(1, mean_csm, 0.5, 'FaceColor', [0.2 0.6 0.8]);
bar(2, mean_nocsm, 0.5, 'FaceColor', [0.8 0.4 0.4]);
errorbar(1, mean_csm, mean_csm - ci_csm(1), ci_csm(2) - mean_csm, 'k', 'LineWidth', 1.5);
errorbar(2, mean_nocsm, mean_nocsm - ci_nocsm(1), ci_nocsm(2) - mean_nocsm, 'k', 'LineWidth', 1.5);
xlim([0.5 2.5]);
xticks([1 2]);
xticklabels({'CSM', 'NoCSM'});
ylabel('Success Rate');
title('Success Rate Comparison: CSM vs NoCSM');
hold off;

%% Fisher's Exact Test
% Pool across animals for a formal test
total_CSM_success = sum(All_CSM_Success);
total_CSM_fail = sum(All_CSM_Total) - total_CSM_success;
total_NoCSM_success = sum(All_NoCSM_Success);
total_NoCSM_fail = sum(All_NoCSM_Total) - total_NoCSM_success;

contingency = [total_CSM_success, total_CSM_fail; total_NoCSM_success, total_NoCSM_fail];

[~, p_fisher] = fishertest(contingency);

fprintf('Fisher''s Exact Test p-value: %.4f\n', p_fisher);
