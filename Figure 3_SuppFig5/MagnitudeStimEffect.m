%Magnitude of Stim effect
AllAnimalsDiff=ExcitAnimalsDiff_CSM;
%ExcitAnimalsDiff=AllAnimalsDiff_CSM;

AllAnimmean=mean(AllAnimalsDiff);
%AllAnimmeanNtsr1=mean(AllAnimalsDiff_ntsr1);
ExcitAnimmean=mean(ExcitAnimalsDiff);

%allcolor=rgb('CadetBlue');
allcolor=rgb('SteelBlue');
excitcolor=rgb('DodgerBlue');
ntsr1color=rgb('LightSteelBlue');
sz=60;

figure()
scatter(AllAnimalsDiff(:,1),AllAnimalsDiff(:,2),sz,'filled','MarkerFaceColor',allcolor,'MarkerEdgeColor','k')
hold on
scatter(ExcitAnimalsDiff(:,1),ExcitAnimalsDiff(:,2),sz,'filled','MarkerFaceColor',excitcolor,'MarkerEdgeColor','k')
% hold on
% scatter(AllAnimalsDiff_ntsr1(:,1),AllAnimalsDiff_ntsr1(:,2),sz,'filled','MarkerFaceColor',ntsr1color,'MarkerEdgeColor','k')
hold on
%quiver(0,0,AllAnimmean(:,1),AllAnimmean(:,2),'LineWidth',5,'Color',[ 0.3711    0.6172    0.6250])
quiver(0,0,AllAnimmean(:,1),AllAnimmean(:,2),'LineWidth',5,'Color',[0.2734    0.5078    0.7031])
hold on
quiver(0,0,ExcitAnimmean(:,1),ExcitAnimmean(:,2),'LineWidth',5,'Color',[ 0.1172    0.5625    1.0000 ])
% hold on
% quiver(0,0,AllAnimmeanNtsr1(:,1),AllAnimmeanNtsr1(:,2),'LineWidth',5,'Color',[0.6875    0.7656    0.8672])

% xlim([-1.5 0.5])
% ylim([-1.5 0.5])
ax=gca;
ax.XAxisLocation='origin';
ax.YAxisLocation='origin';


%%
% --- All Animals ---
mean_all_outward = mean(AllAnimalsDiff(:,1));
sem_all_outward  = std(AllAnimalsDiff(:,1)) / sqrt(size(AllAnimalsDiff,1));

mean_all_upward  = mean(AllAnimalsDiff(:,2));
sem_all_upward   = std(AllAnimalsDiff(:,2)) / sqrt(size(AllAnimalsDiff,1));

% --- Excitatory Animals ---
mean_excit_outward = mean(ExcitAnimalsDiff(:,1));
sem_excit_outward  = std(ExcitAnimalsDiff(:,1)) / sqrt(size(ExcitAnimalsDiff,1));

mean_excit_upward  = mean(ExcitAnimalsDiff(:,2));
sem_excit_upward   = std(ExcitAnimalsDiff(:,2)) / sqrt(size(ExcitAnimalsDiff,1));

% --- Print formatted results to Command Window ---
fprintf('=== Outward (Column 1) ===\n');
fprintf('All Group   : %.3f ± %.3f\n', mean_all_outward, sem_all_outward);
fprintf('Excit Group : %.3f ± %.3f\n', mean_excit_outward, sem_excit_outward);

fprintf('\n=== Upward (Column 2) ===\n');
fprintf('All Group   : %.3f ± %.3f\n', mean_all_upward, sem_all_upward);
fprintf('Excit Group : %.3f ± %.3f\n', mean_excit_upward, sem_excit_upward);

%%
% Example: test if AllMag is normally distributed
AllMag = sqrt(sum(AllAnimalsDiff.^2, 2));
[h,p] = kstest((AllMag - mean(AllMag)) / std(AllMag)); % Kolmogorov-Smirnov test

if h == 1
    disp('Data is NOT normally distributed → use Wilcoxon signed-rank test.');
else
    disp('Data is normally distributed → use one-sample t-test.');
end
%%
AllMag = sqrt(sum(ExcitAnimalsDiff.^2, 2));
[h,p] = kstest((AllMag - mean(AllMag)) / std(AllMag)); % Kolmogorov-Smirnov test

if h == 1
    disp('Data is NOT normally distributed → use Wilcoxon signed-rank test.');
else
    disp('Data is normally distributed → use one-sample t-test.');
end
%% Test Value 
% Compute vector magnitudes
AllMag   = sqrt(sum(AllAnimalsDiff.^2, 2));
ExcitMag = sqrt(sum(ExcitAnimalsDiff.^2, 2));

% % One-sample t-tests against 0
[~, p_all_vs_zero] = ttest(AllMag, 0);
[~, p_excit_vs_zero] = ttest(ExcitMag, 0);

%Alternatively, use sign-rank if non-parametric
% p_all_vs_zero = signrank(AllMag, 0);
% p_excit_vs_zero = signrank(ExcitMag, 0);

%%
% Outward velocity (column 1)
[~, p_outward] = ttest2(AllAnimalsDiff(:,1), ExcitAnimalsDiff(:,1));

% Upward velocity (column 2)
[~, p_upward] = ttest2(AllAnimalsDiff(:,2), ExcitAnimalsDiff(:,2));

% Display results
fprintf('Comparison between groups:\n');
fprintf('Outward velocity difference p = %.4f\n', p_outward);
fprintf('Upward velocity difference p = %.4f\n', p_upward);
%% Magnitude of effect for Intarn AND IntARN HF
HFAnimmean=mean(HFAnimalsDiff);
ExcitAnimmean=mean(ExcitAnimalsDiff);

allcolor=rgb('LightSteelBlue');
excitcolor=rgb('DodgerBlue');
sz=60;

figure()
scatter(HFAnimalsDiff(:,1),HFAnimalsDiff(:,2),sz,'filled','MarkerFaceColor',allcolor,'MarkerEdgeColor','k')
hold on
scatter(ExcitAnimalsDiff(:,1),ExcitAnimalsDiff(:,2),sz,'filled','MarkerFaceColor',excitcolor,'MarkerEdgeColor','k')
hold on
quiver(0,0,HFAnimmean(:,1),HFAnimmean(:,2),'LineWidth',5,'Color',[0.6875    0.7656    0.8672])
hold on
quiver(0,0,ExcitAnimmean(:,1),ExcitAnimmean(:,2),'LineWidth',5,'Color',[ 0.1172    0.5625    1.0000 ])
xlim([-1.5 0.5])
ylim([-1.5 0.5])
ax=gca;
ax.XAxisLocation='origin';
ax.YAxisLocation='origin';


%% Calculate magnitude for  RR and NRR - use data structures in stim magnitude folder
AllAnimmeanNRR=mean(AllAnimalsDiffNRR);
ExcitAnimmeanNRR=mean(ExcitAnimalsDiffNRR);

AllAnimmeanRR=mean(AllAnimalsDiffRR);
ExcitAnimmeanRR=mean(ExcitAnimalsDiffRR);

allRR=rgb('Teal');
allNRR=rgb('MediumAquamarine');

excitcolorRR=rgb('DodgerBlue');
excitcolorNRR=rgb('PowderBlue');
sz=60;

figure()
scatter(AllAnimalsDiffRR(:,1),AllAnimalsDiffRR(:,2),sz,'filled','MarkerFaceColor',allRR,'MarkerEdgeColor','k')
hold on
scatter(ExcitAnimalsDiffRR(:,1),ExcitAnimalsDiffRR(:,2),sz,'filled','MarkerFaceColor',excitcolorRR,'MarkerEdgeColor','k')
hold on
scatter(AllAnimalsDiffNRR(:,1),AllAnimalsDiffNRR(:,2),sz,'filled','MarkerFaceColor',allNRR,'MarkerEdgeColor','k')
hold on
scatter(ExcitAnimalsDiffNRR(:,1),ExcitAnimalsDiffNRR(:,2),sz,'filled','MarkerFaceColor',excitcolorNRR,'MarkerEdgeColor','k')
hold on
quiver(0,0,AllAnimmeanRR(:,1),AllAnimmeanRR(:,2),'LineWidth',5,'Color',[ 0 0.5 0.5])
hold on
quiver(0,0,ExcitAnimmeanRR(:,1),ExcitAnimmeanRR(:,2),'LineWidth',5,'Color',[  0.1172    0.5625    1.0000])
hold on
quiver(0,0,AllAnimmeanNRR(:,1),AllAnimmeanNRR(:,2),'LineWidth',5,'Color',[0.3984    0.8008    0.6641])
hold on
quiver(0,0,ExcitAnimmeanNRR(:,1),ExcitAnimmeanNRR(:,2),'LineWidth',5,'Color',[  0.6875    0.8750    0.8984])
xlim([-1.5 0.5])
ylim([-1.5 0.5])
ax=gca;
ax.XAxisLocation='origin';
ax.YAxisLocation='origin';

%% INdividual Plots 
figure()
scatter(AllAnimalsDiffRR(:,1),AllAnimalsDiffRR(:,2),sz,'filled','MarkerFaceColor',allRR,'MarkerEdgeColor','k')
hold on
scatter(AllAnimalsDiffNRR(:,1),AllAnimalsDiffNRR(:,2),sz,'filled','MarkerFaceColor',allNRR,'MarkerEdgeColor','k')
hold on
quiver(0,0,AllAnimmeanRR(:,1),AllAnimmeanRR(:,2),'LineWidth',5,'Color',[ 0 0.5 0.5])
hold on
quiver(0,0,AllAnimmeanNRR(:,1),AllAnimmeanNRR(:,2),'LineWidth',5,'Color',[0.3984    0.8008    0.6641])

xlim([-1.5 0.5])
ylim([-1.5 0.5])
ax=gca;
ax.XAxisLocation='origin';
ax.YAxisLocation='origin';
title('Change in Velocity (cm/s)')

figure()
scatter(ExcitAnimalsDiffRR(:,1),ExcitAnimalsDiffRR(:,2),sz,'filled','MarkerFaceColor',excitcolorRR,'MarkerEdgeColor','k')
hold on
scatter(ExcitAnimalsDiffNRR(:,1),ExcitAnimalsDiffNRR(:,2),sz,'filled','MarkerFaceColor',excitcolorNRR,'MarkerEdgeColor','k')
hold on
quiver(0,0,ExcitAnimmeanRR(:,1),ExcitAnimmeanRR(:,2),'LineWidth',5,'Color',[  0.1172    0.5625    1.0000])
hold on
quiver(0,0,ExcitAnimmeanNRR(:,1),ExcitAnimmeanNRR(:,2),'LineWidth',5,'Color',[  0.6875    0.8750    0.8984])
xlim([-1.5 0.5])
ylim([-1.5 0.5])
ax=gca;
ax.XAxisLocation='origin';
ax.YAxisLocation='origin';
title('Change in Velocity (cm/s)')