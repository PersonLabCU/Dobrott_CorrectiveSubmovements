% response to reviewers
% IntA neuron reliability and spike sorting (Figure 2)
% I understand Reviewer 3's concern about trial to trial variability in firing,
%but I don't think raw voltage traces should be required here, that isn't standard practice in current systems neuroscience,
% and requiring it would impose a different bar than is typically applied. Standard cluster quality metrics
% (e.g., isolation distance, L ratio, or ISI violation rates) are an appropriate and sufficient substitute to demonstrate
% that the observed variability reflects genuine neural activity rather than sorting error. I'd ask the authors to include one of these rather than raw traces.
% This part is only for the ISI violation
clear
addpath E:/Projects/slanCM/
win=100;
ax=-win:win;
udir='E:\classifiedCells\opto';
fs= dir([udir '\*_opto*']);
c=0;
for i=1:length(fs)
    load([udir '\' fs(i).name]);
    for ii=1:length(cellData)
        c=c+1;
        n=round(length(cellData(ii).spike_times)/1);
        spk_times = cellData(ii).spike_times(1:n); %in s
        l=max(spk_times);
        raster = zeros(1,ceil(max(spk_times)*1000));
        raster(ceil(spk_times.*1000))=1;

        ix= find(raster==1);
        index= repmat(ix(:),1,length(ax))+repmat(ax,length(ix),1);
        index=index(find(index(:,1)>0 & index(:,end)<length(raster)),:);
        tmp= raster(index);
        tmp(:,find(ax==0))=0;
        all_accg(c,:)=nanmean(tmp,1);
        mat(:,2)=raster;

        mat(:,1)=1:size(mat,1);
        mat(:,1)=mat(:,1)./1000;
        accg_3d= get_3d_ACCG(mat,10,win);
        all_acc3d(c,:,:)=accg_3d;
        [Violation_rate(c),nviol] = calculate_isi_violations(spk_times);
        violation_fraction(c)=nviol(2)/length(spk_times);

        all_mFR(c)=cellData(ii).fr;
        all_corr(c)=cellData(ii).Corr;
        all_mode(c)=1;
        clear mat
    end
end
%%
% c=0;
udir='E:\classifiedCells\stim';
fs= dir([udir '\*_CbN*']);
for i=1:length(fs)
    load([udir '\' fs(i).name]);
    for ii=1:length(cellData)
        c=c+1;
        n=round(length(cellData(ii).spike_times)/2);
        spk_times = cellData(ii).spike_times(1:n); %in s
        l=max(spk_times);
        raster = zeros(1,ceil(max(spk_times)*1000));
        raster(ceil(spk_times.*1000))=1;

        ix= find(raster==1);
        index= repmat(ix(:),1,length(ax))+repmat(ax,length(ix),1);
        index=index(find(index(:,1)>0 & index(:,end)<length(raster)),:);
        tmp= raster(index);
        tmp(:,find(ax==0))=0;
        all_accg(c,:)=nanmean(tmp,1);
        mat(:,2)=raster;

        mat(:,1)=1:size(mat,1);
        mat(:,1)=mat(:,1)./1000;
        accg_3d= get_3d_ACCG(mat,10,win);
        all_acc3d(c,:,:)=accg_3d;
        [Violation_rate(c),nviol] = calculate_isi_violations(spk_times);
        violation_fraction(c)=nviol(2)/length(spk_times);
        all_mFR(c)=cellData(ii).fr;
        all_corr(c)=cellData(ii).Corr;
        all_mode(c)=2;
        clear mat
    end
end
Violation_rate(find(Violation_rate==0))=.00001;
%%
for i=1:size(all_accg,1)
    norm_accg(i,:)=all_accg(i,:)./max(all_accg(i,:));
    norm_all_acc3d(i,:,:)=all_acc3d(i,:,:)./max(max(all_acc3d(i,:,:)));
end
%%
all_accg=all_accg;
all_acc3d=all_acc3d;

%%
subplot(3,3,1)
errorbar(ax,mean(all_accg,1,'omitmissing'),std(all_accg,[],1,'omitmissing')./sqrt(size(all_accg,1)),'CapSize',0,'Color','k')
xlabel('lags (ms)')
ylabel('firing rate (sp/s)')
set(gca,'TickDir','out'); box off

subplot(3,3,2)
plot(ax,squeeze(mean(all_acc3d,1,'omitmissing')),'LineWidth',1.5)
xlabel('lags (ms)')
ylabel('firing rate (sp/s)')
set(gca,'TickDir','out'); box off

subplot(3,3,3)
pcolor(ax,1:10,squeeze(mean(all_acc3d,1,'omitmissing')))
shading flat
colormap(flipud(slanCM('rdbu',265)))
xlabel('lags (ms)')
ylabel('firing rate decile')
clim([0 100])
set(gca,'TickDir','out'); box off

subplot(3,3,4)
errorbar(ax,mean(norm_accg,1,'omitmissing'),std(norm_accg,[],1,'omitmissing')./sqrt(size(norm_accg,1)),'CapSize',0,'Color','k')
xlabel('lags (ms)')
ylabel('norm. firing rate (from peak)')
set(gca,'TickDir','out'); box off

subplot(3,3,5)
plot(ax,squeeze(mean(norm_all_acc3d,1,'omitmissing')),'LineWidth',1.5)
xlabel('lags (ms)')
ylabel('norm. firing rate (from peak)')
set(gca,'TickDir','out'); box off

subplot(3,3,6)
pcolor(ax,1:10,squeeze(mean(norm_all_acc3d,1,'omitmissing')))
shading flat
colormap(flipud(slanCM('rdbu',265)))
xlabel('lags (ms)')
ylabel('firing rate decile')
set(gca,'TickDir','out'); box off

violation_fraction(find(violation_fraction==0))=1e-5;
Violation_rate(find(Violation_rate==0))=1e-5;
subplot(3,3,7)
h=histogram(log10(violation_fraction),'BinWidth',.1,'Normalization','probability','FaceColor','k','EdgeColor','none')
% bar(h.BinEdges(1:end-1),cumsum(h.Values),'FaceColor','k')
a=axis;
hold on
plot([log10(.1) log10(.1)],a(3:4),'--r')
hold off
ylabel('fraction')
xlabel('log(Violation fraction)')
set(gca,'TickDir','out'); box off

subplot(3,3,8)
h=histogram(log10(Violation_rate),'BinWidth',.2,'Normalization','probability','FaceColor','k','EdgeColor','none')
% bar(h.BinEdges(1:end-1),cumsum(h.Values),'FaceColor','k')
a=axis;
hold on
plot([log10(.2) log10(.2)],a(3:4),'--r')
hold off
ylabel('fraction')
xlabel('log(Violation rate) ~ contanimation')
set(gca,'TickDir','out'); box off

subplot(3,3,9)
scatter(all_mFR,log10(Violation_rate),50,'ok','fill')
xlabel('average firing rate (sp/s)')
ylabel('log(contanimation)')
set(gca,'TickDir','out'); box off
%% only pos/neg correlated cells
units2take=find(all_corr<2);
subplot(3,3,1)
errorbar(ax,mean(all_accg(units2take,:),1,'omitmissing'),std(all_accg(units2take,:),[],1,'omitmissing')./sqrt(size(all_accg(units2take,:),1)),'CapSize',0,'Color','k')
xlabel('lags (ms)')
ylabel('firing rate (sp/s)')
set(gca,'TickDir','out'); box off

subplot(3,3,2)
plot(ax,squeeze(mean(all_acc3d(units2take,:,:),1,'omitmissing')),'LineWidth',1.5)
xlabel('lags (ms)')
ylabel('firing rate (sp/s)')
set(gca,'TickDir','out'); box off

subplot(3,3,3)
imagesc(ax,1:10,squeeze(mean(all_acc3d(units2take,:,:),1,'omitmissing')))
set(gca,'YDir','normal');
shading flat
colormap(flipud(slanCM('rdbu',265)))
xlabel('lags (ms)')
ylabel('firing rate decile')
clim([0 100])
set(gca,'TickDir','out'); box off

subplot(3,3,4)
errorbar(ax,mean(norm_accg(units2take,:),1,'omitmissing'),std(norm_accg(units2take,:),[],1,'omitmissing')./sqrt(size(norm_accg(units2take,:),1)),'CapSize',0,'Color','k')
xlabel('lags (ms)')
ylabel('norm. firing rate (from peak)')
set(gca,'TickDir','out'); box off

subplot(3,3,5)
plot(ax,squeeze(mean(norm_all_acc3d(units2take,:,:),1,'omitmissing')),'LineWidth',1.5)
xlabel('lags (ms)')
ylabel('norm. firing rate (from peak)')
set(gca,'TickDir','out'); box off

subplot(3,3,6)
imagesc(ax,1:10,squeeze(mean(norm_all_acc3d(units2take,:,:),1,'omitmissing')))
set(gca,'YDir','normal');
shading flat
colormap(flipud(slanCM('rdbu',265)))
xlabel('lags (ms)')
ylabel('firing rate decile')
set(gca,'TickDir','out'); box off

violation_fraction(find(violation_fraction==0))=1e-5;
Violation_rate(find(Violation_rate==0))=1e-5;
subplot(3,3,7)
h=histogram(log10(violation_fraction(units2take)),'BinWidth',.1,'Normalization','probability','FaceColor','k','EdgeColor','none')
% bar(h.BinEdges(1:end-1),cumsum(h.Values),'FaceColor','k')
a=axis;
hold on
plot([log10(.1) log10(.1)],a(3:4),'--r')
hold off
ylabel('fraction')
xlabel('log(Violation fraction)')
set(gca,'TickDir','out'); box off

subplot(3,3,8)
h=histogram(log10(Violation_rate(units2take)),'BinWidth',.2,'Normalization','probability','FaceColor','k','EdgeColor','none')
% bar(h.BinEdges(1:end-1),cumsum(h.Values),'FaceColor','k')
a=axis;
hold on
plot([log10(.2) log10(.2)],a(3:4),'--r')
hold off
ylabel('fraction')
xlabel('log(Violation rate) ~ contanimation')
set(gca,'TickDir','out'); box off

subplot(3,3,9)
scatter(all_mFR(units2take),log10(Violation_rate(units2take)),50,'ok','fill')
xlabel('average firing rate (sp/s)')
ylabel('log(contanimation)')
set(gca,'TickDir','out'); box off

%% only pos/neg correlated cells with CDF
units2take=find(all_corr<2);
subplot(3,3,1)
errorbar(ax,mean(all_accg(units2take,:),1,'omitmissing'),std(all_accg(units2take,:),[],1,'omitmissing')./sqrt(size(all_accg(units2take,:),1)),'CapSize',0,'Color','k')
hold on
a=axis;
patch([-1.5 1.5 1.5 -1.5],[0 0 a(4) a(4)],[0.8 0.80 0.8],'EdgeColor','none','FaceAlpha',.5)
hold off
xlabel('lags (ms)')
ylabel('firing probability')
set(gca,'TickDir','out'); box off
xlim([-15 15])

subplot(3,3,2)
plot(ax,squeeze(mean(all_acc3d(units2take,:,:),1,'omitmissing')),'LineWidth',1.5)
hold on
a=axis;
patch([-1.5 1.5 1.5 -1.5],[0 0 a(4) a(4)],[0.8 0.80 0.8],'EdgeColor','none','FaceAlpha',.5)
hold off
xlabel('lags (ms)')
ylabel('firing probability')
set(gca,'TickDir','out'); box off
xlim([-15 15])

subplot(3,3,3)
imagesc(ax,1:10,squeeze(mean(all_acc3d(units2take,:,:),1,'omitmissing')))
set(gca,'YDir','normal');
shading flat
colormap(flipud(slanCM('rdbu',265)))
xlabel('lags (ms)')
ylabel('firing rate decile')
clim([0 .1])
set(gca,'TickDir','out'); box off
xlim([-15 15])

subplot(3,3,4)
errorbar(ax,mean(norm_accg(units2take,:),1,'omitmissing'),std(norm_accg(units2take,:),[],1,'omitmissing')./sqrt(size(norm_accg(units2take,:),1)),'CapSize',0,'Color','k')
xlabel('lags (ms)')
ylabel('norm. firing rate (from peak)')
set(gca,'TickDir','out'); box off
hold on
a=axis;
patch([-1.5 1.5 1.5 -1.5],[0 0 a(4) a(4)],[0.8 0.80 0.8],'EdgeColor','none','FaceAlpha',.5)
hold off
xlim([-15 15])

subplot(3,3,5)
plot(ax,squeeze(mean(norm_all_acc3d(units2take,:,:),1,'omitmissing')),'LineWidth',1.5)
xlabel('lags (ms)')
ylabel('norm. firing rate (from peak)')
set(gca,'TickDir','out'); box off
hold on
a=axis;
patch([-1.5 1.5 1.5 -1.5],[0 0 a(4) a(4)],[0.8 0.80 0.8],'EdgeColor','none','FaceAlpha',.5)
hold off
xlim([-15 15])

subplot(3,3,6)
imagesc(ax,1:10,squeeze(mean(norm_all_acc3d(units2take,:,:),1,'omitmissing')))
set(gca,'YDir','normal');
shading flat
colormap(flipud(slanCM('rdbu',265)))
xlabel('lags (ms)')
ylabel('firing rate decile')
set(gca,'TickDir','out'); box off
xlim([-15 15])

violation_fraction(find(violation_fraction==0))=1e-5;
Violation_rate(find(Violation_rate==0))=1e-5;
subplot(3,3,7)
h=histogram(log10(violation_fraction(units2take)),'BinWidth',.1,'Normalization','probability','FaceColor','k','EdgeColor','none')
bar(h.BinEdges(1:end-1),cumsum(h.Values),'FaceColor','k')
a=axis;
hold on
plot([log10(.1) log10(.1)],a(3:4),'--r')
hold off
ylabel('fraction')
xlabel('log(Violation fraction)')
set(gca,'TickDir','out'); box off

subplot(3,3,8)
h=histogram(log10(Violation_rate(units2take)),'BinWidth',.2,'Normalization','probability','FaceColor','k','EdgeColor','none')
bar(h.BinEdges(1:end-1),cumsum(h.Values),'FaceColor','k')
a=axis;
hold on
plot([log10(.2) log10(.2)],a(3:4),'--r')
hold off
ylabel('fraction')
xlabel('log(Violation rate) ~ contanimation')
set(gca,'TickDir','out'); box off

subplot(3,3,9)
scatter(all_mFR(units2take),log10(Violation_rate(units2take)),50,'ok','fill')
xlabel('average firing rate (sp/s)')
ylabel('log(contanimation)')
set(gca,'TickDir','out'); box off