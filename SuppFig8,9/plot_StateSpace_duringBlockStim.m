clear
fs= dir('data4StateSpace\*_CbN.mat');
c=0;
gwin=gausswin(151);
gwin=gwin./sum(gwin);
gwin=gwin';
gwin2=gausswin(5);
gwin2=gwin2./sum(gwin2);
gwin2=gwin2';
for i=1:length(fs)
    load(['data4StateSpace\' fs(i).name]);
    for ii=1:length(cellData)
        c=c+1;
        ctrl_fr=conv2(cellData(ii).Ctrl_max_vx.PSTH_smoothbin,gwin2,'same');
        t=cellData(ii).Ctrl_max_vx.time_bins; t=t(2:end-1); t=t./1000;
        [~,tmp]= max(ctrl_fr);
        if t(tmp)<=0 & t(tmp)>=-.25
            timing(c)=1; %change this if you want to compare between early/late units
        elseif t(tmp)<=.25 & t(tmp)>=0
            timing(c)=1; %change this if you want to compare between early/late units
        else
            timing(c)=2; %change this if you want to compare between early/late units
        end
        noCSM_posThresh(c,:)=conv2(cellData(ii).ADNoCSM_pos_thresh.PSTH_raw,gwin,'same');
        CSM_posThresh(c,:)=conv2(cellData(ii).ADCSM_pos_thresh.PSTH_raw,gwin,'same');


        AD_CSM_CSM0(c,:)=conv2(cellData(ii).ADCSM_csmT0.PSTH_raw,gwin,'same'); % stimulation block around T0 data
        AD_CSM_CSM2(c,:)=conv2(cellData(ii).ADCSM_csmT2.PSTH_raw,gwin,'same');% stimulation block around T2 data

        BL_CSM_CSM0(c,:)=conv2(cellData(ii).BLCSM_csmT0.PSTH_raw,gwin,'same');% baseline block around T0 data
        BL_CSM_CSM2(c,:)=conv2(cellData(ii).BLCSM_csmT2.PSTH_raw,gwin,'same');% baseline block around T2 data

        WO_CSM_CSM0(c,:)=conv2(cellData(ii).WOECSM_csmT0.PSTH_raw,gwin,'same'); % early washout block around T0 data
        WO_CSM_CSM2(c,:)=conv2(cellData(ii).WOECSM_csmT2.PSTH_raw,gwin,'same'); % early washout block around T2 data


        AD_CSM_thresh(c,:)=conv2(cellData(ii).ADCSM_pos_thresh.PSTH_raw,gwin,'same');% stimulation block around positional thresh data
        BL_CSM_thresh(c,:)=conv2(cellData(ii).BLCSM_pos_thresh.PSTH_raw,gwin,'same');% baseline block around positional thresh data
        WO_CSM_thresh(c,:)=conv2(cellData(ii).WOECSM_pos_thresh.PSTH_raw,gwin,'same');% early washout block around positional thresh data

        % the next part is to find the relative tming of the differnt
        % behavioral landmarks to each other, for the figure
        for tt=1:size(cellData(1).BLCSM_pos_thresh,1)
            [x,lags] = xcorr(cellData(1).BLCSM_pos_thresh.Raster(tt,:),cellData(1).BLCSM_csmT0.Raster(tt,:));
            [~,ix]= max(x);
            tmp1(tt)=lags(ix);
            [x,lags] = xcorr(cellData(1).BLCSM_pos_thresh.Raster(tt,:),cellData(1).BLCSM_csmT2.Raster(tt,:));
            [~,ix]= max(x);
            tmp2(tt)=lags(ix);
            [x,lags] = xcorr(cellData(1).BLCSM_csmT0.Raster(tt,:),cellData(1).BLCSM_csmT2.Raster(tt,:));
            [~,ix]= max(x);
            tmp3(tt)=lags(ix);
        end
        time_thresh2CSMT0(c)=median(tmp1);
        time_thresh2CSMT2(c)=median(tmp2);
        time_CSM02CSMT2(c)=median(tmp3);

    end
end
%% z score- signal processing
for i=1:size(noCSM_posThresh,1)

    tmp=BL_CSM_CSM0(i,:);
    tmp=tmp-mean(tmp(find(ax<-.5)));
    tmp=tmp./std(tmp(find(ax<-.5)));
    ZBL_CSM_CSM0(i,:)=tmp;
    V(i,:)=~isnan(nanmean(tmp)) & isempty(find(tmp==inf));;

    tmp=BL_CSM_CSM2(i,:);
    tmp=tmp-mean(tmp(find(ax<-.5)));
    tmp=tmp./std(tmp(find(ax<-.5)));
    ZBL_CSM_CSM2(i,:)=tmp;
    v=~isnan(nanmean(tmp)) & isempty(find(tmp==inf));
    V(i,:)=V(i,:)*v;

    tmp=AD_CSM_CSM0(i,:);
    tmp=tmp-mean(tmp(find(ax<-.5)));
    tmp=tmp./std(tmp(find(ax<-.5)));
    ZAD_CSM_CSM0(i,:)=tmp;
    v=~isnan(nanmean(tmp)) & isempty(find(tmp==inf));
    V(i,:)=V(i,:)*v;

    tmp=AD_CSM_CSM2(i,:);
    tmp=tmp-mean(tmp(find(ax<-.5)));
    tmp=tmp./std(tmp(find(ax<-.5)));
    ZAD_CSM_CSM2(i,:)=tmp;
    v=~isnan(nanmean(tmp)) & isempty(find(tmp==inf));
    V(i,:)=V(i,:)*v;

    tmp=WO_CSM_CSM0(i,:);
    tmp=tmp-mean(tmp(find(ax<-.5)));
    tmp=tmp./std(tmp(find(ax<-.5)));
    ZWO_CSM_CSM0(i,:)=tmp;
    v=~isnan(nanmean(tmp)) & isempty(find(tmp==inf));
    V(i,:)=V(i,:)*v;

    tmp=WO_CSM_CSM2(i,:);
    tmp=tmp-mean(tmp(find(ax<-.5)));
    tmp=tmp./std(tmp(find(ax<-.5)));
    ZWO_CSM_CSM2(i,:)=tmp;
    v=~isnan(nanmean(tmp)) & isempty(find(tmp==inf));
    V(i,:)=V(i,:)*v;

    tmp=BL_CSM_thresh(i,:);
    tmp=tmp-mean(tmp(find(ax<-.5)));
    tmp=tmp./std(tmp(find(ax<-.5)));
    ZBL_CSM_thresh(i,:)=tmp;
    v=~isnan(nanmean(tmp)) & isempty(find(tmp==inf));
    V(i,:)=V(i,:)*v;

    tmp=AD_CSM_thresh(i,:);
    tmp=tmp-mean(tmp(find(ax<-.5)));
    tmp=tmp./std(tmp(find(ax<-.5)));
    ZAD_CSM_thresh(i,:)=tmp;
    v=~isnan(nanmean(tmp)) & isempty(find(tmp==inf));
    V(i,:)=V(i,:)*v;

    tmp=WO_CSM_thresh(i,:);
    tmp=tmp-mean(tmp(find(ax<-.5)));
    tmp=tmp./std(tmp(find(ax<-.5)));
    ZWO_CSM_thresh(i,:)=tmp;
        v=~isnan(nanmean(tmp)) & isempty(find(tmp==inf));
    V(i,:)=V(i,:)*v;
end

%% remove any nan/inf signals= this might happen if a unit have 0 std at a specific condition, which might be a result of insufficient data
ix= find(V);
ZBL_CSM_CSM2=ZBL_CSM_CSM2(ix,:);
ZAD_CSM_CSM2=ZAD_CSM_CSM2(ix,:);
ZWO_CSM_CSM2=ZWO_CSM_CSM2(ix,:);
ix1=ix;

ZBL_CSM_CSM0=ZBL_CSM_CSM0(ix,:);
ZAD_CSM_CSM0=ZAD_CSM_CSM0(ix,:);
ZWO_CSM_CSM0=ZWO_CSM_CSM0(ix,:);
ix2=ix;

ZBL_CSM_thresh=ZBL_CSM_thresh(ix,:);
ZAD_CSM_thresh=ZAD_CSM_thresh(ix,:);
ZWO_CSM_thresh=ZWO_CSM_thresh(ix,:);
ix3=ix;

%% plot data
figure
iterations=5000;
k =1:3;
tix=find(ax>-.25 & ax<.5);
new_ax= ax(tix);
ix=find(new_ax==0);
index= find(timing(ix2)==1);
%concatenate the  signals
FR2PC= [ZBL_CSM_thresh(index,tix) ZAD_CSM_thresh(index,tix) ZWO_CSM_thresh(index,tix)];
data1=ZBL_CSM_thresh(index,tix);
data2=ZAD_CSM_thresh(index,tix);
data3=ZWO_CSM_thresh(index,tix);

%PCA
[coeff,score,latent] = pca(FR2PC');

%plot
subplot(2,2,1)
tmp= linspace(0,1,size(FR2PC,2)/3);
C_BL= [tmp(:) tmp(:) tmp(:)];
C_AD= [zeros(length(tmp),1) zeros(length(tmp),1) tmp(:)];
C_WO= [tmp(:) zeros(length(tmp),1) tmp(:)];
l=size(FR2PC,2)/3;
scatter3(score(1:l,1),score(1:l,2),score(1:l,3),20,C_BL)
hold on
scatter3(score(1+l:2*l,1),score(1+l:2*l,2),score(1+l:2*l,3),20,C_AD)
scatter3(score(2*l+1:end,1),score(2*l+1:end,2),score(2*l+1:end,3),20,C_WO)

scatter3(score(ix,1),score(ix,2),score(ix,3),80,'or','fill')
scatter3(score(l+ix,1),score(l+ix,2),score(l+ix,3),80,'or','fill')
scatter3(score(2*l+ix,1),score(2*l+ix,2),score(2*l+ix,3),80,'or','fill')

scatter3(score(ix+median(time_thresh2CSMT0),1),score(ix+median(time_thresh2CSMT0),2),score(ix+median(time_thresh2CSMT0),3),80,'sr','fill')
scatter3(score(l+ix+median(time_thresh2CSMT0),1),score(l+ix+median(time_thresh2CSMT0),2),score(l+ix+median(time_thresh2CSMT0),3),80,'sr','fill')
scatter3(score(2*l+ix+median(time_thresh2CSMT0),1),score(2*l+ix+median(time_thresh2CSMT0),2),score(2*l+ix+median(time_thresh2CSMT0),3),80,'sr','fill')

scatter3(score(ix+median(time_thresh2CSMT2),1),score(ix+median(time_thresh2CSMT2),2),score(ix+median(time_thresh2CSMT2),3),80,'^r','fill')
scatter3(score(l+ix+median(time_thresh2CSMT2),1),score(l+ix+median(time_thresh2CSMT2),2),score(l+ix+median(time_thresh2CSMT2),3),80,'^r','fill')
scatter3(score(2*l+ix+median(time_thresh2CSMT2),1),score(2*l+ix+median(time_thresh2CSMT2),2),score(2*l+ix+median(time_thresh2CSMT2),3),80,'^r','fill')
hold off
title(['Early; around thresh;' num2str(size(FR2PC,1))])
xlabel(['dim1; ' num2str(latent(1)/sum(latent))]); 
ylabel(['dim2; ' num2str(latent(2)/sum(latent))]);
zlabel(['dim3; ' num2str(latent(3)/sum(latent))]);


x1=data1'*coeff(:,k);
x2=data2'*coeff(:,k);
x3=data3'*coeff(:,k);

%calculate the inter-condition distances
DIST_ctrl2stim=sqrt(sum((x1'-x2').^2));
DIST_ctrl2wash=sqrt(sum((x1'-x3').^2));
DIST_stim2wash=sqrt(sum((x3'-x2').^2));
clear DIST_ctrl2stim_shuffle DIST_ctrl2wash_shuffle DIST_stim2wash_shuffle
for i=1:iterations
    X=[ZBL_CSM_thresh(index,tix); ZAD_CSM_thresh(index,tix); ZWO_CSM_thresh(index,tix)];
    [~,ix]= sort(rand(1,size(X,1)));
    n1=size(X,1)/3;
    X=X(ix,:);

    data1=X(1:n1,:);
    data2=X(1+n1:2*n1,:);
    data3=X(1+2*n1:3*n1,:);
    FR2PC=[data1 data2 data3];

    rix=1:size(coeff,1);
    x1=data1'*coeff(rix,k);
    x2=data2'*coeff(rix,k);
    x3=data3'*coeff(rix,k);

    DIST_ctrl2stim_shuffle(i,:)=sqrt(sum((x1'-x2').^2));
    DIST_ctrl2wash_shuffle(i,:)=sqrt(sum((x1'-x3').^2));
    DIST_stim2wash_shuffle(i,:)=sqrt(sum((x3'-x2').^2));

end


subplot(2,2,3)
plot(ax(tix),DIST_ctrl2stim,'b','LineWidth',2)
hold on
plot(ax(tix),DIST_ctrl2wash,'m','LineWidth',2)
plot(ax(tix),DIST_stim2wash,'r','LineWidth',2)
plot(ax(tix),nanmean(DIST_stim2wash_shuffle,1),'-k')
plot(ax(tix),nanmean(DIST_stim2wash_shuffle,1)-(nanstd(DIST_stim2wash_shuffle,1)./sqrt(iterations)).*1.96,'--k')
plot(ax(tix),nanmean(DIST_stim2wash_shuffle,1)+(nanstd(DIST_stim2wash_shuffle,1)./sqrt(iterations)).*1.96,'--k')
a=axis;
plot([median(time_thresh2CSMT2)/1000 median(time_thresh2CSMT2)/1000],a(3:4),'--k');
plot([0 0],a(3:4),'--k');
plot([median(time_thresh2CSMT0)/1000 median(time_thresh2CSMT0)/1000],a(3:4),'--k');

hold off
set(gca,'TickDir','out'); box off
set(gca,'TickDir','out'); box off
xlabel('time to thresh (s)')
ylabel('distance')


for i=1:size(DIST_ctrl2stim_shuffle,2)
    [~,P_ctrl2stim(i)] =ttest(DIST_ctrl2stim_shuffle(:,i),DIST_ctrl2stim(i));
    [~,P_ctrl2wash(i)] =ttest(DIST_ctrl2wash_shuffle(:,i),DIST_ctrl2wash(i));
    [~,P_stim2wash(i)] =ttest(DIST_stim2wash_shuffle(:,i),DIST_stim2wash(i));
end
%%
k =1:3;
tix=find(ax>-.25 & ax<.5);
new_ax= ax(tix);
ix=find(new_ax==0);
index= find(timing(ix2)==1);
FR2PC= [ZBL_CSM_CSM2(index,tix) ZAD_CSM_CSM2(index,tix) ZWO_CSM_CSM2(index,tix)];
[coeff,score,latent] = pca(FR2PC');
subplot(2,2,1)
tmp= linspace(0,1,size(FR2PC,2)/3);
C_BL= [tmp(:) tmp(:) tmp(:)];
C_AD= [zeros(length(tmp),1) zeros(length(tmp),1) tmp(:)];
C_WO= [tmp(:) zeros(length(tmp),1) tmp(:)];
l=size(FR2PC,2)/3;
scatter3(score(1:l,1),score(1:l,2),score(1:l,3),20,C_BL)
hold on
scatter3(score(1+l:2*l,1),score(1+l:2*l,2),score(1+l:2*l,3),20,C_AD)
scatter3(score(2*l+1:end,1),score(2*l+1:end,2),score(2*l+1:end,3),20,C_WO)
scatter3(score(ix,1),score(ix,2),score(ix,3),70,'^r','fill')
scatter3(score(l+ix,1),score(l+ix,2),score(l+ix,3),70,'^r','fill')
scatter3(score(2*l+ix,1),score(2*l+ix,2),score(2*l+ix,3),70,'^r','fill')

scatter3(score(ix-median(time_thresh2CSMT2),1),score(ix-median(time_thresh2CSMT2),2),score(ix-median(time_thresh2CSMT2),3),70,'or','fill')
scatter3(score(l+ix-median(time_thresh2CSMT2),1),score(l+ix-median(time_thresh2CSMT2),2),score(l+ix-median(time_thresh2CSMT2),3),70,'or','fill')
scatter3(score(2*l+ix-median(time_thresh2CSMT2),1),score(2*l+ix-median(time_thresh2CSMT2),2),score(2*l+ix-median(time_thresh2CSMT2),3),70,'or','fill')

scatter3(score(ix-median(time_CSM02CSMT2),1),score(ix-median(time_CSM02CSMT2),2),score(ix-median(time_CSM02CSMT2),3),70,'sr','fill')
scatter3(score(l+ix-median(time_CSM02CSMT2),1),score(l+ix-median(time_CSM02CSMT2),2),score(l+ix-median(time_CSM02CSMT2),3),70,'sr','fill')
scatter3(score(2*l+ix-median(time_CSM02CSMT2),1),score(2*l+ix-median(time_CSM02CSMT2),2),score(2*l+ix-median(time_CSM02CSMT2),3),70,'sr','fill')

hold off
title(['Early; around CSM2;' num2str(size(FR2PC,1))])
xlabel(['dim1; ' num2str(latent(1)/sum(latent))]); 
ylabel(['dim2; ' num2str(latent(2)/sum(latent))]);
zlabel(['dim3; ' num2str(latent(3)/sum(latent))]);


for i=1:3
    subplot(2,3,i+3)
    plot(ax(tix),score(1:l,i),'k')
    hold on
    plot(ax(tix),score(1+l:2*l,i),'b')
    plot(ax(tix),score(2*l+1:end,i),'m')
    hold off
    ylabel(['PC' num2str(i) '(a.u.)'])
end

data1=ZBL_CSM_CSM2(index,tix);
data2=ZAD_CSM_CSM2(index,tix);
data3=ZWO_CSM_CSM2(index,tix);
x1=data1'*coeff(:,k);
x2=data2'*coeff(:,k);
x3=data3'*coeff(:,k);

DIST_ctrl2stim=sqrt(sum((x1'-x2').^2));
DIST_ctrl2wash=sqrt(sum((x1'-x3').^2));
DIST_stim2wash=sqrt(sum((x3'-x2').^2));
clear DIST_ctrl2stim_shuffle DIST_ctrl2wash_shuffle DIST_stim2wash_shuffle
for i=1:iterations
    X=[ZBL_CSM_thresh(index,tix); ZAD_CSM_thresh(index,tix); ZWO_CSM_thresh(index,tix)];
    [~,ix]= sort(rand(1,size(X,1)));
    n1=size(X,1)/3;
    X=X(ix,:);

    data1=X(1:n1,:);
    data2=X(1+n1:2*n1,:);
    data3=X(1+2*n1:3*n1,:);
    FR2PC=[data1 data2 data3];

    rix=1:size(coeff,1);
    x1=data1'*coeff(rix,k);
    x2=data2'*coeff(rix,k);
    x3=data3'*coeff(rix,k);

    DIST_ctrl2stim_shuffle(i,:)=sqrt(sum((x1'-x2').^2));
    DIST_ctrl2wash_shuffle(i,:)=sqrt(sum((x1'-x3').^2));
    DIST_stim2wash_shuffle(i,:)=sqrt(sum((x3'-x2').^2));

end


subplot(2,2,2)
plot(ax(tix),DIST_ctrl2stim,'b','LineWidth',2)
hold on
plot(ax(tix),DIST_ctrl2wash,'m','LineWidth',2)
plot(ax(tix),DIST_stim2wash,'r','LineWidth',2)
plot(ax(tix),nanmean(DIST_stim2wash_shuffle,1),'-k')
plot(ax(tix),nanmean(DIST_stim2wash_shuffle,1)-(nanstd(DIST_stim2wash_shuffle,1)./sqrt(iterations)).*1.96,'--k')
plot(ax(tix),nanmean(DIST_stim2wash_shuffle,1)+(nanstd(DIST_stim2wash_shuffle,1)./sqrt(iterations)).*1.96,'--k')
a=axis;
plot([-median(time_thresh2CSMT2)/1000 -median(time_thresh2CSMT2)/1000],a(3:4),'--k');
plot([0 0],a(3:4),'--k');
plot([-median(time_CSM02CSMT2)/1000 -median(time_CSM02CSMT2)/1000],a(3:4),'--k');
hold off
set(gca,'TickDir','out'); box off
xlabel('time to CSM2 (s)')
ylabel('distance')



%%
figure
iterations=100;
k =1:3;
tix=find(ax>-.25 & ax<.5);
new_ax= ax(tix);
ix=find(new_ax==0);
index= find(timing(ix2)==1);
FR2PC= [ZBL_CSM_thresh(index,tix) ZAD_CSM_thresh(index,tix) ZWO_CSM_thresh(index,tix)];
data1=ZBL_CSM_thresh(index,tix);
data2=ZAD_CSM_thresh(index,tix);
data3=ZWO_CSM_thresh(index,tix);

[coeff,score,latent] = pca(FR2PC');

subplot(2,2,1)
tmp= linspace(0,1,size(FR2PC,2)/3);
C_BL= [tmp(:) tmp(:) tmp(:)];
C_AD= [zeros(length(tmp),1) zeros(length(tmp),1) tmp(:)];
C_WO= [tmp(:) zeros(length(tmp),1) tmp(:)];
l=size(FR2PC,2)/3;
scatter3(score(1:l,1),score(1:l,2),score(1:l,3),20,C_BL)
hold on
scatter3(score(1+l:2*l,1),score(1+l:2*l,2),score(1+l:2*l,3),20,C_AD)
scatter3(score(2*l+1:end,1),score(2*l+1:end,2),score(2*l+1:end,3),20,C_WO)

scatter3(score(ix,1),score(ix,2),score(ix,3),70,'or','fill')
scatter3(score(l+ix,1),score(l+ix,2),score(l+ix,3),70,'or','fill')
scatter3(score(2*l+ix,1),score(2*l+ix,2),score(2*l+ix,3),70,'or','fill')

scatter3(score(ix+median(time_thresh2CSMT0),1),score(ix+median(time_thresh2CSMT0),2),score(ix+median(time_thresh2CSMT0),3),70,'sr','fill')
scatter3(score(l+ix+median(time_thresh2CSMT0),1),score(l+ix+median(time_thresh2CSMT0),2),score(l+ix+median(time_thresh2CSMT0),3),70,'sr','fill')
scatter3(score(2*l+ix+median(time_thresh2CSMT0),1),score(2*l+ix+median(time_thresh2CSMT0),2),score(2*l+ix+median(time_thresh2CSMT0),3),70,'sr','fill')

scatter3(score(ix+median(time_thresh2CSMT2),1),score(ix+median(time_thresh2CSMT2),2),score(ix+median(time_thresh2CSMT2),3),70,'^r','fill')
scatter3(score(l+ix+median(time_thresh2CSMT2),1),score(l+ix+median(time_thresh2CSMT2),2),score(l+ix+median(time_thresh2CSMT2),3),70,'^r','fill')
scatter3(score(2*l+ix+median(time_thresh2CSMT2),1),score(2*l+ix+median(time_thresh2CSMT2),2),score(2*l+ix+median(time_thresh2CSMT2),3),70,'^r','fill')
hold off
title(['Early; around thresh;' num2str(size(FR2PC,1))])
xlabel(['dim1; ' num2str(latent(1)/sum(latent))]); 
ylabel(['dim2; ' num2str(latent(2)/sum(latent))]);
zlabel(['dim3; ' num2str(latent(3)/sum(latent))]);


x1=data1'*coeff(:,k);
x2=data2'*coeff(:,k);
x3=data3'*coeff(:,k);

DIST_ctrl2stim=sqrt(sum((x1'-x2').^2));
DIST_ctrl2wash=sqrt(sum((x1'-x3').^2));
DIST_stim2wash=sqrt(sum((x3'-x2').^2));
clear DIST_ctrl2stim_shuffle DIST_ctrl2wash_shuffle DIST_stim2wash_shuffle
for i=1:iterations
    X=[ZBL_CSM_thresh(index,tix); ZAD_CSM_thresh(index,tix); ZWO_CSM_thresh(index,tix)];
    [~,ix]= sort(rand(1,size(X,1)));
    n1=size(X,1)/3;
    X=X(ix,:);

    data1=X(1:n1,:);
    data2=X(1+n1:2*n1,:);
    data3=X(1+2*n1:3*n1,:);
    FR2PC=[data1 data2 data3];

    rix=1:size(coeff,1);
    x1=data1'*coeff(rix,k);
    x2=data2'*coeff(rix,k);
    x3=data3'*coeff(rix,k);

    DIST_ctrl2stim_shuffle(i,:)=sqrt(sum((x1'-x2').^2));
    DIST_ctrl2wash_shuffle(i,:)=sqrt(sum((x1'-x3').^2));
    DIST_stim2wash_shuffle(i,:)=sqrt(sum((x3'-x2').^2));

end


subplot(6,2,7)
plot(ax(tix),DIST_ctrl2stim,'b','LineWidth',2)
hold on
plot(ax(tix),nanmean(DIST_ctrl2stim_shuffle,1)-(nanstd(DIST_ctrl2stim_shuffle,1)./sqrt(iterations)).*1.96,'--k','LineWidth',1)
plot(ax(tix),nanmean(DIST_ctrl2stim_shuffle,1)+(nanstd(DIST_ctrl2stim_shuffle,1)./sqrt(iterations)).*1.96,'--k','LineWidth',1)
plot(ax(tix),nanmean(DIST_ctrl2stim_shuffle,1),'-k','LineWidth',2)
a=axis;
plot([median(time_thresh2CSMT2)/1000 median(time_thresh2CSMT2)/1000],a(3:4),'--k');
plot([0 0],a(3:4),'--k');
plot([median(time_thresh2CSMT0)/1000 median(time_thresh2CSMT0)/1000],a(3:4),'--k');
hold off
set(gca,'TickDir','out'); box off
set(gca,'TickDir','out'); box off
xlabel('time to thresh (s)')
ylabel('distance')
title('Baseline-stim')

subplot(6,2,9)
plot(ax(tix),DIST_ctrl2wash,'m','LineWidth',2)
hold on
plot(ax(tix),nanmean(DIST_ctrl2wash_shuffle,1)-(nanstd(DIST_ctrl2wash_shuffle,1)./sqrt(iterations)).*1.96,'--k','LineWidth',1)
plot(ax(tix),nanmean(DIST_ctrl2wash_shuffle,1)+(nanstd(DIST_ctrl2wash_shuffle,1)./sqrt(iterations)).*1.96,'--k','LineWidth',1)
plot(ax(tix),nanmean(DIST_ctrl2wash_shuffle,1),'-k','LineWidth',2)
a=axis;
plot([median(time_thresh2CSMT2)/1000 median(time_thresh2CSMT2)/1000],a(3:4),'--k');
plot([0 0],a(3:4),'--k');
plot([median(time_thresh2CSMT0)/1000 median(time_thresh2CSMT0)/1000],a(3:4),'--k');
hold off
set(gca,'TickDir','out'); box off
set(gca,'TickDir','out'); box off
xlabel('time to thresh (s)')
ylabel('distance')
title('Baseline-wash')

subplot(6,2,11)
plot(ax(tix),DIST_stim2wash,'r','LineWidth',2)
hold on
plot(ax(tix),nanmean(DIST_stim2wash_shuffle,1)-(nanstd(DIST_stim2wash_shuffle,1)./sqrt(iterations)).*1.96,'--k','LineWidth',1)
plot(ax(tix),nanmean(DIST_stim2wash_shuffle,1)+(nanstd(DIST_stim2wash_shuffle,1)./sqrt(iterations)).*1.96,'--k','LineWidth',1)
plot(ax(tix),nanmean(DIST_stim2wash_shuffle,1),'-k','LineWidth',2)
a=axis;
plot([median(time_thresh2CSMT2)/1000 median(time_thresh2CSMT2)/1000],a(3:4),'--k');
plot([0 0],a(3:4),'--k');
plot([median(time_thresh2CSMT0)/1000 median(time_thresh2CSMT0)/1000],a(3:4),'--k');
hold off
set(gca,'TickDir','out'); box off
set(gca,'TickDir','out'); box off
xlabel('time to thresh (s)')
ylabel('distance')
title('stim-wash')
%%
k =1:3;
tix=find(ax>-.25 & ax<.5);
new_ax= ax(tix);
ix=find(new_ax==0);
index= find(timing(ix2)==1);
FR2PC= [ZBL_CSM_CSM2(index,tix) ZAD_CSM_CSM2(index,tix) ZWO_CSM_CSM2(index,tix)];
[coeff,score,latent] = pca(FR2PC');
subplot(2,2,2)
tmp= linspace(0,1,size(FR2PC,2)/3);
C_BL= [tmp(:) tmp(:) tmp(:)];
C_AD= [zeros(length(tmp),1) zeros(length(tmp),1) tmp(:)];
C_WO= [tmp(:) zeros(length(tmp),1) tmp(:)];
l=size(FR2PC,2)/3;
scatter3(score(1:l,1),score(1:l,2),score(1:l,3),20,C_BL)
hold on
scatter3(score(1+l:2*l,1),score(1+l:2*l,2),score(1+l:2*l,3),20,C_AD)
scatter3(score(2*l+1:end,1),score(2*l+1:end,2),score(2*l+1:end,3),20,C_WO)
scatter3(score(ix,1),score(ix,2),score(ix,3),70,'^r','fill')
scatter3(score(l+ix,1),score(l+ix,2),score(l+ix,3),70,'^r','fill')
scatter3(score(2*l+ix,1),score(2*l+ix,2),score(2*l+ix,3),70,'^r','fill')

scatter3(score(ix-median(time_thresh2CSMT2),1),score(ix-median(time_thresh2CSMT2),2),score(ix-median(time_thresh2CSMT2),3),70,'or','fill')
scatter3(score(l+ix-median(time_thresh2CSMT2),1),score(l+ix-median(time_thresh2CSMT2),2),score(l+ix-median(time_thresh2CSMT2),3),70,'or','fill')
scatter3(score(2*l+ix-median(time_thresh2CSMT2),1),score(2*l+ix-median(time_thresh2CSMT2),2),score(2*l+ix-median(time_thresh2CSMT2),3),70,'or','fill')

scatter3(score(ix-median(time_CSM02CSMT2),1),score(ix-median(time_CSM02CSMT2),2),score(ix-median(time_CSM02CSMT2),3),70,'sr','fill')
scatter3(score(l+ix-median(time_CSM02CSMT2),1),score(l+ix-median(time_CSM02CSMT2),2),score(l+ix-median(time_CSM02CSMT2),3),70,'sr','fill')
scatter3(score(2*l+ix-median(time_CSM02CSMT2),1),score(2*l+ix-median(time_CSM02CSMT2),2),score(2*l+ix-median(time_CSM02CSMT2),3),70,'sr','fill')

hold off
title(['Early; around CSM2;' num2str(size(FR2PC,1))])
xlabel(['dim1; ' num2str(latent(1)/sum(latent))]); 
ylabel(['dim2; ' num2str(latent(2)/sum(latent))]);
zlabel(['dim3; ' num2str(latent(3)/sum(latent))]);

data1=ZBL_CSM_CSM2(index,tix);
data2=ZAD_CSM_CSM2(index,tix);
data3=ZWO_CSM_CSM2(index,tix);
x1=data1'*coeff(:,k);
x2=data2'*coeff(:,k);
x3=data3'*coeff(:,k);

DIST_ctrl2stim=sqrt(sum((x1'-x2').^2));
DIST_ctrl2wash=sqrt(sum((x1'-x3').^2));
DIST_stim2wash=sqrt(sum((x3'-x2').^2));
clear DIST_ctrl2stim_shuffle DIST_ctrl2wash_shuffle DIST_stim2wash_shuffle
for i=1:iterations
    X=[ZBL_CSM_thresh(index,tix); ZAD_CSM_thresh(index,tix); ZWO_CSM_thresh(index,tix)];
    [~,ix]= sort(rand(1,size(X,1)));
    n1=size(X,1)/3;
    X=X(ix,:);

    data1=X(1:n1,:);
    data2=X(1+n1:2*n1,:);
    data3=X(1+2*n1:3*n1,:);
    FR2PC=[data1 data2 data3];

    rix=1:size(coeff,1);
    x1=data1'*coeff(rix,k);
    x2=data2'*coeff(rix,k);
    x3=data3'*coeff(rix,k);

    DIST_ctrl2stim_shuffle(i,:)=sqrt(sum((x1'-x2').^2));
    DIST_ctrl2wash_shuffle(i,:)=sqrt(sum((x1'-x3').^2));
    DIST_stim2wash_shuffle(i,:)=sqrt(sum((x3'-x2').^2));

end


subplot(6,2,8)
plot(ax(tix),DIST_ctrl2stim,'b','LineWidth',2)
hold on
plot(ax(tix),nanmean(DIST_ctrl2stim_shuffle,1)-(nanstd(DIST_ctrl2stim_shuffle,1)./sqrt(iterations)).*1.96,'--k','LineWidth',1)
plot(ax(tix),nanmean(DIST_ctrl2stim_shuffle,1)+(nanstd(DIST_ctrl2stim_shuffle,1)./sqrt(iterations)).*1.96,'--k','LineWidth',1)
plot(ax(tix),nanmean(DIST_ctrl2stim_shuffle,1),'-k','LineWidth',2)
a=axis;
plot([-median(time_thresh2CSMT2)/1000 -median(time_thresh2CSMT2)/1000],a(3:4),'--k');
plot([0 0],a(3:4),'--k');
plot([-median(time_CSM02CSMT2)/1000 -median(time_CSM02CSMT2)/1000],a(3:4),'--k');
hold off
set(gca,'TickDir','out'); box off
set(gca,'TickDir','out'); box off
xlabel('time to CSM2 (s)')
ylabel('distance')
title('Baseline-stim')

subplot(6,2,10)
plot(ax(tix),DIST_ctrl2wash,'m','LineWidth',2)
hold on
plot(ax(tix),nanmean(DIST_ctrl2wash_shuffle,1)-(nanstd(DIST_ctrl2wash_shuffle,1)./sqrt(iterations)).*1.96,'--k','LineWidth',1)
plot(ax(tix),nanmean(DIST_ctrl2wash_shuffle,1)+(nanstd(DIST_ctrl2wash_shuffle,1)./sqrt(iterations)).*1.96,'--k','LineWidth',1)
plot(ax(tix),nanmean(DIST_ctrl2wash_shuffle,1),'-k','LineWidth',2)
a=axis;
plot([-median(time_thresh2CSMT2)/1000 -median(time_thresh2CSMT2)/1000],a(3:4),'--k');
plot([0 0],a(3:4),'--k');
plot([-median(time_CSM02CSMT2)/1000 -median(time_CSM02CSMT2)/1000],a(3:4),'--k');
hold off
set(gca,'TickDir','out'); box off
set(gca,'TickDir','out'); box off
xlabel('time to CSM2 (s)')
ylabel('distance')
title('Baseline-wash')

subplot(6,2,12)
plot(ax(tix),DIST_stim2wash,'r','LineWidth',2)
hold on
plot(ax(tix),nanmean(DIST_stim2wash_shuffle,1)-(nanstd(DIST_stim2wash_shuffle,1)./sqrt(iterations)).*1.96,'--k','LineWidth',1)
plot(ax(tix),nanmean(DIST_stim2wash_shuffle,1)+(nanstd(DIST_stim2wash_shuffle,1)./sqrt(iterations)).*1.96,'--k','LineWidth',1)
plot(ax(tix),nanmean(DIST_stim2wash_shuffle,1),'-k','LineWidth',2)
a=axis;
plot([-median(time_thresh2CSMT2)/1000 -median(time_thresh2CSMT2)/1000],a(3:4),'--k');
plot([0 0],a(3:4),'--k');
plot([-median(time_CSM02CSMT2)/1000 -median(time_CSM02CSMT2)/1000],a(3:4),'--k');
hold off
set(gca,'TickDir','out'); box off
set(gca,'TickDir','out'); box off
xlabel('time to CSM2 (s)')
ylabel('distance')
title('stim-wash')


