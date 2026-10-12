clear
fs= dir('..\data4StateSpace\*_CbN.mat');
c=0;
gwin=gausswin(15);
gwin=gwin./sum(gwin);
gwin=gwin';
gwin2=gausswin(3);
gwin2=gwin2./sum(gwin2);
gwin2=gwin2';
iteration=20;
win=10;
for i=1:length(fs)
    load(['..\..\data_around_thresh\' fs(i).name]);
    for ii=1:length(cellData)
        c=c+1;
        % ctrl_fr=conv2(cellData(ii).max_vx.PSTH_rawbin,gwin2,'same');
        ctrl_fr=cellData(ii).max_vx.PSTH_smoothbin;
        t=cellData(ii).max_vx.time_bins; t=t(2:end-1); t=t./1000;
        [~,tmp]= max(ctrl_fr);
        if t(tmp)<=0 & t(tmp)>=-.25
            timing(c)=1;
        elseif t(tmp)<=.25 & t(tmp)>=0
            timing(c)=1;
        else
            timing(c)=2;
        end
        if ismember(ii,neg_corr_cells)
            timing(c)=1;
        else
            timing(c)=2;
        end
        if cellData(ii).Corr==1
            timing(c)=1;
        elseif cellData(ii).Corr==-1
            timing(c)=2;
        else
            timing(c)=0;
        end

        noCSM_posThresh(c,:)=conv2(nanmean(bin_data(cellData(ii).NoCSM_pos_thresh.Raster,win),1).*(1000/win),gwin,'same');
        CSM_posThresh(c,:)=conv2(nanmean(bin_data(cellData(ii).AllCSM_pos_thresh.Raster,win),1).*(1000/win),gwin,'same');
        n1=size(cellData(ii).NoCSM_pos_thresh.Raster,1);
        n2=size(cellData(ii).AllCSM_pos_thresh.Raster,1);
        if n1<n2
            for iter=1:iteration
                [~,ix]= sort(rand(1,n2)); ix= ix(1:n1);
                TMP(iter,:)=nanmean(bin_data(cellData(ii).AllCSM_pos_thresh.Raster(ix,:),win),1).*(1000/win);
            end
            CSM_posThresh_matched(c,:)=conv2(nanmean(TMP,1),gwin,'same');
        else
            CSM_posThresh_matched(c,:)=conv2(nanmean(bin_data(cellData(ii).AllCSM_pos_thresh.raster,win),1).*(1000/win),gwin,'same');
        end

        for tt=1:size(cellData(1).AllCSM_csmT0.Raster,1)
            [x,lags] = xcorr(cellData(ii).AllCSM_pos_thresh.Raster(tt,:),cellData(ii).AllCSM_csmT0.Raster(tt,:));
            [~,ix]= max(x);
            tmp1(tt)=lags(ix);
            [x,lags] = xcorr(cellData(ii).AllCSM_pos_thresh.Raster(tt,:),cellData(ii).AllCSM_csmT2.Raster(tt,:));
            [~,ix]= max(x);
            tmp2(tt)=lags(ix);
            [x,lags] = xcorr(cellData(ii).AllCSM_csmT0.Raster(tt,:),cellData(ii).AllCSM_csmT2.Raster(tt,:));
            [~,ix]= max(x);
            tmp3(tt)=lags(ix);
        end
        time_thresh2CSMT0(c)=median(tmp1);
        time_thresh2CSMT2(c)=median(tmp2);
        time_CSM02CSMT2(c)=median(tmp3);
    end
end
%% z score
bix=-.5;
long_ax=ax;
%%
ax=ax(1:win:end);
for i=1:size(noCSM_posThresh,1)
    tmp=noCSM_posThresh(i,:);
    tmp=tmp-mean(tmp(find(ax<bix)));
    % tmp=tmp./range(tmp);
    tmp=tmp./std(tmp(find(ax<bix)));
    ZnoCSM_posThresh(i,:)=tmp;


    tmp=CSM_posThresh(i,:);
    tmp=tmp-mean(tmp(find(ax<bix)));
    % tmp=tmp./range(tmp);
    tmp=tmp./std(tmp(find(ax<bix)));
    ZCSM_posThresh(i,:)=tmp;


    tmp=CSM_posThresh_matched(i,:);
    tmp=tmp-mean(tmp(find(ax<bix)));
    % tmp=tmp./range(tmp);
    tmp=tmp./std(tmp(find(ax<bix)));
    ZCSM_posThresh_matched(i,:)=tmp;

end
%%
iterations=100;
tix=find(ax>-.5 & ax<.5);
new_ax= ax(tix);
ix=find(new_ax==0);
FR2PC= [ZnoCSM_posThresh(find(timing==1),tix) ZCSM_posThresh(find(timing==1),tix)];
[coeff,score,latent] = pca(FR2PC');
subplot(3,2,1)
tmp= linspace(0,1,size(FR2PC,2)/2);
C_noCSM= [tmp(:) tmp(:) tmp(:)];
C_CSM= [tmp(:) zeros(length(tmp),1) tmp(:)];
l=size(FR2PC,2)/2;
scatter3(score(1:l,1),score(1:l,2),score(1:l,3),20,C_noCSM,'filled')
hold on
scatter3(score(ix,1),score(ix,2),score(ix,3),40,'k')
scatter3(score(l,1),score(l,2),score(l,3),60,'^k','fill')
scatter3(score(l+1:end,1),score(l+1:end,2),score(l+1:end,3),20,C_CSM,'filled')
scatter3(score(ix+l+1,1),score(ix+l+1,2),score(ix+l+1,3),40,'m')
scatter3(score(end,1),score(end,2),score(end,3),60,'^m','fill')
hold off
title(['pos; around thresh;' num2str(length(find(timing==1)))])
xlabel('dim1'); ylabel('dim2'); zlabel('dim3');

x1=[score(1:l,1),score(1:l,2),score(1:l,3)];
x2=[score(1+l:end,1),score(1+l:end,2),score(1+l:end,3)];
Dist_real=sqrt(sum((x1-x2).^2,2));

n1=length(find(timing==1));
X= [ZnoCSM_posThresh(find(timing==1),:); ZCSM_posThresh(find(timing==1),:)];
clear Dist_shuffle
for i=1:iterations
    [~,ix1]=sort(rand(1,size(X,1))); ix1=ix1(1:length(find(timing==1)));
    [~,ix2]=sort(rand(1,size(X,1))); ix2=ix2(1:length(find(timing==1)));
    % FR2PC= [X(ix1,tix); X(ix2,tix)];
    data1=X(ix1,tix); data2=X(ix2,tix);
    [coeff,score,latent] = pca(FR2PC');
    % x1=[score(1:l,1),score(1:l,2),score(1:l,3)];
    % x2=[score(1+l:end,1),score(1+l:end,2),score(1+l:end,3)];

    loadings1=mean(coeff(1:n1,1:5));
    loadings2=mean(coeff(n1+1:end,1:5));
    x1=data1'*coeff(:,1:3);
    x2=data2'*coeff(:,1:3);
    % x1=[data.*loadings1(1),score(:,2).*loadings1(2),score(:,3).*loadings1(3)];
    % x2=[score(:,1).*loadings2(1),score(:,2).*loadings2(2),score(:,3).*loadings2(3)];

    Dist_shuffle(i,:)=sqrt(sum((x1-x2).^2,2));
end
%
% subplot(3,2,2)
% tmp= linspace(0,1,size(FR2PC,2)/2);
% C_noCSM= [tmp(:) tmp(:) tmp(:)];
% C_CSM= [tmp(:) zeros(length(tmp),1) tmp(:)];
% l=size(FR2PC,2)/2;
% scatter3(score(1:l,1),score(1:l,2),score(1:l,3),20,C_noCSM,'filled')
% hold on
% scatter3(score(ix,1),score(ix,2),score(ix,3),40,'k')
% scatter3(score(l,1),score(l,2),score(l,3),60,'^k','fill')
% scatter3(score(l+1:end,1),score(l+1:end,2),score(l+1:end,3),20,C_CSM,'filled')
% scatter3(score(ix+l+1,1),score(ix+l+1,2),score(ix+l+1,3),40,'m')
% scatter3(score(end,1),score(end,2),score(end,3),60,'^m','fill')
% hold off
% xlabel('dim1'); ylabel('dim2'); zlabel('dim3');
% title('shuffle')
% xlabel('dim1'); ylabel('dim2');
subplot(3,2,2)
plot(ax(tix),Dist_real,'k','LineWidth',2)
hold on
errorbar(ax(tix),mean(Dist_shuffle,1),std(Dist_shuffle,[],1),'r','LineWidth',0.5,'CapSize',0)
hold off
set(gca,'TickDir','out'); box off
xlabel('Time to thresh')
ylabel('distance')

%%
FR2PC= [ZnoCSM_posThresh(find(timing==1),tix) ZCSM_posThresh(find(timing==1),tix)];
[coeff,score,latent] = pca(FR2PC');
subplot(2,3,1)
tmp= linspace(0,1,size(FR2PC,2)/2);
C_noCSM= [tmp(:) tmp(:) tmp(:)];
C_CSM= [tmp(:) zeros(length(tmp),1) tmp(:)];
l=size(FR2PC,2)/2;
scatter3(score(1:l,1),score(1:l,2),score(1:l,3),20,C_noCSM,'filled')
hold on
scatter3(score(ix,1),score(ix,2),score(ix,3),40,'k')
scatter3(score(l,1),score(l,2),score(l,3),60,'^k','fill')
scatter3(score(l+1:end,1),score(l+1:end,2),score(l+1:end,3),20,C_CSM,'filled')
scatter3(score(ix+l+1,1),score(ix+l+1,2),score(ix+l+1,3),40,'m')
scatter3(score(end,1),score(end,2),score(end,3),60,'^m','fill')
hold off
title(['pos; around thresh;' num2str(length(find(timing==1)))])
xlabel('dim1'); ylabel('dim2'); zlabel('dim3');

x1=[score(1:l,1),score(1:l,2),score(1:l,3)];
x2=[score(1+l:end,1),score(1+l:end,2),score(1+l:end,3)];
Dist_real=sqrt(sum((x1-x2).^2,2));

n1=length(find(timing==1));
X= [ZnoCSM_posThresh(find(timing==1),:); ZCSM_posThresh(find(timing==1),:)];
clear Dist_shuffle
for i=1:iterations
    [~,ix1]=sort(rand(1,size(X,1))); ix1=ix1(1:length(find(timing==1)));
    [~,ix2]=sort(rand(1,size(X,1))); ix2=ix2(1:length(find(timing==1)));
    data1=X(ix1,tix); data2=X(ix2,tix);
    [coeff,score,latent] = pca(FR2PC');

    loadings1=mean(coeff(1:n1,1:5));
    loadings2=mean(coeff(n1+1:end,1:5));
    x1=data1'*coeff(:,1:3);
    x2=data2'*coeff(:,1:3);

    Dist_shuffle(i,:)=sqrt(sum((x1-x2).^2,2));
end
%

subplot(4,3,3)
plot(ax(tix),Dist_real,'k','LineWidth',2)
hold on
errorbar(ax(tix),mean(Dist_shuffle,1),std(Dist_shuffle,[],1),'r','LineWidth',0.5,'CapSize',0)
hold off
set(gca,'TickDir','out'); box off
xlabel('Time to thresh')
ylabel('distance')

%%
FR2PC= [ZnoCSM_posThresh(find(timing==2),tix) ZCSM_posThresh(find(timing==2),tix)];
[coeff,score,latent] = pca(FR2PC');
subplot(2,3,2)
tmp= linspace(0,1,size(FR2PC,2)/2);
C_noCSM= [tmp(:) tmp(:) tmp(:)];
C_CSM= [tmp(:) zeros(length(tmp),1) tmp(:)];
l=size(FR2PC,2)/2;
scatter3(score(1:l,1),score(1:l,2),score(1:l,3),20,C_noCSM,'filled')
hold on
scatter3(score(ix,1),score(ix,2),score(ix,3),40,'k')
scatter3(score(l,1),score(l,2),score(l,3),60,'^k','fill')
scatter3(score(l+1:end,1),score(l+1:end,2),score(l+1:end,3),20,C_CSM,'filled')
scatter3(score(ix+l+1,1),score(ix+l+1,2),score(ix+l+1,3),40,'m')
scatter3(score(end,1),score(end,2),score(end,3),60,'^m','fill')
hold off
title(['neg; around thresh;' num2str(length(find(timing==2)))])
xlabel('dim1'); ylabel('dim2'); zlabel('dim3');

x1=[score(1:l,1),score(1:l,2),score(1:l,3)];
x2=[score(1+l:end,1),score(1+l:end,2),score(1+l:end,3)];
Dist_real=sqrt(sum((x1-x2).^2,2));

n1=length(find(timing==2));
X= [ZnoCSM_posThresh(find(timing==2),:); ZCSM_posThresh(find(timing==2),:)];
clear Dist_shuffle
for i=1:iterations
    [~,ix1]=sort(rand(1,size(X,1))); ix1=ix1(1:length(find(timing==2)));
    [~,ix2]=sort(rand(1,size(X,1))); ix2=ix2(1:length(find(timing==2)));
    data1=X(ix1,tix); data2=X(ix2,tix);
    [coeff,score,latent] = pca(FR2PC');

    loadings1=mean(coeff(1:n1,1:5));
    loadings2=mean(coeff(n1+1:end,1:5));
    x1=data1'*coeff(:,1:3);
    x2=data2'*coeff(:,1:3);

    Dist_shuffle(i,:)=sqrt(sum((x1-x2).^2,2));
end
%

subplot(4,3,6)
plot(ax(tix),Dist_real,'k','LineWidth',2)
hold on
errorbar(ax(tix),mean(Dist_shuffle,1),std(Dist_shuffle,[],1),'r','LineWidth',0.5,'CapSize',0)
hold off
set(gca,'TickDir','out'); box off
xlabel('Time to thresh')
ylabel('distance')

%%
ix0= find(ax(tix)==0);
FR2PC= [ZCSM_posThresh(find(timing==1),tix); ZCSM_posThresh(find(timing==2),tix)];

data1=ZCSM_posThresh(find(timing==1),tix);
data2=ZCSM_posThresh(find(timing==2),tix);
l=length(find(find(timing==1)));

[coeff,score,latent] = pca(FR2PC');

x1=data1'*coeff(1:l,1:3);
x2=data2'*coeff(l+1:end,1:3);

tmp= linspace(0,1,size(FR2PC,2));
c_neg= [zeros(length(tmp),1)+0.3  tmp(:) zeros(length(tmp),1)+0.3];
c_pos= [zeros(length(tmp),1)+0.3 0.3+zeros(length(tmp),1) tmp(:)];


subplot(2,2,3)
scatter3(x1(:,1),x1(:,2),x1(:,3),20,c_pos)
hold on
scatter3(x2(:,1),x2(:,2),x2(:,3),20,c_neg)
scatter3(x1(ix0,1),x1(ix0,2),x1(ix0,3),50,'or','fill')
scatter3(x2(ix0,1),x2(ix0,2),x2(ix0,3),50,'or','fill')
hold off
xlabel('dim1'); ylabel('dim2'); zlabel('dim3')
set(gca,'TickDir','out'); box off



FR2PC= [ZnoCSM_posThresh(find(timing==1),tix) ZCSM_posThresh(find(timing==1),tix); ZnoCSM_posThresh(find(timing==2),tix) ZCSM_posThresh(find(timing==2),tix)];
data1=ZnoCSM_posThresh(find(timing==1),tix);
data11=ZCSM_posThresh(find(timing==1),tix);
data2=ZnoCSM_posThresh(find(timing==2),tix);
data21=ZCSM_posThresh(find(timing==2),tix);
l=length(find(find(timing==1)));

[coeff,score,latent] = pca(FR2PC');

x1=data1'*coeff(1:l,1:3);
x11=data11'*coeff(1:l,1:3);
x2=data2'*coeff(l+1:end,1:3);
x21=data21'*coeff(l+1:end,1:3);

subplot(2,2,4)
plot(ax(tix),sqrt(sum((x11'-x1').^2)),'b','LineWidth',2)
hold on
plot(ax(tix),sqrt(sum((x21'-x2').^2)),'g','LineWidth',2)
a=axis;
plot([0 0],a(3:4),'--k')
plot([median(time_thresh2CSMT0) median(time_thresh2CSMT0)]./1000,a(3:4),'--k')
plot([median(time_thresh2CSMT2) median(time_thresh2CSMT2)]./1000,a(3:4),'--k')
hold off
xlabel('Time to threshold crossing (s)')
ylabel('distance')
set(gca,'TickDir','out'); box off



%% plot for figure
iterations=5000;
tix= find(ax>=-.5 & ax<=.5);
FR2PC= [ZnoCSM_posThresh(find(timing>=1),tix) ZCSM_posThresh(find(timing>=1),tix)];
[coeff,score,latent] = pca(FR2PC');
subplot(2,2,1)
tmp= linspace(0,1,size(FR2PC,2)/2);
C_noCSM= [tmp(:) tmp(:) tmp(:)];
C_CSM= [tmp(:) zeros(length(tmp),1) tmp(:)];
l=size(FR2PC,2)/2;
plot3(score(1:l,1),score(1:l,2),score(1:l,3),'k','LineWidth',2)
hold on
scatter3(score(ix,1),score(ix,2),score(ix,3),40,'k')
scatter3(score(l,1),score(l,2),score(l,3),60,'^k','fill')
plot3(score(l+1:end,1),score(l+1:end,2),score(l+1:end,3),'m','LineWidth',2)
scatter3(score(ix+l,1),score(ix+l,2),score(ix+l,3),40,'m')
scatter3(score(end,1),score(end,2),score(end,3),60,'^m','fill')


tmp=median(time_thresh2CSMT0);
tmp=round(tmp/win);
scatter3(score(ix+tmp,1),score(ix+tmp,2),score(ix+tmp,3),40,'k')
scatter3(score(tmp+ix+l+1,1),score(tmp+ix+l+1,2),score(tmp+ix+l+1,3),40,'m')

tmp=median(time_thresh2CSMT2);
tmp=round(tmp/win);
scatter3(score(ix+tmp,1),score(ix+tmp,2),score(ix+tmp,3),40,'k')
scatter3(score(tmp+ix+l+1,1),score(tmp+ix+l+1,2),score(tmp+ix+l+1,3),40,'m')

hold off
title(['pos; around thresh;' num2str(length(find(timing>=1)))])
xlabel('dim1'); ylabel('dim2'); zlabel('dim3');

x1=[score(1:l,1),score(1:l,2),score(1:l,3)];
x2=[score(1+l:end,1),score(1+l:end,2),score(1+l:end,3)];
Dist_real=sqrt(sum((x1-x2).^2,2));

% angles=acos(sum(X.*Y,2)./(sqrt(sum(X.^2,2)).*sqrt(sum(Y.^2,2))));

n1=length(find(timing>=1));
X= [ZnoCSM_posThresh(find(timing>=1),:); ZCSM_posThresh(find(timing>=1),:)];
clear Dist_shuffle
for i=1:iterations
    [~,ix1]=sort(rand(1,size(X,1))); ix1=ix1(1:length(find(timing>=1)));
    [~,ix2]=sort(rand(1,size(X,1))); ix2=ix2(1:length(find(timing>=1)));
    data1=X(ix1,tix); data2=X(ix2,tix);
    [coeff,score,latent] = pca(FR2PC');

    loadings1=mean(coeff(1:n1,1:5));
    loadings2=mean(coeff(n1+1:end,1:5));
    x1=data1'*coeff(:,1:3);
    x2=data2'*coeff(:,1:3);

    Dist_shuffle(i,:)=sqrt(sum((x1-x2).^2,2));
end
%

subplot(2,2,2)
plot(ax(tix),Dist_real,'r','LineWidth',2)
hold on
plot(ax(tix),mean(Dist_shuffle,1),'k','LineWidth',2)
plot(ax(tix),mean(Dist_shuffle,1)-(std(Dist_shuffle,[],1)./sqrt(iterations)).*1.96,'k','LineWidth',0.5)
plot(ax(tix),mean(Dist_shuffle,1)+(std(Dist_shuffle,[],1)./sqrt(iterations)).*1.96,'k','LineWidth',0.5)
plot([0 0],a(3:4),'--k')
plot([median(time_thresh2CSMT0) median(time_thresh2CSMT0)]./1000,a(3:4),'--k')
plot([median(time_thresh2CSMT2) median(time_thresh2CSMT2)]./1000,a(3:4),'--k')
hold off
set(gca,'TickDir','out'); box off
xlabel('Time to thresh')
ylabel('distance')


ix0= find(ax(tix)==0);
FR2PC= [ZCSM_posThresh(find(timing==1),tix); ZCSM_posThresh(find(timing==2),tix)];

data1=ZCSM_posThresh(find(timing==1),tix);
data2=ZCSM_posThresh(find(timing==2),tix);
l=length(find(find(timing==1)));

[coeff,score,latent] = pca(FR2PC');

x1=data1'*coeff(1:l,1:3);
x2=data2'*coeff(l+1:end,1:3);

tmp= linspace(0,1,size(FR2PC,2));
c_neg= [zeros(length(tmp),1)+0.3  tmp(:) zeros(length(tmp),1)+0.3];
c_pos= [zeros(length(tmp),1)+0.3 0.3+zeros(length(tmp),1) tmp(:)];


subplot(2,2,3)
plot3(x1(:,1),x1(:,2),x1(:,3),'b','LineWidth',2)
hold on
plot3(x2(:,1),x2(:,2),x2(:,3),'g','LineWidth',2)
scatter3(x1(ix0,1),x1(ix0,2),x1(ix0,3),50,'or','fill')
scatter3(x2(ix0,1),x2(ix0,2),x2(ix0,3),50,'or','fill')

scatter3(x1(end,1),x1(end,2),x1(end,3),50,'^r','fill')
scatter3(x2(end,1),x2(end,2),x2(end,3),50,'^r','fill')


tmp=median(time_thresh2CSMT0);
tmp=round(tmp/win);
scatter3(x1(ix0+tmp,1),x1(ix0+tmp,2),x1(ix0+tmp,3),50,'or','fill')
scatter3(x2(ix0+tmp,1),x2(ix0+tmp,2),x2(ix0+tmp,3),50,'or','fill')

tmp=median(time_thresh2CSMT2);
tmp=round(tmp/win);
scatter3(x1(ix0+tmp,1),x1(ix0+tmp,2),x1(ix0+tmp,3),50,'or','fill')
scatter3(x2(ix0+tmp,1),x2(ix0+tmp,2),x2(ix0+tmp,3),50,'or','fill')

hold off
xlabel('dim1'); ylabel('dim2'); zlabel('dim3')
set(gca,'TickDir','out'); box off



FR2PC= [ZnoCSM_posThresh(find(timing==1),tix) ZCSM_posThresh(find(timing==1),tix); ZnoCSM_posThresh(find(timing==2),tix) ZCSM_posThresh(find(timing==2),tix)];
data1=ZnoCSM_posThresh(find(timing==1),tix);
data11=ZCSM_posThresh(find(timing==1),tix);
data2=ZnoCSM_posThresh(find(timing==2),tix);
data21=ZCSM_posThresh(find(timing==2),tix);
l=length(find(find(timing==1)));

[coeff,score,latent] = pca(FR2PC');

x1=data1'*coeff(1:l,1:3);
x11=data11'*coeff(1:l,1:3);
x2=data2'*coeff(l+1:end,1:3);
x21=data21'*coeff(l+1:end,1:3);

subplot(2,2,4)
plot(ax(tix),sqrt(sum((x11'-x1').^2)),'b','LineWidth',2)
hold on
plot(ax(tix),sqrt(sum((x21'-x2').^2)),'g','LineWidth',2)
a=axis;
plot([0 0],a(3:4),'--k')
plot([median(time_thresh2CSMT0) median(time_thresh2CSMT0)]./1000,a(3:4),'--k')
plot([median(time_thresh2CSMT2) median(time_thresh2CSMT2)]./1000,a(3:4),'--k')
hold off
xlabel('Time to threshold crossing (s)')
ylabel('distance')
set(gca,'TickDir','out'); box off

for i=1:size(Dist_shuffle,2)
    P(i) =length(find(Dist_shuffle(:,i)>Dist_real(i)));
end

%%
%% plot for figure matched trial
tix= find(ax>=-.5 & ax<=.5);
FR2PC= [ZnoCSM_posThresh(find(timing>=1),tix) ZCSM_posThresh_matched(find(timing>=1),tix)];
[coeff,score,latent] = pca(FR2PC');
subplot(2,2,1)
tmp= linspace(0,1,size(FR2PC,2)/2);
C_noCSM= [tmp(:) tmp(:) tmp(:)];
C_CSM= [tmp(:) zeros(length(tmp),1) tmp(:)];
l=size(FR2PC,2)/2;
plot3(score(1:l,1),score(1:l,2),score(1:l,3),'k','LineWidth',2)
hold on
scatter3(score(ix,1),score(ix,2),score(ix,3),40,'k')
scatter3(score(l,1),score(l,2),score(l,3),60,'^k','fill')
plot3(score(l+1:end,1),score(l+1:end,2),score(l+1:end,3),'m','LineWidth',2)
scatter3(score(ix+l,1),score(ix+l,2),score(ix+l,3),40,'m')
scatter3(score(end,1),score(end,2),score(end,3),60,'^m','fill')


tmp=median(time_thresh2CSMT0);
tmp=round(tmp/win);
scatter3(score(ix+tmp,1),score(ix+tmp,2),score(ix+tmp,3),40,'k')
scatter3(score(tmp+ix+l+1,1),score(tmp+ix+l+1,2),score(tmp+ix+l+1,3),40,'m')

tmp=median(time_thresh2CSMT2);
tmp=round(tmp/win);
scatter3(score(ix+tmp,1),score(ix+tmp,2),score(ix+tmp,3),40,'k')
scatter3(score(tmp+ix+l+1,1),score(tmp+ix+l+1,2),score(tmp+ix+l+1,3),40,'m')

hold off
title(['pos; around thresh;' num2str(length(find(timing>=1)))])
xlabel('dim1'); ylabel('dim2'); zlabel('dim3');

x1=[score(1:l,1),score(1:l,2),score(1:l,3)];
x2=[score(1+l:end,1),score(1+l:end,2),score(1+l:end,3)];
Dist_real=sqrt(sum((x1-x2).^2,2));
% angles=acos(sum(X.*Y,2)./(sqrt(sum(X.^2,2)).*sqrt(sum(Y.^2,2))));
n1=length(find(timing>=1));
X= [ZnoCSM_posThresh(find(timing>=1),:); ZCSM_posThresh_matched(find(timing>=1),:)];
clear Dist_shuffle
for i=1:iterations
    [~,ix1]=sort(rand(1,size(X,1))); ix1=ix1(1:length(find(timing>=1)));
    [~,ix2]=sort(rand(1,size(X,1))); ix2=ix2(1:length(find(timing>=1)));
    data1=X(ix1,tix); data2=X(ix2,tix);
    [coeff,score,latent] = pca(FR2PC');

    loadings1=mean(coeff(1:n1,1:5));
    loadings2=mean(coeff(n1+1:end,1:5));
    x1=data1'*coeff(:,1:3);
    x2=data2'*coeff(:,1:3);

    Dist_shuffle(i,:)=sqrt(sum((x1-x2).^2,2));
end
%

subplot(2,2,2)
plot(ax(tix),Dist_real,'r','LineWidth',2)
hold on
plot(ax(tix),mean(Dist_shuffle,1),'k','LineWidth',2)
plot(ax(tix),mean(Dist_shuffle,1)-(std(Dist_shuffle,[],1)./sqrt(iterations)).*3.291,'k','LineWidth',0.5)
plot(ax(tix),mean(Dist_shuffle,1)+(std(Dist_shuffle,[],1)./sqrt(iterations)).*3.291,'k','LineWidth',0.5)
plot([0 0],a(3:4),'--k')
plot([median(time_thresh2CSMT0) median(time_thresh2CSMT0)]./1000,a(3:4),'--k')
plot([median(time_thresh2CSMT2) median(time_thresh2CSMT2)]./1000,a(3:4),'--k')
hold off
set(gca,'TickDir','out'); box off
xlabel('Time to thresh')
ylabel('distance')


ix0= find(ax(tix)==0);
FR2PC= [ZCSM_posThresh_matched(find(timing==1),tix); ZCSM_posThresh_matched(find(timing==2),tix)];

data1=ZCSM_posThresh_matched(find(timing==1),tix);
data2=ZCSM_posThresh_matched(find(timing==2),tix);
l=length(find(find(timing==1)));

[coeff,score,latent] = pca(FR2PC');

x1=data1'*coeff(1:l,1:3);
x2=data2'*coeff(l+1:end,1:3);

tmp= linspace(0,1,size(FR2PC,2));
c_neg= [zeros(length(tmp),1)+0.3  tmp(:) zeros(length(tmp),1)+0.3];
c_pos= [zeros(length(tmp),1)+0.3 0.3+zeros(length(tmp),1) tmp(:)];


subplot(2,2,3)
plot3(x1(:,1),x1(:,2),x1(:,3),'b','LineWidth',2)
hold on
plot3(x2(:,1),x2(:,2),x2(:,3),'g','LineWidth',2)
scatter3(x1(ix0,1),x1(ix0,2),x1(ix0,3),50,'or','fill')
scatter3(x2(ix0,1),x2(ix0,2),x2(ix0,3),50,'or','fill')

scatter3(x1(end,1),x1(end,2),x1(end,3),50,'^r','fill')
scatter3(x2(end,1),x2(end,2),x2(end,3),50,'^r','fill')



tmp=median(time_thresh2CSMT0);
tmp=round(tmp/win);
scatter3(x1(ix0+tmp,1),x1(ix0+tmp,2),x1(ix0+tmp,3),50,'or','fill')
scatter3(x2(ix0+tmp,1),x2(ix0+tmp,2),x2(ix0+tmp,3),50,'or','fill')

tmp=median(time_thresh2CSMT2);
tmp=round(tmp/win);
scatter3(x1(ix0+tmp,1),x1(ix0+tmp,2),x1(ix0+tmp,3),50,'or','fill')
scatter3(x2(ix0+tmp,1),x2(ix0+tmp,2),x2(ix0+tmp,3),50,'or','fill')


hold off
xlabel('dim1'); ylabel('dim2'); zlabel('dim3')
set(gca,'TickDir','out'); box off



FR2PC= [ZnoCSM_posThresh(find(timing==1),tix) ZCSM_posThresh_matched(find(timing==1),tix); ZnoCSM_posThresh(find(timing==2),tix) ZCSM_posThresh_matched(find(timing==2),tix)];
data1=ZnoCSM_posThresh(find(timing==1),tix);
data11=ZCSM_posThresh_matched(find(timing==1),tix);
data2=ZnoCSM_posThresh(find(timing==2),tix);
data21=ZCSM_posThresh_matched(find(timing==2),tix);
l=length(find(find(timing==1)));

[coeff,score,latent] = pca(FR2PC');

x1=data1'*coeff(1:l,1:3);
x11=data11'*coeff(1:l,1:3);
x2=data2'*coeff(l+1:end,1:3);
x21=data21'*coeff(l+1:end,1:3);

subplot(2,2,4)
plot(ax(tix),sqrt(sum((x11'-x1').^2)),'b','LineWidth',2)
hold on
plot(ax(tix),sqrt(sum((x21'-x2').^2)),'g','LineWidth',2)
a=axis;
plot([0 0],a(3:4),'--k')
plot([median(time_thresh2CSMT0) median(time_thresh2CSMT0)]./1000,a(3:4),'--k')
plot([median(time_thresh2CSMT2) median(time_thresh2CSMT2)]./1000,a(3:4),'--k')
hold off
xlabel('Time to threshold crossing (s)')
ylabel('distance')
set(gca,'TickDir','out'); box off