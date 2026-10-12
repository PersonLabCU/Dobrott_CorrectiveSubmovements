%alleaches bellshape analysis

%stimvel
%aligned to peak v
reaches_stimvel = ExcitRR_Ctrl;


a_s = fieldnames(reaches_stimvel);

for i = 1:length(a_s)
    
    %if you want to filter out animals with less than n reaches
%     r_count(i) = length(fieldnames(reaches_stimvel.(a_s{i}).Reaches));
%     
%     if r_count(i) < 50
%         continue
%     end
    
    C = AlignPeakV(reaches_stimvel.(a_s{i}).CtrlReaches,7,200);
    
    %NOTE if using this , dim 5 becomes speed instead of vel
    C = AlignPeakSp(reaches_stimvel.(a_s{i}).CtrlReaches,7,200);

    C_a = stimvel_accel(C);
    
    reaches_stimvel.(a_s{i}).stimvel = C_a;
    
    reaches_stimvel.(a_s{i}).BKins = BasicKinematics(reaches_stimvel.(a_s{i}).CtrlReaches);
end


%%
%plotting
figure; hold on
bmap = brewermap(length(a_s),'Set3');

All_big = struct;
j = 1;
for i = 1: length(a_s)
    
    if isfield(reaches_stimvel.(a_s{i}),'stimvel')
        col = bmap(1,:);
%         col = [rand rand rand .4];
        StimvelPlot2(reaches_stimvel.(a_s{i}).stimvel,5,[col .3],0)
        
        All_big.stimvel(:,:,j) = reaches_stimvel.(a_s{i}).stimvel.Average;
        j = j+1;
    end
    drawnow
    
end

ylim auto
xlim auto

All_big.Average = nanmean(All_big.stimvel,3);
All_big.Std = nanstd(All_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(All_big.Average),1);
for nnum = 1:length(All_big.Average)
    samp = squeeze(All_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
All_big.N = N;

StimvelPlot2(All_big,5,col-.15,0,1)
ylim auto
xlim([-.2 .2])
% figset
%%
figure; hold on
All_big2 = struct;
j = 1;

for i = 1: length(a_s)
    
    if isfield(reaches_stimvel.(a_s{i}),'stimvel')
        %duplicate for normalization
        reaches_stimvel.(a_s{i}).stimvel_peaknorm = reaches_stimvel.(a_s{i}).stimvel;
        maxv = max(reaches_stimvel.(a_s{i}).stimvel_peaknorm.Average(:,5));
        reaches_stimvel.(a_s{i}).stimvel_peaknorm.Average(:,5) = reaches_stimvel.(a_s{i}).stimvel_peaknorm.Average(:,5) / maxv;
        
        col = bmap(1,:);
%         col = [rand rand rand .4];
        StimvelPlot2(reaches_stimvel.(a_s{i}).stimvel_peaknorm,5,[col .3],0)
        
        All_big2.stimvel(:,:,j) = reaches_stimvel.(a_s{i}).stimvel_peaknorm.Average;
        j = j+1;
    end
    drawnow
    
end

ylim auto
xlim ([-.2 .2])

All_big2.Average = nanmean(All_big2.stimvel,3);
All_big2.Std = nanstd(All_big2.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(All_big2.Average),1);
for nnum = 1:length(All_big2.Average)
    samp = squeeze(All_big2.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
All_big2.N = N;

StimvelPlot2(All_big2,5,[0 0 0],0,1)

ylim auto
xlim ([-.2 .2])




%%
%sort reaches based on distance. look at asymmetry.
dists = [];
asymmetry = [];

Sym.J1 = struct;
Sym.J2 = struct;
Sym.J3 = struct;
Sym.J4 = struct;

jn = fieldnames(Sym);

for i = 1:length(a_s)
    
    rdists = reaches_stimvel.(a_s{i}).BKins.dist;
    %do sorting
    [~,sort_inx] = sort(rdists);
    
    segs = 4;
    
    for j = 1:segs
        
        

        
        %get top and bottom of segments 
        min_inx = max([1 floor(((j-1)/segs)*length(sort_inx))]);
        max_inx = floor((j/segs)*length(sort_inx));
        
        test_inx = sort_inx(min_inx:max_inx);
        
        dists(i,j) = mean(rdists(test_inx));
        
        %get reaches
        temp_reaches = Extract(reaches_stimvel.(a_s{i}).Reaches,test_inx);
        
        %align to peakv
        C = AlignPeakSp(temp_reaches,7,200);
        
        %get asymmetry statistic
        firsthalf = trapz(0:.01:.2,C.Average(1:21,5));
        secondhalf = trapz(0:.01:.2,C.Average(21:41,5));
        asymmetry(i,j) = firsthalf / (firsthalf + secondhalf);
        
        
        %save stimvel for plotting
        Sym.(jn{j}).stimvel(:,:,i) = C.Average;

    end


end

cmap = brewermap(segs,'*RdBu');
%plotting
figure;
hold on
for s = 1:segs
    scatter(ones(length(a_s),1)*s,asymmetry(:,s),30,cmap(s,:),...
        'filled','jitter','on','jitterAmount',.1,...
        'MarkerFaceAlpha',.33)
    scatter(s,nanmean(asymmetry(:,s)),170,cmap(s,:),'filled')

end
xlim([.5 s+.5])
errorbar((1:segs),nanmean(asymmetry),nanstd(asymmetry),'LineStyle','none','Color','k','CapSize',14)
set(gcf,'color','w');


[h,p] = ttest(asymmetry(:,1),asymmetry(:,4))
nanmean(asymmetry(:,1))
nanmean(asymmetry(:,4))
std(asymmetry(:,1))
std(asymmetry(:,4))

%%
%plot stimvel curves
figure; hold on
for j = 1:length(jn)
    Sym.(jn{j}).Average = nanmean(Sym.(jn{j}).stimvel,3);
    Sym.(jn{j}).Std = nanstd(Sym.(jn{j}).stimvel,0,3);
    
    %loop through to count the number of reaches contributing to each time pt
    N = zeros(length(Sym.(jn{j}).Average),1);
    for nnum = 1:length(Sym.(jn{j}).Average)
        samp = squeeze(Sym.(jn{j}).stimvel(nnum,2,:));
        N(nnum) = length(samp(~isnan(samp)));
    end
    Sym.(jn{j}).N = N;
    
%     subplot(2,2,j)
    hold on
    StimvelPlot2(Sym.(jn{j}),5,cmap(j,:),0,0)
    xlim([-.1 .1])
    ylim([0 40])
end



%%
%schematic figure attempt
a = 15;
C = AlignPeakSp(reaches_stimvel.(a_s{i}).Reaches,7,200);

bmap = brewermap(3,'Set2');


figure; hold on
subplot(2,1,1); hold on
scatter3(C.Average(:,1),C.Average(:,2),C.Average(:,3),30,bmap(1,:),'filled')
xlim([-.1 .1])
subplot(2,1,2); hold on
scatter3(C.Average(:,1),C.Average(:,5),C.Average(:,6),30,bmap(1,:),'filled')
xlim([-.1 .1])




%%
%sort reaches based on distance. look at asymmetry in normalized profiles
dists = [];
asymmetry2 = [];

Sym.J1 = struct;
Sym.J2 = struct;
Sym.J3 = struct;
Sym.J4 = struct;

jn = fieldnames(Sym);

for i = 1:length(a_s)
    
    rdists = reaches_stimvel.(a_s{i}).BKins.dist;
    %do sorting
    [~,sort_inx] = sort(rdists);
    
    segs = 4;
    
    for j = 1:segs
        
        

        
        %get top and bottom of segments 
        min_inx = max([1 floor(((j-1)/segs)*length(sort_inx))]);
        max_inx = floor((j/segs)*length(sort_inx));
        
        test_inx = sort_inx(min_inx:max_inx);
        
        dists(i,j) = mean(rdists(test_inx));
        
        %get reaches
        temp_reaches = Extract(reaches_stimvel.(a_s{i}).Reaches,test_inx);
        
        %align to peakv
        C = AlignPeakSp(temp_reaches,7,200);
        
        %normalize to peak
        maxv = max(C.Average(:,5));
        C_peaknorm = C.Average(:,5) / maxv;
       
        
        %get asymmetry statistic
        firsthalf = trapz(0:.01:.2,C.Average(1:21,5));
        secondhalf = trapz(0:.01:.2,C.Average(21:41,5));
        asymmetry2(i,j) = firsthalf / (firsthalf + secondhalf);
        
        
        %save stimvel for plotting
        Sym.(jn{j}).stimvel(:,:,i) = C_peaknorm;

    end


end

cmap = brewermap(segs,'*RdBu');
%plotting
figure;
hold on
for s = 1:segs
    scatter(ones(length(a_s),1)*s,asymmetry2(:,s),30,cmap(s,:),...
        'filled','jitter','on','jitterAmount',.1,...
        'MarkerFaceAlpha',.33)
    scatter(s,nanmean(asymmetry2(:,s)),170,cmap(s,:),'filled')

end
xlim([.5 s+.5])
errorbar((1:segs),nanmean(asymmetry2),nanstd(asymmetry2),'LineStyle','none','Color','k','CapSize',14)
set(gcf,'color','w');


[h,p] = ttest(asymmetry2(:,1),asymmetry2(:,4))
nanmean(asymmetry2(:,1))
nanmean(asymmetry2(:,4))
nanstd(asymmetry2(:,1))
nanstd(asymmetry2(:,4))


%%
%check asymm of min jerk

%data from dylan 20200416

load('M:\analysis2\FB_norm_curves.mat')
load('M:\analysis2\HF_norm_curves.mat')

asymmetry3 = [];
for a = 1:size(norm_emp,2)
    %get asymmetry statistic
        firsthalf = trapz(0:.01:.2,norm_emp(1:21,a));
        secondhalf = trapz(0:.01:.2,norm_emp(21:41,a));
        asymmetry3(a,1) = firsthalf / (firsthalf + secondhalf);
        
        firsthalf = trapz(0:.01:.2,norm_mod(1:21,a));
        secondhalf = trapz(0:.01:.2,norm_mod(21:41,a));
        asymmetry3(a,2) = firsthalf / (firsthalf + secondhalf);
end
    
ma = mean(asymmetry3(:,2));
mad = std(asymmetry3(:,2));

%code to plot mean and std lines for MJ on asymmetry by dist plot
plot([.5 4.5], [ma ma],'k')
plot([.5 4.5], [ma-mad ma-mad],'--k')
plot([.5 4.5],[ma+mad ma+mad],'--k')



%%
%plot dylans data (fb and hf speed profiles with mj model)
emp_mean = mean(norm_emp,2);
emp_ste = std(norm_emp,0,2)/sqrt(size(norm_emp,2));

mod_mean = mean(norm_mod,2);
mod_ste = std(norm_mod,0,2)/sqrt(size(norm_mod,2));

ts = -.2:.01:.2;
col = [bmap(1,:)];
col = bmap(11,:); %for HF

figure; hold on


%model individs
for a = 1:size(norm_mod,2)
plot(ts,norm_mod(:,a),'Color',[.7 .7 .7 .3],'LineWidth',3)
end
%empirical individs
for a = 1:size(norm_emp,2)
    plot(ts,norm_emp(:,a),'Color',[col .3],'LineWidth',3)
end

%mod mean
plot(ts,mod_mean,'Color',[.5 .5 .5],'LineWidth',3)
plot(ts,mod_mean-mod_ste,'--','Color',[.5 .5 .5],'LineWidth',1)
plot(ts,mod_mean+mod_ste,'--','Color',[.5 .5 .5],'LineWidth',1)

%emp mean
plot(ts,emp_mean,'Color',col-.15,'LineWidth',3)
plot(ts,emp_mean-emp_ste,'--','Color',col-.15,'LineWidth',1)
plot(ts,emp_mean+emp_ste,'--','Color',col-.15,'LineWidth',1)

figset



%%
%measure peak delta from MJ curve for each animal, before and after peak
peak_diff_1 = []; peak_diff_2 = [];

for a = 1:size(norm_emp,2)
    delta_first = norm_emp(1:21,a) - norm_mod(1:21,a);
    peak_diff_1(a) = max(delta_first);
    
    delta_last = norm_emp(21:41,a) - norm_mod(21:41,a);
    peak_diff_2(a) = max(delta_last);

end


[h,p] = signrank(peak_diff_1,peak_diff_2)
nanmean(peak_diff_1)
nanmean(peak_diff_2)
nanstd(peak_diff_1)
nanstd(peak_diff_2)



%%
%dylans code test

tst = HFR_local;

% names of animals
fn = fieldnames(tst);
% parameters
tails = 40; % number of 10ms points on each side of peak
tot = (2*tails)+1; % total length
colnum = 5;
% wh
across_animals = nan(tot,length(fn));
for i = 1:length(fn)
    % names of reaches
    re = fieldnames(tst.(fn{i}).Reaches);
    % allocate empty array for reaches and minjerk points
    animal(i).reach = nan(tot,length(re));
    points(i).animal = nan(length(re),8);
    for ii = 1:length(re)
        % interpolation at 10 ms
        Segment = tst.(fn{i}).Reaches.(re{ii}).kin;
        % get points for min jerk
        points(i).animal(ii,:) = [Segment(1,2) Segment(end,2) Segment(1,3) Segment(end,3) Segment(1,4) Segment(end,4) Segment(end,1)-Segment(1,1) size(Segment,1)];
        % find peak for speed alignment
        [val,idx] = max(Segment(:,6));
        new_t = [Segment(idx,1)-(tails/100):0.01:Segment(idx,1)+(tails/100)]';
        t = Segment(:,1);
        %interpolate
        pos = interp1(t,Segment(:,2:4),new_t);
        vel = interp1(t,Segment(:,7:9),new_t);
        v = sqrt(sum(vel.^2,2));
        Reach = [new_t pos v vel];
        % populate reach array with variable of interest
        animal(i).reach(:,ii) = Reach(:,colnum);
        clear Reach Segment
    end
    % average by animal
    across_animals(:,i) = nanmean(animal(i).reach,2);
end



%%
%dylan code matt edit


tst = HFR_local;

% names of animals
fn = fieldnames(tst);
% parameters
tails = 40; % number of 10ms points on each side of peak
tot = (2*tails)+1; % total length
colnum = 5;
% wh
across_animals = nan(tot,length(fn));
for i = 1:length(fn)
    % names of reaches
    re = fieldnames(tst.(fn{i}).Reaches);
    % allocate empty array for reaches and minjerk points
    animal(i).reach = nan(tot,length(re));
    points(i).animal = nan(length(re),8);
    for ii = 1:length(re)
        % interpolation at 10 ms
        Segment = tst.(fn{i}).Reaches.(re{ii}).kin;
        % get points for min jerk
        points(i).animal(ii,:) = [Segment(1,2) Segment(end,2) Segment(1,3) Segment(end,3) Segment(1,4) Segment(end,4) Segment(end,1)-Segment(1,1) size(Segment,1)];
        % find peak for speed alignment
        [val,idx] = max(Segment(:,6));
        new_t = [Segment(idx,1)-(tails/100):0.01:Segment(idx,1)+(tails/100)]';
        t = Segment(:,1);
        %interpolate
        pos = interp1(t,Segment(:,2:4),new_t);
        vel = interp1(t,Segment(:,6:9),new_t);
        Reach = [new_t pos vel];
        % populate reach array with variable of interest
        animal(i).reach(:,ii) = Reach(:,colnum);
        clear Reach Segment
    end
    % average by animal
    across_animals(:,i) = nanmean(animal(i).reach,2);
end
