%alleaches bellshape analysis

%stimvel
%aligned to peak v
%reaches_stimvel = ChR2;
reaches_stimvel=ExcitRR_Ctrl;
%reaches_stimvel=ExcitRR_AllSuc;


a_s = fieldnames(reaches_stimvel);

for i = 1:length(a_s)
    
    %if you want to filter out animals with less than n reaches
%     r_count(i) = length(fieldnames(reaches_stimvel.(a_s{i}).Reaches));
%     
%     if r_count(i) < 50
%         continue
%     end
  
    T0Reaches=reaches_stimvel.(a_s{i}).RRT0;
    T0= AlignPeakV(reaches_stimvel.(a_s{i}).RRT0,7,500)
    
     T1Reaches=reaches_stimvel.(a_s{i}).RRT1;
     T1=AlignPeakV(reaches_stimvel.(a_s{i}).RRT1,7,500);
   
    T2Reaches=reaches_stimvel.(a_s{i}).RRT2;
    T2= AlignPeakV(reaches_stimvel.(a_s{i}).RRT2,7,500)
    
     T3Reaches=reaches_stimvel.(a_s{i}).RRT3;
     T3=AlignPeakV(reaches_stimvel.(a_s{i}).RRT3,7,500);
   

    
    %NOTE if using this , dim 5 becomes speed instead of vel
    %C = AlignPeakSp(reaches_stimvel.(a_s{i}).Reaches,7,200);

    T0_a = stimvel_accel(T0);% CS_a = stimvel_accel(CS); CF_a = stimvel_accel(CF);
    T1_a=stimvel_accel(T1);%SS_a=stimvel_accel(SS);
    T2_a=stimvel_accel(T2);
    T3_a=stimvel_accel(T3);%PS_a=stimvel_accel(PS);
 
    
    reaches_stimvel.(a_s{i}).stimvelT0 = T0_a;% reaches_stimvel.(a_s{i}).stimvelCS = CS_a;    reaches_stimvel.(a_s{i}).stimvelCF = CF_a;
    reaches_stimvel.(a_s{i}).stimvelT1 = T1_a;% reaches_stimvel.(a_s{i}).stimvelSS = SS_a;
    reaches_stimvel.(a_s{i}).stimvelT2 = T2_a;
    reaches_stimvel.(a_s{i}).stimvelT3 = T3_a; %reaches_stimvel.(a_s{i}).stimvelPS = PS_a; 

    
    
    reaches_stimvel.(a_s{i}).T0Kins = BasicKinematics(reaches_stimvel.(a_s{i}).RRT0);
    reaches_stimvel.(a_s{i}).T1Kins = BasicKinematics(reaches_stimvel.(a_s{i}).RRT1);
    reaches_stimvel.(a_s{i}).T2Kins = BasicKinematics(reaches_stimvel.(a_s{i}).RRT2);
    reaches_stimvel.(a_s{i}).T3Kins = BasicKinematics(reaches_stimvel.(a_s{i}).RRT3);
    
end


%%
%plotting


bmap = brewermap(length(a_s),'Set3');

All_bigT0 = struct;
All_bigT1=struct;
All_bigT2=struct;
All_bigT3=struct;

for i = 1:length(a_s)
   
    All_bigT0.stimvel(:,:,i) =reaches_stimvel.(a_s{i}).stimvelT0.Average;
    All_bigT1.stimvel(:,:,i) = reaches_stimvel.(a_s{i}).stimvelT1.Average;
    All_bigT2.stimvel(:,:,i) = reaches_stimvel.(a_s{i}).stimvelT2.Average;
    All_bigT3.stimvel(:,:,i) = reaches_stimvel.(a_s{i}).stimvelT3.Average;

    
end

% 
% j = 1;
% for i = 1: length(a_s)
%     
%     if isfield(reaches_stimvel.(a_s{i}),'stimvel')
%         col = bmap(1,:);
% %         col = [rand rand rand .4];
%         %StimvelPlot2(reaches_stimvel.(a_s{i}).stimvelC,5,0,0)
%         StimvelPlot2(reaches_stimvel.(a_s{i}).stimvelP,5,0,0)
%         
%         All_bigC.stimvel(:,:,j) = reaches_stimvel.(a_s{i}).stimvelC.Average;
%         
%         j = j+1;
%     end
%     drawnow
%     
% end
% 
% ylim auto
% xlim auto

All_bigT0.Average = nanmean(All_bigT0.stimvel,3);
All_bigT0.Std = nanstd(All_bigT0.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(All_bigT0.Average),1);
for nnum = 1:length(All_bigT0.Average)
    samp = squeeze(All_bigT0.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
All_bigT0.N = N;

%STIM
All_bigT1.Average = nanmean(All_bigT1.stimvel,3);
All_bigT1.Std = nanstd(All_bigT1.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(All_bigT1.Average),1);
for nnum = 1:length(All_bigT1.Average)
    samp = squeeze(All_bigT1.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
All_bigT1.N = N;

%STIM
All_bigT2.Average = nanmean(All_bigT2.stimvel,3);
All_bigT2.Std = nanstd(All_bigT2.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(All_bigT2.Average),1);
for nnum = 1:length(All_bigT2.Average)
    samp = squeeze(All_bigT2.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
All_bigT2.N = N;

%STIM
All_bigT3.Average = nanmean(All_bigT3.stimvel,3);
All_bigT3.Std = nanstd(All_bigT3.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(All_bigT3.Average),1);
for nnum = 1:length(All_bigT3.Average)
    samp = squeeze(All_bigT3.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
All_bigT3.N = N;


%Add acceleration column to stimvel
T0color=rgb('Thistle');
T1color=rgb('MediumOrchid');
T2color=rgb('SlateBlue');
T3color=rgb('DarkSlateBlue');


linecolor=rgb('DarkGray');


dim=5;

StimvelPlot(All_bigT0,dim,T0color,1)
StimvelPlot(All_bigT1,dim,T1color,0)
StimvelPlot(All_bigT2,dim,T2color,0)
StimvelPlot(All_bigT3,dim,T3color,0)

xlabel('Time from max V')


ylim auto
xlim([-.2 .2])
% figset
%%
All_big2T0 = struct;
All_big2T1 = struct;
All_big2T2 = struct;
All_big2T3 = struct;



for i = 1:length(a_s)
   
    reaches_stimvel.(a_s{i}).stimvel_peaknormT0 =reaches_stimvel.(a_s{i}).stimvelT0;
    maxvT0 = max(reaches_stimvel.(a_s{i}).stimvel_peaknormT0.Average(:,5));
    reaches_stimvel.(a_s{i}).stimvel_peaknormT0.Average(:,5) = reaches_stimvel.(a_s{i}).stimvel_peaknormT0.Average(:,5) / maxvT0;
    All_big2T0.stimvel(:,:,i) = reaches_stimvel.(a_s{i}).stimvel_peaknormT0.Average;  
    
    reaches_stimvel.(a_s{i}).stimvel_peaknormT1 =reaches_stimvel.(a_s{i}).stimvelT1;
    maxvT1 = max(reaches_stimvel.(a_s{i}).stimvel_peaknormT1.Average(:,5));
    reaches_stimvel.(a_s{i}).stimvel_peaknormT1.Average(:,5) = reaches_stimvel.(a_s{i}).stimvel_peaknormT1.Average(:,5) / maxvT1;
    All_big2T1.stimvel(:,:,i) = reaches_stimvel.(a_s{i}).stimvel_peaknormT1.Average; 

    reaches_stimvel.(a_s{i}).stimvel_peaknormT2 =reaches_stimvel.(a_s{i}).stimvelT2;
    maxvT2 = max(reaches_stimvel.(a_s{i}).stimvel_peaknormT2.Average(:,5));
    reaches_stimvel.(a_s{i}).stimvel_peaknormT2.Average(:,5) = reaches_stimvel.(a_s{i}).stimvel_peaknormT2.Average(:,5) / maxvT2;
    All_big2T2.stimvel(:,:,i) = reaches_stimvel.(a_s{i}).stimvel_peaknormT2.Average; 
    
    reaches_stimvel.(a_s{i}).stimvel_peaknormT3 =reaches_stimvel.(a_s{i}).stimvelT3;
    maxvT3 = max(reaches_stimvel.(a_s{i}).stimvel_peaknormT3.Average(:,5));
    reaches_stimvel.(a_s{i}).stimvel_peaknormT3.Average(:,5) = reaches_stimvel.(a_s{i}).stimvel_peaknormT3.Average(:,5) / maxvT3;
    All_big2T3.stimvel(:,:,i) = reaches_stimvel.(a_s{i}).stimvel_peaknormT3.Average; 
    
end



All_big2T0.Average = nanmean(All_big2T0.stimvel,3);
All_big2T0.Std = nanstd(All_big2T0.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(All_big2T0.Average),1);
for nnum = 1:length(All_big2T0.Average)
    samp = squeeze(All_big2T0.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
All_big2T0.N = N;

All_big2T1.Average = nanmean(All_big2T1.stimvel,3);
All_big2T1.Std = nanstd(All_big2T1.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(All_big2T1.Average),1);
for nnum = 1:length(All_big2T1.Average)
    samp = squeeze(All_big2T1.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
All_big2T1.N = N;

All_big2T2.Average = nanmean(All_big2T2.stimvel,3);
All_big2T2.Std = nanstd(All_big2T2.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(All_big2T2.Average),1);
for nnum = 1:length(All_big2T2.Average)
    samp = squeeze(All_big2T2.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
All_big2T2.N = N;

All_big2T3.Average = nanmean(All_big2T3.stimvel,3);
All_big2T3.Std = nanstd(All_big2T3.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(All_big2T3.Average),1);
for nnum = 1:length(All_big2T3.Average)
    samp = squeeze(All_big2T3.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
All_big2T3.N = N;




dim=5;

StimvelPlot(All_big2T0,dim,T0color,1)
StimvelPlot(All_big2T1,dim,T1color,0)
StimvelPlot(All_big2T2,dim,T2color,0)
StimvelPlot(All_big2T3,dim,T3color,0)


ylim auto
xlim ([-.2 .2])



%% Look at asymmmetry per bin 
asymmetryT0 = [];
asymmetryT1 = [];
asymmetryT2 = [];
asymmetryT3 = [];
SymT0=[];SymT1=[];SymT2=[];SymT3=[];
for i = 1:length(a_s)

        temp_reachesT0 = reaches_stimvel.(a_s{i}).RRT0;
        temp_reachesT1 = reaches_stimvel.(a_s{i}).RRT1;
        temp_reachesT2 = reaches_stimvel.(a_s{i}).RRT2;
        temp_reachesT3 = reaches_stimvel.(a_s{i}).RRT3;
        %align to peakv
        T0sp = AlignPeakSp(temp_reachesT0,7,200);
        T1sp = AlignPeakSp(temp_reachesT1,7,200);
        T2sp = AlignPeakSp(temp_reachesT2,7,200);
        T3sp = AlignPeakSp(temp_reachesT3,7,200);
        
        %get asymmetry statistic
        firsthalfT0 = trapz(0:.01:.2,T0sp.Average(1:21,5));
        secondhalfT0 = trapz(0:.01:.2,T0sp.Average(21:41,5));
        asymmetryT0(i) = firsthalfT0 / (firsthalfT0 + secondhalfT0);

        %get asymmetry statistic
        firsthalfT1 = trapz(0:.01:.2,T1sp.Average(1:21,5));
        secondhalfT1 = trapz(0:.01:.2,T1sp.Average(21:41,5));
        asymmetryT1(i) = firsthalfT1 / (firsthalfT1 + secondhalfT1);

        %get asymmetry statistic
        firsthalfT2 = trapz(0:.01:.2,T2sp.Average(1:21,5));
        secondhalfT2 = trapz(0:.01:.2,T2sp.Average(21:41,5));
        asymmetryT2(i) = firsthalfT2 / (firsthalfT2 + secondhalfT2);

        %get asymmetry statistic
        firsthalfT3 = trapz(0:.01:.2,T3sp.Average(1:21,5));
        secondhalfT3 = trapz(0:.01:.2,T3sp.Average(21:41,5));
        asymmetryT3(i) = firsthalfT3 / (firsthalfT3 + secondhalfT3);
        
     
        
        %save stimvel for plotting
        SymT0.stimvel(:,:,i) = T0sp.Average;
        SymT1.stimvel(:,:,i) = T1sp.Average;
        SymT2.stimvel(:,:,i) = T2sp.Average;
        SymT3.stimvel(:,:,i) = T3sp.Average;
end


cmap = brewermap(segs,'*RdBu');
pmap=brewermap(segs,'PiYg');
%pmap=brewermap(segs,'PuOr');
%plotting
figure;
hold on
for s = 1:segs
    scatter(ones(length(a_s),1)*s,asymmetryC(:,s),30,cmap(s,:),...
        'filled','jitter','on','jitterAmount',.1,...
        'MarkerFaceAlpha',.33)
    scatter(s,nanmean(asymmetryC(:,s)),170,cmap(s,:),'filled')
    
  
end
xlim([.5 s+.5])
ylim([0.35 0.55])
yline(0.5,'--k')
errorbar((1:segs),nanmean(asymmetryC),nanstd(asymmetryC),'LineStyle','none','Color','k','CapSize',14)
set(gcf,'color','w');
%%
%sort reaches based on distance. look at asymmetry.
distsT0 = [];
asymmetryT0 = [];
distsT3 = [];
asymmetryT3 = [];

SymT0.J1 = struct;
SymT0.J2 = struct;
SymT0.J3 = struct;
SymT0.J4 = struct;


SymT3.J1 = struct;
SymT3.J2 = struct;
SymT3.J3 = struct;
SymT3.J4 = struct;

jn = fieldnames(SymT0);

for i = 1:length(a_s)
    
  

    rdistsT0 = reaches_stimvel.(a_s{i}).T0Kins.dist;
    rdistsT3 = reaches_stimvel.(a_s{i}).T3Kins.dist;
    %do sorting
    [~,sort_inxT0] = sort(rdistsT0);
    [~,sort_inxT3] = sort(rdistsT3);
    
    segs = 4;
    
    for j = 1:segs
        
        

        
        %get top and bottom of segments 
        min_inxT0 = max([1 floor(((j-1)/segs)*length(sort_inxT0))]);
        max_inxT0 = floor((j/segs)*length(sort_inxT0));
        
        min_inxT3 = max([1 floor(((j-1)/segs)*length(sort_inxT3))]);
        max_inxT3 = floor((j/segs)*length(sort_inxT3));
        
        test_inxT0 = sort_inxT0(min_inxT0:max_inxT0);
        test_inxT3 = sort_inxT3(min_inxT3:max_inxT3);
        
        distsC(i,j) = mean(rdistsT0(test_inxT0));
        distsP(i,j) = mean(rdistsT3(test_inxT3));
        
        %get reaches
        %temp_reachesC = Extract(reaches_stimvel.(a_s{i}).ControlAll,test_inxC);
        %temp_reachesP = Extract(reaches_stimvel.(a_s{i}).Plus1All,test_inxP);
        temp_reachesC = Extract(reaches_stimvel.(a_s{i}).RRT0,test_inxT0);
        temp_reachesP = Extract(reaches_stimvel.(a_s{i}).RRT3,test_inxT3)
        %align to peakv
        C = AlignPeakSp(temp_reachesC,7,200);
        P = AlignPeakSp(temp_reachesP,7,200);
        
        %get asymmetry statistic
        firsthalfC = trapz(0:.01:.2,C.Average(1:21,5));
        secondhalfC = trapz(0:.01:.2,C.Average(21:41,5));
        asymmetryC(i,j) = firsthalfC / (firsthalfC + secondhalfC);
        
        firsthalfP = trapz(0:.01:.2,P.Average(1:21,5));
        secondhalfP = trapz(0:.01:.2,P.Average(21:41,5));
        asymmetryP(i,j) = firsthalfP / (firsthalfP + secondhalfP);
        
        %save stimvel for plotting
        SymT0.(jn{j}).stimvel(:,:,i) = C.Average;
        SymT3.(jn{j}).stimvel(:,:,i) = P.Average;

    end


end

cmap = brewermap(segs,'*RdBu');
pmap=brewermap(segs,'PiYg');
%pmap=brewermap(segs,'PuOr');
%plotting
figure;
hold on
for s = 1:segs
    scatter(ones(length(a_s),1)*s,asymmetryC(:,s),30,cmap(s,:),...
        'filled','jitter','on','jitterAmount',.1,...
        'MarkerFaceAlpha',.33)
    scatter(s,nanmean(asymmetryC(:,s)),170,cmap(s,:),'filled')
    
  
end
xlim([.5 s+.5])
ylim([0.35 0.55])
yline(0.5,'--k')
errorbar((1:segs),nanmean(asymmetryC),nanstd(asymmetryC),'LineStyle','none','Color','k','CapSize',14)
set(gcf,'color','w');


[h,p] = ttest(asymmetryC(:,1),asymmetryC(:,4))
nanmean(asymmetryC(:,1))
nanmean(asymmetryC(:,4))
std(asymmetryC(:,1))
std(asymmetryC(:,4))
xlabel('Distance Quartile')
ylabel('Asymmetry')
specific_ticks=[1,2,3,4];
xticks(specific_ticks)
title('Control Fail Reaches')

figure;
hold on
for s = 1:segs

    scatter(ones(length(a_s),1)*s,asymmetryP(:,s),30,pmap(s,:),...
        'filled','jitter','on','jitterAmount',.1,...
        'MarkerFaceAlpha',.33)
    scatter(s,nanmean(asymmetryP(:,s)),170,pmap(s,:),'filled')

end
xlim([.5 s+.5])
errorbar((1:segs),nanmean(asymmetryP),nanstd(asymmetryP),'LineStyle','none','Color','k','CapSize',14)
set(gcf,'color','w');


[hP,pP] = ttest(asymmetryP(:,1),asymmetryP(:,4))
nanmean(asymmetryP(:,1))
nanmean(asymmetryP(:,4))
std(asymmetryP(:,1))
std(asymmetryP(:,4))
%ylim([0.35 0.55])
yline(0.5,'--k')
xlabel('Distance Quartile')
ylabel('Asymmetry')
xticks(specific_ticks)
title('Plus1 fail. reaches')

%% Get asymmetry for every reach 
%for N+1 after suc stim
%reaches_stimvel=ConFon;
a_s=fieldnames(reaches_stimvel)

for i = 1:length(a_s)
    
    Asymmetryc=[];
    Asymmetryp=[];
    Asymmetrypf=[];
    Asymmetryss=[];
  

    %get reaches
     temp_reachesC = reaches_stimvel.(a_s{i}).RRT0;
     temp_reachesSS = reaches_stimvel.(a_s{i}).RRT1;
     temp_reachesSF = reaches_stimvel.(a_s{i}).RRT2;
    temp_reachesP = reaches_stimvel.(a_s{i}).RRT3;
 

        %align to peakv
    C = AlignPeakSp(temp_reachesC,7,200);
    SS = AlignPeakSp(temp_reachesSS,7,200);
    SF = AlignPeakSp(temp_reachesSF,7,200);
    P = AlignPeakSp(temp_reachesP,7,200);

    
    ctrlR=fieldnames(temp_reachesC);
    stimRS=fieldnames(temp_reachesSS);
    stimRF=fieldnames(temp_reachesSF);
    plusR=fieldnames(temp_reachesP);
 

        %ctrl
    for a=1:length(ctrlR)
        %get asymmetry statistic
        firsthalfC = trapz(0:.01:.2,C.stimvel(1:21,5,a));
        secondhalfC = trapz(0:.01:.2,C.stimvel(21:41,5,a));
        asymmetryC = firsthalfC / (firsthalfC + secondhalfC);
        Asymmetryc(a)=asymmetryC;
    end
        % Create a logical mask for NaN values
    nanMask = isnan(Asymmetryc);
    
    % Use the mask to extract non-NaN values
    AsymmetryC = Asymmetryc(~nanMask);
    reaches_stimvel.(a_s{i}).AsymmetryC=AsymmetryC';
   
    % stim suc 
    for p=1:length(stimRS)
        %get asymmetry statistic
        firsthalfSS = trapz(0:.01:.2,SS.stimvel(1:21,5,p));
        secondhalfSS = trapz(0:.01:.2,SS.stimvel(21:41,5,p));
        asymmetrySS = firsthalfSS / (firsthalfSS + secondhalfSS);
        AsymmetrySs(p)=asymmetrySS;
    end
        % Create a logical mask for NaN values
    nanMaskSS = isnan(AsymmetrySs);
    
    % Use the mask to extract non-NaN values
    AsymmetrySS = AsymmetrySs(~nanMaskSS);
    reaches_stimvel.(a_s{i}).AsymmetrySS=AsymmetrySS';

      % stim fail 
    for d=1:length(stimRF)
        %get asymmetry statistic
        firsthalfSF = trapz(0:.01:.2,SF.stimvel(1:21,5,d));
        secondhalfSF = trapz(0:.01:.2,SF.stimvel(21:41,5,d));
        asymmetrySF = firsthalfSF / (firsthalfSF + secondhalfSF);
        AsymmetrySf(d)=asymmetrySF;
    end
        % Create a logical mask for NaN values
    nanMaskSF = isnan(AsymmetrySf);
    
    % Use the mask to extract non-NaN values
    AsymmetrySF = AsymmetrySf(~nanMaskSF);
    reaches_stimvel.(a_s{i}).AsymmetrySF=AsymmetrySF';


    %plus1
    for b=1:length(plusR)
        firsthalfP = trapz(0:.01:.2,P.stimvel(1:21,5,b));
        secondhalfP = trapz(0:.01:.2,P.stimvel(21:41,5,b));
        asymmetryP = firsthalfP / (firsthalfP + secondhalfP);
        Asymmetryp(b)=asymmetryP;
    end
    % Create a logical mask for NaN values
    nanMaskp = isnan(Asymmetryp);
    
    % Use the mask to extract non-NaN values
    AsymmetryP = Asymmetryp(~nanMaskp);
        
    reaches_stimvel.(a_s{i}).AsymmetryP=AsymmetryP';

  

end
%% Now plot CDFS for each animal and condition 
% Different way to calculate CDF 
ctrlcolor=rgb('DarkGray');
% p1color=rgb('DarkSeaGreen');
% p1fcolor=rgb('DarkGreen');
% sscolor=rgb('LightSkyBlue');
% sfcolor=rgb('MediumBlue');

p1color=rgb('Orange');
p1fcolor=rgb('OrangeRed');
sscolor=rgb('LightSeaGreen');
sfcolor=rgb('Teal');

% Initialize variables to store all CDFs
all_cdfs_ctrl = [];
all_cdfs_ss = []; all_cdfs_sf = [];
all_cdfs_p1 = []; all_cdfs_p1f = [];

an=fieldnames(reaches_stimvel);

global_min=inf; global_max=-inf; 
% Loop through each animal to determine the global min and max
for a = 1:length(an)
   

    c = reaches_stimvel.(an{a}).AsymmetryC;
    p1 = reaches_stimvel.(an{a}).AsymmetryP;
    p1f = reaches_stimvel.(an{a}).AsymmetryPF;

    % Update global min and max
    global_min = min(global_min, min([c; p1; p1f]));
    global_max = max(global_max, max([c; p1; p1f]));
end

% Define common x_values based on global min and max
x_values = linspace(global_min, global_max, 1000);

% Plot CDFs for each animal and interpolate to common x_values
figure;
hold on;

% Loop through each animal
for a = 1:length(an)
 
    ctrl = reaches_stimvel.(an{a}).AsymmetryC;
    plus1 = reaches_stimvel.(an{a}).AsymmetryP;
    plus1f = reaches_stimvel.(an{a}).AsymmetryPF;

    % Sort data and calculate CDF
    sorted_data_ctrl = sort(ctrl);
    cdf_values_ctrl = (1:length(sorted_data_ctrl)) / length(sorted_data_ctrl);

    % plus 1 after suc
    sorted_data_p1 = sort(plus1);
    cdf_values_p1 = (1:length(sorted_data_p1)) / length(sorted_data_p1);

      % plus 1 after fail 
    sorted_data_p1f = sort(plus1f);
    cdf_values_p1f = (1:length(sorted_data_p1f)) / length(sorted_data_p1f);

    
    % Plot the CDF for the current animal
   % Use the same figure for all CDFs
    
    plot(sorted_data_p1, cdf_values_p1, 'LineWidth', 1.5,'Color',p1color);
    hold on
    plot(sorted_data_p1f, cdf_values_p1f, 'LineWidth', 1.5,'Color',p1fcolor);
    hold on
    plot(sorted_data_ctrl, cdf_values_ctrl, 'LineWidth', 1.5,'Color',ctrlcolor);
    


     % Interpolate CDFs to the common x_values
    interpolated_cdf_ctrl = interp1(sorted_data_ctrl, cdf_values_ctrl, x_values, 'linear', 'extrap');
   
    interpolated_cdf_p1 = interp1(sorted_data_p1, cdf_values_p1, x_values, 'linear', 'extrap');
    interpolated_cdf_p1f = interp1(sorted_data_p1f, cdf_values_p1f, x_values, 'linear', 'extrap');

    % Ensure CDF values are between 0 and 1
    interpolated_cdf_ctrl = min(max(interpolated_cdf_ctrl, 0), 1);

    interpolated_cdf_p1 = min(max(interpolated_cdf_p1, 0), 1);
    interpolated_cdf_p1f = min(max(interpolated_cdf_p1f, 0), 1);

    % Store interpolated CDFs
    all_cdfs_ctrl = [all_cdfs_ctrl; interpolated_cdf_ctrl];

    all_cdfs_p1 = [all_cdfs_p1; interpolated_cdf_p1];
    all_cdfs_p1f = [all_cdfs_p1f; interpolated_cdf_p1f];
end


% Calculate the mean CDF across all animals
mean_cdf_ctrl = mean(all_cdfs_ctrl, 1);

mean_cdf_p1 = mean(all_cdfs_p1, 1);
mean_cdf_p1f = mean(all_cdfs_p1f, 1);

% Plot the mean CDF on the same plot

plot(x_values, mean_cdf_p1, 'LineWidth', 5,'Color',p1color); 
hold on
plot(x_values, mean_cdf_p1f, 'LineWidth', 5,'Color',p1fcolor); 
hold on
plot(x_values, mean_cdf_ctrl, 'LineWidth', 5,'Color','k'); % Plot in black with a thicker line


% Customize the plot
%ylim([0 1])
%xlim([-20 1])
%set(gca, 'XDir', 'reverse'); % Reverse the x-axis direction
xlabel('Asymmetry');
ylabel('CDF');
title('CDFs for All Animals and Mean CDF');

% [h, p, k] = kstest2(Ctrl_VmaxP,Plus1S_VmaxP);
% [h1, p1, k1] = kstest2(Ctrl_VmaxP,Plus1F_VmaxP);
% legend( ['p-value (ctrl,Plus1S) = ' num2str(p, '%.4f');'p-value (ctrl,Plus1F) = ' num2str(p1, '%.4f')], 'Location', 'SouthEast');


% Perform K-S tests between different conditions
[data_ctrl, data_ss, data_sf, data_p1, data_p1f] = deal([]);
for a = 1:length(an)
 

    data_ctrl = [data_ctrl; reaches_stimvel.(an{a}).AsymmetryC];
    
    data_p1 = [data_p1; reaches_stimvel.(an{a}).AsymmetryP];
    data_p1f = [data_p1f;reaches_stimvel.(an{a}).AsymmetryPF];
end

% Run K-S tests between different conditions

[~, p_ctrl_p1] = kstest2(data_ctrl, data_p1);
[~, p_ctrl_p1f] = kstest2(data_ctrl, data_p1f);


% Display results
fprintf('K-S Test Results:\n');

fprintf('Control vs Plus1: p = %.4f\n', p_ctrl_p1);
fprintf('Control vs Plus1Fail: p = %.4f\n', p_ctrl_p1f);


legend( ['p-value (ctrl,Plus1S) = ' num2str(p_ctrl_p1, '%.4f');'p-value (ctrl,Plus1F) = ' num2str( p_ctrl_p1f, '%.4f')], 'Location', 'SouthEast')
%legend( ['p-value (ctrl,StimS) = ' num2str(p_ctrl_ss, '%.4f');'p-value (ctrl,StimF) = ' num2str( p_ctrl_sf, '%.4f')], 'Location', 'NorthWest')
hold off; % Release the plot

%% all categories CDF plots 
% Different way to calculate CDF 
ctrlcolor=rgb('Thistle');
p1color=rgb('DarkSlateBlue');

sscolor=rgb('MediumOrchid');
sfcolor=rgb('SlateBlue');

% Initialize variables to store all CDFs
all_cdfs_ctrl = [];
all_cdfs_ss = []; all_cdfs_sf = [];
all_cdfs_p1 = []; all_cdfs_p1f = [];

global_min=inf; global_max=-inf; 
% Loop through each animal to determine the global min and max
for a = 1:length(an)
    % Extract data for the current animal
    c = reaches_stimvel.(an{a}).AsymmetryC;
    ss = reaches_stimvel.(an{a}).AsymmetrySS;
    sf = reaches_stimvel.(an{a}).AsymmetrySF;
    p1 =  reaches_stimvel.(an{a}).AsymmetryP;


    % Update global min and max
    global_min = min(global_min, min([c; ss; sf; p1]));
    global_max = max(global_max, max([c; ss; sf; p1]));
end

% Define common x_values based on global min and max
x_values = linspace(global_min, global_max, 1000);

% Plot CDFs for each animal and interpolate to common x_values
figure;
hold on;

% Loop through each animal
for a = 1:length(an)
    % Extract data for the current animal
 

    ctrl = reaches_stimvel.(an{a}).AsymmetryC;
    stims = reaches_stimvel.(an{a}).AsymmetrySS;
    stimf = reaches_stimvel.(an{a}).AsymmetrySF;
    plus1 = reaches_stimvel.(an{a}).AsymmetryP;


    % Sort data and calculate CDF
    sorted_data_ctrl = sort(ctrl);
    cdf_values_ctrl = (1:length(sorted_data_ctrl)) / length(sorted_data_ctrl);
    
    %stim suc
    sorted_data_ss = sort(stims);
    cdf_values_ss = (1:length(sorted_data_ss)) / length(sorted_data_ss);

    %stim fail
    sorted_data_sf = sort(stimf);
    cdf_values_sf = (1:length(sorted_data_sf)) / length(sorted_data_sf);


    % plus 1 after suc
    sorted_data_p1 = sort(plus1);
    cdf_values_p1 = (1:length(sorted_data_p1)) / length(sorted_data_p1);

    
    
    % Plot the CDF for the current animal
   % Use the same figure for all CDFs
    plot(sorted_data_ctrl, cdf_values_ctrl, 'LineWidth', 1.5,'Color',ctrlcolor);
    hold on; % Hold the plot for the next animal's CDF
    plot(sorted_data_ss, cdf_values_ss, 'LineWidth', 1.5,'Color',sscolor);
    hold on
    plot(sorted_data_sf, cdf_values_sf, 'LineWidth', 1.5,'Color',sfcolor);
    hold on
    plot(sorted_data_p1, cdf_values_p1, 'LineWidth', 1.5,'Color',p1color);

    


     % Interpolate CDFs to the common x_values
    interpolated_cdf_ctrl = interp1(sorted_data_ctrl, cdf_values_ctrl, x_values, 'linear', 'extrap');
    interpolated_cdf_ss = interp1(sorted_data_ss, cdf_values_ss, x_values, 'linear', 'extrap');
    interpolated_cdf_sf = interp1(sorted_data_sf, cdf_values_sf, x_values, 'linear', 'extrap');
    interpolated_cdf_p1 = interp1(sorted_data_p1, cdf_values_p1, x_values, 'linear', 'extrap');


    % Ensure CDF values are between 0 and 1
    interpolated_cdf_ctrl = min(max(interpolated_cdf_ctrl, 0), 1);
    interpolated_cdf_ss = min(max(interpolated_cdf_ss, 0), 1);
    interpolated_cdf_sf = min(max(interpolated_cdf_sf, 0), 1);
    interpolated_cdf_p1 = min(max(interpolated_cdf_p1, 0), 1);


    % Store interpolated CDFs
    all_cdfs_ctrl = [all_cdfs_ctrl; interpolated_cdf_ctrl];
    all_cdfs_ss = [all_cdfs_ss; interpolated_cdf_ss];
    all_cdfs_sf = [all_cdfs_sf; interpolated_cdf_sf];
    all_cdfs_p1 = [all_cdfs_p1; interpolated_cdf_p1];

end


% Calculate the mean CDF across all animals
mean_cdf_ctrl = mean(all_cdfs_ctrl, 1);
mean_cdf_ss = mean(all_cdfs_ss, 1);
mean_cdf_sf = mean(all_cdfs_sf, 1);
mean_cdf_p1 = mean(all_cdfs_p1, 1);


% Plot the mean CDF on the same plot
plot(x_values, mean_cdf_ctrl, 'LineWidth', 5,'Color',ctrlcolor); % Plot in black with a thicker line
hold on
plot(x_values, mean_cdf_ss, 'LineWidth', 5,'Color',sscolor); 
hold on
plot(x_values, mean_cdf_sf, 'LineWidth', 5,'Color',sfcolor); 
hold on
plot(x_values, mean_cdf_p1, 'LineWidth', 5,'Color',p1color); 


% Customize the plot
%ylim([0 1])
%xlim([-20 1])
%set(gca, 'XDir', 'reverse'); % Reverse the x-axis direction
xlabel('Asymmetry');
ylabel('CDF');
title('CDFs for All Animals and Mean CDF');


% Perform K-S tests between different conditions
[data_ctrl, data_ss, data_sf, data_p1, data_p1f] = deal([]);
for a = 1:length(an)
    
    data_ctrl = reaches_stimvel.(an{a}).AsymmetryC;
    data_ss = reaches_stimvel.(an{a}).AsymmetrySS;
    data_sf = reaches_stimvel.(an{a}).AsymmetrySF;
    data_p1 = reaches_stimvel.(an{a}).AsymmetryP;
  
  
end

% % Run K-S tests between different conditions
% [~, p_ctrl_ss] = kstest2(data_ctrl, data_ss);
% [~, p_ctrl_sf] = kstest2(data_ctrl, data_sf);
% [~, p_ctrl_p1] = kstest2(data_ctrl, data_p1);
% [~, p_ctrl_p1f] = kstest2(data_ctrl, data_p1f);
% [~, p_ss_sf] = kstest2(data_ss, data_sf);
% [~, p_ss_p1] = kstest2(data_ss, data_p1);
% [~, p_ss_p1f] = kstest2(data_ss, data_p1f);
% [~, p_sf_p1] = kstest2(data_sf, data_p1);
% [~, p_sf_p1f] = kstest2(data_sf, data_p1f);
% [~, p_p1_p1f] = kstest2(data_p1, data_p1f);
% 
% % Display results
% fprintf('K-S Test Results:\n');
% fprintf('Control vs StimSuc: p = %.4f\n', p_ctrl_ss);
% fprintf('Control vs StimFail: p = %.4f\n', p_ctrl_sf);
% fprintf('Control vs Plus1: p = %.4f\n', p_ctrl_p1);
% fprintf('Control vs Plus1Fail: p = %.4f\n', p_ctrl_p1f);
% fprintf('StimSuc vs StimFail: p = %.4f\n', p_ss_sf);
% fprintf('StimSuc vs Plus1: p = %.4f\n', p_ss_p1);
% fprintf('StimSuc vs Plus1Fail: p = %.4f\n', p_ss_p1f);
% fprintf('StimFail vs Plus1: p = %.4f\n', p_sf_p1);
% fprintf('StimFail vs Plus1Fail: p = %.4f\n', p_sf_p1f);
% fprintf('Plus1 vs Plus1Fail: p = %.4f\n', p_p1_p1f);
% 
% %legend( ['p-value (ctrl,Plus1S) = ' num2str(p_ctrl_p1, '%.4f');'p-value (ctrl,Plus1F) = ' num2str( p_ctrl_p1f, '%.4f')], 'Location', 'NorthWest')
% legend( ['p-value (ctrl,StimS) = ' num2str(p_ctrl_ss, '%.4f');'p-value (ctrl,StimF) = ' num2str( p_ctrl_sf, '%.4f')], 'Location', 'SouthEast')
% hold off; % Release the plot


%%
cmap = brewermap(segs,'*RdBu');
pmap=brewermap(segs,'PiYg');
%pmap=brewermap(segs,'PuOr');
%plotting
figure;
hold on
for s = 1:segs
    scatter(ones(length(a_s),1)*s,asymmetryC(:,s),30,cmap(s,:),...
        'filled','jitter','on','jitterAmount',.1,...
        'MarkerFaceAlpha',.33)
    scatter(s,nanmean(asymmetryC(:,s)),170,cmap(s,:),'filled')
    
  
end
xlim([.5 s+.5])
ylim([0.35 0.55])
yline(0.5,'--k')
errorbar((1:segs),nanmean(asymmetryC),nanstd(asymmetryC),'LineStyle','none','Color','k','CapSize',14)
set(gcf,'color','w');


[h,p] = ttest(asymmetryC(:,1),asymmetryC(:,4))
nanmean(asymmetryC(:,1))
nanmean(asymmetryC(:,4))
std(asymmetryC(:,1))
std(asymmetryC(:,4))
xlabel('Distance Quartile')
ylabel('Asymmetry')
specific_ticks=[1,2,3,4];
xticks(specific_ticks)
title('Control Fail Reaches')

figure;
hold on
for s = 1:segs

    scatter(ones(length(a_s),1)*s,asymmetryP(:,s),30,pmap(s,:),...
        'filled','jitter','on','jitterAmount',.1,...
        'MarkerFaceAlpha',.33)
    scatter(s,nanmean(asymmetryP(:,s)),170,pmap(s,:),'filled')

end
xlim([.5 s+.5])
errorbar((1:segs),nanmean(asymmetryP),nanstd(asymmetryP),'LineStyle','none','Color','k','CapSize',14)
set(gcf,'color','w');


[hP,pP] = ttest(asymmetryP(:,1),asymmetryP(:,4))
nanmean(asymmetryP(:,1))
nanmean(asymmetryP(:,4))
std(asymmetryP(:,1))
std(asymmetryP(:,4))
ylim([0.25 0.6])
yline(0.5,'--k')
xlabel('Distance Quartile')
ylabel('Asymmetry')
xticks(specific_ticks)
title('Plus1 fail. reaches')

%% cdf plot 
colors=cool(4);
colorsP=autumn(4);
Cmean=mean(asymmetryC,2);
Pmean=mean(asymmetryP,2);
[cdfC,xC]=ecdf(Cmean);
   [cdfP,xP]=ecdf(Pmean);

   figure;
    plot(xC,cdfC,'Linewidth',4,'Color','k')

%     hold on
    %plot (xp,cdfp,'Linewidth',2,'Color',colorsP(i,:));
%     plot (xp,cdfp,'Linewidth',2,'Color','g');
    hold on
    plot(xP,cdfP,'Linewidth',4,'Color','g')
    
    [h,p,kstat]=kstest2(Cmean,Pmean)

%%
% go through asymmetr
figure; hold on;
for i=1:4
    data=asymmetryC(:,i);
    dataP=asymmetryP(:,i)
    [cdf,x]=ecdf(data);
    [cdfC,xC]=ecdf(Cmean);
   
    
    [cdfp,xp]=ecdf(dataP);
    [cdfP,xP]=ecdf(Pmean);
    
   % plot (x,cdf,'Linewidth',2,'Color',colors(i,:));
    plot (x,cdf,'Linewidth',2,'Color','k');
    hold on
    plot(xC,cdfC,'Linewidth',4,'Color','k')
    hold on 
    plot(x,error,'b-','LineWidth',2)
    hold on
   % plot (xp,cdfp,'Linewidth',2,'Color',colorsP(i,:));
    plot (xp,cdfp,'Linewidth',2,'Color','g');
    hold on
    plot(xP,cdfP,'Linewidth',4,'Color','g')
end
xlim([min(x), max(x)]); % Adjust the limits based on your data
ylim([0, 1]); % Y-axis limit for CDF
xlabel('X-axis Label');
ylabel('Cumulative Probability'); % Y-axis label
xlim([0.445 0.48])

%%
[cdf1,x1]=ecdf(asymmetryC);
% [cdf2,x2]=ecdf(asymmetryC(:,2));
% [cdf3,x3]=ecdf(asymmetryC(:,3));
% [cdf4,x4]=ecdf(asymmetryC(:,4));

figure()
plot (x1,cdf1,'Linewidth',2,'k');
% hold on
% plot (x2,cdf2,'Linewidth',2,'g');
% hold on
% plot (x3,cdf3,'Linewidth',2,'r');
% hold on
% plot (x4,cdf4,'Linewidth',2,'b');


%%
%plot stimvel curves
figure; hold on
for j = 1:length(jn)
    SymC.(jn{j}).Average = nanmean(SymP.(jn{j}).stimvel,3);
    SymC.(jn{j}).Std = nanstd(SymP.(jn{j}).stimvel,0,3);
    
    SymP.(jn{j}).Average = nanmean(SymP.(jn{j}).stimvel,3);
    SymP.(jn{j}).Std = nanstd(SymP.(jn{j}).stimvel,0,3);
    
    %loop through to count the number of reaches contributing to each time pt
    N = zeros(length(SymC.(jn{j}).Average),1);
    for nnum = 1:length(SymC.(jn{j}).Average)
        samp = squeeze(SymC.(jn{j}).stimvel(nnum,2,:));
        N(nnum) = length(samp(~isnan(samp)));
    end
    SymC.(jn{j}).N = N;
    
    %PLus
     N = zeros(length(SymC.(jn{j}).Average),1);
    for nnum = 1:length(SymC.(jn{j}).Average)
        samp = squeeze(SymP.(jn{j}).stimvel(nnum,2,:));
        N(nnum) = length(samp(~isnan(samp)));
    end
    SymP.(jn{j}).N = N;
    
%     subplot(2,2,j)
    hold on
    StimvelPlot2(SymC.(jn{j}),5,cmap(j,:),1)
    StimvelPlot2(SymP.(jn{j}),5,cmap(j,:),0)
    xlim([-.1 .1])
    ylim([0 40])
end



%%
%schematic figure attempt
a = 15;
C = AlignPeakSp(reaches_stimvel.(a_s{i}).T0,7,200);

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
distsC = [];
asymmetry2C = [];

distsP = [];
asymmetry2P = [];

SymC.J1 = struct;
SymC.J2 = struct;
SymC.J3 = struct;
SymC.J4 = struct;


SymP.J1 = struct;
SymP.J2 = struct;
SymP.J3 = struct;
SymP.J4 = struct;


jn = fieldnames(SymC);

for i = 1:length(a_s)
    
    rdistsC = reaches_stimvel.(a_s{i}).CKins.dist;
    %do sorting
    [~,sort_inxC] = sort(rdistsC);
    
    rdistsP = reaches_stimvel.(a_s{i}).PKins.dist;
    %do sorting
    [~,sort_inxP] = sort(rdistsP);
    
    segs = 4;
    
    for j = 1:segs
        

        %get top and bottom of segments 
        min_inxC = max([1 floor(((j-1)/segs)*length(sort_inxC))]);
        max_inxC = floor((j/segs)*length(sort_inxC));
        min_inxP = max([1 floor(((j-1)/segs)*length(sort_inxP))]);
        max_inxP = floor((j/segs)*length(sort_inxP));
        
        test_inxC = sort_inxC(min_inxC:max_inxC);
        test_inxP = sort_inxC(min_inxP:max_inxP);
        
        distsC(i,j) = mean(rdistsC(test_inxC));
        
        distsP(i,j) = mean(rdistsP(test_inxP));
        
        %get reaches
        temp_reachesC = Extract(reaches_stimvel.(a_s{i}).Control,test_inxC);
        temp_reachesP = Extract(reaches_stimvel.(a_s{i}).Plus1Suc,test_inxP);
        %align to peakv
        C = AlignPeakSp(temp_reachesC,7,200);
         P = AlignPeakSp(temp_reachesP,7,200);
        %normalize to peak
        maxvC = max(C.Average(:,5));
        C_peaknorm = C.Average(:,5) / maxvC;
       
        maxvP = max(P.Average(:,5));
        P_peaknorm = P.Average(:,5) / maxvP;
       
        
        %get asymmetry statistic
        firsthalfC = trapz(0:.01:.2,C.Average(1:21,5));
        secondhalfC = trapz(0:.01:.2,C.Average(21:41,5));
        asymmetry2C(i,j) = firsthalfC/ (firsthalfC + secondhalfC);
        
        firsthalfP = trapz(0:.01:.2,P.Average(1:21,5));
        secondhalfP = trapz(0:.01:.2,P.Average(21:41,5));
        asymmetry2P(i,j) = firsthalfP / (firsthalfP + secondhalfP);
        
        
        
        
        %save stimvel for plotting
        SymC.(jn{j}).stimvel(:,:,i) = C_peaknorm;
        SymP.(jn{j}).stimvel(:,:,i) = P_peaknorm;
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

pmap=brewermap(segs,'PiYg');

%plotting
figure;
hold on
for s = 1:segs
    scatter(ones(length(a_s),1)*s,asymmetry2P(:,s),30,Pmap(s,:),...
        'filled','jitter','on','jitterAmount',.1,...
        'MarkerFaceAlpha',.33)
    scatter(s,nanmean(asymmetry2P(:,s)),170,pmap(s,:),'filled')

end
xlim([.5 s+.5])
errorbar((1:segs),nanmean(asymmetry2P),nanstd(asymmetry2P),'LineStyle','none','Color','k','CapSize',14)
set(gcf,'color','w');


[h,p] = ttest(asymmetry2P(:,1),asymmetry2P(:,4))
nanmean(asymmetry2P(:,1))
nanmean(asymmetry2P(:,4))
nanstd(asymmetry2P(:,1))
nanstd(asymmetry2P(:,4))


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
