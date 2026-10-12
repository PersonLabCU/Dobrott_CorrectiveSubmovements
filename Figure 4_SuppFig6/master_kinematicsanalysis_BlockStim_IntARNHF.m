%name your master struct
% This code is not clipping reaches 

ConFon = struct;


fname = uigetdir;
direc = dir(strcat(fname,'\*.mat'));
%%


for i = 1:length(direc)
    %load reaches
    load(fullfile(fname,direc(i).name));
    
    %get name
    fullname = direc(i).name;
    mouse_name = regexprep(fullname,'.mat','');
    ConFon.(mouse_name).anim_data = anim_data;
    dn=fieldnames(anim_data);
    
    %setup variables 
    ag = 6;
    filt_state = 1;
    pct = .8;
    leng = 5;
    %leng = 10; % changed for Elena's data
    %thresh = .16;
  
    dirdir = 1;

    % Set up structures 
    BL = struct;
    AD = struct;
    AD_E = struct;
    AD_L = struct;
    WO_ALL = struct;
    WO = struct;
    WO_L=struct;
    
   

    
    for d=[1:length(dn)]
    %save reaches
    Reaches=anim_data.(dn{d}).Reaches;
    rn = fieldnames(Reaches);
    blk_array = [];
    for r = 1:length(rn)
        blk_array(r) = Reaches.(rn{r}).block;
    end
    %assume at least three blocks (baseline, adaptation, washout). separate
    %blocks out 
    bl_inx = find(blk_array == 1);
    BASE = Extract(Reaches,bl_inx);
    ad_inx = find(blk_array == 2);
    ADAPT = Extract(Reaches,ad_inx);
    wo_inx = find(blk_array == 3);
    WASH = Extract(Reaches,wo_inx);
      
    % % % basline % % %
    [BASE_S,BASE_C] = StimFilter(BASE);
    BASE_L= rfilts_standard(BASE_C,thresh,leng,pct);

    BL = Combine(BL,BASE_L);
    BL=accel_toreaches2(BL); % add accel to the column
     
    ConFon.(mouse_name).BL = BL; 
    
    % % % adaptation % % %
    [ADAPT_S,ADAPT_C] = StimFilter(ADAPT);
    
    if filt_state == 1
        [ADAPT_S,bada] = rfilts_standard(ADAPT_S,thresh,leng,pct);
    end

    AD = Combine(AD,ADAPT_S);
    AD=accel_toreaches2(AD); % add accel
    ConFon.(mouse_name).AD = AD;
    
    rn = fieldnames(ADAPT_S);
    %ADAPT_E = Extract(ADAPT_S,[1:10]);
    ADAPT_E = Extract(ADAPT_S,[1:5]); % extract early adapt only 5 
    %ADAPT_L = Extract(ADAPT_S,[(length(rn)-10) : length(rn)]);
    ADAPT_L = Extract(ADAPT_S,[(length(rn)-5) : length(rn)]); %extract
    %late adapte only 5
    
    if filt_state == 1
        ADAPT_E = rfilts_standard(ADAPT_E,thresh,leng,pct);
        ADAPT_L = rfilts_standard(ADAPT_L,thresh,leng,pct);
    end
    
    AD_E = Combine(AD_E,ADAPT_E);
    AD_E=accel_toreaches2(AD_E); % adding accel
    AD_L = Combine(AD_L,ADAPT_L);
    AD_L=accel_toreaches2(AD_L);
    
    ConFon.(mouse_name).AD_E = AD_E;
    ConFon.(mouse_name).AD_L = AD_L;
    
    
    % % % washout % % %
    [WASH_S,WASH_C] = StimFilter(WASH);
    WASH_Filt=rfilts_standard(WASH_C,thresh,leng,pct);
    %     WO_all = Combine(WO_all,WASH_C);
    
    rn = fieldnames(WASH_C);
    WASH_E = Extract(WASH_C,1:5); % first 5 reaches 
    %WASH_E = Extract(WASH_C,[1:2]);% first 2 reaches 
    
    %WASH_E = Extract(WASH_C,[1]);
    WASH_L = Extract(WASH_C,[(length(rn)-5) : length(rn)]);
    
    if filt_state == 1
       [ WASH_Efilt,badw] = rfilts_standard(WASH_E,thresh,leng,pct);
       [ WASH_Lfilt,badw1] = rfilts_standard(WASH_L,thresh,leng,pct);
    end
    
    WO_ALL=Combine(WO_ALL,WASH_Filt);
    WO_ALL=accel_toreaches2(WO_ALL);
    WO = Combine(WO,WASH_Efilt);
    WO=accel_toreaches2(WO);
    WO_L = Combine(WO_L,WASH_Lfilt);
    WO_L=accel_toreaches2(WO_L);
    
    ConFon.(mouse_name).WO_ALL = WO_ALL;
    ConFon.(mouse_name).WO = WO;
    ConFon.(mouse_name).WO_L = WO_L;
    
    % % % do stimvel
    % get stim good for averaging directly off of stim point instead of thresh point
    %additional filtering steps 
    %all adapt reaches 
    [AD2,bad] = stimfilt_thresh(AD,thresh);
    [AD3] = stimfilt_firstreach(AD2, thresh,1);
    
    [AD_E2,bad2] = stimfilt_thresh(AD_E,thresh);
    [AD_E3] = stimfilt_firstreach(AD_E2, thresh,1);
    
    [AD_L2,bad3] = stimfilt_thresh(AD_L,thresh);
    [AD_L3] = stimfilt_firstreach(AD_L2, thresh,1);

    
    BLN = AlignThresh(BL,thresh,1,500); %non stim trials get aligned to thresh
    BLN=stimvel_accel(BLN);
    BLN_CSM=AlignThreshCSM(BL,thresh,1,500);
    BLN_CSM=stimvel_accel(BLN_CSM);
    
    ADP = AlignStim(AD3,500); %use align stim 
    ADP=stimvel_accel(ADP);
    ADP_CSM = AlignCSM(AD3,500); %use align stim 
    ADP_CSM=stimvel_accel(ADP_CSM);
    
    ADP_E = AlignStim(AD_E3,500);
    ADP_E=stimvel_accel(ADP_E);
    ADP_E_CSM = AlignCSM(AD_E3,500);
    ADP_E_CSM=stimvel_accel(ADP_E_CSM);
    
    ADP_L = AlignStim(AD_L3,500);
    ADP_L=stimvel_accel(ADP_L);
    ADP_L_CSM = AlignCSM(AD_L3,500); %use align stim 
    ADP_L_CSM=stimvel_accel(ADP_L_CSM);
    
    WSH_ALL = AlignThresh(WO_ALL,thresh,1,500);
    %WSH_ALL = stimvel_accel(WSH_ALL);

    WSH_ALLCSM=AlignThreshCSM(WO_ALL,thresh,1,500);
    %WSH_ALLCSM = stimvel_accel(WSH_ALLCSM);

    WSH = AlignThresh(WO,thresh,1,500);
    %WSH = stimvel_accel(WSH);
    WSH_CSM=AlignThreshCSM(WO,thresh,1,500);
    %WSH_CSM = stimvel_accel(WSH_CSM);
    %WSH_CSM=stimvel_accel(WSH_CSM);
    % WSH=stimvel_accel(WSH);
    % if isempty(WSH)
    %     disp('Structure is empty');
    % end 
    
    WSH_L = AlignThresh(WO_L,thresh,1,500);
   % WSH_L = stimvel_accel(WSH_L);
    WSH_L_CSM=AlignThreshCSM(WO_L,thresh,1,500);
   % WSH_L_CSM= stimvel_accel(WSH_L_CSM);
    %WSH_L_CSM=stimvel_accel(WSH_L_CSM);
    %WSH_L=stimvel_accel(WSH_L);
%     CTCH = AlignThresh(CH,thresh,1,200);
    
%Save stimvel curves in ConFon Structure 
    ConFon.(mouse_name).BLN = BLN;
    ConFon.(mouse_name).BLN_CSM = BLN_CSM;
    ConFon.(mouse_name).ADP = ADP;
    ConFon.(mouse_name).ADP_CSM = ADP_CSM;
    ConFon.(mouse_name).ADP_E = ADP_E;
    ConFon.(mouse_name).ADP_E_CSM = ADP_E_CSM;
    ConFon.(mouse_name).ADP_L = ADP_L;
    ConFon.(mouse_name).ADP_L_CSM = ADP_L_CSM;
    ConFon.(mouse_name).WSH_ALL = WSH_ALL;
    ConFon.(mouse_name).WSH_ALLCSM = WSH_ALLCSM;
    ConFon.(mouse_name).WSH = WSH; % will need to add acceleration later- per animal basis 
    ConFon.(mouse_name).WSH_CSM = WSH_CSM;
    ConFon.(mouse_name).WSH_L = WSH_L;
    ConFon.(mouse_name).WSH_L_CSM= WSH_L_CSM;
    
    % Get Kinematic Data for blocks 
   % BL_clipped=localmaxclip(BL); % make sure this is clipped - clipped reaches 
    %BL_Kin=BasicKinematics(BL_clipped);
    BL_Kin=BasicKinematics(BL);
    ConFon.(mouse_name).BL_Kin=BL_Kin;
    
    %AD_clipped=localmaxclip(AD);
    %AD_Kin=BasicKinematics(AD_clipped);
    AD_Kin=BasicKinematics(AD);
    ConFon.(mouse_name).AD_Kin=AD_Kin;
    
    %AD_Eclipped=localmaxclip(AD_E);
    %ADE_Kin=BasicKinematics(AD_Eclipped);
    ADE_Kin=BasicKinematics(AD_E);
    ConFon.(mouse_name).ADE_Kin=ADE_Kin;
    
    %AD_Lclipped=localmaxclip(AD_L);
    %ADL_Kin=BasicKinematics(AD_Lclipped);
    ADL_Kin=BasicKinematics(AD_L);
    ConFon.(mouse_name).ADL_Kin=ADL_Kin; 
    
    %WO_clipped=localmaxclip(WO);
    %WO_Kin=BasicKinematics(WO_clipped);
    WO_Kin=BasicKinematics(WO);
    ConFon.(mouse_name).WO_Kin=WO_Kin;
    
    %WO_Lclipped=localmaxclip(WO_L);
   % WOL_Kin=BasicKinematics(WO_Lclipped);
    WOL_Kin=BasicKinematics(WO_L);
    ConFon.(mouse_name).WOL_Kin=WOL_Kin;
    
    end 
        clearvars -except ConFon fname direc %need to do this in order to get each animal/session separate 
end 
    
 
clearvars -except ConFon

%% add stimvel if needed - IntARN cells
an=fieldnames(ConFon);

for a=1:length(an)
    washall=ConFon.(an{a}).WSH_ALL;
    washall=stimvel_accel(washall);
    
    washallcsm=ConFon.(an{a}).WSH_ALLCSM;
    washallcsm=stimvel_accel(washallcsm);

    wash=ConFon.(an{a}).WSH;

    if isempty(wash)
        wash=[];
    else
        wash=stimvel_accel(wash);
    end

    washcsm=ConFon.(an{a}).WSH_CSM;
    washcsm=stimvel_accel(washcsm);

    washl=ConFon.(an{a}).WSH_L;
    washl=stimvel_accel(washl);
    
    washlcsm=ConFon.(an{a}).WSH_L_CSM;
    washlcsm=stimvel_accel(washlcsm);

    ConFon.(an{a}).WSH_ALL=washall;
    ConFon.(an{a}).WSH_ALLCSM=washallcsm;
    ConFon.(an{a}).WSH=wash;
    ConFon.(an{a}).WSH_CSM=washcsm;
    ConFon.(an{a}).WSH_L=washl;
    ConFon.(an{a}).WSH_L_CSM=washlcsm;


end 

%% Get diff in stimvel 
an=fieldnames(ConFon);
stimveldiffX_early=[];stimveldiffY_early=[];
stimveldiffX_late=[];stimveldiffY_late=[];

for a=1:length(an)
     % row 56 for stimvel 500; row 26 for stimvel 200
     % get magnitude for early trials 
    %CavgX_early=ConFon.(an{a}).BLN.Average(56,5);
    CavgX_early=ConFon.(an{a}).BLN.Average(58,5);
    %SavgX_early=ConFon.(an{a}).ADP_E.Average(56,5);
    SavgX_early=ConFon.(an{a}).WSH.Average(58,5);

    %CavgY_early=ConFon.(an{a}).BLN.Average(56,6);
    CavgY_early=ConFon.(an{a}).BLN.Average(58,6);
    %SavgY_early=ConFon.(an{a}).ADP_E.Average(56,6);
    SavgY_early=ConFon.(an{a}).WSH.Average(58,6);

    DiffX_early=SavgX_early-CavgX_early;
    DiffY_early=SavgY_early-CavgY_early;
    stimveldiffX_early(a,1)=DiffX_early;
    stimveldiffY_early(a,1)=DiffY_early;

     % get magnitude for late trials 
    %CavgX_late=ConFon.(an{a}).BLN.Average(56,5);
    CavgX_late=ConFon.(an{a}).BLN.Average(58,5);
    %SavgX_late=ConFon.(an{a}).ADP_L.Average(56,5);
    SavgX_late=ConFon.(an{a}).WSH_L.Average(58,5);

    %CavgY_late=ConFon.(an{a}).BLN.Average(56,6);
    CavgY_late=ConFon.(an{a}).BLN.Average(58,6);
   % SavgY_late=ConFon.(an{a}).ADP_L.Average(56,6);
    SavgY_late=ConFon.(an{a}).WSH_L.Average(58,6);

    DiffX_late=SavgX_late-CavgX_late;
    DiffY_late=SavgY_late-CavgY_late;
    stimveldiffX_late(a,1)=DiffX_late;
    stimveldiffY_late(a,1)=DiffY_late;
end 

ExcitAnimalsDiff_early=[stimveldiffX_early,stimveldiffY_early];
ExcitAnimalsDiff_late=[stimveldiffX_late,stimveldiffY_late];

%%
%get average combined stimvel plot for 'out significant' animals

BL_big = struct;
ADP_big = struct;
ADP_E_big=struct;
ADP_L_big=struct;
WO_big= struct;
WO_E_big= struct;
WO_L_big=struct;

% could add late WO to this but we don't really see an effect as it is 

an=fieldnames(ConFon);

for a = 1:length(an)
   
    BL_big.stimvel(:,:,a) = ConFon.(an{a}).BLN.Average;
    ADP_big.stimvel(:,:,a) = ConFon.(an{a}).ADP.Average;
    ADP_E_big.stimvel(:,:,a) = ConFon.(an{a}).ADP_E.Average;
    ADP_L_big.stimvel(:,:,a) = ConFon.(an{a}).ADP_L.Average;
    WO_big.stimvel(:,:,a) = ConFon.(an{a}).WSH_ALL.Average;
    WO_E_big.stimvel(:,:,a) = ConFon.(an{a}).WSH.Average;
    WO_L_big.stimvel(:,:,a) = ConFon.(an{a}).WSH_L.Average;
    
end




BL_big.Average = nanmean(BL_big.stimvel,3);
BL_big.Std = nanstd(BL_big.stimvel,0,3);
% 
% %loop through to count the number of reaches contributing to each time pt
N = zeros(length(BL_big.Average),1); 
for nnum = 1:length(BL_big.Average)
    samp = squeeze(BL_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
BL_big.N = N;



ADP_big.Average = nanmean(ADP_big.stimvel,3);
ADP_big.Std = nanstd(ADP_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(BL_big.Average),1); 
for nnum = 1:length(BL_big.Average)
    samp = squeeze(ADP_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
ADP_big.N = N;



ADP_E_big.Average = nanmean(ADP_E_big.stimvel,3);
ADP_E_big.Std = nanstd(ADP_E_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(BL_big.Average),1); 
for nnum = 1:length(BL_big.Average)
    samp = squeeze(ADP_E_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
ADP_E_big.N = N;


ADP_L_big.Average = nanmean(ADP_L_big.stimvel,3);
ADP_L_big.Std = nanstd(ADP_L_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(BL_big.Average),1); 
for nnum = 1:length(BL_big.Average)
    samp = squeeze(ADP_L_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
ADP_L_big.N = N;


WO_big.Average = nanmean(WO_big.stimvel,3);
WO_big.Std = nanstd(WO_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(BL_big.Average),1); 
for nnum = 1:length(BL_big.Average)
    samp = squeeze(WO_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
WO_big.N = N;

WO_E_big.Average = nanmean(WO_E_big.stimvel,3);
WO_E_big.Std = nanstd(WO_E_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(WO_E_big.Average),1); 
for nnum = 1:length(WO_E_big.Average)
    samp = squeeze(WO_E_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
WO_E_big.N = N;

WO_L_big.Average = nanmean(WO_L_big.stimvel,3);
WO_L_big.Std = nanstd(WO_L_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(BL_big.Average),1); 
for nnum = 1:length(BL_big.Average)
    samp = squeeze(WO_L_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
WO_L_big.N = N;

%% CSM
%get average combined stimvel plot for 'out significant' animals

BL_CSM_big = struct;
ADP_CSM_big = struct;
ADP_E_CSM_big=struct;
ADP_L_CSM_big=struct;
WO_ALLCSM_big= struct;
%WO_L_CSM_big=struct;

% could add late WO to this but we don't really see an effect as it is 

an=fieldnames(ConFon);

for a = 1:length(an)
   
    BL_CSM_big.stimvel(:,:,a) = ConFon.(an{a}).BLN_CSM.Average;
    ADP_CSM_big.stimvel(:,:,a) = ConFon.(an{a}).ADP_CSM.Average;
    ADP_E_CSM_big.stimvel(:,:,a) = ConFon.(an{a}).ADP_E_CSM.Average;
    ADP_L_CSM_big.stimvel(:,:,a) = ConFon.(an{a}).ADP_L_CSM.Average;
    WO_ALLCSM_big.stimvel(:,:,a) = ConFon.(an{a}).WSH_ALLCSM.Average;
    %WO_L_CSM_big.stimvel(:,:,a) = ConFon.(an{a}).WSH_L_CSM.Average;
    
end

%For structures aligned to CSM
BL_CSM_big.Average = nanmean(BL_CSM_big.stimvel,3);
BL_CSM_big.Std = nanstd(BL_CSM_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(BL_CSM_big.Average),1); 
for nnum = 1:length(BL_CSM_big.Average)
    samp = squeeze(BL_CSM_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
BL_CSM_big.N = N;


%all stim
ADP_CSM_big.Average = nanmean(ADP_CSM_big.stimvel,3);
ADP_CSM_big.Std = nanstd(ADP_CSM_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(BL_CSM_big.Average),1); 
for nnum = 1:length(BL_CSM_big.Average)
    samp = squeeze(ADP_CSM_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
ADP_CSM_big.N = N;

%early stim only
ADP_E_CSM_big.Average = nanmean(ADP_E_CSM_big.stimvel,3);
ADP_E_CSM_big.Std = nanstd(ADP_E_CSM_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(BL_CSM_big.Average),1); 
for nnum = 1:length(BL_CSM_big.Average)
    samp = squeeze(ADP_E_CSM_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
ADP_E_CSM_big.N = N;

%late stim only
ADP_L_CSM_big.Average = nanmean(ADP_L_CSM_big.stimvel,3);
ADP_L_CSM_big.Std = nanstd(ADP_L_CSM_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(BL_CSM_big.Average),1); 
for nnum = 1:length(BL_CSM_big.Average)
    samp = squeeze(ADP_L_CSM_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
ADP_L_CSM_big.N = N;

%WO
WO_ALLCSM_big.Average = nanmean(WO_ALLCSM_big.stimvel,3);
WO_ALLCSM_big.Std = nanstd(WO_ALLCSM_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(BL_CSM_big.Average),1); 
for nnum = 1:length(BL_CSM_big.Average)
    samp = squeeze(WO_ALLCSM_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
WO_ALLCSM_big.N = N;

% %WO late
% WO_L_CSM_big.Average = nanmean(WO_L_CSM_big.stimvel,3);
% WO_L_CSM_big.Std = nanstd(WO_L_CSM_big.stimvel,0,3);
% 
% %loop through to count the number of reaches contributing to each time pt
% N = zeros(length(BL_CSM_big.Average),1); 
% for nnum = 1:length(BL_CSM_big.Average)
%     samp = squeeze(WO_L_CSM_big.stimvel(nnum,2,:));
%     N(nnum) = length(samp(~isnan(samp)));
% end
% WO_L_CSM_big.N = N;
%% Add accel
BL_big=stimvel_accel(BL_big);
ADP_big=stimvel_accel(ADP_big);
ADP_E_big=stimvel_accel(ADP_E_big);
ADP_L_big=stimvel_accel(ADP_L_big);
WO_E_big= stimvel_accel(WO_E_big);
WO_L_big=stimvel_accel(WO_L_big);
%%
BL_CSM_big=stimvel_accel(BL_CSM_big);
ADP_CSM_big=stimvel_accel(ADP_CSM_big);
ADP_E_CSM_big=stimvel_accel(ADP_E_CSM_big);
ADP_L_CSM_big=stimvel_accel(ADP_L_CSM_big);
WO_CSM_big= stimvel_accel(WO_CSM_big);
WO_L_CSM_big=stimvel_accel(WO_L_CSM_big);

%% plot stimvel
stimLcolor= rgb('RoyalBlue');
stimEcolor= rgb('LightSteelBlue');
adaptEcolor=rgb('Violet');
adaptLcolor=rgb('Magenta');


dim = 2;

x1=0 ;x2=0.05;
y_upper=0.27; y_lower=0;
StimvelPlot(BL_big,dim,'k',1);
StimvelPlot(ADP_E_big,dim,stimEcolor,0);
StimvelPlot(ADP_L_big,dim,stimLcolor,0);
%StimvelPlot(WO_big,dim,adaptEcolor,0);
StimvelPlot(WO_E_big,dim,adaptEcolor,0);
StimvelPlot(WO_L_big,dim,adaptLcolor,0);

hold on 
% Fill the area between the vertical lines
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
%xlim([-.2 .2])
%ylim('auto')
%ylim([-.5,1.5])
ylim([0 0.27])
%ylim([-20 20])
yline(0,'--')
xlabel('Time from Stim.(s)')
ylabel('Outward Position (cm)')


%%
dim = 8;

x1=0 ;x2=0.05;
y_upper=0.25; y_lower=0;

StimvelPlot(BL_big,dim,[.5 .5 .5],1);
StimvelPlot(ADP_E_big,dim,stimEcolor,0);
StimvelPlot(ADP_L_big,dim,stimLcolor,0);
StimvelPlot(WO_E_big,dim,adaptEcolor,0);
StimvelPlot(WO_L_big,dim,adaptLcolor,0);


hold on
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
xlim([-.2 .2])
ylim('auto')
%ylim([-.5,1.5])
yline(0,'--')
xline(0.1)
%pbaspect([1 2 1])

xlabel('Time from Stim.(s)')
ylabel('Outward Velocity (dm/s)')
%% Get average endpoint to plot eumetria line on position plot 
BL_avgEndpt=[];
for a=1:length(an)
    BLendpt=ConFon.(an{a}).BL_Kin.endpt(:,1);
    BLmean=mean(BLendpt);
    BL_avgEndpt=[BL_avgEndpt;BLmean];
end

avg_endpt=mean(BL_avgEndpt);
%% plot position data across animals 

dim = 2; 

stimLcolor= rgb('RoyalBlue');
stimEcolor= rgb('LightSteelBlue');
adaptEcolor=rgb('Violet');
adaptLcolor=rgb('Magenta');

x1=0 ;x2=0.05;
y_upper=0.27; y_lower=0;

subplot(1,4,1); hold on 
StimvelPlot(BL_big,dim,[.5 .5 .5],0);
%shadedErrorBars(BL_big.Average(:,1),BL_big.Average(:,dim),BL_big.Std(:,dim)./sqrt(BL_big.N),[.5 .5 .5])
StimvelPlot(ADP_E_big,dim,stimEcolor,0);
%shadedErrorBars(ADP_E_big.Average(:,1),ADP_E_big.Average(:,dim),ADP_E_big.Std(:,dim)./sqrt(ADP_E_big.N),[.2 .4 .8])
hold on
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
xlim([-.2 .2])
ylim([0 0.27])

xlabel('Time from Stim.(s)')
ylabel('Outward Position (cm)')
pbaspect([1 2 1]) % for position plots 

subplot(1,4,2); hold on 
StimvelPlot(BL_big,dim,[.5 .5 .5],0);
%shadedErrorBars(BL_big.Average(:,1),BL_big.Average(:,dim),BL_big.Std(:,dim)./sqrt(BL_big.N),[.5 .5 .5])
StimvelPlot(ADP_L_big,dim,stimLcolor,0);
%shadedErrorBars(ADP_L_big.Average(:,1),ADP_L_big.Average(:,dim),ADP_L_big.Std(:,dim)./sqrt(ADP_L_big.N),[.2 .4 .8])
hold on
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
xlim([-.2 .2])
ylim([0 0.27])

xlabel('Time from Stim.(s)')
pbaspect([1 2 1]) % for position plots 

subplot(1,4,3); hold on
StimvelPlot(BL_big,dim,[.5 .5 .5],0);
%shadedErrorBars(BL_big.Average(:,1),BL_big.Average(:,dim),BL_big.Std(:,dim)./sqrt(BL_big.N),[.5 .5 .5])
StimvelPlot(WO_E_big,dim,adaptEcolor,0);
%shadedErrorBars(WO_E_big.Average(:,1),WO_E_big.Average(:,dim),WO_E_big.Std(:,dim)./sqrt(WO_E_big.N),adaptEcolor)
hold on
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
xlim([-.2 .2])
ylim([0 0.27])

xlabel('Time from Stim.(s)')
pbaspect([1 2 1]) 

subplot(1,4,4); hold on
StimvelPlot(BL_big,dim,[.5 .5 .5],0);
%shadedErrorBars(BL_big.Average(:,1),BL_big.Average(:,dim),BL_big.Std(:,dim)./sqrt(BL_big.N),[.5 .5 .5])
StimvelPlot(WO_L_big,dim,adaptLcolor,0);
%shadedErrorBars(WO_L_big.Average(:,1),WO_L_big.Average(:,dim),WO_L_big.Std(:,dim)./sqrt(WO_L_big.N),adaptLcolor)
hold on
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
xlim([-.2 .2])
ylim([0. 0.27])

xlabel('Time from Stim.(s)')
pbaspect([1 2 1]) 

%% Plot Velocity 
dim = 6; 


stimLcolor= rgb('RoyalBlue');
stimEcolor= rgb('LightSteelBlue');
adaptEcolor=rgb('Violet');
adaptLcolor=rgb('Magenta');

subplot(2,2,1);hold on; 
StimvelPlot(BL_big,dim,[.5 .5 .5],0);
shadedErrorBars(BL_big.Average(:,1),BL_big.Average(:,dim),BL_big.Std(:,dim)./sqrt(BL_big.N),[.5 .5 .5])
StimvelPlot(ADP_E_big,dim,stimEcolor,0);
shadedErrorBars(ADP_E_big.Average(:,1),ADP_E_big.Average(:,dim),ADP_E_big.Std(:,dim)./sqrt(ADP_E_big.N),stimEcolor)
yline(0,'--')
xlim([-.2 .2])
ylim([-.5,1.5])

subplot (2,2,2); hold on; 
StimvelPlot(BL_big,dim,[.5 .5 .5],0);
shadedErrorBars(BL_big.Average(:,1),BL_big.Average(:,dim),BL_big.Std(:,dim)./sqrt(BL_big.N),[.5 .5 .5])
StimvelPlot(ADP_L_big,dim,stimLcolor,0);
shadedErrorBars(ADP_L_big.Average(:,1),ADP_L_big.Average(:,dim),ADP_L_big.Std(:,dim)./sqrt(ADP_L_big.N),stimLcolor)
yline(0,'--')
xlim([-.2 .2])
ylim([-.5,1.5])

subplot(2,2,3); hold on; 
StimvelPlot(BL_big,dim,[.5 .5 .5],0);
shadedErrorBars(BL_big.Average(:,1),BL_big.Average(:,dim),BL_big.Std(:,dim)./sqrt(BL_big.N),[.5 .5 .5])
StimvelPlot(WO_E_big,dim,adaptEcolor,0);
shadedErrorBars(WO_E_big.Average(:,1),WO_E_big.Average(:,dim),WO_E_big.Std(:,dim)./sqrt(WO_E_big.N),adaptEcolor)
yline(0,'--')
xlim([-.2 .2])
ylim([-.5,1.5])

%pbaspect([2 1 1]) 

subplot(2,2,4); hold on; 
StimvelPlot(BL_big,dim,[.5 .5 .5],0);
shadedErrorBars(BL_big.Average(:,1),BL_big.Average(:,dim),BL_big.Std(:,dim)./sqrt(BL_big.N),[.5 .5 .5])
StimvelPlot(WO_L_big,dim,adaptLcolor,0);
shadedErrorBars(WO_L_big.Average(:,1),WO_L_big.Average(:,dim),WO_L_big.Std(:,dim)./sqrt(WO_L_big.N),adaptLcolor)
yline(0,'--')
xlim([-.2 .2])
ylim([-.5,1.5])

%pbaspect([1 2 1]) 

%% PLotting CSM for each block 
dim = 5; 

stimLcolor= rgb('RoyalBlue');
stimEcolor= rgb('LightSteelBlue');
adaptEcolor=rgb('Violet');
adaptLcolor=rgb('Magenta');

figure; hold on
StimvelPlot(BL_CSM_big,dim,[0 0 0],0);
StimvelPlot(ADP_E_CSM_big,dim,stimEcolor,0);
StimvelPlot(ADP_L_CSM_big,dim,stimLcolor,0);
StimvelPlot(WO_ALLCSM_big,dim,[1 0 1],0);


yline(0,'--k')
xlim([-.2 .2])
%ylim([-0.5,2])
%ylim([0 0.2])
%ylim([-50 30])
ylim('auto')
xlabel('Time from CSM (s)')
ylabel('Outward Velocity (dm/s)')
%pbaspect([1 2 1]) % for position plots 

%% Now get ReReach Info for each Block 
%% Use ReReach Filter - Control Reaches 
thresh=0.16;
dirdir=1;
an=fieldnames(ConFon);
ExcitRR_BL=[]; % Savethis structure to use for graphing rereach timing in sep. script
for a=1:length(an)

    Creaches=ConFon.(an{a}).BL;

    % ReReach
    [ReReaches,StimRR,MaxVel]=ReReachFilter_accel_ctrlreaches_edit(Creaches); % currently looking within 200ms 
    RRKin=BasicKinematics(ReReaches);

    C = AlignThresh(ReReaches,thresh,dirdir,500);
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
    
    ExcitRR_BL.(an{a}).successrate.ReReach = length(fieldnames(S_C)) / length(fieldnames(Creaches));
    ExcitRR_BL.(an{a}).successrate.NoReReach = length(fieldnames(S_S)) / length(fieldnames(Creaches));
    ExcitRR_BL.(an{a}).BL=Creaches;
    ExcitRR_BL.(an{a}).BLN=C;
    ExcitRR_BL.(an{a}).ReReaches=ReReaches;
    ExcitRR_BL.(an{a}).RRKin=RRKin;
    ExcitRR_BL.(an{a}).NoReReach=NoReach;
    ExcitRR_BL.(an{a}).NRKin=NRKin;
    ExcitRR_BL.(an{a}).BLN_N=CN;
    ExcitRR_BL.(an{a}).StimRR=StimRR;
    ExcitRR_BL.(an{a}).MaxVel=MaxVel;
    % get % of reaches with a rereach
    cnum=length(fields1);
    rrnum=length(fields2);
    percentrr=rrnum/cnum;
    ExcitRR_BL.(an{a}).Cnum=cnum;
    ExcitRR_BL.(an{a}).RRnum=rrnum;
    ExcitRR_BL.(an{a}).PercentRR=percentrr;
  
end 

%% Use ReReach Filter- STIM reaches 
ExcitRR_AD=[]; % Savethis structure to use for graphing rereach timing in sep. script
for a=1:length(an)
    Sreaches=ConFon.(an{a}).AD;
  

    % Stim Re
    [ReReaches,StimRR,MaxVel]=ReReachFilter_accel(Sreaches);
    RRKin=BasicKinematics(ReReaches);
    
    S = AlignStim(ReReaches,500);
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
    
    ExcitRR_AD.(an{a}).successrate.ReReach = length(fieldnames(S_C)) / length(fieldnames(Sreaches));
    ExcitRR_AD.(an{a}).successrate.NoReReach = length(fieldnames(S_S)) / length(fieldnames(Sreaches));
    ExcitRR_AD.(an{a}).AD=Sreaches;
    ExcitRR_AD.(an{a}).ADP=S;
    ExcitRR_AD.(an{a}).ReReaches=ReReaches;
    ExcitRR_AD.(an{a}).RRKin=RRKin;
    ExcitRR_AD.(an{a}).NoReReach=NoReach;
    ExcitRR_AD.(an{a}).NRKin=NRKin;
    ExcitRR_AD.(an{a}).ADP_N=SN;
    ExcitRR_AD.(an{a}).StimRR=StimRR;
    ExcitRR_AD.(an{a}).MaxVel=MaxVel;
    snum=length(fields1);
    rrnum=length(fields2);
    percentrr=rrnum/snum;
    ExcitRR_AD.(an{a}).Snum=snum;
    ExcitRR_AD.(an{a}).RRnum=rrnum;
    ExcitRR_AD.(an{a}).PercentRR=percentrr;
end 
%% Use ReReach Filter- STIM reaches - early 
ExcitRR_AD_E=[]; % Savethis structure to use for graphing rereach timing in sep. script
for a=1:length(an)
    Sreaches=ConFon.(an{a}).AD_E;
  

    % Stim Re
    [ReReaches,StimRR,MaxVel]=ReReachFilter_accel(Sreaches);
    RRKin=BasicKinematics(ReReaches);
    
    S = AlignStim(ReReaches,500);
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
    
    ExcitRR_AD_E.(an{a}).successrate.ReReach = length(fieldnames(S_C)) / length(fieldnames(Sreaches));
    ExcitRR_AD_E.(an{a}).successrate.NoReReach = length(fieldnames(S_S)) / length(fieldnames(Sreaches));
    ExcitRR_AD_E.(an{a}).AD_E=Sreaches;
    ExcitRR_AD_E.(an{a}).ADP_E=S;
    ExcitRR_AD_E.(an{a}).ReReaches=ReReaches;
    ExcitRR_AD_E.(an{a}).RRKin=RRKin;
    ExcitRR_AD_E.(an{a}).NoReReach=NoReach;
    ExcitRR_AD_E.(an{a}).NRKin=NRKin;
    ExcitRR_AD_E.(an{a}).ADP_EN=SN;
    ExcitRR_AD_E.(an{a}).StimRR=StimRR;
    ExcitRR_AD_E.(an{a}).MaxVel=MaxVel;
    snum=length(fields1);
    rrnum=length(fields2);
    percentrr=rrnum/snum;
    ExcitRR_AD_E.(an{a}).Snum=snum;
    ExcitRR_AD_E.(an{a}).RRnum=rrnum;
    ExcitRR_AD_E.(an{a}).PercentRR=percentrr;
end 

%% Use ReReach Filter- STIM reaches 
ExcitRR_AD_L=[]; % Savethis structure to use for graphing rereach timing in sep. script
for a=1:length(an)
    Sreaches=ConFon.(an{a}).AD_L;
  

    % Stim Re
    [ReReaches,StimRR,MaxVel]=ReReachFilter_accel(Sreaches);
    RRKin=BasicKinematics(ReReaches);
    
    S = AlignStim(ReReaches,500);
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
    
    ExcitRR_AD_L.(an{a}).successrate.ReReach = length(fieldnames(S_C)) / length(fieldnames(Sreaches));
    ExcitRR_AD_L.(an{a}).successrate.NoReReach = length(fieldnames(S_S)) / length(fieldnames(Sreaches));
    ExcitRR_AD_L.(an{a}).AD_L=Sreaches;
    ExcitRR_AD_L.(an{a}).ADP_L=S;
    ExcitRR_AD_L.(an{a}).ReReaches=ReReaches;
    ExcitRR_AD_L.(an{a}).RRKin=RRKin;
    ExcitRR_AD_L.(an{a}).NoReReach=NoReach;
    ExcitRR_AD_L.(an{a}).NRKin=NRKin;
    ExcitRR_AD_L.(an{a}).ADP_N=SN;
    ExcitRR_AD_L.(an{a}).StimRR=StimRR;
    ExcitRR_AD_L.(an{a}).MaxVel=MaxVel;
    snum=length(fields1);
    rrnum=length(fields2);
    percentrr=rrnum/snum;
    ExcitRR_AD_L.(an{a}).Snum=snum;
    ExcitRR_AD_L.(an{a}).RRnum=rrnum;
    ExcitRR_AD_L.(an{a}).PercentRR=percentrr;
end 
%% Use ReReach Filter -Washout Reaches 

ExcitRR_WO=[]; % Savethis structure to use for graphing rereach timing in sep. script
for a=1:length(an)

    Creaches=ConFon.(an{a}).WO_ALL;

 
    % ReReach
    [ReReaches,StimRR,MaxVel]=ReReachFilter_accel_ctrlreaches_edit(Creaches); % currently looking within 200ms 
    RRKin=BasicKinematics(ReReaches);

    C = AlignThresh(ReReaches,thresh,dirdir,500);
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
    
    ExcitRR_WO.(an{a}).successrate.ReReach = length(fieldnames(S_C)) / length(fieldnames(Creaches));
    ExcitRR_WO.(an{a}).successrate.NoReReach = length(fieldnames(S_S)) / length(fieldnames(Creaches));
    ExcitRR_WO.(an{a}).WO=Creaches;
    ExcitRR_WO.(an{a}).WSH=C;
    ExcitRR_WO.(an{a}).ReReaches=ReReaches;
    ExcitRR_WO.(an{a}).RRKin=RRKin;
    ExcitRR_WO.(an{a}).NoReReach=NoReach;
    ExcitRR_WO.(an{a}).NRKin=NRKin;
    ExcitRR_WO.(an{a}).WSH_N=CN;
    ExcitRR_WO.(an{a}).StimRR=StimRR;
    ExcitRR_WO.(an{a}).MaxVel=MaxVel;
    % get % of reaches with a rereach
    cnum=length(fields1);
    rrnum=length(fields2);
    percentrr=rrnum/cnum;
    ExcitRR_WO.(an{a}).Cnum=cnum;
    ExcitRR_WO.(an{a}).RRnum=rrnum;
    ExcitRR_WO.(an{a}).PercentRR=percentrr;
  
end 
%% Use ReReach Filter -Washout Reaches 

ExcitRR_WO_E=[]; % Savethis structure to use for graphing rereach timing in sep. script
for a=1:length(an)

    Creaches=ConFon.(an{a}).WO;

 
    % ReReach
    [ReReaches,StimRR,MaxVel]=ReReachFilter_accel_ctrlreaches_edit(Creaches); % currently looking within 200ms 
    RRKin=BasicKinematics(ReReaches);

    C = AlignThresh(ReReaches,thresh,dirdir,500);
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
        %CN=stimvel_accel(CN); %causing issues in the code 
    end 

    [S_C, F_C] = SuccessFilter(ReReaches);
    [S_S, F_S] = SuccessFilter(NoReach);
    
    ExcitRR_WO_E.(an{a}).successrate.ReReach = length(fieldnames(S_C)) / length(fieldnames(Creaches));
    ExcitRR_WO_E.(an{a}).successrate.NoReReach = length(fieldnames(S_S)) / length(fieldnames(Creaches));
    ExcitRR_WO_E.(an{a}).WO_E=Creaches;
    ExcitRR_WO_E.(an{a}).WSH_E=C;
    ExcitRR_WO_E.(an{a}).ReReaches=ReReaches;
    ExcitRR_WO_E.(an{a}).RRKin=RRKin;
    ExcitRR_WO_E.(an{a}).NoReReach=NoReach;
    ExcitRR_WO_E.(an{a}).NRKin=NRKin;
    ExcitRR_WO_E.(an{a}).WSH_N=CN;
    ExcitRR_WO_E.(an{a}).StimRR=StimRR;
    ExcitRR_WO_E.(an{a}).MaxVel=MaxVel;
    % get % of reaches with a rereach
    cnum=length(fields1);
    rrnum=length(fields2);
    percentrr=rrnum/cnum;
    ExcitRR_WO_E.(an{a}).Cnum=cnum;
    ExcitRR_WO_E.(an{a}).RRnum=rrnum;
    ExcitRR_WO_E.(an{a}).PercentRR=percentrr;
  
end 

%% Use ReReach Filter -Washout Reaches 

ExcitRR_WO_L=[]; % Savethis structure to use for graphing rereach timing in sep. script
for a=1:length(an)

    Creaches=ConFon.(an{a}).WO_L;

 
    % ReReach
    [ReReaches,StimRR,MaxVel]=ReReachFilter_accel_ctrlreaches_edit(Creaches); % currently looking within 200ms 
    RRKin=BasicKinematics(ReReaches);

    C = AlignThresh(ReReaches,thresh,dirdir,500);
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
        %CN=stimvel_accel(CN); %causing issues in the code 
    end 

    [S_C, F_C] = SuccessFilter(ReReaches);
    [S_S, F_S] = SuccessFilter(NoReach);
    
    ExcitRR_WO_L.(an{a}).successrate.ReReach = length(fieldnames(S_C)) / length(fieldnames(Creaches));
    ExcitRR_WO_L.(an{a}).successrate.NoReReach = length(fieldnames(S_S)) / length(fieldnames(Creaches));
    ExcitRR_WO_L.(an{a}).WO=Creaches;
    ExcitRR_WO_L.(an{a}).WSH=C;
    ExcitRR_WO_L.(an{a}).ReReaches=ReReaches;
    ExcitRR_WO_L.(an{a}).RRKin=RRKin;
    ExcitRR_WO_L.(an{a}).NoReReach=NoReach;
    ExcitRR_WO_L.(an{a}).NRKin=NRKin;
    ExcitRR_WO_L.(an{a}).WSH_N=CN;
    ExcitRR_WO_L.(an{a}).StimRR=StimRR;
    ExcitRR_WO_L.(an{a}).MaxVel=MaxVel;
    % get % of reaches with a rereach
    cnum=length(fields1);
    rrnum=length(fields2);
    percentrr=rrnum/cnum;
    ExcitRR_WO_L.(an{a}).Cnum=cnum;
    ExcitRR_WO_L.(an{a}).RRnum=rrnum;
    ExcitRR_WO_L.(an{a}).PercentRR=percentrr;
  
end 