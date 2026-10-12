%% NP data org per session for population analysis


%load NP data
%establish directory: navigate to 'recording' subfolder
pathy = uigetdir;
pathx = string(pathy);

%Make sure to change npholder for the 
data=np_holder(pathx); % view 'np_holder_usage' for more info

data.get_event_times(.05);


%%
%load behavior data
[fz,pz] = uigetfile;
load(fullfile(pz,fz));

%process behav file into reaches
[Reaches,ided,Kins] = HFKinsWrap_blk(dta,suc,r_act_time,blk_track); % adding blk_track so that blocks can bbe added to reaches 

%%
%unify timebase
NP_events = data.event_times;
time_dif = TimeOffset_NP(NP_events,ided);

%correct reach timebase
ided.m1(:,1) = ided.m1(:,1) + time_dif;
Kins(:,1) = Kins(:,1) + time_dif;
rn = fieldnames(Reaches);
for r = 1:length(rn)
    Reaches.(rn{r}).kin(:,1) = Reaches.(rn{r}).kin(:,1) + time_dif;
end
thresh=0.1; leng=6; pct=0.5;
[Reaches,bad1] = shortfilter(Reaches,leng);
[Reaches,bad2] = highstartfilter(Reaches,thresh);
[Reaches,bad3] = threshcrossfilter(Reaches,thresh);
[Reaches,bad4] = coveragefilter(Reaches,pct);

[Stim,Control]=StimFilter(Reaches);

Reaches=accel_toreaches_jerk(Reaches);
Stim=accel_toreaches_jerk(Stim);
Control=accel_toreaches_jerk(Control);
% 
% %% looking at kins 
% dirdir=1;
% 
% goodr=AlignThresh(Reaches,thresh,dirdir,500);
% goodr=stimvel_accel(goodr);
% badr=AlignThresh(bad,thresh,dirdir,500);
% badr=stimvel_accel(badr);
% 
% dim=8;
% StimvelPlot(goodr,dim,'k',1)
% StimvelPlot(badr,dim,'r',0)
% xlim([-0.2 0.2])
% ylim('auto')

%% Use ReReach Filter - Control Reaches 
thresh=0.1;
dirdir=1;

ExcitRR=[]; 
%Creaches=Reaches;
Creaches=Control;

% ReReach
[ReReaches,StimRR,MaxVel]=ReReachFilter_accel_ctrlreachesHF_edit(Creaches); % currently looking within 200ms and pos accel has to be sustained for 50ms
RRKin=BasicKinematics(ReReaches);
C = AlignThresh(ReReaches,thresh,dirdir,500);
C=stimvel_accel(C);
C=stimvel_accel_jerk(C);

% find difference in fieldnames to get stim re-reaches without RR
fields1 = fieldnames(Creaches);
fields2 = fieldnames(ReReaches);

% Find field names present in structure1 but not in structure2
fields_only_in_structure1 = setdiff(fields1, fields2);
NoCSM=struct;
for g=1:length(fields_only_in_structure1)
    if isempty(fields_only_in_structure1)
        NoCSM=[];
    else
    NoCSM.(fields_only_in_structure1{g})=Creaches.(fields_only_in_structure1{g});
    end
end

if isempty(NoCSM)
    NRkin=[];
    CN=[];
else
    NRKin=BasicKinematics(NoCSM);
    CN = AlignThresh(NoCSM,thresh,dirdir,500);
    CN=stimvel_accel(CN); 
    CN=stimvel_accel_jerk(CN);
end 

[S_CSM, F_CSM] = SuccessFilter(ReReaches);
[S_NoCSM, F_NoCSM] = SuccessFilter(NoCSM);

ExcitRR.successrate.ReReach = length(fieldnames(S_CSM)) / length(fieldnames(Creaches));
ExcitRR.successrate.NoReReach = length(fieldnames(S_NoCSM)) / length(fieldnames(Creaches));
ExcitRR.Creaches=Creaches;
ExcitRR.C=C;
ExcitRR.ReReaches=ReReaches;
ExcitRR.RRKin=RRKin;
ExcitRR.NoCSM=NoCSM;
ExcitRR.NRKin=NRKin;
ExcitRR.BLN_N=CN;
ExcitRR.StimRR=StimRR;
ExcitRR.MaxVel=MaxVel;
% get % of reaches with a rereach
cnum=length(fields1);
rrnum=length(fields2);
percentrr=rrnum/cnum;
ExcitRR.Cnum=cnum;
ExcitRR.RRnum=rrnum;
ExcitRR.PercentRR=percentrr;

%% plot histo of times 
bs=ExcitRR.StimRR(:,1);

blcolor=rgb('Black');

figure;
histogram(bs,'BinWidth',0.02,'Normalization','Count','EdgeColor','none','FaceColor',blcolor)
xlim([0 0.2])
%% Separate out into different types of reaches 
thresh=0.1; dirdir= 1;

Rreaches=ReReaches;
rrt=ExcitRR.StimRR; rrt=rrt(:,1);

% Timings taken from histo from all animals 
T0=[];
T1=[];
T2=[];
T3=[];

 % % Bins covering all - for headfixed reaches only 
T0bin=find(rrt==0);
T1bin=find(rrt>0 & rrt<=0.06);
T2bin=find(rrt>0.06 & rrt<=0.12);
T3bin=find(rrt>0.12 & rrt<=0.2);

rnum=fieldnames(Rreaches);

% for all reaches with a 0 time 
for d=1:length(T0bin)
    idx0=T0bin(d);
    T0.(rnum{idx0})= ExcitRR.ReReaches.(rnum{idx0});  
end 

if isempty(T0)
    T0Kin=[];
    ST0=[];
    ST0_CSM=[];
else
     T0Kin=BasicKinematics(T0);
     ST0 = AlignThresh(T0,thresh,dirdir,500);
     ST0=stimvel_accel(ST0); ST0=stimvel_accel_jerk(ST0);
     ST0_CSM=AlignThreshCSM_T0HF(T0,thresh,dirdir,500); ST0_CSM=stimvel_accel(ST0_CSM);
end 


for f=1:length(T1bin)
    idx1=T1bin(f);
    T1.(rnum{idx1})= ExcitRR.ReReaches.(rnum{idx1});  
end 

if isempty(T1)
    T1Kin=[];
    ST1=[];
    ST1_CSM=[];
else
     T1Kin=BasicKinematics(T1);
     ST1 = AlignThresh(T1,thresh,dirdir,500);
     ST1=stimvel_accel(ST1); ST1=stimvel_accel_jerk(ST1);
     ST1_CSM=AlignThreshCSM_T1HF(T1,thresh,dirdir,500); ST1_CSM=stimvel_accel(ST1_CSM);
end 


for h=1:length(T2bin)
    idx2=T2bin(h);
    T2.(rnum{idx2})= ExcitRR.ReReaches.(rnum{idx2});  
end 
T2Kin=BasicKinematics(T2);
ST2 = AlignThresh(T2,thresh,dirdir,500);
ST2=stimvel_accel(ST2); ST2=stimvel_accel_jerk(ST2);

ST2_CSM=AlignThreshCSM_T2HF(T2,thresh,dirdir,500);
ST2_CSM=stimvel_accel(ST2_CSM);

for k=1:length(T3bin)
    idx3=T3bin(k);
    T3.(rnum{idx3})= ExcitRR.ReReaches.(rnum{idx3});  
end 

if isempty(T3)
    T3Kin=[];
    ST3=[];
    ST3_CSM=[];
else
   
     T3Kin=BasicKinematics(T3);
     ST3 = AlignThresh(T3,thresh,dirdir,500);
     ST3=stimvel_accel(ST3); ST3=stimvel_accel_jerk(ST3);
  
     ST3_CSM=AlignThreshCSM_T3HF(T3,thresh,dirdir,500);
     ST3_CSM=stimvel_accel(ST3_CSM);
end 


if isempty(T0)
    t0fields=[];
    rrnum0=[];
else
    t0fields=fieldnames(T0);
    rrnum0=length(t0fields);
end

if isempty(T1)
    t1fields=[];
    rrnum1=[];
else
    t1fields=fieldnames(T1);
    rrnum1=length(t1fields);
end

t2fields=fieldnames(T2);
rrnum2=length(t2fields);

if isempty(T3)
    t3fields=[];
    rrnum3=[];
else
    t3fields=fieldnames(T3);
    rrnum3=length(t3fields);
end



ExcitRR.rrnum0=rrnum0;
ExcitRR.rrnum1=rrnum1;
ExcitRR.rrnum2=rrnum2;
ExcitRR.rrnum3=rrnum3;

 % Kinematics
 ExcitRR.T0Kin=T0Kin;
 ExcitRR.T1Kin=T1Kin;
 ExcitRR.T2Kin=T2Kin;
 ExcitRR.T3Kin=T3Kin;

 %Stimvel
 ExcitRR.ST0=ST0;
 ExcitRR.ST1=ST1;
 ExcitRR.ST2=ST2;
 ExcitRR.ST3=ST3;

 ExcitRR.ST0_CSM=ST0_CSM;
 ExcitRR.ST1_CSM=ST1_CSM;
 ExcitRR.ST2_CSM=ST2_CSM;
 ExcitRR.ST3_CSM=ST3_CSM;

 %% Get data aligned to CSM
thresh = .1; dirdir = 1;
% 
T0R=AlignThresh(T0,thresh,dirdir,500);
T0R=stimvel_accel(T0R);
T0R=stimvel_accel_jerk(T0R);

T1R=AlignThresh(T1,thresh,dirdir,500);
T1R=stimvel_accel(T1R);
T1R=stimvel_accel_jerk(T1R);

T2R=AlignThresh(T2,thresh,dirdir,500);
T2R=stimvel_accel(T2R);
T2R=stimvel_accel_jerk(T2R);

T3R=AlignThresh(T3,thresh,dirdir,500);
T3R=stimvel_accel(T3R);
T3R=stimvel_accel_jerk(T3R);
%% Plot
dim = 8; 

%Add acceleration column to stimvel
T0color=rgb('Thistle');
T1color=rgb('MediumOrchid');
T2color=rgb('SlateBlue');
T3color=rgb('DarkSlateBlue');
NoCSMcolor=rgb('OrangeRed');

StimvelPlot(ST0,dim,T0color,1)
StimvelPlot(ST1,dim,T1color,0)
StimvelPlot(ST2,dim,T2color,0)
StimvelPlot(ST3,dim,T3color,0)
StimvelPlot(CN,dim,NoCSMcolor,0)

yline(0,'--k')
xlim([-.2 .2])
%ylim([-0.5,2])
%ylim([0 0.2])
%ylim([-50 30])
ylim('auto')
xlabel('Time from threshold crossing (s)')
ylabel('Outward Velocity (dm/s)')
%pbaspect([1 2 1]) % for position plots 
%% Now get TOreach StimRR timing

thresh=0.1;

[ReReachT0, StimRRT0, MaxVelT0,ReachesWithoutSustainedJerkT0] = CSMEarlyThresh_jerk_edit5(T0,thresh);
ExcitRR.StimRRT0=StimRRT0;
%% Take out0s from StimRRy

RR_Excit=ExcitRR.StimRR;
nonZeroRows_ex = [];
for i = 1:size(RR_Excit, 1)
    if any(RR_Excit(i,:) == 0)
        continue;
    end
    nonZeroRows_ex = [nonZeroRows_ex; RR_Excit(i,:)];
end
disp(nonZeroRows_ex);
ExcitRR.StimRREdit=nonZeroRows_ex;
%% Combine all timings to look at distribution

excitall = ExcitRR.StimRREdit(:,1);
t0bins=ExcitRR.StimRRT0(:,1);
allbins=[excitall;t0bins];
ExcitRR.AllBins=allbins;

%% plot histo of times 
bs=ExcitRR.AllBins(:,1);

blcolor=rgb('Black');

figure;
histogram(bs,'BinWidth',0.01,'Normalization','Count','EdgeColor','none','FaceColor',blcolor)
%xlim([0 0.2])
%% Align NoRR with biggest bin the T2
NoRR_T2aligned=AlignThreshCSM_T2HF(NoCSM,thresh,dirdir,500);
NoRR_T2aligned=stimvel_accel(NoRR_T2aligned);

% Align to T0
NoRR_T0aligned=AlignThreshCSM_T0HF(NoCSM,thresh,dirdir,500);
NoRR_T0aligned=stimvel_accel(NoRR_T0aligned);
%% PLot aligned to CSM
NoCSMcolorT0=rgb('Orange');
dim=8;
StimvelPlot(NoRR_T2aligned,dim,NoCSMcolor,1)
StimvelPlot(NoRR_T0aligned,dim,NoCSMcolorT0,0)
StimvelPlot(ST2_CSM,dim,T2color,0)
StimvelPlot(ST0_CSM,dim,T0color,0)


yline(0,'--k')
xlim([-.2 .2])
%ylim([-0.5,2])
%ylim([0 0.2])
%ylim([-50 50])
ylim('auto')
xlabel('Time from CSM (s)')
ylabel('Outward Velocity (dm/s)')
%% Get Events times for each set of reaches 
% Base
Events_All = GetEvents_NP(Reaches);
Events_NoCSM=GetEvents_NP(NoCSM);
Events_T0= GetEvents_NP(T0);
Events_T1= GetEvents_NP(T1);
Events_T2= GetEvents_NP(T2); % seems to be the most reaches in this group
Events_T3= GetEvents_NP(T3);
%% add average CSM crossing time for NoCSM

reach_ind = fieldnames(NoCSM);

thresh=0.1; % For headfixed reaches 

for reachnum = 1:length(reach_ind)
    reach = NoCSM.(reach_ind{reachnum}).kin;
    stim_inx = find(reach(:,2)>=thresh,1); 

    windowstart=reach(stim_inx,1); windowend=windowstart+0.087; % CSMT2 
    
    differences = abs(reach(:,1) -windowend); % find value that is closest to 100 ms after stim 

    % Find the index of the minimum absolute difference
    [minDifference, index] = min(differences);
    
    Events_NoCSM.csmT2(reachnum) = (reach(index,1));
   
end 
%% add average CSM crossing time for NoCSM aligned to T0

reach_ind = fieldnames(NoCSM);

thresh=0.1; % For headfixed reaches 

for reachnum = 1:length(reach_ind)
    reach = NoCSM.(reach_ind{reachnum}).kin;
    stim_inx = find(reach(:,2)>=thresh,1); 

    windowstart=reach(stim_inx,1); windowend=windowstart+(-0.055); % CSMT0
    
    differences = abs(reach(:,1) -windowend); % find value that is closest to 100 ms after stim 

    % Find the index of the minimum absolute difference
    [minDifference, index] = min(differences);
    
    Events_NoCSM.csmT0(reachnum) = (reach(index,1));
   
end 
%% add average CSM crossing time for T0

reach_ind = fieldnames(T0);

thresh=0.1; % For headfixed reaches 

for reachnum = 1:length(reach_ind)
    reach = T0.(reach_ind{reachnum}).kin;
    stim_inx = find(reach(:,2)>=thresh,1); 

    windowstart=reach(stim_inx,1); windowend=windowstart+(-0.055); % CSMT0 
    
    differences = abs(reach(:,1) -windowend); % find value that is closest to 100 ms after stim 

    % Find the index of the minimum absolute difference
    [minDifference, index] = min(differences);
    
    Events_T0.csmT0(reachnum) = (reach(index,1));
   
end 
%% add average CSM crossing time for T2

reach_ind = fieldnames(T2);

thresh=0.1; % For headfixed reaches 

for reachnum = 1:length(reach_ind)
    reach = T2.(reach_ind{reachnum}).kin;
    stim_inx = find(reach(:,2)>=thresh,1); 

    windowstart=reach(stim_inx,1); windowend=windowstart+0.087; % CSMT2 
    
    differences = abs(reach(:,1) -windowend); % find value that is closest to 100 ms after stim 

    % Find the index of the minimum absolute difference
    [minDifference, index] = min(differences);
    
    Events_T2.csmT2(reachnum) = (reach(index,1));
   
end 
%% Get Events for Baseline Trials -T0
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["pos_thresh"]; 
%event_name = ["end_point","pos_thresh","reach_start","max_vx"]; 
enames = fieldnames(Events_T0);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_T0.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    % 
    data.units(c).T0_pos_thresh(enum).Referenced = Ref;
    data.units(c).T0_pos_thresh(enum).EventTrial = EvTr;
    data.units(c).T0_pos_thresh(enum).Raster = BinMast;

    [psth_raw,psth_smooth] = Psth1(BinMast,5); % was 5 before 
    data.units(c).T0_pos_thresh(enum).PSTH_raw = psth_raw;
    data.units(c).T0_pos_thresh(enum).PSTH_smooth = psth_smooth;   

  
end
end

%% Get Events for Baseline Trials -T1
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["pos_thresh"]; 
%event_name = ["end_point","pos_thresh","reach_start","max_vx"]; 
enames = fieldnames(Events_T1);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_T1.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    % 
    data.units(c).T1_pos_thresh(enum).Referenced = Ref;
    data.units(c).T1_pos_thresh(enum).EventTrial = EvTr;
    data.units(c).T1_pos_thresh(enum).Raster = BinMast;

    [psth_raw,psth_smooth] = Psth1(BinMast,5); % was 5 before 
    data.units(c).T1_pos_thresh(enum).PSTH_raw = psth_raw;
    data.units(c).T1_pos_thresh(enum).PSTH_smooth = psth_smooth;   
end
end

%% Get Events for Baseline Trials -T2
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["pos_thresh"]; 
%event_name = ["end_point","pos_thresh","reach_start","max_vx"]; 
enames = fieldnames(Events_T2);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_T2.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    % 
    data.units(c).T2_pos_thresh(enum).Referenced = Ref;
    data.units(c).T2_pos_thresh(enum).EventTrial = EvTr;
    data.units(c).T2_pos_thresh(enum).Raster = BinMast;

    [psth_raw,psth_smooth] = Psth1(BinMast,5); % was 5 before 
    data.units(c).T2_pos_thresh(enum).PSTH_raw = psth_raw;
    data.units(c).T2_pos_thresh(enum).PSTH_smooth = psth_smooth;   
end
end


%% Get Events for Baseline Trials -T3
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["pos_thresh"]; 
%event_name = ["end_point","pos_thresh","reach_start","max_vx"]; 
enames = fieldnames(Events_T3);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_T3.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    % 
    data.units(c).T3_pos_thresh(enum).Referenced = Ref;
    data.units(c).T3_pos_thresh(enum).EventTrial = EvTr;
    data.units(c).T3_pos_thresh(enum).Raster = BinMast;

    [psth_raw,psth_smooth] = Psth1(BinMast,5); % was 5 before 
    data.units(c).T3_pos_thresh(enum).PSTH_raw = psth_raw;
    data.units(c).T3_pos_thresh(enum).PSTH_smooth = psth_smooth;   
end
end
%% Get Events for Baseline Trials - No CSM
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["pos_thresh"]; 
%event_name = ["end_point","pos_thresh","reach_start","max_vx"]; 
enames = fieldnames(Events_NoCSM);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_NoCSM.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    % 
    data.units(c).NoCSM_pos_thresh(enum).Referenced = Ref;
    data.units(c).NoCSM_pos_thresh(enum).EventTrial = EvTr;
    data.units(c).NoCSM_pos_thresh(enum).Raster = BinMast;

    [psth_raw,psth_smooth] = Psth1(BinMast,5); % was 5 before 
    data.units(c).NoCSM_pos_thresh(enum).PSTH_raw = psth_raw;
    data.units(c).NoCSM_pos_thresh(enum).PSTH_smooth = psth_smooth;   
end
end
%% NoCSM align CSM
% Testing out new code to plot a different type of PSTH 
baseline_epoch = 1000;
event_epoch = 1000;

% Choose event to align to
event_name = ["csmT2"]; 
enames = fieldnames(Events_NoCSM);
evnt = find(ismember(enames, event_name));

% Loop through all units
for c = 1:length(data.units)
    spike_times = data.units(c).spike_times * 1000; % Convert to ms
    
    % Loop through each event type
    for enum = 1:length(evnt)
        event_times = Events_NoCSM.(enames{evnt(enum)}) * 1000;
        window = [baseline_epoch; event_epoch];

        % Compute referenced spikes and raster binning
        [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
        data.units(c).NoCSM_csmT2(enum).Referenced = Ref;
        data.units(c).NoCSM_csmT2(enum).EventTrial = EvTr;
        data.units(c).NoCSM_csmT2(enum).Raster = BinMast;

        % Compute PSTH using original method
        [psth_raw, psth_smooth] = Psth1(BinMast, 5); % Was 5 before
        data.units(c).NoCSM_csmT2(enum).PSTH_raw = psth_raw;
        data.units(c).NoCSM_csmT2(enum).PSTH_smooth = psth_smooth;

        % Compute PSTH using new computePSTH method (10ms bins + smoothing)
        bin_size = 10; % 10ms bin size
        smoothing_bins = 3; % Moving average over 3 bins

       [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);
        
        % Store results for the new method
        data.units(c).NoCSM_csmT2(enum).PSTH_rawbin = psth_rawbin;
        data.units(c).NoCSM_csmT2(enum).PSTH_smoothbin = psth_smoothbin;
        data.units(c).NoCSM_csmT2(enum).time_bins = time_bins;

    end
end


%% NoCSM align CSM T0
% Testing out new code to plot a different type of PSTH 
baseline_epoch = 1000;
event_epoch = 1000;

% Choose event to align to
event_name = ["csmT0"]; 
enames = fieldnames(Events_NoCSM);
evnt = find(ismember(enames, event_name));

% Loop through all units
for c = 1:length(data.units)
    spike_times = data.units(c).spike_times * 1000; % Convert to ms
    
    % Loop through each event type
    for enum = 1:length(evnt)
        event_times = Events_NoCSM.(enames{evnt(enum)}) * 1000;
        window = [baseline_epoch; event_epoch];

        % Compute referenced spikes and raster binning
        [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
        data.units(c).NoCSM_csmT0(enum).Referenced = Ref;
        data.units(c).NoCSM_csmT0(enum).EventTrial = EvTr;
        data.units(c).NoCSM_csmT0(enum).Raster = BinMast;

        % Compute PSTH using original method
        [psth_raw, psth_smooth] = Psth1(BinMast, 5); % Was 5 before
        data.units(c).NoCSM_csmT0(enum).PSTH_raw = psth_raw;
        data.units(c).NoCSM_csmT0(enum).PSTH_smooth = psth_smooth;

        % Compute PSTH using new computePSTH method (10ms bins + smoothing)
        bin_size = 10; % 10ms bin size
        smoothing_bins = 3; % Moving average over 3 bins

       [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);
        
        % Store results for the new method
        data.units(c).NoCSM_csmT0(enum).PSTH_rawbin = psth_rawbin;
        data.units(c).NoCSM_csmT0(enum).PSTH_smoothbin = psth_smoothbin;
        data.units(c).NoCSM_csmT0(enum).time_bins = time_bins;

    end
end
%% T0 align CSM
%
%make rasters
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["csmT0"]; 
enames = fieldnames(Events_T0);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_T0.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    
    data.units(c).T0_csm(enum).Referenced = Ref;
    data.units(c).T0_csm(enum).EventTrial = EvTr;
    data.units(c).T0_csm(enum).Raster = BinMast;
    
    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).T0_csm(enum).PSTH_raw = psth_raw;
    data.units(c).T0_csm(enum).PSTH_smooth = psth_smooth;

    % Compute PSTH using new computePSTH method (10ms bins + smoothing)
    bin_size = 10; % 10ms bin size
    smoothing_bins = 3; % Moving average over 3 bins

   [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);
    
    % Store results for the new method
    data.units(c).T0_csm(enum).PSTH_rawbin = psth_rawbin;
    data.units(c).T0_csm(enum).PSTH_smoothbin = psth_smoothbin;
    data.units(c).T0_csm(enum).time_bins = time_bins;

    
end
end

%% T2 align CSM
%
%make rasters
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["csmT2"]; 
enames = fieldnames(Events_T2);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_T2.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    
    data.units(c).T2_csm(enum).Referenced = Ref;
    data.units(c).T2_csm(enum).EventTrial = EvTr;
    data.units(c).T2_csm(enum).Raster = BinMast;
    
    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).T2_csm(enum).PSTH_raw = psth_raw;
    data.units(c).T2_csm(enum).PSTH_smooth = psth_smooth;

    % Compute PSTH using new computePSTH method (10ms bins + smoothing)
    bin_size = 10; % 10ms bin size
    smoothing_bins = 3; % Moving average over 3 bins

   [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);
    
    % Store results for the new method
    data.units(c).T2_csm(enum).PSTH_rawbin = psth_rawbin;
    data.units(c).T2_csm(enum).PSTH_smoothbin = psth_smoothbin;
    data.units(c).T2_csm(enum).time_bins = time_bins;

    
end
end

%% all trials aligned to endpoint
%
%make rasters
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["end_point"]; 
enames = fieldnames(Events_All);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_All.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    
    data.units(c).(enames{evnt(enum)}).Referenced = Ref;
    data.units(c).(enames{evnt(enum)}).EventTrial = EvTr;
    data.units(c).(enames{evnt(enum)}).Raster = BinMast;
    
    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).(enames{evnt(enum)}).PSTH_raw = psth_raw;
    data.units(c).(enames{evnt(enum)}).PSTH_smooth = psth_smooth;
    
    
end
end
%% all trials aligned to pos thresh 
%
%make rasters
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["pos_thresh"]; 
enames = fieldnames(Events_All);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_All.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    
    data.units(c).(enames{evnt(enum)}).Referenced = Ref;
    data.units(c).(enames{evnt(enum)}).EventTrial = EvTr;
    data.units(c).(enames{evnt(enum)}).Raster = BinMast;
    
    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).(enames{evnt(enum)}).PSTH_raw = psth_raw;
    data.units(c).(enames{evnt(enum)}).PSTH_smooth = psth_smooth;
    
    
end
end

%%
% optotagging analysis - get event times for opto 
NP_opto=data.event_times;
diffNP_opto=diff(NP_opto);
idx=find(diffNP_opto>=2.95 & diffNP_opto<=3.1);%find first stim given during 
%eventopto=NP_opto(87:end); %
eventopto=NP_opto(idx(1,:):end)';
Events.eventopto=eventopto; %add opto event times to structure 

%% Get rasters for opto
%think about changing this - shorten to 500
baseline_epoch = 500;
event_epoch = 500;


event_name = ["eventopto"]; 
enames = fieldnames(Events);
evnt= find(ismember(enames, event_name));

%choose optoevents to align to 

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000;
for enum = 1:length(evnt)
    event_times = Events.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    
    data.units(c).eventopto(enum).Referenced = Ref;
    data.units(c).eventopto(enum).EventTrial = EvTr;
    data.units(c).eventopto(enum).Raster = BinMast;
    
    [psth_raw,psth_smooth] = Psth1opto(BinMast,2); % was 5 before 
    data.units(c).eventopto(enum).PSTH_raw = psth_raw;
    data.units(c).eventopto(enum).PSTH_smooth = psth_smooth;
    
    
end
end

%% Optotagging script 

cnames = length(data.units);

pre_fr = zeros(1,length(cnames));
post_fr = zeros(1,length(cnames));
pval= zeros(1,length(cnames));
Stim_zscore=zeros(1,length(cnames));

pre_epoch = 100;
post_epoch = 10;

for c = 1:length(data.units)

    pre_stim = []; post_stim = [];
    for i  = 1:size(data.units(c).eventopto.Raster,1)
        pre_stim(i) = nnz(data.units(c).eventopto.Raster(i,(500-pre_epoch):500));
        post_stim(i) = nnz(data.units(c).eventopto.Raster(i,500:(500+ post_epoch)));
    end
        
    pre_fr(c) = mean(pre_stim) / (pre_epoch/1000);
    post_fr(c) = mean(post_stim) / (post_epoch/1000);
   
    
    pre_stim = pre_stim / (pre_epoch/1000);
    post_stim = post_stim / (post_epoch/1000);
    [~,pval(c)] = ttest(pre_stim,post_stim);
   
   
    pre_mean = mean(pre_stim);
    post_mean=mean(post_stim);
    pre_std = std(pre_stim);
    stim_z = (post_mean - pre_mean) / pre_std;
    Stim_zscore(c)=stim_z;
    data.units(c).zscore=stim_z;
end

 diff_fr = post_fr - pre_fr;

%  %% unity plot
% figure; hold on
% scatter(pre_fr,post_fr,50,'filled')
% plot([1:150],[1:150],'k')
% xlabel('pre FR(100ms)')
% ylabel('post _FR (10ms)')
% 
% figure; hold on
% scatter(Stim_zscore,diff_fr,50,'filled')
% xlabel('Z Score')
% ylabel('Diff in FR (post-pre)')

%% Take cells with sig. change in FR 
opto_index=find(Stim_zscore>=1.96)
%% Now use salt to find cells that are significantly higher 

dt=0.001; % spike raster is in ms (1ms=0.001s)
wn=0.01;% size window - looking at 10ms windows 
%wn=0.005;
Pvalue=[];

for c=1:length(data.units)
    spt_baseline= data.units(c).eventopto.Raster(:,250:499); %250ms baseline period before stim  
    %spt_test=data.units(c).eventopto.Raster(:,500:549); % 50 ms stim train
    spt_test=data.units(c).eventopto.Raster(:,500:509); % 10 ms of stim- hopefully rule out any orthodromically activated cells 
    [p, I]=salt(spt_baseline,spt_test,dt,wn);
    data.units(c).pvalue=p;
    data.units(c).I=I;
    Pvalue{c,1}=p;
end 

%% Find index using salt 
salt_index=cell2mat(Pvalue);
Salt_index = find(salt_index<=0.05);
Salt_index=Salt_index'

%% Now need to plot optotagged cells for each type of trial 

%plotting

%epo=500; % for reaching endpoint 
epo = 250;
%epo=100;


evnt="eventopto";
%folder_path="C:\Users\Courtney\Documents\NeuropixelsResults\Cell Catalog\051022_SAM012\opto_raster";

%for c = opto_index
for c= Salt_index

    
    
    
    figure; hold on
    
    subplot(2,1,1); hold on
    %h = title(['cluster' ' ' num2str(data.units(c).cluster_ID) ' ' 'depth' ' ' num2str(data.units(c).depth) ' ' 'zscore' ' ' num2str(data.units(c).zscore)]) ;
    h = title(['cluster' ' ' num2str(data.units(c).cluster_ID) ' ' 'depth' ' ' num2str(data.units(c).depth) ' ' 'pvalue' ' ' num2str(data.units(c).pvalue) ' ' 'zscore ' ' ' num2str(data.units(c).zscore)]) ;
    set(h,'interpreter','none')
    
    % Plots Raster from binary matrix
    tick_length = 1; % length of tick marks for spikes
    spike_tick_width = 0.5; % width of lines marking spikes
    event_tick_width = 2; % width of lines marking secondary events
    
    %rastr = data.units(c).(evnt).Raster(:,(1000-epo):(2000-(1000-epo)));
    rastr = data.units(c).(evnt).Raster(:,(500-epo):(1000-(500-epo)));
    %rastr = data.units(c).(evnt).Raster(:,(200-epo):(400-(200-epo)));
    [y_co,x_co] = find(rastr);
    x_co = x_co - epo;
    line([x_co';x_co'],[y_co'-(tick_length/2);y_co'+(tick_length/2)],'Color', [0 0 0.5], 'LineWidth', spike_tick_width);
    line([zeros(1,size(rastr,1));zeros(1,size(rastr,1))],[(1:1:size(rastr,1))-.5;(1:1:size(rastr,1))+.5],'Color', [0 0 0.5], 'LineWidth', event_tick_width);
    axis ij;
    ylim([0,size(rastr,1)])
    xlim([-epo,epo])
    ylabel('Trial')
    
    
    subplot(2,1,2); hold on
    title('peth')
    
    psth_smooth = data.units(c).(evnt).PSTH_smooth;
    psth_raw = data.units(c).(evnt).PSTH_raw;
    td1 = length(psth_raw)-length(psth_smooth)+1;
    td2 = td1/2;
    %plot((-epo:epo),psth_smooth(((1000-epo)-td2):((2000-(1000-epo))-td2)),'Color',[.2 .4 .8],'LineWidth',3)
    plot((-epo:epo),psth_smooth(((500-epo)-td2):((1000-(500-epo))-td2)),'Color',[0 0 0.5],'LineWidth',3)
    %plot((-epo:epo),psth_smooth(((200-epo)-td2):((400-(200-epo))-td2)),'Color',[.2 .4 .8],'LineWidth',3)
    xlim([-epo,epo])
    xlabel('Time from Stim')
    ylabel('FR (Hz)')
%     ylim([0,100])
    
%      filename=sprintf('opto_raster%d.svg',c);
%     full_path=fullfile(folder_path,filename);
%     saveas(gcf,full_path);
%     
%     close;
if mod(c,10) == 0
    drawnow
    input("press enter when done looking at figures")
    close all
end

end

%% Get NC index 
% Convert struct array to cell array for easier processing
neurontypes = {data.units.neurontype}; % Extract neuron types as a cell array

% Convert to string if needed (sometimes stored as char arrays)
neurontypes = string(neurontypes)'; 

% Find the indices where neurontype is exactly "NC" (ignore extra spaces)
nc_rows = find(strtrim(neurontypes) == "NC"); 

% Extract corresponding cluster_IDs (assuming that's the index you need)
nc_indices = [data.units(nc_rows).cluster_ID]; 

% Display results
disp(nc_indices);
nc_rows=nc_rows'
length(nc_rows)

%% Go through units and identify those that are modulated by reaching 
% Go through units and identify those that are modulated by reaching 
Q_frs = struct;
Q = struct;

for d = 1:length(data.units)
%for d = nc_rows
   cname = sprintf('C%d', d);
   spike_times = data.units(d).spike_times * 1000; 
   spike_times = round(spike_times);
   spiketimes = int32((spike_times));
   spk_times = double(spiketimes);
   [inst_smooth, inst] = spktimes_instFR_NP(spk_times, 20);

   Q_frs.(cname).inst_smooth = inst_smooth;
   Q_frs.(cname).inst = inst;
   Q.(cname).spike_times = spk_times;

   kins = Kins;
   inst_zscore = zscore_instFR_NP(inst_smooth, kins, [100 100], Reaches);
   Q_frs.(cname).z = inst_zscore;
   Q.(cname).Kinematics = Kins;

   [reach_ts, nonreach_ts, control_ts, stim_ts] = GetReachTimes_NP(inst_smooth, [100 100], Reaches);
   Q_frs.(cname).reach_ts = reach_ts;
   Q_frs.(cname).nonreach_ts = nonreach_ts;
   Q_frs.(cname).control_ts = control_ts;
   Q_frs.(cname).stim_ts = stim_ts;

   Events_Control = GetEvents_NP(Control);
   Events_Stim = GetEvents_NP(Stim);
   Q_frs.(cname).Events_Control = Events_Control;
   Q_frs.(cname).Events_Stim = Events_Stim;

   Q_reach.(cname).Control = Control;
   Q_reach.(cname).Stim = Stim;
end


%% Analyze Firing Rates and Modulation
offset = -45;  % Offset for alignment
dim = 5;       % Dimensionality parameter (not used explicitly in this part)

edges = [0 20; 20 40; 40 60; 60 80; 80 100;...33
         100 120; 120 140; 140 1500];

cn = fieldnames(Q);
FRBins = struct;

for c = 1:length(cn)
    fr = Q_frs.(cn{c}).inst_smooth;  % Get smoothed firing rates
    kins = Q.(cn{c}).Kinematics;     % Get kinematics
    time = round(kins(:, 1) * 1000); % Convert to ms for matching

    % Retrieve timestamps
    Ts.reach_ts = Q_frs.(cn{c}).reach_ts;      % Reach-related timestamps
    Ts.nonreach_ts = Q_frs.(cn{c}).nonreach_ts; % Non-reach timestamps

    % Ensure FRBins has the necessary fields
    FRBins.(cn{c}).control_ts.frs = [];
    FRBins.(cn{c}).nonreach_ts.frs = [];

    [x, index] = unique(time); % Get unique timestamps for interpolation
    pos = interp1(x, kins(index, 2:4), fr(:, 1));  % Interpolate position
    vel = interp1(x, kins(index, 7:9), fr(:, 1));  % Interpolate velocity
    kins_int = [fr(:, 1) pos vel];

    tn = fieldnames(Ts);
    for t = 1:length(tn)
        times = Ts.(tn{t});
        times_o = Ts.(tn{t}) + offset;

        % Find the closest available firing rate for each timestamp
        fr_values = interp1(fr(:, 1), fr(:, 2), times, 'nearest', NaN);
        t_kin_values = interp1(kins_int(:, 1), kins_int, times, 'nearest', NaN);

        % Store in FRBins
        FRBins.(cn{c}).(tn{t}).frs = fr_values;
        FRBins.(cn{c}).(tn{t}).kin_subt = t_kin_values;
    end

    % Assign firing rate data correctly
    FRBins.(cn{c}).control_ts.frs = FRBins.(cn{c}).reach_ts.frs(:);
    FRBins.(cn{c}).nonreach_ts.frs = FRBins.(cn{c}).nonreach_ts.frs(:);

    % Ensure data is non-empty before computing statistics
    if ~isempty(FRBins.(cn{c}).control_ts.frs) && ~isempty(FRBins.(cn{c}).nonreach_ts.frs)
        FR_reach = nanmean(FRBins.(cn{c}).control_ts.frs);
        FR_nonreach = nanmean(FRBins.(cn{c}).nonreach_ts.frs);

        % Compute modulation index (handle zero rates properly)
        if (FR_reach + FR_nonreach) > 0
            modulation_index = (FR_reach - FR_nonreach) / (FR_reach + FR_nonreach);
        else
            modulation_index = NaN; % Prevent division by zero
        end

        FRBins.(cn{c}).modulation_index = modulation_index;
        FRBins.(cn{c}).FR_reach = FR_reach;
        FRBins.(cn{c}).FR_nonreach = FR_nonreach;

        % Perform Wilcoxon rank-sum test if enough data exists
        if numel(FRBins.(cn{c}).control_ts.frs) > 1 && numel(FRBins.(cn{c}).nonreach_ts.frs) > 1
            p_value = ranksum(FRBins.(cn{c}).control_ts.frs, FRBins.(cn{c}).nonreach_ts.frs);
        else
            p_value = NaN;
        end
        FRBins.(cn{c}).p_value = p_value;

        % Set modulation flag with NaN check
        %FRBins.(cn{c}).modulated = ~isnan(modulation_index) && modulation_index >= 0.2;
        FRBins.(cn{c}).modulated = ~isnan(modulation_index) && modulation_index >= 0.1;

        % Compute effect size (Cohen's d)
        n_reach = numel(FRBins.(cn{c}).control_ts.frs);
        n_nonreach = numel(FRBins.(cn{c}).nonreach_ts.frs);
        pooled_std = sqrt(((n_reach - 1) * nanvar(FRBins.(cn{c}).control_ts.frs) + ...
                           (n_nonreach - 1) * nanvar(FRBins.(cn{c}).nonreach_ts.frs)) / ...
                           (n_reach + n_nonreach - 2));

        if pooled_std > 0
            effect_size = (FR_reach - FR_nonreach) / pooled_std; % Cohen's d
        else
            effect_size = NaN; % Prevent division by zero
        end
        FRBins.(cn{c}).effect_size = effect_size;
    else
        warning(['Missing data for cell: ', cn{c}, ' - Skipping computation.']);
    end
end

% Identify Significantly Modulated Cells
modulated_cells = {};
modulated_indices = [];
for c = 1:length(cn)
    if isfield(FRBins.(cn{c}), 'modulated') && FRBins.(cn{c}).modulated
        modulated_cells{end + 1} = cn{c};
        modulated_indices(end + 1) = c;
    end
end

% Display modulated cells
disp('Cells modulated by reaching:');
disp(modulated_cells);

length(modulated_cells)

%%
% Get NC indices
% Convert struct array to cell array for easier processing
neurontypes = {data.units.neurontype}; % Extract neuron types as a cell array

% Convert to string if needed (sometimes stored as char arrays)
neurontypes = string(neurontypes)'; 

% Find the indices where neurontype is exactly "NC" (ignore extra spaces)
nc_rows = find(strtrim(neurontypes) == "NC"); 

% Extract corresponding cluster_IDs (assuming that's the index you need)
nc_indices = [data.units(nc_rows).cluster_ID]; 

% Initialize array to hold filtered modulated cells
filtered_modulated_cells = {};
filtered_modulated_indices = [];

% Compare modulated indices with NC indices
for i = 1:length(modulated_indices)
    if ismember(data.units(modulated_indices(i)).cluster_ID, nc_indices)
        filtered_modulated_cells{end + 1} = modulated_cells{i}; % Store the name of the modulated cell
        filtered_modulated_indices(end + 1) = modulated_indices(i); % Store the index
    end
end

% Display results
disp('Modulated cells with NC label:');
disp(filtered_modulated_cells);
disp('Indices of modulated cells with NC label:');
disp(filtered_modulated_indices);
length(filtered_modulated_indices)
%% Plot all units together- looking at all good units 

epo = 500;
%epo = 250;

evnt = "pos_thresh";

good_units = false(1, length(data.units)); 

nc_rows=nc_rows';
%folder_path="C:\Users\PersonLab\Documents\NeuropixelsResults\ShamAnimals\050922_SAM009\mod cells";
%folder_path= "C:\Users\Courtney\Documents\NeuropixelsResults\Cell Catalog\051022_SAM012\Vmax_raster";
%for c = 1:length(data.units)
%for c = nc_rows
%for c=modulated_indices
for c=filtered_modulated_indices


    figure;hold on
   
    subplot(3,1,1); hold on
    dim=8;
    StimvelPlot(CN,dim,NoCSMcolor,0)
    StimvelPlot(ST0,dim,T0color,0)
    %StimvelPlot(ST1,dim,T1color,0)
    StimvelPlot(ST2,dim,T2color,0)
   % StimvelPlot(ST3,dim,T3color,0)
    %ylim([-1 2])
    ylim('auto')
    xlim([-.5,.5])
    yline(0,'--')
    % %xlim([-0.5 1])

  
    subplot(3,1,2); hold on
    h = title(['cluster' ' ' num2str(data.units(c).neurontype) ' ' 'depth' ' ' num2str(data.units(c).depth)]) ;
    set(h,'interpreter','none')

    % Plots Raster from binary matrix
    tick_length = 1; % length of tick marks for spikes
    spike_tick_width = 0.5; % width of lines marking spikes
    event_tick_width = 2; % width of lines marking secondary events

    %baseline NoRR\

    rastr = data.units(c).pos_thresh.Raster(:,(1000-epo):(2000-(1000-epo)));

    [y_co,x_co] = find(rastr);
    x_co = x_co - epo;
    line([x_co';x_co'],[y_co'-(tick_length/2);y_co'+(tick_length/2)],'Color', [0.5 0.5 0.5], 'LineWidth', spike_tick_width);
    line([zeros(1,size(rastr,1));zeros(1,size(rastr,1))],[(1:1:size(rastr,1))-.5;(1:1:size(rastr,1))+.5],'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);
    axis ij;
    ylim([0,size(rastr,1)])
    xlim([-epo,epo])


    subplot(3,1,3); hold on
    title('peth')

    %ctrl
    psth_smooth = data.units(c).pos_thresh.PSTH_smooth;
    psth_raw = data.units(c).pos_thresh.PSTH_raw;


    td1 = length(psth_raw)-length(psth_smooth)+1;
    td2 = td1/2;

    plot((-epo:epo),psth_smooth(((1000-epo)-td2):((2000-(1000-epo))-td2)),'Color',[0.5 0.5 0.5],'LineWidth',1.5)  
    xlim([-epo,epo])
%     ylim([0,100])

    xlabel('Time from Threshold Crossing')
    ylabel('FR (Hz)')
    % filename=sprintf('posthresh_raster%d.svg',c);
    % %filename=sprintf('vmax_raster%d.svg',c);
    % full_path=fullfile(folder_path,filename);
    % saveas(gcf,full_path);
    % 
    % close;
% % 
if mod(c,10) == 0
    drawnow
    input("press enter when done looking at figures")
    close all
 end
end
%   drawnow
%     response = input(['Mark unit ' num2str(c) ' as good? (y/n): '], 's');
% 
%     if strcmpi(response, 'y')
%         good_units(c) = true; % Mark this unit as good
%     end
% 
%     close all % Close the figure after input
% end

% Extract the indices of good units
% good_unit_indices = find(good_units);
% 
% % Display the selected good units
% disp("Good Units:");
% disp(good_unit_indices);


%% Now need to get spike counts 
% Define parameters for baseline and event duration
baseline_epoch = 500; % milliseconds before the event
event_epoch = 500;     % milliseconds after the event
bin_size = 1;         % size of time bins in milliseconds

% Choose event to align to
event_name = ["csmT2"];
enames = fieldnames(Events_NoCSM);
evnt = find(ismember(enames, event_name));

% Loop through each unit
for c = 1:length(data.units)
    spike_times = data.units(c).spike_times * 1000; % Convert to ms
    
    for enum = 1:length(evnt)
        event_times = Events_NoCSM.(enames{evnt(enum)}) * 1000; % Get event times
        Trial_spike_times = []; % Initialize trial spike times structure
        
        % Initialize spike counts array
        num_trials = length(event_times);
        spike_counts = zeros(num_trials, ceil((baseline_epoch + event_epoch) / bin_size));
        
        % Get spike times aligned to the events
        reach_names = 'R%d';
        reach_count = 1;
        
        for event = 1:length(event_times)
            ref_time = event_times(event);
            spk_times_inds = (spike_times >= (ref_time - baseline_epoch)) & (spike_times <= (ref_time + event_epoch));
            trial_spike_times = spike_times(spk_times_inds) - ref_time; 
            reach_name = sprintf(reach_names, reach_count);
            reach_count = reach_count + 1;
            Trial_spike_times.(reach_name) = trial_spike_times; % Store trial spike times
            
            % Calculate spike counts in time bins
            for b = 1:ceil((baseline_epoch + event_epoch) / bin_size)
                bin_start = (b - 1) * bin_size - baseline_epoch; % Calculate bin start time
                bin_end = bin_start + bin_size; % Calculate bin end time
                spike_counts(event, b) = sum(trial_spike_times >= bin_start & trial_spike_times < bin_end);
            end
        end
        
        % Store trial spike times and spike counts
        data.units(c).SpikeTimes_NoCSMT2(enum).trial_spktimes = Trial_spike_times;
        data.units(c).SpikeTimes_NoCSMT2(enum).spike_counts = spike_counts; % Store spike counts for this event
    end
end

%% Now need to get spike counts 
% Define parameters for baseline and event duration
baseline_epoch = 500; % milliseconds before the event
event_epoch = 500;     % milliseconds after the event
bin_size = 1;         % size of time bins in milliseconds

% Choose event to align to
event_name = ["csmT0"];
enames = fieldnames(Events_NoCSM);
evnt = find(ismember(enames, event_name));

% Loop through each unit
for c = 1:length(data.units)
    spike_times = data.units(c).spike_times * 1000; % Convert to ms
    
    for enum = 1:length(evnt)
        event_times = Events_NoCSM.(enames{evnt(enum)}) * 1000; % Get event times
        Trial_spike_times = []; % Initialize trial spike times structure
        
        % Initialize spike counts array
        num_trials = length(event_times);
        spike_counts = zeros(num_trials, ceil((baseline_epoch + event_epoch) / bin_size));
        
        % Get spike times aligned to the events
        reach_names = 'R%d';
        reach_count = 1;
        
        for event = 1:length(event_times)
            ref_time = event_times(event);
            spk_times_inds = (spike_times >= (ref_time - baseline_epoch)) & (spike_times <= (ref_time + event_epoch));
            trial_spike_times = spike_times(spk_times_inds) - ref_time; 
            reach_name = sprintf(reach_names, reach_count);
            reach_count = reach_count + 1;
            Trial_spike_times.(reach_name) = trial_spike_times; % Store trial spike times
            
            % Calculate spike counts in time bins
            for b = 1:ceil((baseline_epoch + event_epoch) / bin_size)
                bin_start = (b - 1) * bin_size - baseline_epoch; % Calculate bin start time
                bin_end = bin_start + bin_size; % Calculate bin end time
                spike_counts(event, b) = sum(trial_spike_times >= bin_start & trial_spike_times < bin_end);
            end
        end
        
        % Store trial spike times and spike counts
        data.units(c).SpikeTimes_NoCSMT0(enum).trial_spktimes = Trial_spike_times;
        data.units(c).SpikeTimes_NoCSMT0(enum).spike_counts = spike_counts; % Store spike counts for this event
    end
end
%% Now need to get spike counts 
% Define parameters for baseline and event duration
baseline_epoch = 500; % milliseconds before the event
event_epoch = 500;     % milliseconds after the event
bin_size = 1;         % size of time bins in milliseconds

% Choose event to align to
event_name = ["csmT2"];
enames = fieldnames(Events_T2);
evnt = find(ismember(enames, event_name));

% Loop through each unit
for c = 1:length(data.units)
    spike_times = data.units(c).spike_times * 1000; % Convert to ms
    
    for enum = 1:length(evnt)
        event_times = Events_T2.(enames{evnt(enum)}) * 1000; % Get event times
        Trial_spike_times = []; % Initialize trial spike times structure
        
        % Initialize spike counts array
        num_trials = length(event_times);
        spike_counts = zeros(num_trials, ceil((baseline_epoch + event_epoch) / bin_size));
        
        % Get spike times aligned to the events
        reach_names = 'R%d';
        reach_count = 1;
        
        for event = 1:length(event_times)
            ref_time = event_times(event);
            spk_times_inds = (spike_times >= (ref_time - baseline_epoch)) & (spike_times <= (ref_time + event_epoch));
            trial_spike_times = spike_times(spk_times_inds) - ref_time; 
            reach_name = sprintf(reach_names, reach_count);
            reach_count = reach_count + 1;
            Trial_spike_times.(reach_name) = trial_spike_times; % Store trial spike times
            
            % Calculate spike counts in time bins
            for b = 1:ceil((baseline_epoch + event_epoch) / bin_size)
                bin_start = (b - 1) * bin_size - baseline_epoch; % Calculate bin start time
                bin_end = bin_start + bin_size; % Calculate bin end time
                spike_counts(event, b) = sum(trial_spike_times >= bin_start & trial_spike_times < bin_end);
            end
        end
        
        % Store trial spike times and spike counts
        data.units(c).SpikeTimes_T2(enum).trial_spktimes = Trial_spike_times;
        data.units(c).SpikeTimes_T2(enum).spike_counts = spike_counts; % Store spike counts for this event
    end
end

%% Now need to get spike counts 
% Define parameters for baseline and event duration
baseline_epoch = 500; % milliseconds before the event
event_epoch = 500;     % milliseconds after the event
bin_size = 1;         % size of time bins in milliseconds

% Choose event to align to
event_name = ["csmT0"];
enames = fieldnames(Events_T0);
evnt = find(ismember(enames, event_name));

% Loop through each unit
for c = 1:length(data.units)
    spike_times = data.units(c).spike_times * 1000; % Convert to ms
    
    for enum = 1:length(evnt)
        event_times = Events_T0.(enames{evnt(enum)}) * 1000; % Get event times
        Trial_spike_times = []; % Initialize trial spike times structure
        
        % Initialize spike counts array
        num_trials = length(event_times);
        spike_counts = zeros(num_trials, ceil((baseline_epoch + event_epoch) / bin_size));
        
        % Get spike times aligned to the events
        reach_names = 'R%d';
        reach_count = 1;
        
        for event = 1:length(event_times)
            ref_time = event_times(event);
            spk_times_inds = (spike_times >= (ref_time - baseline_epoch)) & (spike_times <= (ref_time + event_epoch));
            trial_spike_times = spike_times(spk_times_inds) - ref_time; 
            reach_name = sprintf(reach_names, reach_count);
            reach_count = reach_count + 1;
            Trial_spike_times.(reach_name) = trial_spike_times; % Store trial spike times
            
            % Calculate spike counts in time bins
            for b = 1:ceil((baseline_epoch + event_epoch) / bin_size)
                bin_start = (b - 1) * bin_size - baseline_epoch; % Calculate bin start time
                bin_end = bin_start + bin_size; % Calculate bin end time
                spike_counts(event, b) = sum(trial_spike_times >= bin_start & trial_spike_times < bin_end);
            end
        end
        
        % Store trial spike times and spike counts
        data.units(c).SpikeTimes_T0(enum).trial_spktimes = Trial_spike_times;
        data.units(c).SpikeTimes_T0(enum).spike_counts = spike_counts; % Store spike counts for this event
    end
end
%% Use stat test to see if cells are modulated during CSM -looking at difference between NoCSM and T2
% Define post-CSM window
% post_csm_start = 0;   % 0 ms
% post_csm_end = 200;   % 200 ms
post_csm_start = -50;   % 0 ms
post_csm_end = 150;   % 200 ms
% Define time vector for 25 ms bins from -500 to 500 ms
%time_vector = -500:10:500;  
time_vector = -500:1:500; 
post_csm_indices = find(time_vector >= post_csm_start & time_vector <= post_csm_end);

% Initialize array to store indices of significantly modulated cells
csm_cells = [];

% Loop through all units
for c = filtered_modulated_indices
%for c=20
%for c =1:length(data.units)

    % Get spike counts for NoCSM and CSM trials in the post-CSM window
    post_csm_spike_counts_NoCSM = sum(data.units(c).SpikeTimes_NoCSMT2.spike_counts(:, post_csm_indices), 2);
    post_csm_spike_counts_CSM = sum(data.units(c).SpikeTimes_T2.spike_counts(:, post_csm_indices), 2);

    % Check for normality using Lilliefors test (alternative to swtest)
    [h_norm_NCSM, p_norm_NCSM] = lillietest(post_csm_spike_counts_NoCSM);
    [h_norm_CSM, p_norm_CSM] = lillietest(post_csm_spike_counts_CSM);
    
    % Non-normal distributions -> Use Wilcoxon rank-sum test
    [p, h] = ranksum(post_csm_spike_counts_CSM, post_csm_spike_counts_NoCSM);
    test_type = 'Wilcoxon rank-sum test';

    % Store results
    data.units(c).PostCSMT2_Stats.p_value = p;
    data.units(c).PostCSMT2_Stats.significant = h;
    data.units(c).PostCSMT2_Stats.test_type = test_type;

    % Store index if significant
    if p < 0.05
        csm_cells = [csm_cells, c]; 
    end

    % Display results
    disp(['Unit ', num2str(c), ': ', test_type, ' p = ', num2str(p)]);
end

% Display the indices of significant cells
disp('Indices of significantly modulated cells:');
disp(csm_cells);
%% other way to look for csm cells - GLM which seems to work better
% Use GLM to identify significantly modulated cells during CSM
% Define post-CSM window
post_csm_start = -50;   
post_csm_end = 150;   

% Define time vector for 1 ms bins from -500 to 500 ms
time_vector = -500:1:500;  
post_csm_indices = find(time_vector >= post_csm_start & time_vector <= post_csm_end);

% Initialize array to store indices of significantly modulated cells
csm_cells = [];

% Loop through all units
for c = filtered_modulated_indices

    % Get spike counts for NoCSM and CSM trials in the post-CSM window
    post_csm_spike_counts_NoCSM = sum(data.units(c).SpikeTimes_NoCSMT2.spike_counts(:, post_csm_indices), 2);
    post_csm_spike_counts_CSM = sum(data.units(c).SpikeTimes_T2.spike_counts(:, post_csm_indices), 2);

    % Combine spike counts into one vector
    spike_counts = [post_csm_spike_counts_CSM; post_csm_spike_counts_NoCSM];

    % Create predictor variable (1 = CSM, 0 = NoCSM)
    conditions = [ones(size(post_csm_spike_counts_CSM)); zeros(size(post_csm_spike_counts_NoCSM))];

    % Fit GLM (Poisson regression)
    glm_model = fitglm(conditions, spike_counts, 'linear', 'Distribution', 'poisson');

    % Extract p-value for condition effect
    p_value = glm_model.Coefficients.pValue(2);  % Second coefficient is for condition (CSM vs NoCSM)

    % Store results
    data.units(c).PostCSMT2_Stats.p_value_glm = p_value;
    data.units(c).PostCSMT2_Stats.significant_glm = p_value < 0.05;
    data.units(c).PostCSMT2_Stats.test_type_glm = 'Poisson GLM';

    % Store index if significant
    if p_value < 0.05
        csm_cells = [csm_cells, c]; 
    end

    % Display results
    disp(['Unit ', num2str(c), ': Poisson GLM p = ', num2str(p_value)]);
end

% Display the indices of significant cells
disp('Indices of significantly modulated cells (GLM):');
disp(csm_cells);

%% other way to look for csm cells - T0 vs NOCSM
% Use GLM to identify significantly modulated cells during CSM
% Define post-CSM window
post_csm_start = -50;   
post_csm_end = 150;   

% Define time vector for 1 ms bins from -500 to 500 ms
time_vector = -500:1:500;  
post_csm_indices = find(time_vector >= post_csm_start & time_vector <= post_csm_end);

% Initialize array to store indices of significantly modulated cells
csm_cells = [];

% Loop through all units
for c = filtered_modulated_indices

    % Get spike counts for NoCSM and CSM trials in the post-CSM window
    post_csm_spike_counts_NoCSM = sum(data.units(c).SpikeTimes_NoCSMT0.spike_counts(:, post_csm_indices), 2);
    post_csm_spike_counts_CSM = sum(data.units(c).SpikeTimes_T0.spike_counts(:, post_csm_indices), 2);

    % Combine spike counts into one vector
    spike_counts = [post_csm_spike_counts_CSM; post_csm_spike_counts_NoCSM];

    % Create predictor variable (1 = CSM, 0 = NoCSM)
    conditions = [ones(size(post_csm_spike_counts_CSM)); zeros(size(post_csm_spike_counts_NoCSM))];

    % Fit GLM (Poisson regression)
    glm_model = fitglm(conditions, spike_counts, 'linear', 'Distribution', 'poisson');

    % Extract p-value for condition effect
    p_value = glm_model.Coefficients.pValue(2);  % Second coefficient is for condition (CSM vs NoCSM)

    % Store results
    data.units(c).PostCSMT0_Stats.p_value_glm = p_value;
    data.units(c).PostCSMT0_Stats.significant_glm = p_value < 0.05;
    data.units(c).PostCSMT0_Stats.test_type_glm = 'Poisson GLM';

    % Store index if significant
    if p_value < 0.05
        csm_cells = [csm_cells, c]; 
    end

    % Display results
    disp(['Unit ', num2str(c), ': Poisson GLM p = ', num2str(p_value)]);
end

% Display the indices of significant cells
disp('Indices of significantly modulated cells (GLM):');
disp(csm_cells);



%% T2 vs NoCSM 
% Loop through each unit to calculate selectivity and plot
folder_path="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\ShamAnimals\032322_SAM005\NewAnalysis\T0CSMaligned";

for c = filtered_modulated_indices
%for c = csm_cells
%for c=15:23
    for enum = 1:length(evnt)
        % Extract PSTH for CSM and NoCSM conditions
        psth_NoCSM_smooth = data.units(c).NoCSM_csmT0(enum).PSTH_smoothbin; % Assuming this is already smoothed
        psth_CSM_smooth = data.units(c).T0_csm(enum).PSTH_smoothbin; % Assuming this is already smoothed

        % Calculate the duration based on the PSTH length
        total_duration = length(psth_NoCSM_smooth) * 10; % Assuming 10 ms bin size
        time_vector = linspace(-1000, 1000, length(psth_NoCSM_smooth)); % Adjust based on window size

        % Create a Gaussian filter
        smoother = 10; % Standard deviation of the Gaussian
        gauss_filter = fspecial('gaussian', [1, 6 * smoother], smoother);

        % Smooth the PSTHs using the Gaussian filter
        psth_NoCSM_gaussian = conv(psth_NoCSM_smooth, gauss_filter, 'same');
        psth_CSM_gaussian = conv(psth_CSM_smooth, gauss_filter, 'same');

        % Calculate selectivity
        selectivity = abs(psth_CSM_gaussian - psth_NoCSM_gaussian) ./ max(abs(psth_CSM_gaussian), abs(psth_NoCSM_gaussian));
        selectivity(isnan(selectivity)) = 0; % Handle potential NaNs if both are zero

        % Smooth the selectivity using Gaussian filter
        selectivity_smooth = conv(selectivity, gauss_filter, 'same'); % Using Gaussian smoothing

        % Plot PSTHs
        figure('Position', [100, 100, 800, 900]); % Adjusts the figure size (width x height)
       
        subplot(3,1,1); hold on
        dim=8;
      
        StimvelPlot(NoRR_T0aligned,dim,NoCSMcolor,0)
        StimvelPlot(ST0_CSM,dim,T0color,0)
    
        ylim('auto')
        xlim([-.5,.5])
        yline(0,'--')

      
        subplot(3, 1, 2);
        plot(time_vector, psth_NoCSM_gaussian,'Color',NoCSMcolor, 'LineWidth', 1.5); % Smooth PSTH
        hold on;
        plot(time_vector, psth_CSM_gaussian,'Color',T0color, 'LineWidth', 1.5); % Smooth PSTH
        xlabel('Time (ms)');
        ylabel('Firing Rate (Hz)');
        title(['PSTH for Unit ', num2str(c), ' - Event ', num2str(enum)]);
        %legend('No CSM', 'CSM');
        xline(0,'--k','LineWidth',2)
        xline(200,'--k','LineWidth',2)
        xlim([-500 500]); % Change to view from -500 to 500 ms
        grid off;

        % Plot Selectivity
        subplot(3, 1, 3);
        plot(time_vector, selectivity_smooth, 'k', 'LineWidth', 1.5); % Smooth selectivity
        xlabel('Time (ms)');
        ylabel('Selectivity');
        title(['Selectivity for Unit ', num2str(c), ' - Event ', num2str(enum)]);
        xlim([-500 500]); % Change to view from -500 to 500 ms
        xline(0,'--k','LineWidth',2)
        xline(200,'--k','LineWidth',2)
        grid off;

        % Adjust figure layout
        sgtitle(['Unit ', num2str(c), ' - Event ', num2str(enum)]);
        % 
          filename=sprintf('psth_selective%d.svg',c);
        full_path=fullfile(folder_path,filename);
        saveas(gcf,full_path);

        close;
    end
end
%% T2 vs NoCSM
% Loop through each unit to calculate selectivity and plot
folder_path="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\ShamAnimals\032322_SAM005\NewAnalysis\T2CSMaligned";

for c = filtered_modulated_indices
%for c = csm_cells
%for c=15:23
    for enum = 1:length(evnt)
        % Extract PSTH for CSM and NoCSM conditions
        psth_NoCSM_smooth = data.units(c).NoCSM_csmT2(enum).PSTH_smoothbin; % Assuming this is already smoothed
        psth_CSM_smooth = data.units(c).T2_csm(enum).PSTH_smoothbin; % Assuming this is already smoothed

        % Calculate the duration based on the PSTH length
        total_duration = length(psth_NoCSM_smooth) * 10; % Assuming 10 ms bin size
        time_vector = linspace(-1000, 1000, length(psth_NoCSM_smooth)); % Adjust based on window size

        % Create a Gaussian filter
        smoother = 10; % Standard deviation of the Gaussian
        gauss_filter = fspecial('gaussian', [1, 6 * smoother], smoother);

        % Smooth the PSTHs using the Gaussian filter
        psth_NoCSM_gaussian = conv(psth_NoCSM_smooth, gauss_filter, 'same');
        psth_CSM_gaussian = conv(psth_CSM_smooth, gauss_filter, 'same');

        % Calculate selectivity
        selectivity = abs(psth_CSM_gaussian - psth_NoCSM_gaussian) ./ max(abs(psth_CSM_gaussian), abs(psth_NoCSM_gaussian));
        selectivity(isnan(selectivity)) = 0; % Handle potential NaNs if both are zero

        % Smooth the selectivity using Gaussian filter
        selectivity_smooth = conv(selectivity, gauss_filter, 'same'); % Using Gaussian smoothing

        % Plot PSTHs
        figure('Position', [100, 100, 800, 900]); % Adjusts the figure size (width x height)
       
        subplot(3,1,1); hold on
        dim=8;
      
        StimvelPlot(NoRR_T2aligned,dim,NoCSMcolor,0)
        StimvelPlot(ST2_CSM,dim,T2color,0)
    
        ylim('auto')
        xlim([-.5,.5])
        yline(0,'--')

      
        subplot(3, 1, 2);
        plot(time_vector, psth_NoCSM_gaussian,'Color',NoCSMcolor, 'LineWidth', 1.5); % Smooth PSTH
        hold on;
        plot(time_vector, psth_CSM_gaussian,'Color',T2color, 'LineWidth', 1.5); % Smooth PSTH
        xlabel('Time (ms)');
        ylabel('Firing Rate (Hz)');
        title(['PSTH for Unit ', num2str(c), ' - Event ', num2str(enum)]);
        %legend('No CSM', 'CSM');
        xline(0,'--k','LineWidth',2)
        xline(200,'--k','LineWidth',2)
        xlim([-500 500]); % Change to view from -500 to 500 ms
        grid off;

        % Plot Selectivity
        subplot(3, 1, 3);
        plot(time_vector, selectivity_smooth, 'k', 'LineWidth', 1.5); % Smooth selectivity
        xlabel('Time (ms)');
        ylabel('Selectivity');
        title(['Selectivity for Unit ', num2str(c), ' - Event ', num2str(enum)]);
        xlim([-500 500]); % Change to view from -500 to 500 ms
        xline(0,'--k','LineWidth',2)
        xline(200,'--k','LineWidth',2)
        grid off;

        % Adjust figure layout
        sgtitle(['Unit ', num2str(c), ' - Event ', num2str(enum)]);
        % 
          filename=sprintf('psth_selective%d.svg',c);
        full_path=fullfile(folder_path,filename);
        saveas(gcf,full_path);

        close;
    end
end
 
%% Another way to calculate 
% Loop through each unit to calculate selectivity and plot

for c = filtered_modulated_indices
    for enum = 1:length(evnt)
        % Extract Binary Master PSTHs for standard error calculation
        BinaryMaster_NoCSM = data.units(c).NoCSM_csmT2(enum).Raster;
        BinaryMaster_CSM = data.units(c).T2_csm(enum).Raster;
        
        % Compute mean PSTH
        psth_NoCSM_mean = mean(BinaryMaster_NoCSM, 1) * 1000; % Convert to Hz
        psth_CSM_mean = mean(BinaryMaster_CSM, 1) * 1000; % Convert to Hz
        
        % Compute standard error (SE = std/sqrt(N))
        SE_NoCSM = std(BinaryMaster_NoCSM, 0, 1) * 1000 / sqrt(size(BinaryMaster_NoCSM, 1));
        SE_CSM = std(BinaryMaster_CSM, 0, 1) * 1000 / sqrt(size(BinaryMaster_CSM, 1));

        % Generate Gaussian filter for smoothing
        smoothing_bins = 5;
        gauss_filter = fspecial('gaussian', [1, smoothing_bins * 6], smoothing_bins);

        % Smooth PSTH using Gaussian filter
        psth_NoCSM_smooth = conv2(psth_NoCSM_mean, gauss_filter, 'same');
        psth_CSM_smooth = conv2(psth_CSM_mean, gauss_filter, 'same');
        
        % Smooth standard error
        SE_NoCSM_smooth = conv2(SE_NoCSM, gauss_filter, 'same');
        SE_CSM_smooth = conv2(SE_CSM, gauss_filter, 'same');

        % Compute selectivity
        selectivity = abs(psth_CSM_smooth - psth_NoCSM_smooth) ./ max(abs(psth_CSM_smooth), abs(psth_NoCSM_smooth));
        selectivity(isnan(selectivity)) = 0;

        % Smooth selectivity
        selectivity_smooth = conv2(selectivity, gauss_filter, 'same');

        % Time vector
        time_vector = linspace(-1000, 1000, length(psth_NoCSM_smooth));

        % Plot PSTHs with standard error
        figure;
        subplot(2, 1, 1);
        hold on;
        fill([time_vector, fliplr(time_vector)], ...
             [psth_NoCSM_smooth + SE_NoCSM_smooth, fliplr(psth_NoCSM_smooth - SE_NoCSM_smooth)], ...
             'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Shaded error for NoCSM
        fill([time_vector, fliplr(time_vector)], ...
             [psth_CSM_smooth + SE_CSM_smooth, fliplr(psth_CSM_smooth - SE_CSM_smooth)], ...
             'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Shaded error for CSM
        plot(time_vector, psth_NoCSM_smooth, 'b', 'LineWidth', 1.5);
        plot(time_vector, psth_CSM_smooth, 'r', 'LineWidth', 1.5);
        xlabel('Time (ms)');
        ylabel('Firing Rate (Hz)');
        title(['PSTH for Unit ', num2str(c), ' - Event ', num2str(enum)]);
        legend('No CSM', 'CSM');
        xlim([-500 500]); % View from -500 to 500 ms
        grid on;

        % Plot Selectivity
        subplot(2, 1, 2);
        plot(time_vector, selectivity_smooth, 'k', 'LineWidth', 1.5);
        xlabel('Time (ms)');
        ylabel('Selectivity');
        title(['Selectivity for Unit ', num2str(c), ' - Event ', num2str(enum)]);
        xlim([-500 500]);
        grid on;

        % Adjust figure layout
        sgtitle(['Unit ', num2str(c), ' - Event ', num2str(enum)]);
    end
  
end

%% another way to calculate 
% Loop through each unit to calculate selectivity and plot
for c = filtered_modulated_indices
    for enum = 1:length(evnt)
        % Extract Binary Master PSTHs for standard error calculation
        BinaryMaster_NoCSM = data.units(c).NoCSM_csmT2(enum).Raster;
        BinaryMaster_CSM = data.units(c).T2_csm(enum).Raster;
        
        % Compute mean PSTH (in Hz)
        psth_NoCSM_mean = mean(BinaryMaster_NoCSM, 1) * 1000;
        psth_CSM_mean = mean(BinaryMaster_CSM, 1) * 1000;
        
        % Compute SEM = std/sqrt(N)
        SEM_NoCSM = std(BinaryMaster_NoCSM, 0, 1) * 1000 / sqrt(size(BinaryMaster_NoCSM, 1));
        SEM_CSM = std(BinaryMaster_CSM, 0, 1) * 1000 / sqrt(size(BinaryMaster_CSM, 1));

        % Generate Gaussian filter for smoothing
        smoothing_bins = 5; % Adjust for optimal smoothing
        gauss_filter = fspecial('gaussian', [1, smoothing_bins * 3], smoothing_bins);

        % Smooth PSTH using Gaussian filter
        psth_NoCSM_smooth = conv2(psth_NoCSM_mean, gauss_filter, 'same');
        psth_CSM_smooth = conv2(psth_CSM_mean, gauss_filter, 'same');

        % Lightly smooth SEM (use smaller kernel for SEM)
        SEM_NoCSM_smooth = conv2(SEM_NoCSM, fspecial('gaussian', [1, smoothing_bins], smoothing_bins / 2), 'same');
        SEM_CSM_smooth = conv2(SEM_CSM, fspecial('gaussian', [1, smoothing_bins], smoothing_bins / 2), 'same');

        % Compute selectivity
        selectivity = abs(psth_CSM_smooth - psth_NoCSM_smooth) ./ max(abs(psth_CSM_smooth), abs(psth_NoCSM_smooth));
        selectivity(isnan(selectivity)) = 0;

        % Smooth selectivity
        selectivity_smooth = conv2(selectivity, gauss_filter, 'same');

        % Time vector
        time_vector = linspace(-1000, 1000, length(psth_NoCSM_smooth));

        % Plot PSTHs with SEM as shaded error
        figure;
        subplot(2, 1, 1);
        hold on;
        fill([time_vector, fliplr(time_vector)], ...
             [psth_NoCSM_smooth + SEM_NoCSM_smooth, fliplr(psth_NoCSM_smooth - SEM_NoCSM_smooth)], ...
             'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Lighter SEM shading
        fill([time_vector, fliplr(time_vector)], ...
             [psth_CSM_smooth + SEM_CSM_smooth, fliplr(psth_CSM_smooth - SEM_CSM_smooth)], ...
             'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Lighter SEM shading
        plot(time_vector, psth_NoCSM_smooth, 'b', 'LineWidth', 2);
        plot(time_vector, psth_CSM_smooth, 'r', 'LineWidth', 2);
        xlabel('Time (ms)');
        ylabel('Firing Rate (Hz)');
        title(['PSTH for Unit ', num2str(c), ' - Event ', num2str(enum)]);
        legend('No CSM', 'CSM');
        xlim([-500 500]); % Focus on -500 to 500 ms
        grid on;

        % Plot Selectivity
        subplot(2, 1, 2);
        plot(time_vector, selectivity_smooth, 'k', 'LineWidth', 2);
        xlabel('Time (ms)');
        ylabel('Selectivity');
        title(['Selectivity for Unit ', num2str(c), ' - Event ', num2str(enum)]);
        xlim([-500 500]);
        grid on;

        % Adjust figure layout
        sgtitle(['Unit ', num2str(c), ' - Event ', num2str(enum)]);
    end
end

%% now plot a heat map
% Define parameters
num_units = length(csm_cells);
%num_units = length(filtered_modulated_indices);
time_vector = linspace(-1000, 1000, size(data.units(1).NoCSM_csmT2(1).Raster, 2));

% Preallocate matrix to store selectivity for each unit
selectivity_matrix = zeros(num_units, length(time_vector));

% Smoothing parameters
smoothing_bins = 5;
gauss_filter = fspecial('gaussian', [1, smoothing_bins * 3], smoothing_bins);

% Loop through each unit
for unit_idx = 1:num_units
    c = csm_cells(unit_idx); % Get unit ID
    %c = filtered_modulated_indices(unit_idx); % Get unit ID
    
    % Extract Binary Master PSTHs
    BinaryMaster_NoCSM = data.units(c).NoCSM_csmT2(1).Raster; % Trials x Time
    BinaryMaster_CSM = data.units(c).T2_csm(1).Raster; % Trials x Time

    % Compute mean PSTH (convert to Hz)
    psth_NoCSM_mean = mean(BinaryMaster_NoCSM, 1) * 1000;
    psth_CSM_mean = mean(BinaryMaster_CSM, 1) * 1000;

    % Apply Gaussian smoothing
    psth_NoCSM_smooth = conv2(psth_NoCSM_mean, gauss_filter, 'same');
    psth_CSM_smooth = conv2(psth_CSM_mean, gauss_filter, 'same');

    % Compute selectivity
    selectivity = abs(psth_CSM_smooth - psth_NoCSM_smooth) ./ max(abs(psth_CSM_smooth), abs(psth_NoCSM_smooth));
    selectivity(isnan(selectivity)) = 0; % Handle NaNs

    % Normalize selectivity (Z-score)
    selectivity_z = (selectivity - mean(selectivity)) / std(selectivity);

    % Store in matrix (each row = one unit)
    selectivity_matrix(unit_idx, :) = selectivity_z;
end

% Plot Selectivity Heatmap
figure;
imagesc(time_vector, 1:num_units, selectivity_matrix);
colormap(parula); % Heatmap colormap
colorbar;
xlabel('Time (ms)');
ylabel('Unit Index');
title('Normalized Selectivity Heatmap (Z-score)');
xlim([-500 500]);
xline(0,'--r','LineWidth',2)
xline(200,'--r','LineWidth',2)
set(gca, 'YDir', 'normal'); % Ensure proper axis orientation

%% other ways 
% Define parameters
num_units = length(csm_cells);
time_vector = linspace(-1000, 1000, size(data.units(1).NoCSM_csmT2(1).Raster, 2));

% Preallocate matrix to store selectivity for each unit
selectivity_matrix = zeros(num_units, length(time_vector));

% Smoothing parameters
smoothing_bins = 5;
gauss_filter = fspecial('gaussian', [1, smoothing_bins * 3], smoothing_bins);

% Loop through each unit
for unit_idx = 1:num_units
    c = csm_cells(unit_idx); % Get unit ID
    
    % Extract Binary Master PSTHs
    BinaryMaster_NoCSM = data.units(c).NoCSM_csmT2(1).Raster; % Trials x Time
    BinaryMaster_CSM = data.units(c).T2_csm(1).Raster; % Trials x Time

    % Compute mean PSTH (convert to Hz)
    psth_NoCSM_mean = mean(BinaryMaster_NoCSM, 1) * 1000;
    psth_CSM_mean = mean(BinaryMaster_CSM, 1) * 1000;

    % Apply Gaussian smoothing
    psth_NoCSM_smooth = conv2(psth_NoCSM_mean, gauss_filter, 'same');
    psth_CSM_smooth = conv2(psth_CSM_mean, gauss_filter, 'same');

    % Compute selectivity
    selectivity = abs(psth_CSM_smooth - psth_NoCSM_smooth) ./ max(abs(psth_CSM_smooth), abs(psth_NoCSM_smooth));
    selectivity(isnan(selectivity)) = 0; % Handle NaNs

    % Store in matrix (each row = one unit)
    selectivity_matrix(unit_idx, :) = selectivity;
end

% Plot Selectivity Heatmap
figure;
imagesc(time_vector, 1:num_units, selectivity_matrix);
colormap(parula); % Heatmap colormap
colorbar;
xlabel('Time (ms)');
ylabel('Unit Index');
title('Selectivity Heatmap (Each Row = One Unit)');
xlim([-500 500]);
set(gca, 'YDir', 'normal'); % Ensure proper axis orientation

%% Generate example raster trace 
% Plot rasters aligned to each one 


epo = 500;
%epo = 250;

evnt = "pos_thresh";


%folder_path="C:\Users\PersonLab\Documents\NeuropixelsResults\ShamAnimals\050922_SAM009\mod cells";

for c=filtered_modulated_indices
%for c=modulated_indices


    figure;hold on
    figure('Position', [100, 100, 800, 900]); hold on
   
    subplot(3,1,1); hold on
    h = title(['PSTH for Unit ', num2str(c), ' ' 'depth' ' ' num2str(data.units(c).depth)]) ;
    set(h,'interpreter','none')

    dim=8;
    StimvelPlot(NoRR_T0aligned,dim,NoCSMcolor,0)
    StimvelPlot(ST0_CSM,dim,T0color,0)
    %StimvelPlot(ST1,dim,T1color,0)
    %StimvelPlot(ST2,dim,T2color,0)
   % StimvelPlot(ST3,dim,T3color,0)
    %ylim([-1 2])
    xlabel('CSM Aligned-T0')
    ylim('auto')
    xlim([-.5,.5])
    yline(0,'--')
    % %xlim([-0.5 1])

  
    subplot(3,1,2); hold on
    title('NoCSM Aligned to T0')

    % Plots Raster from binary matrix
    tick_length = 1; % length of tick marks for spikes
    spike_tick_width = 0.5; % width of lines marking spikes
    event_tick_width = 2; % width of lines marking secondary events

    %baseline NoRR\

    rastr = data.units(c).NoCSM_csmT0.Raster(:,(1000-epo):(2000-(1000-epo)));

    [y_co,x_co] = find(rastr);
    x_co = x_co - epo;
    line([x_co';x_co'],[y_co'-(tick_length/2);y_co'+(tick_length/2)],'Color', [0.5 0.5 0.5], 'LineWidth', spike_tick_width);
    line([zeros(1,size(rastr,1));zeros(1,size(rastr,1))],[(1:1:size(rastr,1))-.5;(1:1:size(rastr,1))+.5],'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);
    axis ij;
    ylim([0,size(rastr,1)])
    xlim([-epo,epo])

    subplot(3,1,3); hold on
 
    title('T0_csm')
    % Plots Raster from binary matrix
    tick_length = 1; % length of tick marks for spikes
    spike_tick_width = 0.5; % width of lines marking spikes
    event_tick_width = 2; % width of lines marking secondary events

    %baseline NoRR\

    rastr = data.units(c).T0_csm.Raster(:,(1000-epo):(2000-(1000-epo)));

    [y_co,x_co] = find(rastr);
    x_co = x_co - epo;
    line([x_co';x_co'],[y_co'-(tick_length/2);y_co'+(tick_length/2)],'Color', [0.5 0.5 0.5], 'LineWidth', spike_tick_width);
    line([zeros(1,size(rastr,1));zeros(1,size(rastr,1))],[(1:1:size(rastr,1))-.5;(1:1:size(rastr,1))+.5],'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);
    axis ij;
    ylim([0,size(rastr,1)])
    xlim([-epo,epo])



if mod(c,10) == 0
    drawnow
    input("press enter when done looking at figures")
    close all
 end
end
%   drawnow
%     response = input(['Mark unit ' num2str(c) ' as good? (y/n): '], 's');
% 
%     if strcmpi(response, 'y')
%         good_units(c) = true; % Mark this unit as good
%     end
% 
%     close all % Close the figure after input
% end

% Extract the indices of good units
% good_unit_indices = find(good_units);
% 
% % Display the selected good units
% disp("Good Units:");
% disp(good_unit_indices);
%%
epo = 500; % Epoch length
evnt = "pos_thresh"; % Event type
% folder_path = "C:\Users\PersonLab\Documents\NeuropixelsResults\ShamAnimals\050922_SAM009\mod cells";

for c = filtered_modulated_indices
    % Create a new figure for each unit
    figure('Position', [100, 100, 800, 900]); hold on
    
    % Plot Stimulus aligned to T0 (subplot 1)
    subplot(4, 1, 1); hold on
    h = title(['PSTH for Unit ', num2str(c), ' ' 'depth' ' ' num2str(data.units(c).depth)]);
    set(h, 'interpreter', 'none')

    dim = 8;
    StimvelPlot(NoRR_T0aligned, dim, NoCSMcolor, 0)
    StimvelPlot(ST0_CSM, dim, T0color, 0)

    xlabel('CSM Aligned-T0')
    ylim('auto')
    xlim([-.5, .5])
    yline(0, '--')
    
    % Plot raster aligned to NoCSM condition (subplot 2)
    subplot(4, 1, 2); hold on
    title('NoCSM Aligned to T0')

    tick_length = 1; % length of tick marks for spikes
    spike_tick_width = 0.5; % width of lines marking spikes
    event_tick_width = 2; % width of lines marking secondary events

    % Extract raster for NoCSM condition
    rastr = data.units(c).NoCSM_csmT0.Raster(:, (1000 - epo):(2000 - (1000 - epo)));
    [y_co, x_co] = find(rastr);
    x_co = x_co - epo; % Align to the event
    line([x_co'; x_co'], [y_co' - (tick_length / 2); y_co' + (tick_length / 2)], 'Color', [0.5 0.5 0.5], 'LineWidth', spike_tick_width);
    line([zeros(1, size(rastr, 1)); zeros(1, size(rastr, 1))], [(1:1:size(rastr, 1)) - 0.5; (1:1:size(rastr, 1)) + 0.5], 'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);
    axis ij;
    %ylim([0, size(rastr, 1)])
    ylim([0, size(rastr, 1)] +0.5) % Adjust the ylim based on the height scaling
    xlim([-epo, epo])
    
    % Plot raster aligned to T0_csm condition (subplot 3)
    subplot(4, 1, 3); hold on
    title('T0_csm')

    % Extract raster for T0_csm condition
    rastr = data.units(c).T0_csm.Raster(:, (1000 - epo):(2000 - (1000 - epo)));
    [y_co, x_co] = find(rastr);
    x_co = x_co - epo; % Align to the event
    line([x_co'; x_co'], [y_co' - (tick_length / 2); y_co' + (tick_length / 2)], 'Color', [0.5 0.5 0.5], 'LineWidth', spike_tick_width);
    line([zeros(1, size(rastr, 1)); zeros(1, size(rastr, 1))], [(1:1:size(rastr, 1)) - 0.5; (1:1:size(rastr, 1)) + 0.5], 'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);
    axis ij;
    ylim([0, size(rastr, 1)]+ 0.5)
    xlim([-epo, epo])

    % --- PSTH Overlay (subplot 4) ---
    subplot(4, 1, 4); hold on
    % Extract PSTH data (already smoothed)
    psth_NoCSM_smooth = data.units(c).NoCSM_csmT0.PSTH_smoothbin;
    psth_CSM_smooth = data.units(c).T0_csm.PSTH_smoothbin;
    
    % Time vector based on bin size (assuming 10ms bins)
    time_vector = linspace(-1000, 1000, length(psth_NoCSM_smooth)); % Adjust based on window size
    
    % Gaussian filter for smoothing
    smoother = 10; % Standard deviation of the Gaussian
    gauss_filter = fspecial('gaussian', [1, 6 * smoother], smoother);
    
    % Smooth the PSTHs using the Gaussian filter
    psth_NoCSM_gaussian = conv(psth_NoCSM_smooth, gauss_filter, 'same');
    psth_CSM_gaussian = conv(psth_CSM_smooth, gauss_filter, 'same');
    
    % Plot PSTH for NoCSM and CSM conditions
    plot(time_vector, psth_NoCSM_gaussian, 'Color', NoCSMcolor, 'LineWidth', 1.5); % Smooth PSTH for NoCSM
    hold on;
    plot(time_vector, psth_CSM_gaussian, 'Color', T0color, 'LineWidth', 1.5); % Smooth PSTH for CSM
    xlabel('Time (ms)');
    ylabel('Firing Rate (Hz)');
    title(['PSTH for Unit ', num2str(c)]);
    xline(0, '--k', 'LineWidth', 2);
    xlim([-500 500]); % Adjust to view from -500 to 500 ms
    grid off;

    % Optionally save the figure
    % filename = sprintf('PSTH_Overlay_Unit_%d.svg', c);
    % full_path = fullfile(folder_path, filename);
    % saveas(gcf, full_path);

    % Optionally close the figure after viewing
    if mod(c, 10) == 0
        drawnow
        input("Press Enter when done looking at figures");
        close all
    end
end


%% Bootstrapping PSTH - to help with low sample size number 

% Parameters
epo = 500; % 500 ms before and after event
n_boot = 1000; % number of bootstrap iterations

% Loop over units
for c = filtered_modulated_indices
    figure('Position', [100, 100, 800, 900]); hold on

    % Subplot 1: Stimulus Velocity Plots
    subplot(4,1,1); hold on
    h = title(['PSTH for Unit ', num2str(c), ' depth ' num2str(data.units(c).depth)]);
    set(h, 'interpreter', 'none')

    dim = 8;
    StimvelPlot(NoRR_T0aligned, dim, NoCSMcolor, 0)
    StimvelPlot(ST0_CSM, dim, T0color, 0)
    xlabel('CSM Aligned-T0')
    ylim('auto')
    xlim([-.5, .5])
    yline(0, '--')

    % Subplot 2: Raster NoCSM aligned to T0
    subplot(4,1,2); hold on
    title('NoCSM Aligned to T0')
    tick_length = 1; spike_tick_width = 0.5; event_tick_width = 2;

    rastr = data.units(c).NoCSM_csmT0.Raster(:, (1000-epo):(1000+epo-1));
    [y_co, x_co] = find(rastr);
    x_co = x_co - epo;
    line([x_co'; x_co'], [y_co' - (tick_length/2); y_co' + (tick_length/2)], 'Color', [0.5 0.5 0.5], 'LineWidth', spike_tick_width);
    line([zeros(1,size(rastr,1)); zeros(1,size(rastr,1))], [(1:size(rastr,1)) - 0.5; (1:size(rastr,1)) + 0.5], 'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);
    axis ij;
    ylim([0, size(rastr,1)]); xlim([-epo, epo])

    % Subplot 3: Raster T0_csm
    subplot(4,1,3); hold on
    title('T0_csm')
    rastr2 = data.units(c).T0_csm.Raster(:, (1000-epo):(1000+epo-1));
    [y_co2, x_co2] = find(rastr2);
    x_co2 = x_co2 - epo;
    line([x_co2'; x_co2'], [y_co2' - (tick_length/2); y_co2' + (tick_length/2)], 'Color', [0.5 0.5 0.5], 'LineWidth', spike_tick_width);
    line([zeros(1,size(rastr2,1)); zeros(1,size(rastr2,1))], [(1:size(rastr2,1)) - 0.5; (1:size(rastr2,1)) + 0.5], 'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);
    axis ij;
    ylim([0, size(rastr2,1)]); xlim([-epo, epo])

    % Subplot 4: PSTHs with bootstrap CI
    subplot(4,1,4); hold on
    title('PSTHs with Bootstrap 95% CI')

    time_vector = linspace(-0.5, 0.5, size(rastr, 2)); % robust length

    [mean_psth1, ci1] = BootstrapPSTH(rastr, n_boot);
    [mean_psth2, ci2] = BootstrapPSTH(rastr2, n_boot);

    fill([time_vector fliplr(time_vector)], [ci1(1,:) fliplr(ci1(2,:))], [0.5 0.5 0.5], 'FaceAlpha', 0.3, 'EdgeColor', 'none')
    plot(time_vector, mean_psth1, 'Color', NoCSMcolor, 'LineWidth', 1.5)

    fill([time_vector fliplr(time_vector)], [ci2(1,:) fliplr(ci2(2,:))], T0color, 'FaceAlpha', 0.3, 'EdgeColor', 'none')
    plot(time_vector, mean_psth2, 'Color', T0color, 'LineWidth', 1.5)

    xlabel('Time (s)'); ylabel('Firing Rate (Hz)')
    xlim([-0.5 0.5])
    yline(0, '--')
    legend('NoCSM CI', 'NoCSM Mean', 'T0_csm CI', 'T0_csm Mean')

    if mod(c,10) == 0
        drawnow
        input("press enter when done looking at figures")
        close all
    end
end



%%
epo = 500;
evnt = "pos_thresh";

for c = filtered_modulated_indices

    figure; hold on

    % ---- Subplot 1: Stimulus-aligned velocity traces ----
    subplot(4,1,1); hold on
    h = title(['PSTH for Unit ', num2str(c), ' ' 'depth' ' ' num2str(data.units(c).depth)]) ;
    set(h,'interpreter','none')

    dim = 8;
    StimvelPlot(NoRR_T0aligned, dim, NoCSMcolor, 0)
    StimvelPlot(ST0_CSM, dim, T0color, 0)
    xlabel('CSM Aligned-T0')
    ylim('auto')
    xlim([-.5, .5])
    yline(0, '--')

    % ---- Subplot 2: Raster - NoCSM ----
    subplot(4,1,2); hold on
    title('NoCSM Aligned to T0')

    tick_length = 1;
    spike_tick_width = 0.5;
    event_tick_width = 2;

    rastr = data.units(c).NoCSM_csmT0.Raster(:, (1000-epo):(2000-(1000-epo)));
    [y_co, x_co] = find(rastr);
    x_co = x_co - epo;

    line([x_co'; x_co'], [y_co'-(tick_length/2); y_co'+(tick_length/2)], 'Color', [0.5 0.5 0.5], 'LineWidth', spike_tick_width);
    line([zeros(1,size(rastr,1)); zeros(1,size(rastr,1))], [(1:size(rastr,1))-.5; (1:size(rastr,1))+.5], 'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);
    axis ij;
    ylim([0, size(rastr,1)])
    xlim([-epo, epo])

    % ---- Subplot 3: Raster - T0_csm ----
    subplot(4,1,3); hold on
    title('T0_csm')

    rastr = data.units(c).T0_csm.Raster(:, (1000-epo):(2000-(1000-epo)));
    [y_co, x_co] = find(rastr);
    x_co = x_co - epo;

    line([x_co'; x_co'], [y_co'-(tick_length/2); y_co'+(tick_length/2)], 'Color', [0.5 0.5 0.5], 'LineWidth', spike_tick_width);
    line([zeros(1,size(rastr,1)); zeros(1,size(rastr,1))], [(1:size(rastr,1))-.5; (1:size(rastr,1))+.5], 'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);
    axis ij;
    ylim([0, size(rastr,1)])
    xlim([-epo, epo])

    % ---- Subplot 4: PSTHs (NoCSM and CSM) ----
    subplot(4,1,4); hold on
    title('PSTH - NoCSM vs CSM')
    xlabel('Time (s)')
    ylabel('Firing Rate (Hz)')

    psth_nocsm_smooth = data.units(c).NoCSM_csmT0.PSTH_smooth;
    psth_nocsm_raw = data.units(c).NoCSM_csmT0.PSTH_raw;

    psth_t0_smooth = data.units(c).T0_csm.PSTH_smooth;
    psth_t0_raw = data.units(c).T0_csm.PSTH_raw;

    td_nocsm = length(psth_nocsm_raw) - length(psth_nocsm_smooth);
    td2_nocsm = td_nocsm / 2;

    td_t0 = length(psth_t0_raw) - length(psth_t0_smooth);
    td2_t0 = td_t0 / 2;

    time = (-epo:epo) / 1000;

    plot(time, psth_nocsm_smooth((1000-epo)-td2_nocsm:(2000-(1000-epo))-td2_nocsm), 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5)
    plot(time, psth_t0_smooth((1000-epo)-td2_t0:(2000-(1000-epo))-td2_t0), 'Color', [0 0 0], 'LineWidth', 1.5)

    legend({'NoCSM', 'CSM'}, 'Location', 'northeast')
    xlim([-epo, epo]/1000)

    % ---- Optional: Wait every 10 figures ----
    if mod(c,10) == 0
        drawnow
        input("press enter when done looking at figures")
        close all
    end
end

 %% trial wise comparison T0
win_size = 50;   % Window size in ms (or samples, assuming 1 kHz)
step_size = 10;  % Step size in ms (or samples)
epo = 500;       % Half-width of the peri-event window
fs = 1000;       % Sampling rate (Hz)

win_edges = 1:step_size:(2*epo - win_size + 1); % Sliding window indices
win_centers = win_edges + win_size/2;           % For plotting (relative to event)

for c = filtered_modulated_indices
    % Extract rasters
    ras_nocsm = data.units(c).NoCSM_csmT0.Raster(:,(1000-epo):(2000-(1000-epo)));
    ras_csm   = data.units(c).T0_csm.Raster(:,(1000-epo):(2000-(1000-epo)));

    % Ensure both rasters have the same number of trials
    num_trials = min(size(ras_nocsm, 1), size(ras_csm, 1));
    ras_nocsm = ras_nocsm(1:num_trials, :);  % Trim to min number of trials
    ras_csm = ras_csm(1:num_trials, :);      % Trim to min number of trials

    % Preallocate
    fr_diff = zeros(size(win_edges));

    for i = 1:length(win_edges)
        t1 = win_edges(i);
        t2 = t1 + win_size - 1;

        % Calculate firing rate for each trial within the window
        fr_nocsm_trial = sum(ras_nocsm(:,t1:t2), 2) / (win_size/fs);  % Firing rate for each trial
        fr_csm_trial   = sum(ras_csm(:,t1:t2), 2) / (win_size/fs);    % Firing rate for each trial

        % Compute the difference for each trial and then average
        fr_diff(i) = mean(fr_csm_trial - fr_nocsm_trial);
    end

    % Plot firing rate difference
    figure;
    plot(win_centers-epo, fr_diff, 'k', 'LineWidth', 2);
    xlabel('Time (ms) relative to T0');
    ylabel('Firing Rate Diff (CSM - NoCSM)');
    title(['Unit ', num2str(c), ' | Depth: ', num2str(data.units(c).depth)]);
    yline(0, '--');
    xlim([-epo epo]);
    grid on;
end



%% FR chanes and kinematic differences 
% Sliding window parameters
win_size = 25; step_size = 10; epo = 500; fs = 1000;
win_edges = 1:step_size:(2*epo - win_size + 1);
win_centers = win_edges + win_size/2;

% Kinematic differences
T0CSM_diff = ST0_CSM.Average(:,5) - NoRR_T0aligned.Average(:,5); % outward accel
T2CSM_diff = ST2_CSM.Average(:,5) - NoRR_T2aligned.Average(:,5);
t = ST0_CSM.Average(:,1);  % Time vector for kinematics (should be -0.5 to 0.5 s)

for c = filtered_modulated_indices
    % Extract rasters
    ras_nocsm = data.units(c).NoCSM_csmT0.Raster(:,(1000-epo):(2000-(1000-epo)));
    ras_csm   = data.units(c).T0_csm.Raster(:,(1000-epo):(2000-(1000-epo)));
    ras_t2    = data.units(c).T2_csm.Raster(:,(1000-epo):(2000-(1000-epo)));

    % Match trial count
    num_trials = min([size(ras_nocsm, 1), size(ras_csm, 1), size(ras_t2, 1)]);
    ras_nocsm = ras_nocsm(1:num_trials, :);
    ras_csm   = ras_csm(1:num_trials, :);
    ras_t2    = ras_t2(1:num_trials, :);

    % Preallocate FR diffs
    fr_diff_csm = zeros(size(win_edges));
    fr_diff_t2  = zeros(size(win_edges));

    for i = 1:length(win_edges)
        t1 = win_edges(i); t2win = t1 + win_size - 1;

        fr_nocsm_trial = sum(ras_nocsm(:,t1:t2win), 2) / (win_size/fs);
        fr_csm_trial   = sum(ras_csm(:,t1:t2win), 2) / (win_size/fs);
        fr_t2_trial    = sum(ras_t2(:,t1:t2win), 2) / (win_size/fs);

        fr_diff_csm(i) = mean(fr_csm_trial - fr_nocsm_trial);
        fr_diff_t2(i)  = mean(fr_t2_trial - fr_nocsm_trial);
    end

    % ==== Plot ====
    figure;

    % Subplot 1: Kinematic difference (TOP)
    subplot(2,1,1);
    hold on;
    plot(t, T0CSM_diff,'Color', T0color, 'LineWidth', 2);
    plot(t, T2CSM_diff,'Color',T2color, 'LineWidth', 2);
    yline(0, '--k');
    xlabel('Time (s) relative to event');
    ylabel('Outward Accel. Diff');
    legend({'T0-CSM vs NoCSM', 'T2-CSM vs NoCSM'}, 'Location', 'best');
    title(['Unit ', num2str(c), ' | Depth: ', num2str(data.units(c).depth)]);
    xlim([-0.5 0.5]);
    grid on;

    % Subplot 2: Firing rate diff (BOTTOM)
    subplot(2,1,2);
    hold on;
    plot(win_centers - epo, fr_diff_csm, 'Color',T0color, 'LineWidth', 2);
    plot(win_centers - epo, fr_diff_t2, 'Color',T2color, 'LineWidth', 2);
    yline(0, '--');
    xlabel('Time (ms) relative to T0');
    ylabel('FR Diff (Hz)');
    legend({'T0-CSM - NoCSM', 'T2-CSM - NoCSM'});
    xlim([-epo epo]);
    grid on;
end
%% More smoothing version 
% FR changes and kinematic differences
% Sliding window parameters
win_size = 25; step_size = 10; epo = 500; fs = 1000;
win_edges = 1:step_size:(2*epo - win_size + 1);
win_centers = win_edges + win_size/2;

% Gaussian smoothing kernel
sigma = 1.5;  % Controls smoothness
kernel_size = 2*ceil(3*sigma) + 1;  % Ensures odd length
gauss_kernel = fspecial('gaussian', [1 kernel_size], sigma);

% Kinematic differences (raw)
T0CSM_diff = ST0_CSM.Average(:,5) - NoRR_T0aligned.Average(:,5); % outward accel
T2CSM_diff = ST2_CSM.Average(:,5) - NoRR_T2aligned.Average(:,5);
t = ST0_CSM.Average(:,1);  % Time vector for kinematics (should be -0.5 to 0.5 s)

% Optional: smooth kinematics
smooth_kin = false;  % Set to true if you want to smooth kinematics too
if smooth_kin
    T0CSM_diff = conv(T0CSM_diff, gauss_kernel, 'same');
    T2CSM_diff = conv(T2CSM_diff, gauss_kernel, 'same');
end

for c = filtered_modulated_indices
%for c = modulated_indices
    % Extract rasters
    ras_nocsm = data.units(c).NoCSM_csmT0.Raster(:,(1000-epo):(2000-(1000-epo)));
    ras_csm   = data.units(c).T0_csm.Raster(:,(1000-epo):(2000-(1000-epo)));
    ras_t2    = data.units(c).T2_csm.Raster(:,(1000-epo):(2000-(1000-epo)));

    % Match trial count
    num_trials = min([size(ras_nocsm, 1), size(ras_csm, 1), size(ras_t2, 1)]);
    ras_nocsm = ras_nocsm(1:num_trials, :);
    ras_csm   = ras_csm(1:num_trials, :);
    ras_t2    = ras_t2(1:num_trials, :);

    % Preallocate FR diffs
    fr_diff_csm = zeros(size(win_edges));
    fr_diff_t2  = zeros(size(win_edges));

    for i = 1:length(win_edges)
        t1 = win_edges(i); t2win = t1 + win_size - 1;

        fr_nocsm_trial = sum(ras_nocsm(:,t1:t2win), 2) / (win_size/fs);
        fr_csm_trial   = sum(ras_csm(:,t1:t2win), 2) / (win_size/fs);
        fr_t2_trial    = sum(ras_t2(:,t1:t2win), 2) / (win_size/fs);

        fr_diff_csm(i) = mean(fr_csm_trial - fr_nocsm_trial);
        fr_diff_t2(i)  = mean(fr_t2_trial - fr_nocsm_trial);
    end

    % Smooth FR difference traces
    fr_diff_csm_smooth = conv(fr_diff_csm, gauss_kernel, 'same');
    fr_diff_t2_smooth  = conv(fr_diff_t2, gauss_kernel, 'same');

    % ==== Plot ====
    figure;

    % Subplot 1: Kinematic difference (TOP)
    subplot(2,1,1);
    hold on;
    plot(t, T0CSM_diff,'Color', T0color, 'LineWidth', 2);
    plot(t, T2CSM_diff,'Color',T2color, 'LineWidth', 2);
    yline(0, '--k');
    xlabel('Time (s) relative to event');
    ylabel('Outward Accel. Diff');
    legend({'T0-CSM vs NoCSM', 'T2-CSM vs NoCSM'}, 'Location', 'best');
    title(['Unit ', num2str(c), ' | Depth: ', num2str(data.units(c).depth)]);
    xlim([-0.5 0.5]);
    grid on;

    % Subplot 2: Firing rate diff (BOTTOM)
    subplot(2,1,2);
    hold on;
    plot(win_centers - epo, fr_diff_csm_smooth, 'Color',T0color, 'LineWidth', 2);
    plot(win_centers - epo, fr_diff_t2_smooth, 'Color',T2color, 'LineWidth', 2);
    yline(0, '--');
    xlabel('Time (ms) relative to T0');
    ylabel('FR Diff (Hz)');
    legend({'T0-CSM - NoCSM', 'T2-CSM - NoCSM'});
    xlim([-epo epo]);
 
end

%% Same but with more smoothing 
% Parameters
win_size = 25;          % Window size in ms (or samples, assuming 1 kHz)
step_size = 10;         % Step size in ms (or samples)
epo = 500;              % Half-width of the peri-event window
fs = 1000;              % Sampling rate (Hz)
num_boot = 1000;        % Number of bootstrap samples
%sigma = 1.5;  % Smoothing std dev for Gaussian
sigma = 1;

win_edges = 1:step_size:(2*epo - win_size + 1);
win_centers = win_edges + win_size/2;

% Gaussian smoothing kernel
gauss_kernel = exp(-0.5 * ((-2*sigma):(2*sigma)).^2 / sigma^2);
gauss_kernel = gauss_kernel / sum(gauss_kernel);

% Kinematic differences
T0CSM_diff = ST0_CSM.Average(:,5) - NoRR_T0aligned.Average(:,5); % outward accel
T2CSM_diff = ST2_CSM.Average(:,5) - NoRR_T2aligned.Average(:,5);
t = ST0_CSM.Average(:,1);

% Loop through units
for c = filtered_modulated_indices
    % ==== Extract rasters ====
    ras_nocsmT0 = data.units(c).NoCSM_csmT0.Raster(:,(1000-epo):(2000-(1000-epo)));
    ras_csmT0   = data.units(c).T0_csm.Raster(:,(1000-epo):(2000-(1000-epo)));
    ras_nocsmT2 = data.units(c).NoCSM_csmT2.Raster(:,(1000-epo):(2000-(1000-epo)));
    ras_t2      = data.units(c).T2_csm.Raster(:,(1000-epo):(2000-(1000-epo)));

    % Match trial count
    num_trials = min([size(ras_nocsmT0, 1), size(ras_csmT0, 1), size(ras_nocsmT2, 1), size(ras_t2, 1)]);
    ras_nocsmT0 = ras_nocsmT0(1:num_trials, :);
    ras_csmT0   = ras_csmT0(1:num_trials, :);
    ras_nocsmT2 = ras_nocsmT2(1:num_trials, :);
    ras_t2      = ras_t2(1:num_trials, :);

    % Preallocate
    fr_diff_t0 = zeros(size(win_edges));
    fr_diff_t2 = zeros(size(win_edges));
    ci_t0 = zeros(2, length(win_edges));
    ci_t2 = zeros(2, length(win_edges));

    % ==== Compute firing rate differences and bootstrapped CIs ====
    for i = 1:length(win_edges)
        t1 = win_edges(i);
        t2win = t1 + win_size - 1;

        fr_nocsmT0 = sum(ras_nocsmT0(:,t1:t2win), 2) / (win_size/fs);
        fr_csmT0   = sum(ras_csmT0(:,t1:t2win), 2) / (win_size/fs);
        fr_nocsmT2 = sum(ras_nocsmT2(:,t1:t2win), 2) / (win_size/fs);
        fr_csmT2   = sum(ras_t2(:,t1:t2win), 2) / (win_size/fs);

        % Differences
        d_t0 = fr_csmT0 - fr_nocsmT0;
        d_t2 = fr_csmT2 - fr_nocsmT2;

        % Store mean differences
        fr_diff_t0(i) = mean(d_t0);
        fr_diff_t2(i) = mean(d_t2);

        % Bootstrap
        boot_diff_t0 = bootstrp(num_boot, @mean, d_t0);
        boot_diff_t2 = bootstrp(num_boot, @mean, d_t2);
        ci_t0(:,i) = prctile(boot_diff_t0, [2.5, 97.5]);
        ci_t2(:,i) = prctile(boot_diff_t2, [2.5, 97.5]);
    end

    % Smooth FR differences and CI
    fr_diff_t0_smooth = conv(fr_diff_t0, gauss_kernel, 'same');
    fr_diff_t2_smooth = conv(fr_diff_t2, gauss_kernel, 'same');
    ci_t0_smooth = [conv(ci_t0(1,:), gauss_kernel, 'same'); conv(ci_t0(2,:), gauss_kernel, 'same')];
    ci_t2_smooth = [conv(ci_t2(1,:), gauss_kernel, 'same'); conv(ci_t2(2,:), gauss_kernel, 'same')];

    % ==== Plot ====
    figure;

    % === Top: Kinematic differences ===
    subplot(2,1,1);
    hold on;
    plot(t, T0CSM_diff, 'Color', T0color, 'LineWidth', 2);
    plot(t, T2CSM_diff, 'Color', T2color, 'LineWidth', 2);
    yline(0, '--k');
    xlabel('Time (s) relative to event');
    ylabel('Outward Accel. Diff');
    legend({'T0-CSM vs NoCSM', 'T2-CSM vs NoCSM'}, 'Location', 'best');
    title(['Unit ', num2str(c), ' | Depth: ', num2str(data.units(c).depth)]);
    xlim([-0.5 0.5]);
    grid on;

    % === Bottom: FR differences + CIs ===
    subplot(2,1,2);
    hold on;
    % Shaded CIs
    fill([win_centers-epo, fliplr(win_centers-epo)], ...
         [ci_t0_smooth(1,:), fliplr(ci_t0_smooth(2,:))], ...
         T0color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    fill([win_centers-epo, fliplr(win_centers-epo)], ...
         [ci_t2_smooth(1,:), fliplr(ci_t2_smooth(2,:))], ...
         T2color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

    % Smoothed FR differences
    plot(win_centers - epo, fr_diff_t0_smooth, 'Color', T0color, 'LineWidth', 2);
    plot(win_centers - epo, fr_diff_t2_smooth, 'Color', T2color, 'LineWidth', 2);
    yline(0, '--');
    xlabel('Time (ms) relative to T0/T2');
    ylabel('FR Diff (Hz)');
    legend({'T0 diff', 'T2 diff'});
    xlim([-epo epo]);
    grid on;
end
%% shorter window to look at 
% Parameters
win_size = 25;          % Window size in ms (or samples, assuming 1 kHz)
step_size = 10;         % Step size in ms (or samples)
epo = 500;              % Half-width of the peri-event window
fs = 1000;              % Sampling rate (Hz)
num_boot = 1000;        % Number of bootstrap samples
sigma = 1;              % Smoothing std dev for Gaussian

win_edges = 1:step_size:(2*epo - win_size + 1);
win_centers = win_edges + win_size/2;

% Gaussian smoothing kernel
gauss_kernel = exp(-0.5 * ((-2*sigma):(2*sigma)).^2 / sigma^2);
gauss_kernel = gauss_kernel / sum(gauss_kernel);

% % Kinematic differences
% T0CSM_diff = ST0_CSM.Average(:,5) - NoRR_T0aligned.Average(:,5); % outward accel
% T2CSM_diff = ST2_CSM.Average(:,5) - NoRR_T2aligned.Average(:,5);
% t = ST0_CSM.Average(:,1);

% Bootstrapped smooth difference and CI for T0 and T2 kinematics
vel1 = squeeze(ST0_CSM.stimvel(:,5,:));
vel2 = squeeze(NoRR_T0aligned.stimvel(:,5,:));
[mean_t0_diff, ci_t0_kin_smooth] = bootstrap_kinematic_diff(vel1, vel2, sigma, num_boot);

vel3 = squeeze(ST2_CSM.stimvel(:,5,:));
vel4 = squeeze(NoRR_T2aligned.stimvel(:,5,:));
[mean_t2_diff, ci_t2_kin_smooth] = bootstrap_kinematic_diff(vel3, vel4, sigma, num_boot);
t = ST0_CSM.Average(:,1);  % time vector

% Loop through units
for c = filtered_modulated_indices
    % ==== Extract rasters ====
    ras_nocsmT0 = data.units(c).NoCSM_csmT0.Raster(:,(1000-epo):(2000-(1000-epo)));
    ras_csmT0   = data.units(c).T0_csm.Raster(:,(1000-epo):(2000-(1000-epo)));
    ras_nocsmT2 = data.units(c).NoCSM_csmT2.Raster(:,(1000-epo):(2000-(1000-epo)));
    ras_t2      = data.units(c).T2_csm.Raster(:,(1000-epo):(2000-(1000-epo)));

    % Match trial count
    num_trials = min([size(ras_nocsmT0, 1), size(ras_csmT0, 1), size(ras_nocsmT2, 1), size(ras_t2, 1)]);
    ras_nocsmT0 = ras_nocsmT0(1:num_trials, :);
    ras_csmT0   = ras_csmT0(1:num_trials, :);
    ras_nocsmT2 = ras_nocsmT2(1:num_trials, :);
    ras_t2      = ras_t2(1:num_trials, :);

    % Preallocate
    fr_diff_t0 = zeros(size(win_edges));
    fr_diff_t2 = zeros(size(win_edges));
    ci_t0 = zeros(2, length(win_edges));
    ci_t2 = zeros(2, length(win_edges));

    % ==== Compute firing rate differences and bootstrapped CIs ====
    for i = 1:length(win_edges)
        t1 = win_edges(i);
        t2win = t1 + win_size - 1;

        fr_nocsmT0 = sum(ras_nocsmT0(:,t1:t2win), 2) / (win_size/fs);
        fr_csmT0   = sum(ras_csmT0(:,t1:t2win), 2) / (win_size/fs);
        fr_nocsmT2 = sum(ras_nocsmT2(:,t1:t2win), 2) / (win_size/fs);
        fr_csmT2   = sum(ras_t2(:,t1:t2win), 2) / (win_size/fs);

        % Differences
        d_t0 = fr_csmT0 - fr_nocsmT0;
        d_t2 = fr_csmT2 - fr_nocsmT2;

        % Store mean differences
        fr_diff_t0(i) = mean(d_t0);
        fr_diff_t2(i) = mean(d_t2);

        % Bootstrap
        boot_diff_t0 = bootstrp(num_boot, @mean, d_t0);
        boot_diff_t2 = bootstrp(num_boot, @mean, d_t2);
        ci_t0(:,i) = prctile(boot_diff_t0, [2.5, 97.5]);
        ci_t2(:,i) = prctile(boot_diff_t2, [2.5, 97.5]);
    end

    % Smooth FR differences and CI
    fr_diff_t0_smooth = conv(fr_diff_t0, gauss_kernel, 'same');
    fr_diff_t2_smooth = conv(fr_diff_t2, gauss_kernel, 'same');
    ci_t0_smooth = [conv(ci_t0(1,:), gauss_kernel, 'same'); conv(ci_t0(2,:), gauss_kernel, 'same')];
    ci_t2_smooth = [conv(ci_t2(1,:), gauss_kernel, 'same'); conv(ci_t2(2,:), gauss_kernel, 'same')];

    % === Restrict to -200 to 200 ms for visualization ===
    vis_range = (-200 <= (win_centers - epo)) & ((win_centers - epo) <= 200);

    % Trim smoothed signals and CIs
    plot_x = (win_centers - epo);
    plot_x_vis = plot_x(vis_range);
    fr_diff_t0_smooth = fr_diff_t0_smooth(vis_range);
    fr_diff_t2_smooth = fr_diff_t2_smooth(vis_range);
    ci_t0_smooth = ci_t0_smooth(:, vis_range);
    ci_t2_smooth = ci_t2_smooth(:, vis_range);

    % % ==== Plot ====
    % figure;
    % 
    % % === Top: Kinematic differences ===
    % subplot(2,1,1);
    % hold on;
    % plot(t, T0CSM_diff, 'Color', T0color, 'LineWidth', 2);
    % plot(t, T2CSM_diff, 'Color', T2color, 'LineWidth', 2);
    % yline(0, '--k');
    % xlabel('Time (s) relative to event');
    % ylabel('Outward Accel. Diff');
    % legend({'T0-CSM vs NoCSM', 'T2-CSM vs NoCSM'}, 'Location', 'best');
    % title(['Unit ', num2str(c), ' | Depth: ', num2str(data.units(c).depth)]);
    % xlim([-0.2 0.2]);
    % grid on;
    
    % % Plot kinematics
    figure;
    subplot(2,1,1);
    hold on;

    % Trim kinematic CI and mean to match [-200 200 ms] window
    kin_vis_range = (t >= -0.2) & (t <= 0.2);
    t_vis = t(kin_vis_range);
    mean_t0_vis = mean_t0_diff(kin_vis_range);
    mean_t2_vis = mean_t2_diff(kin_vis_range);
    ci_t0_vis = ci_t0_kin_smooth(kin_vis_range, :);
    ci_t2_vis = ci_t2_kin_smooth(kin_vis_range, :);

    % Shaded CI
    fill([t_vis; flipud(t_vis)], ...
        [ci_t0_vis(:,1); flipud(ci_t0_vis(:,2))], ...
        T1color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

    fill([t_vis; flipud(t_vis)], ...
        [ci_t2_vis(:,1); flipud(ci_t2_vis(:,2))], ...
        T2color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

    % Mean difference lines
    plot(t_vis, mean_t0_vis, 'Color', T0color, 'LineWidth', 2);
    plot(t_vis, mean_t2_vis, 'Color', T2color, 'LineWidth', 2);

    yline(0, '--k');
    xlabel('Time (s) relative to event');
    ylabel('Outward Vel. Diff (cm/s)');
    legend({'T0 CI', 'T2 CI', 'T0 Mean', 'T2 Mean'}, 'Location', 'best');
    title(['Unit ', num2str(c), ' | Depth: ', num2str(data.units(c).depth)]);
    xlim([-0.2 0.2]);
    grid on;

    % === Bottom: FR differences + CIs ===
    subplot(2,1,2);
    hold on;
    % Shaded CIs
    fill([plot_x_vis, fliplr(plot_x_vis)], ...
         [ci_t0_smooth(1,:), fliplr(ci_t0_smooth(2,:))], ...
         T0color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    fill([plot_x_vis, fliplr(plot_x_vis)], ...
         [ci_t2_smooth(1,:), fliplr(ci_t2_smooth(2,:))], ...
         T2color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

    % Smoothed FR differences
    plot(plot_x_vis, fr_diff_t0_smooth, 'Color', T0color, 'LineWidth', 2);
    plot(plot_x_vis, fr_diff_t2_smooth, 'Color', T2color, 'LineWidth', 2);
    yline(0, '--');
    xlabel('Time (ms) relative to T0/T2');
    ylabel('FR Diff (Hz)');
    legend({'T0 diff', 'T2 diff'});
    xlim([-200 200]);
    grid on;
end


%%
plot(t, ci_t0_kin_smooth(1,:), 'r');  % Lower CI
plot(t, ci_t0_kin_smooth(2,:), 'b');  % Upper CI
Error using plot
Specify the coordinates as vectors or matrices of the same size, or as a vector and a matrix that share the same length in at
least one dimension.
 

%% Without Smoothing 
% Parameters
win_size = 25;          % Window size in ms
step_size = 10;         % Step size in ms
epo = 500;              % Half-width of the peri-event window
fs = 1000;              % Sampling rate (Hz)
num_boot = 1000;        % Number of bootstrap samples

win_edges = 1:step_size:(2*epo - win_size + 1);
win_centers = win_edges + win_size/2;

% Kinematic differences
T0CSM_diff = ST0_CSM.Average(:,5) - NoRR_T0aligned.Average(:,5); % outward accel
T2CSM_diff = ST2_CSM.Average(:,5) - NoRR_T2aligned.Average(:,5);
t = ST0_CSM.Average(:,1);

% Loop through units
for c = filtered_modulated_indices
    % ==== Extract rasters ====
    ras_nocsmT0 = data.units(c).NoCSM_csmT0.Raster(:,(1000-epo):(2000-(1000-epo)));
    ras_csmT0   = data.units(c).T0_csm.Raster(:,(1000-epo):(2000-(1000-epo)));
    ras_nocsmT2 = data.units(c).NoCSM_csmT2.Raster(:,(1000-epo):(2000-(1000-epo)));
    ras_t2      = data.units(c).T2_csm.Raster(:,(1000-epo):(2000-(1000-epo)));

    % Match trial count
    num_trials = min([size(ras_nocsmT0, 1), size(ras_csmT0, 1), size(ras_nocsmT2, 1), size(ras_t2, 1)]);
    ras_nocsmT0 = ras_nocsmT0(1:num_trials, :);
    ras_csmT0   = ras_csmT0(1:num_trials, :);
    ras_nocsmT2 = ras_nocsmT2(1:num_trials, :);
    ras_t2      = ras_t2(1:num_trials, :);

    % Preallocate
    fr_diff_t0 = zeros(size(win_edges));
    fr_diff_t2 = zeros(size(win_edges));
    ci_t0 = zeros(2, length(win_edges));
    ci_t2 = zeros(2, length(win_edges));

    % ==== Compute firing rate differences and bootstrapped CIs ====
    for i = 1:length(win_edges)
        t1 = win_edges(i);
        t2win = t1 + win_size - 1;

        fr_nocsmT0 = sum(ras_nocsmT0(:,t1:t2win), 2) / (win_size/fs);
        fr_csmT0   = sum(ras_csmT0(:,t1:t2win), 2) / (win_size/fs);
        fr_nocsmT2 = sum(ras_nocsmT2(:,t1:t2win), 2) / (win_size/fs);
        fr_csmT2   = sum(ras_t2(:,t1:t2win), 2) / (win_size/fs);

        % Differences
        d_t0 = fr_csmT0 - fr_nocsmT0;
        d_t2 = fr_csmT2 - fr_nocsmT2;

        % Store mean differences
        fr_diff_t0(i) = mean(d_t0);
        fr_diff_t2(i) = mean(d_t2);

        % Bootstrap
        boot_diff_t0 = bootstrp(num_boot, @mean, d_t0);
        boot_diff_t2 = bootstrp(num_boot, @mean, d_t2);
        ci_t0(:,i) = prctile(boot_diff_t0, [2.5, 97.5]);
        ci_t2(:,i) = prctile(boot_diff_t2, [2.5, 97.5]);
    end

    % ==== Plot ====
    figure;

    % === Top: Kinematic differences ===
    subplot(2,1,1);
    hold on;
    plot(t, T0CSM_diff, 'Color', T0color, 'LineWidth', 2);
    plot(t, T2CSM_diff, 'Color', T2color, 'LineWidth', 2);
    yline(0, '--k');
    xlabel('Time (s) relative to event');
    ylabel('Outward Accel. Diff');
    legend({'T0-CSM vs NoCSM', 'T2-CSM vs NoCSM'}, 'Location', 'best');
    title(['Unit ', num2str(c), ' | Depth: ', num2str(data.units(c).depth)]);
    xlim([-0.5 0.5]);
    grid on;

    % === Bottom: FR differences + CIs ===
    subplot(2,1,2);
    hold on;
    % Shaded CIs
    fill([win_centers-epo, fliplr(win_centers-epo)], ...
         [ci_t0(1,:), fliplr(ci_t0(2,:))], ...
         T0color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    fill([win_centers-epo, fliplr(win_centers-epo)], ...
         [ci_t2(1,:), fliplr(ci_t2(2,:))], ...
         T2color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

    % Raw FR differences
    plot(win_centers - epo, fr_diff_t0, 'Color', T0color, 'LineWidth', 2);
    plot(win_centers - epo, fr_diff_t2, 'Color', T2color, 'LineWidth', 2);
    yline(0, '--');
    xlabel('Time (ms) relative to T0/T2');
    ylabel('FR Diff (Hz)');
    legend({'T0 diff', 'T2 diff'});
    xlim([-epo epo]);
    grid on;
end
%%

%% Generate example raster trace  : T2
% Plot rasters aligned to each one 


epo = 500;
%epo = 250;



%folder_path="C:\Users\PersonLab\Documents\NeuropixelsResults\ShamAnimals\050922_SAM009\mod cells";

for c=filtered_modulated_indices


    figure;hold on
   
    subplot(3,1,1); hold on
    h = title(['PSTH for Unit ', num2str(c), ' ' 'depth' ' ' num2str(data.units(c).depth)]) ;
    set(h,'interpreter','none')

    dim=8;
    StimvelPlot(NoRR_T2aligned,dim,NoCSMcolor,0)
    StimvelPlot(ST2_CSM,dim,T2color,0)
  
    %ylim([-1 2])
    xlabel('CSM Aligned-T2')
    ylim('auto')
    xlim([-.5,.5])
    yline(0,'--')
    % %xlim([-0.5 1])

  
    subplot(3,1,2); hold on
    title('NoCSM Aligned to T2')

    % Plots Raster from binary matrix
    tick_length = 1; % length of tick marks for spikes
    spike_tick_width = 0.5; % width of lines marking spikes
    event_tick_width = 2; % width of lines marking secondary events

    %baseline NoRR\

    rastr = data.units(c).NoCSM_csmT2.Raster(:,(1000-epo):(2000-(1000-epo)));

    [y_co,x_co] = find(rastr);
    x_co = x_co - epo;
    line([x_co';x_co'],[y_co'-(tick_length/2);y_co'+(tick_length/2)],'Color', [0.5 0.5 0.5], 'LineWidth', spike_tick_width);
    line([zeros(1,size(rastr,1));zeros(1,size(rastr,1))],[(1:1:size(rastr,1))-.5;(1:1:size(rastr,1))+.5],'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);
    axis ij;
    ylim([0,size(rastr,1)])
    xlim([-epo,epo])

    subplot(3,1,3); hold on
 
    title('T2_csm')
    % Plots Raster from binary matrix
    tick_length = 1; % length of tick marks for spikes
    spike_tick_width = 0.5; % width of lines marking spikes
    event_tick_width = 2; % width of lines marking secondary events

    %baseline NoRR\

    rastr = data.units(c).T2_csm.Raster(:,(1000-epo):(2000-(1000-epo)));

    [y_co,x_co] = find(rastr);
    x_co = x_co - epo;
    line([x_co';x_co'],[y_co'-(tick_length/2);y_co'+(tick_length/2)],'Color', [0.5 0.5 0.5], 'LineWidth', spike_tick_width);
    line([zeros(1,size(rastr,1));zeros(1,size(rastr,1))],[(1:1:size(rastr,1))-.5;(1:1:size(rastr,1))+.5],'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);
    axis ij;
    ylim([0,size(rastr,1)])
    xlim([-epo,epo])



if mod(c,10) == 0
    drawnow
    input("press enter when done looking at figures")
    close all
 end
end
%   drawnow
%     response = input(['Mark unit ' num2str(c) ' as good? (y/n): '], 's');
% 
%     if strcmpi(response, 'y')
%         good_units(c) = true; % Mark this unit as good
%     end
% 
%     close all % Close the figure after input
% end

% Extract the indices of good units
% good_unit_indices = find(good_units);
% 
% % Display the selected good units
% disp("Good Units:");
% disp(good_unit_indices);


%% Now plot modulated NC RN cells for each RR 


epo = 500;
%epo = 250;

folder_path="C:\Users\PersonLab\Documents\NeuropixelsResults\ShamAnimals\032522_SAM005\CSMalggnedpsths";
for c = filtered_modulated_indices
%for c=rn_unit_indices
    figure;hold on
   
    subplot(2,1,1); hold on
    dim=8;
    StimvelPlot(NoRR_T2aligned,dim,NoCSMcolor,0)
    StimvelPlot(ST2_CSM,dim,T2color,0)

    ylim('auto')
    xlim([-.5,.5])
    yline(0,'--')
    % %xlim([-0.5 1])

    subplot(2,1,2); hold on
    title('peth')

    %T0
    psth_smooth_NoCSM = data.units(c).NoCSM_csm.PSTH_smooth;
    psth_raw_NoCSM = data.units(c).NoCSM_csm.PSTH_raw;
    % psth_smooth_NoCSM = data.units(c).NoCSM_csm.PSTH_smoothbin;
    % psth_raw_NoCSM = data.units(c).NoCSM_csm.PSTH_rawbin;

    psth_smooth_T2 = data.units(c).T2_csm.PSTH_smooth;
    psth_raw_T2 = data.units(c).T2_csm.PSTH_raw;
    % psth_smooth_T2 = data.units(c).T2_csm.PSTH_smoothbin;
    % psth_raw_T2 = data.units(c).T2_csm.PSTH_rawbin;
    
    % Time windows 
    td1_NoCSM = length(psth_raw_NoCSM)-length(psth_smooth_NoCSM)+1;
    td2_NoCSM= td1_NoCSM/2;


    td1_T2 = length(psth_raw_T2)-length(psth_smooth_T2)+1;
    td2_T2= td1_T2/2;


    plot((-epo:epo),psth_smooth_NoCSM(((1000-epo)-td2_NoCSM):((2000-(1000-epo))-td2_NoCSM)),'Color',NoCSMcolor,'LineWidth',2)  
    hold on
    plot((-epo:epo),psth_smooth_T2(((1000-epo)-td2_T2):((2000-(1000-epo))-td2_T2)),'Color',T2color,'LineWidth',2)
 

    xlim([-epo,epo])
%     ylim([0,100])

    xlabel('Time from CSM- Aligned to T2')
    ylabel('FR (Hz)')
    filename=sprintf('csm_raster%d.svg',c);
    %filename=sprintf('vmax_raster%d.svg',c);
    full_path=fullfile(folder_path,filename);
    saveas(gcf,full_path);

    close;
% 
% if mod(c,10) == 0
%     drawnow
%     input("press enter when done looking at figures")
%     close all
% end
end



%%
Cell_numbers=[];
csm_mod=length(csm_cells);
behav_mod=length(filtered_modulated_indices);
perc_mod=(csm_mod/behav_mod)*100;
% rncells=length(rn_unit_indices);
% nocells=length(no_unit_indices);
Cell_numbers.csm_modulated=csm_mod;
Cell_numbers.reach_modulated=behav_mod;
Cell_numbers.perc_modulated=perc_mod;
% Cell_numbers.rn_cells=rncells;
% Cell_numbers.no_cells=nocells;

%% Get corelations for NC cells that are modulated during behavior 

FRCorrStats_edit = struct;

cn = fieldnames(FRBins);
Vd_corr = [];
Vd_corr_p = [];
corrctrl = [];

%for i = 1:length(rn_unit_indices)
%for i=1:length(csm_cells)
for i = 1:length(filtered_modulated_indices)
   % unit_name = cn{rn_unit_indices(i)}; % Get corresponding field name
    unit_name = cn{filtered_modulated_indices(i)};

    allkins = FRBins.(unit_name).reach_ts.kin_subt(:,5);
    allfrs = FRBins.(unit_name).reach_ts.frs;

    % Get correlation from control
    allkins_c = FRBins.(unit_name).reach_ts.kin_subt(:,5);
    allfrs_c = FRBins.(unit_name).reach_ts.frs;
    [vd_corr_ctrl, vd_corr_p_ctrl] = corr(allfrs_c, allkins_c);
    FRCorrStats_edit.(unit_name).corrcoef_ctrl = vd_corr_ctrl;

    corrctrl = [corrctrl; vd_corr_ctrl];

    % Fit a linear model
    mdl = fitlm(allfrs, allkins);
    [vd_corr, vd_corr_p] = corr(allfrs, allkins);
    Vd_corr = [Vd_corr; vd_corr];
    Vd_corr_p = [Vd_corr_p; vd_corr_p];

    % Store stats properly
    FRCorrStats_edit.(unit_name).coeffs = mdl.Coefficients.Estimate(:);
    FRCorrStats_edit.(unit_name).pvals = mdl.Coefficients.pValue(:);
    FRCorrStats_edit.(unit_name).corrcoef = vd_corr;
    FRCorrStats_edit.(unit_name).corrcoef_pvalue = vd_corr_p;

    if any(FRCorrStats_edit.(unit_name).pvals >= .05)
        disp(unit_name) % Display unit name if any p-value >= 0.05
    end

    ft = flipud(FRCorrStats_edit.(unit_name).coeffs);
    all_coeffs(i, :) = ft;

    xs = [min(allfrs) max(allfrs)];
    ys = polyval(ft, xs);
end
%%
figure; hold on
cn = fieldnames(FRBins);
% Plot results
%for i = 1:length(rn_unit_indices)
for i = 1:length(filtered_modulated_indices)
%for i=1:length(csm_cells)
    %unit_name = cn{csm_cells(i)}; 
    %unit_name = cn{rn_unit_indices(i)}; % Get corresponding field name
    unit_name = cn{filtered_modulated_indices(i)};
 
    allkins = FRBins.(unit_name).reach_ts.kin_subt(:,5);
    allfrs = FRBins.(unit_name).reach_ts.frs;

    ft = flipud(FRCorrStats_edit.(unit_name).coeffs);
    xs = [min(allfrs) max(allfrs)];
    ys = polyval(ft, xs);
    all_ps(i) = FRCorrStats_edit.(unit_name).pvals(2);

    % Assign color based on correlation
    if FRCorrStats_edit.(unit_name).corrcoef_ctrl > 0
        corr_col = [1.0000, 0.3867, 0.2773]; % Positive correlation
    else
        corr_col = [0.2734, 0.5078, 0.7031]; % Negative correlation
    end

    plot(xs, ys, 'LineWidth', 2, 'Color', corr_col);
    mid_x = mean(xs); % Find a suitable x position (midpoint)
    mid_y = polyval(ft, mid_x); % Compute corresponding y position
    text(mid_x, mid_y, unit_name, 'Color', corr_col, 'FontSize', 10, 'FontWeight', 'bold');
end

% Compute and plot mean line
allc_mean = mean(all_coeffs);
xs = [0 600];
ys = polyval(allc_mean, xs);
xlim([0 400])
title(['b ' num2str(allc_mean(2)) '; m ' num2str(allc_mean(1))])
xlabel('Firing Rate Hz')
ylabel('V Outward (cm/s)')

 %% Now plot pos and negative cells 
pos_corr_cells = []; % Stores indices of positively correlated cells
neg_corr_cells = []; % Stores indices of negatively correlated cells
cn=fieldnames(FRCorrStats_edit);
for c = 1:length(cn)
    if FRCorrStats_edit.(cn{(c)}).corrcoef_ctrl > 0  
        pos_corr_cells = [pos_corr_cells; c]; % Store positive correlation indices
    else
        neg_corr_cells = [neg_corr_cells; c]; % Store negative correlation indices
    end
end

%% Plot normalized
baseline_epoch = 1000;
event_epoch = 1000;

% Choose event
event_name = "pos_thresh";
enames = fieldnames(Events_T2);
evnt = find(ismember(enames, event_name));


% extract_psth = @(cell_indices) arrayfun(@(c) struct(...
%     'T0', data.units(rn_unit_indices(c)).T0_pos_thresh.PSTH_smooth, ...
%     'T1', data.units(rn_unit_indices(c)).T1_pos_thresh.PSTH_smooth, ...
%     'T2', data.units(rn_unit_indices(c)).T2_pos_thresh.PSTH_smooth, ...
%     'T3', data.units(rn_unit_indices(c)).T3_pos_thresh.PSTH_smooth ...
% ), cell_indices, 'UniformOutput', false);

extract_psth = @(cell_indices) arrayfun(@(c) struct(...
    'T0', data.units(filtered_modulated_indices(c)).NoCSM_csm.PSTH_smooth, ...
    'T2', data.units(filtered_modulated_indices(c)).T2_csm.PSTH_smooth ...
), cell_indices, 'UniformOutput', false);

% extract_psth = @(cell_indices) arrayfun(@(c) struct(...
%     'T0', data.units(rn_unit_indices(c)).NoCSM_csm.PSTH_smooth, ...
%     'T2', data.units(rn_unit_indices(c)).T2_csm.PSTH_smooth ...
% ), cell_indices, 'UniformOutput', false);




% Get PSTH data for positive and negative correlation groups
psth_pos = extract_psth(pos_corr_cells);
psth_neg = extract_psth(neg_corr_cells);

% Normalize firing rates for each group
normalize_psth = @(psth_group) cellfun(@(p) struct( ...
    'T0', p.T0 / max([p.T0(:);p.T2(:)]), ...
    'T2', p.T2 / max([p.T0(:);p.T2(:)]) ...
), psth_group, 'UniformOutput', false);

psth_pos = normalize_psth(psth_pos);
psth_neg = normalize_psth(psth_neg);

% Compute averages and SEM
compute_avg_sem = @(psth_group, field) struct( ...
    'avg', mean(cell2mat(cellfun(@(p) p.(field)(:)', psth_group, 'UniformOutput', false)), 1, 'omitnan'), ...
    'sem', std(cell2mat(cellfun(@(p) p.(field)(:)', psth_group, 'UniformOutput', false)), 0, 1, 'omitnan') / sqrt(length(psth_group)) ...
);

psth_T0_pos = compute_avg_sem(psth_pos, 'T0');
%psth_T1_pos = compute_avg_sem(psth_pos, 'T1');
psth_T2_pos = compute_avg_sem(psth_pos, 'T2');
%psth_T3_pos = compute_avg_sem(psth_pos, 'T3');

psth_T0_neg = compute_avg_sem(psth_neg, 'T0');
%psth_T1_neg = compute_avg_sem(psth_neg, 'T1');
psth_T2_neg = compute_avg_sem(psth_neg, 'T2');
%psth_T3_neg = compute_avg_sem(psth_neg, 'T3');

% Define time axis
time_axis = linspace(-baseline_epoch, event_epoch, length(psth_T0_pos.avg));
%%
% Plot positive and negative correlation groups
figure('Position', [100, 100, 800, 900]); % Adjusts the figure size (width x height)
subplot(3,1,1); hold on
dim = 8;

StimvelPlot(NoRR_T2aligned,dim,NoCSMcolor,0)
StimvelPlot(ST2_CSM,dim,T2color,0)

xlim([-.25, .25])
ylim([-60 60])
yline(0, '--')
yline(0,'--k')
xlim([-.25 .25])

xline(0, '--');
xline(0.05, '--');
%ylim([-0.5,2])
%ylim([0 0.2])
%ylim([-50 30])
ylim('auto')
xlabel('Time from threshold crossing (s)')
ylabel('Outward Accel. (dm/s)')


subplot(3,1,2); hold on;
fill([time_axis, fliplr(time_axis)], [psth_T0_pos.avg + psth_T0_pos.sem, fliplr(psth_T0_pos.avg - psth_T0_pos.sem)],[0.8438    0.7461    0.8438], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
%fill([time_axis, fliplr(time_axis)], [psth_T1_pos.avg + psth_T1_pos.sem, fliplr(psth_T1_pos.avg - psth_T1_pos.sem)],[0.7266    0.3320    0.8242] , 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([time_axis, fliplr(time_axis)], [psth_T2_pos.avg + psth_T2_pos.sem, fliplr(psth_T2_pos.avg - psth_T2_pos.sem)],[0.4141    0.3516    0.8008] , 'FaceAlpha', 0.2, 'EdgeColor', 'none');
%fill([time_axis, fliplr(time_axis)], [psth_T3_pos.avg + psth_T2_pos.sem, fliplr(psth_T3_pos.avg - psth_T3_pos.sem)],[ 0.2812    0.2383    0.5430], 'FaceAlpha', 0.2, 'EdgeColor', 'none');

plot(time_axis, psth_T0_pos.avg,'Color',NoCSMcolor, 'LineWidth', 2);
%plot(time_axis, psth_T1_pos.avg, 'Color',T1color, 'LineWidth', 2);
plot(time_axis, psth_T2_pos.avg, 'Color',T2color, 'LineWidth', 2);
%plot(time_axis, psth_T3_pos.avg, 'Color',T3color, 'LineWidth', 2);
xlabel('Time from threshold crossing (ms)');
ylabel('Normalized Firing Rate');
title('Positively Correlated Cells');
xlim([-250 250]);
ylim([0 1])
xline(0, '--');
xline(50, '--');
grid off;

subplot(3,1,3); hold on;
fill([time_axis, fliplr(time_axis)], [psth_T0_neg.avg + psth_T0_neg.sem, fliplr(psth_T0_neg.avg - psth_T0_neg.sem)],[0.8438    0.7461    0.8438], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
%fill([time_axis, fliplr(time_axis)], [psth_T1_neg.avg + psth_T1_neg.sem, fliplr(psth_T1_neg.avg - psth_T1_neg.sem)],[0.7266    0.3320    0.8242], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([time_axis, fliplr(time_axis)], [psth_T2_neg.avg + psth_T2_neg.sem, fliplr(psth_T2_neg.avg - psth_T2_neg.sem)],[0.4141    0.3516    0.8008], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
%fill([time_axis, fliplr(time_axis)], [psth_T3_neg.avg + psth_T3_neg.sem, fliplr(psth_T3_neg.avg - psth_T3_neg.sem)],[ 0.2812    0.2383    0.5430], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(time_axis, psth_T0_neg.avg,'Color',NoCSMcolor, 'LineWidth', 2);
%plot(time_axis, psth_T1_neg.avg,'Color',T1color, 'LineWidth', 2);
plot(time_axis, psth_T2_neg.avg,'Color',T2color, 'LineWidth', 2);
%plot(time_axis, psth_T3_neg.avg,'Color',T3color, 'LineWidth', 2);
xlabel('Time from threshold crossing (ms)');
ylabel('Normalized Firing Rate');
title('Negatively Correlated Cells');
xlim([-250 250]);
ylim([0 1])
xline(0, '--');
xline(50, '--');
grid off;


%%

cd("C:\Users\PersonLab\Documents\NeuropixelsResults\ShamAnimals\ReachandCSMModulated\GLM Method")
dataname = '220509_SAM006_CSMNC';
save(dataname, 'data','Reaches','ExcitRR', 'blk_track','blk_len','dta','suc','Salt_index','opto_index','Q','Q_frs','Q_reach','Kins','FRBins','FRCorrStats_edit','filtered_modulated_indices','pos_corr_cells','neg_corr_cells','csm_cells','Cell_numbers','-v7.3');

%Reaches includes all reaches 