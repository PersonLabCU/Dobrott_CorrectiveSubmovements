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

%% Get Events times for each set of reaches 

Events_All = GetEvents_NP(Reaches);

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

%
%% Calculate rasters and PSTHs for all trials aligned to endpoint
% Define baseline duration and event duration (time pre and post 'event') in ms
baseline_epoch = 1000;
event_epoch = 1000;

% Choose event to align to
event_name = ["end_point"]; 
enames = fieldnames(Events_All);
evnt = find(ismember(enames, event_name));

% Process spike data
for c = 1:length(data.units)
    spike_times = data.units(c).spike_times * 1000; % convert to ms
    
    for enum = 1:length(evnt)
        event_times = Events_All.(enames{evnt(enum)}) * 1000;
        window = [baseline_epoch; event_epoch];
        
        [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
        
        data.units(c).(enames{evnt(enum)}).Referenced = Ref;
        data.units(c).(enames{evnt(enum)}).EventTrial = EvTr;
        data.units(c).(enames{evnt(enum)}).Raster = BinMast;
        
        [psth_raw, psth_smooth] = Psth1(BinMast, 5);
        data.units(c).(enames{evnt(enum)}).PSTH_raw = psth_raw;
        data.units(c).(enames{evnt(enum)}).PSTH_smooth = psth_smooth;
    end
end


%%
epo = 500;
%epo = 250;

evnt = "end_point";



nc_rows=nc_rows';
%folder_path="C:\Users\PersonLab\Documents\NeuropixelsResults\ShamAnimals\050922_SAM009\mod cells";
%folder_path= "C:\Users\Courtney\Documents\NeuropixelsResults\Cell Catalog\051022_SAM012\Vmax_raster";

%for c=filtered_modulated_indices
for c=1
    %figure('Position', [100, 100, 500, 600]); hold on
    figure('Position', [100, 100, 500, 400]); hold on
    %figure; hold on
    subplot(2,1,1); hold on
    title(['cluster ' num2str(data.units(c).neurontype) ' depth ' num2str(data.units(c).depth)], 'Interpreter', 'none');

    % Binary raster matrix
    rastr = data.units(c).end_point.Raster(:,(1000-epo):(2000-(1000-epo)));

    % Find spike coordinates
    [y_co, x_co] = find(rastr);
    x_co = x_co - epo;

    % Set spike tick parameters
    tick_length = 1;  % height of each tick
    spike_tick_width = 0.5;

    % Vectorize spike ticks for `plot()`
    X = [x_co'; x_co'; nan(size(x_co'))];
    Y = [y_co'-tick_length/2; y_co'+tick_length/2; nan(size(y_co'))];
    plot(X(:), Y(:), 'Color', [0.5 0.5 0.5], 'LineWidth', spike_tick_width);

    % Add alignment event line
    event_tick_width = 2;
    plot([0 0], [0.5 size(rastr, 1) + 0.5], 'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);

    axis ij;
    ylim([0, size(rastr, 1)]);
    xlim([-epo, epo]);

    subplot(2,1,2); hold on
    title('peth')

    %ctrl
    psth_smooth = data.units(c).end_point.PSTH_smooth;
    psth_raw = data.units(c).end_point.PSTH_raw;


    td1 = length(psth_raw)-length(psth_smooth)+1;
    td2 = td1/2;

    plot((-epo:epo),psth_smooth(((1000-epo)-td2):((2000-(1000-epo))-td2)),'Color',[0.5 0.5 0.5],'LineWidth',1.5)  
    xlim([-epo,epo])
  
    xticks(-epo:250:epo)

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
% 
% T1R=AlignThresh(T1,thresh,dirdir,500);
% T1R=stimvel_accel(T1R);
% T1R=stimvel_accel_jerk(T1R);

T2R=AlignThresh(T2,thresh,dirdir,500);
T2R=stimvel_accel(T2R);
T2R=stimvel_accel_jerk(T2R);

T3R=AlignThresh(T3,thresh,dirdir,500);
T3R=stimvel_accel(T3R);
T3R=stimvel_accel_jerk(T3R);

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

%% Get Events times for each set of reaches 
% Base

Events_NoCSM=GetEvents_NP(NoCSM);
 Events_T0= GetEvents_NP(T0);
% Events_T1= GetEvents_NP(T1);
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


%% other way to look for csm cells - GLM which seems to work better
% Use GLM to identify significantly modulated cells during CSM
% Define post-CSM window
post_csm_start = -50;   
post_csm_end = 150;   

% Define time vector for 1 ms bins from -500 to 500 ms
time_vector = -500:1:500;  
post_csm_indices = find(time_vector >= post_csm_start & time_vector <= post_csm_end);

% Initialize array to store indices of significantly modulated cells
csm_cellsT2 = [];

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
        csm_cellsT2 = [csm_cellsT2, c]; 
    end

    % Display results
    disp(['Unit ', num2str(c), ': Poisson GLM p = ', num2str(p_value)]);
end

% Display the indices of significant cells
disp('Indices of significantly modulated cells (GLM):');
disp(csm_cellsT2);


%% other way to look for csm cells - T0 vs NOCSM
% Use GLM to identify significantly modulated cells during CSM
% Define post-CSM window
post_csm_start = -50;   
post_csm_end = 150;   

% Define time vector for 1 ms bins from -500 to 500 ms
time_vector = -500:1:500;  
post_csm_indices = find(time_vector >= post_csm_start & time_vector <= post_csm_end);

% Initialize array to store indices of significantly modulated cells
csm_cellsT0 = [];

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
        csm_cellsT0 = [csm_cellsT0, c]; 
    end

    % Display results
    disp(['Unit ', num2str(c), ': Poisson GLM p = ', num2str(p_value)]);
end

% Display the indices of significant cells
disp('Indices of significantly modulated cells (GLM):');
disp(csm_cellsT0);


%%
epo = 500; % Epoch length

T0color=rgb('Thistle');
T1color=rgb('MediumOrchid');
T2color=rgb('SlateBlue');
T3color=rgb('DarkSlateBlue');
NoCSMcolor=rgb('OrangeRed');

% Define heights
small_h = 0.08;
large_h = 0.18;
gap = 0.02;
width = 0.35;

%folder_path ="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\ShamAnimals\050522_SAM009\t0VSt2csm";
for c = 1
    % Create a new figure for each unit
    figure('Position', [100, 100, 800, 700]); hold on

   
    % Plot Stimulus aligned to T0 (subplot 1)
    subplot(4, 2, 1); hold on
    %axes('Position', [0.1, 0.95 - small_h, width, small_h]);
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
    subplot(4, 2, 3); hold on
    %axes('Position', [0.1, 0.95 - small_h - gap - large_h, width, large_h]);
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
    subplot(4, 2, 5); hold on
    %axes('Position', [0.1, 0.95 - small_h - gap - large_h - gap - large_h, width, large_h]);
    title('T0 csm')

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
    subplot(4, 2, 7); hold on
   % axes('Position', [0.1, 0.95 - small_h - gap - large_h - gap - large_h - gap - large_h, width, large_h]);
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

    % % Aligned to T2
        % Plot Stimulus aligned to T0 (subplot 1)

    subplot(4, 2, 2); hold on
    %axes('Position', [0.55, 0.95 - small_h, width, small_h]);
    h = title(['PSTH for Unit ', num2str(c), ' ' 'depth' ' ' num2str(data.units(c).depth)]);
    set(h, 'interpreter', 'none')

    dim = 8;
    StimvelPlot(NoRR_T2aligned, dim, NoCSMcolor, 0)
    StimvelPlot(ST2_CSM, dim, T2color, 0)

    xlabel('CSM Aligned-T2')
    ylim('auto')
    xlim([-.5, .5])
    yline(0, '--')
    
    % Plot raster aligned to NoCSM condition (subplot 2)
     subplot(4, 2, 4); hold on
    %axes('Position', [0.55, 0.95 - small_h - gap - large_h, width, large_h]);
    title('NoCSM Aligned to T2')

    tick_length = 1; % length of tick marks for spikes
    spike_tick_width = 0.5; % width of lines marking spikes
    event_tick_width = 2; % width of lines marking secondary events

    % Extract raster for NoCSM condition
    rastr = data.units(c).NoCSM_csmT2.Raster(:, (1000 - epo):(2000 - (1000 - epo)));
    [y_co, x_co] = find(rastr);
    x_co = x_co - epo; % Align to the event
    line([x_co'; x_co'], [y_co' - (tick_length / 2); y_co' + (tick_length / 2)], 'Color', [0.5 0.5 0.5], 'LineWidth', spike_tick_width);
    line([zeros(1, size(rastr, 1)); zeros(1, size(rastr, 1))], [(1:1:size(rastr, 1)) - 0.5; (1:1:size(rastr, 1)) + 0.5], 'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);
    axis ij;
    %ylim([0, size(rastr, 1)])
    ylim([0, size(rastr, 1)] +0.5) % Adjust the ylim based on the height scaling
    xlim([-epo, epo])
    
    % Plot raster aligned to T0_csm condition (subplot 3)
    subplot(4, 2, 6); hold on
    %axes('Position', [0.55, 0.95 - small_h - gap - large_h - gap - large_h, width, large_h]);
    title('T2 csm')

    % Extract raster for T0_csm condition
    rastr = data.units(c).T2_csm.Raster(:, (1000 - epo):(2000 - (1000 - epo)));
    [y_co, x_co] = find(rastr);
    x_co = x_co - epo; % Align to the event
    line([x_co'; x_co'], [y_co' - (tick_length / 2); y_co' + (tick_length / 2)], 'Color', [0.5 0.5 0.5], 'LineWidth', spike_tick_width);
    line([zeros(1, size(rastr, 1)); zeros(1, size(rastr, 1))], [(1:1:size(rastr, 1)) - 0.5; (1:1:size(rastr, 1)) + 0.5], 'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);
    axis ij;
    ylim([0, size(rastr, 1)]+ 0.5)
    xlim([-epo, epo])

    % --- PSTH Overlay (subplot 4) ---
    subplot(4, 2, 8); hold on
    %axes('Position', [0.55, 0.95 - small_h - gap - large_h - gap - large_h - gap - large_h, width, large_h]);
    % Extract PSTH data (already smoothed)
    psth_NoCSM_smooth = data.units(c).NoCSM_csmT2.PSTH_smoothbin;
    psth_CSM_smooth = data.units(c).T2_csm.PSTH_smoothbin;
    
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
    plot(time_vector, psth_CSM_gaussian, 'Color', T2color, 'LineWidth', 1.5); % Smooth PSTH for CSM
    xlabel('Time (ms)');
    ylabel('Firing Rate (Hz)');
    title(['PSTH for Unit ', num2str(c)]);
    xline(0, '--k', 'LineWidth', 2);
    xlim([-500 500]); % Adjust to view from -500 to 500 ms
    grid off;


    % Optionally close the figure after viewing
    if mod(c, 10) == 0
        drawnow
        input("Press Enter when done looking at figures");
        close all
    end
end

%% 
epo = 500; % Epoch length

T0color = rgb('Thistle');
T1color = rgb('MediumOrchid');
T2color = rgb('SlateBlue');
T3color = rgb('DarkSlateBlue');
NoCSMcolor = rgb('OrangeRed');

% Define spike tick parameters
tick_length = 1;  % height of each tick
spike_tick_width = 0.5;
event_tick_width = 2;

for c = 1
    % Create figure
    figure('Position', [100, 100, 800, 700]); hold on

    % --- PSTH Left: T0-aligned --- %%
    subplot(4, 2, 1); hold on
    title(['PSTH for Unit ', num2str(c), ' depth ', num2str(data.units(c).depth)], 'Interpreter', 'none');
    StimvelPlot(NoRR_T0aligned, 8, NoCSMcolor, 0)
    StimvelPlot(ST0_CSM, 8, T0color, 0)
    xlabel('CSM Aligned-T0')
    ylim('auto'); xlim([-.5, .5]); yline(0, '--')
    xticks(-epo:250:epo)

    % --- Raster: NoCSM T0 --- %%
    subplot(4, 2, 3); hold on
    title('NoCSM Aligned to T0')
    rastr = data.units(c).NoCSM_csmT0.Raster(:, (1000 - epo):(2000 - (1000 - epo)));
    [y_co, x_co] = find(rastr); x_co = x_co - epo;
    X = [x_co'; x_co'; nan(size(x_co'))];
    Y = [y_co'-tick_length/2; y_co'+tick_length/2; nan(size(y_co'))];
    plot(X(:), Y(:), 'Color', [0.5 0.5 0.5], 'LineWidth', spike_tick_width);
    plot([0 0], [0.5 size(rastr, 1) + 0.5], 'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);
    axis ij; ylim([0, size(rastr, 1)]+0.5); xlim([-epo, epo])
    xticks(-epo:250:epo)

    % --- Raster: T0 CSM --- %%
    subplot(4, 2, 5); hold on
    title('T0 CSM')
    rastr = data.units(c).T0_csm.Raster(:, (1000 - epo):(2000 - (1000 - epo)));
    [y_co, x_co] = find(rastr); x_co = x_co - epo;
    X = [x_co'; x_co'; nan(size(x_co'))];
    Y = [y_co'-tick_length/2; y_co'+tick_length/2; nan(size(y_co'))];
    plot(X(:), Y(:), 'Color', [0.5 0.5 0.5], 'LineWidth', spike_tick_width);
    plot([0 0], [0.5 size(rastr, 1) + 0.5], 'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);
    axis ij; ylim([0, size(rastr, 1)]+0.5); xlim([-epo, epo])
    xticks(-epo:250:epo)

    % --- PSTH: Overlay --- %%
    subplot(4, 2, 7); hold on
    psth_NoCSM_smooth = data.units(c).NoCSM_csmT0.PSTH_smoothbin;
    psth_CSM_smooth = data.units(c).T0_csm.PSTH_smoothbin;
    time_vector = linspace(-1000, 1000, length(psth_NoCSM_smooth));
    smoother = 10;
    gauss_filter = fspecial('gaussian', [1, 6 * smoother], smoother);
    psth_NoCSM_gaussian = conv(psth_NoCSM_smooth, gauss_filter, 'same');
    psth_CSM_gaussian = conv(psth_CSM_smooth, gauss_filter, 'same');
    plot(time_vector, psth_NoCSM_gaussian, 'Color', NoCSMcolor, 'LineWidth', 1.5);
    plot(time_vector, psth_CSM_gaussian, 'Color', T0color, 'LineWidth', 1.5);
    xlabel('Time (ms)'); ylabel('Firing Rate (Hz)');
    title(['PSTH Overlay for Unit ', num2str(c)]);
    xline(0, '--k', 'LineWidth', 2); xlim([-500 500]);
    xticks(-epo:250:epo)

    % --- PSTH Right: T2-aligned --- %%
    subplot(4, 2, 2); hold on
    title(['PSTH for Unit ', num2str(c), ' depth ', num2str(data.units(c).depth)], 'Interpreter', 'none');
    StimvelPlot(NoRR_T2aligned, 8, NoCSMcolor, 0)
    StimvelPlot(ST2_CSM, 8, T2color, 0)
    xlabel('CSM Aligned-T2')
    ylim('auto'); xlim([-.5, .5]); yline(0, '--')
    xticks(-epo:250:epo)
    % --- Raster: NoCSM T2 --- %%
    subplot(4, 2, 4); hold on
    title('NoCSM Aligned to T2')
    rastr = data.units(c).NoCSM_csmT2.Raster(:, (1000 - epo):(2000 - (1000 - epo)));
    [y_co, x_co] = find(rastr); x_co = x_co - epo;
    X = [x_co'; x_co'; nan(size(x_co'))];
    Y = [y_co'-tick_length/2; y_co'+tick_length/2; nan(size(y_co'))];
    plot(X(:), Y(:), 'Color', [0.5 0.5 0.5], 'LineWidth', spike_tick_width);
    plot([0 0], [0.5 size(rastr, 1) + 0.5], 'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);
    axis ij; ylim([0, size(rastr, 1)]+0.5); xlim([-epo, epo])
    xticks(-epo:250:epo)

    % --- Raster: T2 CSM --- %%
    subplot(4, 2, 6); hold on
    title('T2 CSM')
    rastr = data.units(c).T2_csm.Raster(:, (1000 - epo):(2000 - (1000 - epo)));
    [y_co, x_co] = find(rastr); x_co = x_co - epo;
    X = [x_co'; x_co'; nan(size(x_co'))];
    Y = [y_co'-tick_length/2; y_co'+tick_length/2; nan(size(y_co'))];
    plot(X(:), Y(:), 'Color', [0.5 0.5 0.5], 'LineWidth', spike_tick_width);
    plot([0 0], [0.5 size(rastr, 1) + 0.5], 'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);
    axis ij; ylim([0, size(rastr, 1)]+0.5); xlim([-epo, epo])
    xticks(-epo:250:epo)

     % --- PSTH Overlay (subplot 4) ---
    subplot(4, 2, 8); hold on
    %axes('Position', [0.55, 0.95 - small_h - gap - large_h - gap - large_h - gap - large_h, width, large_h]);
    % Extract PSTH data (already smoothed)
    psth_NoCSM_smooth = data.units(c).NoCSM_csmT2.PSTH_smoothbin;
    psth_CSM_smooth = data.units(c).T2_csm.PSTH_smoothbin;
    
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
    plot(time_vector, psth_CSM_gaussian, 'Color', T2color, 'LineWidth', 1.5); % Smooth PSTH for CSM
    xlabel('Time (ms)');
    ylabel('Firing Rate (Hz)');
    title(['PSTH for Unit ', num2str(c)]);
    xline(0, '--k', 'LineWidth', 2);
    xticks(-epo:250:epo)
    xlim([-500 500]); % Adjust to view from -500 to 500 ms
    grid off;


end