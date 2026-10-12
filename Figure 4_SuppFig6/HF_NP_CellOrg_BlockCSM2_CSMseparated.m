%% NP data org per session for population analysis

%%
%load NP data
%establish directory: navigate to 'recording' subfolder
% if using files without neurotype have to go in and change code 
pathy = uigetdir;
pathx = string(pathy);

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
%% Separate Out Reaches 
% adaptation ephys analysis
% Extract blocks and find the index so that it can be used 
rn = fieldnames(Reaches);
blk_array = [];
for r = 1:length(rn)
    blk_array(r) = Reaches.(rn{r}).block;
end
blk_array=blk_array';
%assume at least three blocks (baseline, adaptation, washout). separate
%blocks out 
% some sessions have something weird going on
bl_inx = find(blk_array == 1);
BASE = Extract(Reaches,bl_inx);
ad_inx = find(blk_array == 2);
ADAPT = Extract(Reaches,ad_inx);
wo_inx = find(blk_array == 3);
WASHALL = Extract(Reaches,wo_inx);

% use for session 20220510 SAM010
% % bl_inx = find(blk_array == 1);
% % BASE = Extract(Reaches,bl_inx);
% % ad_inx = find(blk_array == 3); % block 2 doesn't have an stims??
% % ADAPT = Extract(Reaches,ad_inx);
% % wo_inx = find(blk_array == 4);
% % WASHALL = Extract(Reaches,wo_inx);
%process reaches and identify kinematic event times for alignment

[Adapt, ctrl] = StimFilter(ADAPT);

Base=accel_toreaches(BASE);
Adapt=accel_toreaches(Adapt);
Wash=accel_toreaches(WASHALL);
Wash_E = Extract(Wash,1:5);% first 2 reaches 
  
%% Get data aligned to CSM
thresh = .1; dirdir = 1;
S=AlignStim(Adapt,500); %always check to make sure that the window is correct
S=stimvel_accel(S);

BL=AlignThresh(Base,thresh,dirdir,500);
BL=stimvel_accel(BL);

WO=AlignThresh(Wash,thresh,dirdir,500);
WO=stimvel_accel(WO);

WOE=AlignThresh(Wash_E,thresh,dirdir,500);
WOE=stimvel_accel(WOE);
%% Plot
dim = 8; 

BL_color=rgb('Black');
Adaptcolor=rgb('DodgerBlue');
WO_color=rgb('Magenta');
WOE_color=rgb('Orchid');


StimvelPlot(BL,dim,BL_color,1);
StimvelPlot(S,dim,Adaptcolor,0);
StimvelPlot(WO,dim,WO_color,0);
StimvelPlot(WOE,dim,WOE_color,0);

yline(0,'--k')
xlim([-.2 .2])
%ylim([-0.5,2])
%ylim([0 0.2])
%ylim([-50 30])
ylim('auto')
xlabel('Time from threshold (s)')
ylabel('Outward Velocity (dm/s)')
%pbaspect([1 2 1]) % for position plots 

%% Use ReReach Filter - Control Reaches 
thresh=0.1;
dirdir=1;

ExcitRR_Stimoff=[]; 
Creaches=Control;
% ReReach
[ReReaches,StimRR,MaxVel]=ReReachFilter_accel_ctrlreachesHF_edit(Creaches); % currently looking within 200ms 
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

ExcitRR_Stimoff.successrate.ReReach = length(fieldnames(S_C)) / length(fieldnames(Creaches));
ExcitRR_Stimoff.successrate.NoReReach = length(fieldnames(S_S)) / length(fieldnames(Creaches));
ExcitRR_Stimoff.BL=Creaches;
ExcitRR_Stimoff.BLN=C;
ExcitRR_Stimoff.ReReaches=ReReaches;
ExcitRR_Stimoff.RRKin=RRKin;
ExcitRR_Stimoff.NoReReach=NoReach;
ExcitRR_Stimoff.NRKin=NRKin;
ExcitRR_Stimoff.BLN_N=CN;
ExcitRR_Stimoff.StimRR=StimRR;
ExcitRR_Stimoff.MaxVel=MaxVel;
% get % of reaches with a rereach
cnum=length(fields1);
rrnum=length(fields2);
percentrr=rrnum/cnum;
ExcitRR_Stimoff.Cnum=cnum;
ExcitRR_Stimoff.RRnum=rrnum;
ExcitRR_Stimoff.PercentRR=percentrr;
%% Use ReReach Filter - Control Reaches 
thresh=0.1;
dirdir=1;

ExcitRR_BL=[]; 
Creaches=Base;
% ReReach
[ReReaches,StimRR,MaxVel]=ReReachFilter_accel_ctrlreachesHF_edit(Creaches); % currently looking within 200ms 
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

ExcitRR_BL.successrate.ReReach = length(fieldnames(S_C)) / length(fieldnames(Creaches));
ExcitRR_BL.successrate.NoReReach = length(fieldnames(S_S)) / length(fieldnames(Creaches));
ExcitRR_BL.BL=Creaches;
ExcitRR_BL.BLN=C;
ExcitRR_BL.ReReaches=ReReaches;
ExcitRR_BL.RRKin=RRKin;
ExcitRR_BL.NoReReach=NoReach;
ExcitRR_BL.NRKin=NRKin;
ExcitRR_BL.BLN_N=CN;
ExcitRR_BL.StimRR=StimRR;
ExcitRR_BL.MaxVel=MaxVel;
% get % of reaches with a rereach
cnum=length(fields1);
rrnum=length(fields2);
percentrr=rrnum/cnum;
ExcitRR_BL.Cnum=cnum;
ExcitRR_BL.RRnum=rrnum;
ExcitRR_BL.PercentRR=percentrr;


%% Use ReReach Filter- STIM reaches 
ExcitRR_AD=[]; % Savethis structure to use for graphing rereach timing in sep. script

Sreaches=Adapt;
  

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

ExcitRR_AD.successrate.ReReach = length(fieldnames(S_C)) / length(fieldnames(Sreaches));
ExcitRR_AD.successrate.NoReReach = length(fieldnames(S_S)) / length(fieldnames(Sreaches));
ExcitRR_AD.AD=Sreaches;
ExcitRR_AD.ADP=S;
ExcitRR_AD.ReReaches=ReReaches;
% ExcitRR_AD.RRKin=RRKin;
% ExcitRR_AD.NoReReach=NoReach;
% ExcitRR_AD.NRKin=NRKin;
%ExcitRR_AD.ADP_N=SN;
ExcitRR_AD.StimRR=StimRR;
ExcitRR_AD.MaxVel=MaxVel;
snum=length(fields1);
rrnum=length(fields2);
percentrr=rrnum/snum;
ExcitRR_AD.Snum=snum;
ExcitRR_AD.RRnum=rrnum;
ExcitRR_AD.PercentRR=percentrr;


%% Use ReReach Filter -Washout Reaches 

ExcitRR_WO=[]; % Savethis structure to use for graphing rereach timing in sep. script

%Creaches=Wash;
Creaches=Wash;% First 5 reaches only 

% ReReach
[ReReaches,StimRR,MaxVel]=ReReachFilter_accel_ctrlreachesHF_edit(Creaches); % currently looking within 200ms 
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
   % CN=stimvel_accel(CN); %causing issues in the code 
end 

[S_C, F_C] = SuccessFilter(ReReaches);
[S_S, F_S] = SuccessFilter(NoReach);

ExcitRR_WO.successrate.ReReach = length(fieldnames(S_C)) / length(fieldnames(Creaches));
ExcitRR_WO.successrate.NoReReach = length(fieldnames(S_S)) / length(fieldnames(Creaches));
ExcitRR_WO.WO=Creaches;
ExcitRR_WO.WSH=C;
ExcitRR_WO.ReReaches=ReReaches;
ExcitRR_WO.RRKin=RRKin;
ExcitRR_WO.NoReReach=NoReach;
ExcitRR_WO.NRKin=NRKin;
ExcitRR_WO.WSH_N=CN;
ExcitRR_WO.StimRR=StimRR;
ExcitRR_WO.MaxVel=MaxVel;
% get % of reaches with a rereach
cnum=length(fields1);
rrnum=length(fields2);
percentrr=rrnum/cnum;
ExcitRR_WO.Cnum=cnum;
ExcitRR_WO.RRnum=rrnum;
ExcitRR_WO.PercentRR=percentrr;

%% Use ReReach Filter -Washout Early Reaches 

ExcitRR_WOE=[]; % Savethis structure to use for graphing rereach timing in sep. script

%Creaches=Wash;
Creaches=Wash_E;% First 5 reaches only 

% ReReach
[ReReaches,StimRR,MaxVel]=ReReachFilter_accel_ctrlreachesHF_edit(Creaches); % currently looking within 200ms 
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
   % CN=stimvel_accel(CN); %causing issues in the code 
end 

[S_C, F_C] = SuccessFilter(ReReaches);
[S_S, F_S] = SuccessFilter(NoReach);

ExcitRR_WOE.successrate.ReReach = length(fieldnames(S_C)) / length(fieldnames(Creaches));
ExcitRR_WOE.successrate.NoReReach = length(fieldnames(S_S)) / length(fieldnames(Creaches));
ExcitRR_WOE.WO=Creaches;
ExcitRR_WOE.WSH=C;
ExcitRR_WOE.ReReaches=ReReaches;
ExcitRR_WOE.RRKin=RRKin;
ExcitRR_WOE.NoReReach=NoReach;
ExcitRR_WOE.NRKin=NRKin;
ExcitRR_WOE.WSH_N=CN;
ExcitRR_WOE.StimRR=StimRR;
ExcitRR_WOE.MaxVel=MaxVel;
% get % of reaches with a rereach
cnum=length(fields1);
rrnum=length(fields2);
percentrr=rrnum/cnum;
ExcitRR_WOE.Cnum=cnum;
ExcitRR_WOE.RRnum=rrnum;
ExcitRR_WOE.PercentRR=percentrr;
%% plot histo of times 
bs=ExcitRR_BL.StimRR(:,1);
ad=ExcitRR_AD.StimRR(:,1);
wo=ExcitRR_WO.StimRR(:,1);
woe=ExcitRR_WOE.StimRR(:,1);
blcolor=rgb('Black');
adcolor=rgb('DodgerBlue');
wocolor=rgb('Magenta');
woecolor=rgb('OrangeRed');

figure;
histogram(bs,'BinWidth',0.01,'Normalization','Count','EdgeColor','none','FaceColor',blcolor)
%xlim([0 0.2])
hold on
histogram(ad,'BinWidth',0.01,'Normalization','Count','EdgeColor','none','FaceColor',adcolor)
%xlim([0 0.2])
hold on
histogram(wo,'BinWidth',0.01,'Normalization','Count','EdgeColor','none','FaceColor',wocolor)
hold on
histogram(woe,'BinWidth',0.01,'Normalization','Count','EdgeColor','none','FaceColor',woecolor)
xlim([0 0.2])


%% -----------  CSM vs NoCSM  -----------------
Stimoff_CSM   = ExcitRR_Stimoff.ReReaches;
Stimoff_NoCSM = ExcitRR_Stimoff.NoReReach;

BL_CSM   = ExcitRR_BL.ReReaches;
BL_NoCSM = ExcitRR_BL.NoReReach;

AD_CSM = ExcitRR_AD.ReReaches;
if isfield(ExcitRR_AD, 'NoReReach') && isstruct(ExcitRR_AD.NoReReach) && ~isempty(fieldnames(ExcitRR_AD.NoReReach))
    AD_NoCSM = ExcitRR_AD.NoReReach;
else
    AD_NoCSM = [];  % nothing to analyse
end

WO_CSM   = ExcitRR_WO.ReReaches;
WO_NoCSM = ExcitRR_WO.NoReReach;

WOE_CSM = ExcitRR_WOE.ReReaches;
if isfield(ExcitRR_WOE, 'NoReReach') && isstruct(ExcitRR_WOE.NoReReach) && ~isempty(fieldnames(ExcitRR_WOE.NoReReach))
    WOE_NoCSM = ExcitRR_WOE.NoReReach;
else
    WOE_NoCSM = [];  % nothing to analyse
end


%% -----------  Get event-times for each reach set  -----------------
% All non-stimmed trials
Events_Ctrl = GetEvents_NP(Control);
Events_StimoffCSM = GetEvents_NP(Stimoff_CSM);
Events_StimoffNoCSM = GetEvents_NP(Stimoff_NoCSM);

% Base
Events_BLCSM   = GetEvents_NP(BL_CSM);
if isstruct(BL_NoCSM) && ~isempty(fieldnames(BL_NoCSM))
    Events_BLNoCSM = GetEvents_NP(BL_NoCSM);
else
    Events_BLNoCSM = struct(); % Or [] if that's what downstream expects
end

% Stim – only if AD_NoCSM is valid
Events_ADCSM = GetEvents_NP(AD_CSM);
if isstruct(AD_NoCSM) && ~isempty(fieldnames(AD_NoCSM))
    Events_ADNoCSM = GetEvents_NP(AD_NoCSM);
else
    Events_ADNoCSM = struct(); % Or [] if that's what downstream expects
end

% Washout
Events_WOCSM   = GetEvents_NP(WO_CSM);

if isstruct(WO_NoCSM) && ~isempty(fieldnames(WO_NoCSM))
    Events_WONoCSM = GetEvents_NP(WO_NoCSM);
else
    Events_WONoCSM = struct(); % Or [] depending on what your code expects
end

Events_WOECSM = GetEvents_NP(WOE_CSM);
if isstruct(WOE_NoCSM) && ~isempty(fieldnames(WOE_NoCSM))
    Events_WOENoCSM = GetEvents_NP(WOE_NoCSM);
else
    Events_WOENoCSM = struct(); % Or [] depending on what your code expects
end
%% add average CSM crossing time for Baseline -T2- 
% % CSM reaches 

reach_ind = fieldnames(Stimoff_CSM);

thresh=0.1; % For headfixed reaches 

for reachnum = 1:length(reach_ind)
    reach = Stimoff_CSM.(reach_ind{reachnum}).kin;
    stim_inx = find(reach(:,2)>=thresh,1); 

   % windowstart=reach(stim_inx,1); windowend=windowstart+0.087; % CSMT2 
    windowstart=reach(stim_inx,1); windowend=windowstart+0.0815; % CSMT2 
    
    differences = abs(reach(:,1) -windowend); % find value that is closest to 100 ms after stim 

    % Find the index of the minimum absolute difference
    [minDifference, index] = min(differences);
    
    Events_StimoffCSM.csmT2(reachnum) = (reach(index,1));
   
end 
%% add average CSM crossing time for Baseline -T2- 
% % NO CSM reaches 

reach_ind = fieldnames(Stimoff_NoCSM);

thresh=0.1; % For headfixed reaches 

for reachnum = 1:length(reach_ind)
    reach = Stimoff_NoCSM.(reach_ind{reachnum}).kin;
    stim_inx = find(reach(:,2)>=thresh,1); 

    %windowstart=reach(stim_inx,1); windowend=windowstart+0.087; % CSMT2 
    windowstart=reach(stim_inx,1); windowend=windowstart+0.0815; % CSMT2 
    
    differences = abs(reach(:,1) -windowend); % find value that is closest to 100 ms after stim 

    % Find the index of the minimum absolute difference
    [minDifference, index] = min(differences);
    
    Events_StimoffNoCSM.csmT2(reachnum) = (reach(index,1));
   
end 

%% add average CSM crossing time for Baseline -T0- 
% % CSM reaches 

reach_ind = fieldnames(Stimoff_CSM);

thresh=0.1; % For headfixed reaches 

for reachnum = 1:length(reach_ind)
    reach = Stimoff_CSM.(reach_ind{reachnum}).kin;
    stim_inx = find(reach(:,2)>=thresh,1); 

    windowstart=reach(stim_inx,1); windowend=windowstart+(-0.055); % CSMT2 
    
    differences = abs(reach(:,1) -windowend); % find value that is closest to 100 ms after stim 

    % Find the index of the minimum absolute difference
    [minDifference, index] = min(differences);
    
    Events_StimoffCSM.csmT0(reachnum) = (reach(index,1));
   
end 
%% add average CSM crossing time for Baseline -T0- 
% % NO CSM reaches 

reach_ind = fieldnames(Stimoff_NoCSM);

thresh=0.1; % For headfixed reaches 

for reachnum = 1:length(reach_ind)
    reach = Stimoff_NoCSM.(reach_ind{reachnum}).kin;
    stim_inx = find(reach(:,2)>=thresh,1); 

    windowstart=reach(stim_inx,1); windowend=windowstart+(-0.055); %
    
    differences = abs(reach(:,1) -windowend); % find value that is closest to 100 ms after stim 

    % Find the index of the minimum absolute difference
    [minDifference, index] = min(differences);
    
    Events_StimoffNoCSM.csmT0(reachnum) = (reach(index,1));
   
end 

%% add average CSM crossing time for Baseline -T2- 
% % CSM reaches 

reach_ind = fieldnames(BL_CSM);

thresh=0.1; % For headfixed reaches 

for reachnum = 1:length(reach_ind)
    reach = BL_CSM.(reach_ind{reachnum}).kin;
    stim_inx = find(reach(:,2)>=thresh,1); 

    %windowstart=reach(stim_inx,1); windowend=windowstart+0.087; % CSMT2 
    windowstart=reach(stim_inx,1); windowend=windowstart+0.0815; % CSMT2 
    
    differences = abs(reach(:,1) -windowend); % find value that is closest to 100 ms after stim 

    % Find the index of the minimum absolute difference
    [minDifference, index] = min(differences);
    
    Events_BLCSM.csmT2(reachnum) = (reach(index,1));
   
end 

%% add average CSM crossing time for ADAPT

% Guard against the case where AD_NoCSM is missing or empty
if isstruct(BL_NoCSM) && ~isempty(fieldnames(BL_NoCSM))

    reach_ind = fieldnames(BL_NoCSM);
    thresh    = 0.1;                % (head-fixed reaches, velocity-based)

    % Pre-allocate and loop
    Events_BLNoCSM.csmT2 = nan(1, numel(reach_ind));

    for reachnum = 1:numel(reach_ind)
        reach = BL_NoCSM.(reach_ind{reachnum}).kin;  % [time  vel  …]

        stim_inx = find(reach(:,2) >= thresh, 1, 'first');
        if isempty(stim_inx)                       % no crossing → skip
            continue
        end

        windowstart = reach(stim_inx,1);
       % windowend   = windowstart + 0.087;         % CSM-T2 (= +87 ms)
        windowend   = windowstart + 0.0815;   

        [~, idx] = min( abs(reach(:,1) - windowend) );
        Events_BLNoCSM.csmT2(reachnum) = reach(idx,1);
    end

else
    % Structure absent or empty → return safe placeholder
    Events_BLNoCSM.csmT2 = [];
end

%% add average CSM crossing time aligned to T0- BASELINE
% % CSM REACHES
reach_ind = fieldnames(BL_CSM);

thresh=0.1; % For headfixed reaches 

for reachnum = 1:length(reach_ind)
    reach = BL_CSM.(reach_ind{reachnum}).kin;
    stim_inx = find(reach(:,2)>=thresh,1); 

    windowstart=reach(stim_inx,1); windowend=windowstart+(-0.055); % CSMT0
    
    differences = abs(reach(:,1) -windowend); % find value that is closest to 100 ms after stim 

    % Find the index of the minimum absolute difference
    [minDifference, index] = min(differences);
    
    Events_BLCSM.csmT0(reachnum) = (reach(index,1));
   
end 

%% add average CSM crossing time for ADAPT

% Guard against the case where AD_NoCSM is missing or empty
if isstruct(BL_NoCSM) && ~isempty(fieldnames(BL_NoCSM))

    reach_ind = fieldnames(BL_NoCSM);
    thresh    = 0.1;                % (head-fixed reaches, velocity-based)

    % Pre-allocate and loop
    Events_BLNoCSM.csmT0 = nan(1, numel(reach_ind));

    for reachnum = 1:numel(reach_ind)
        reach = BL_NoCSM.(reach_ind{reachnum}).kin;  % [time  vel  …]

        stim_inx = find(reach(:,2) >= thresh, 1, 'first');
        if isempty(stim_inx)                       % no crossing → skip
            continue
        end

        windowstart = reach(stim_inx,1);
        windowend   = windowstart + (-0.055);

        [~, idx] = min( abs(reach(:,1) - windowend) );
        Events_BLNoCSM.csmT0(reachnum) = reach(idx,1);
    end

else
    % Structure absent or empty → return safe placeholder
    Events_BLNoCSM.csmT0 = [];
    Events_BLNoCSM.max_vx = [];
    Events_BLNoCSM.pos_thresh = [];
end
%% add average CSM crossing time for ADAPT
% % CSM REACHES
reach_ind = fieldnames(AD_CSM);

thresh=0.1; % For headfixed reaches 

for reachnum = 1:length(reach_ind)
    reach = AD_CSM.(reach_ind{reachnum}).kin;
    stim_inx = find(reach(:,2)>=thresh,1); 

   % windowstart=reach(stim_inx,1); windowend=windowstart+0.087; % CSMT2 
    windowstart=reach(stim_inx,1); windowend=windowstart+0.0815; % CSMT2
    
    differences = abs(reach(:,1) -windowend); % find value that is closest to 100 ms after stim 

    % Find the index of the minimum absolute difference
    [minDifference, index] = min(differences);
    
    Events_ADCSM.csmT2(reachnum) = (reach(index,1));
   
end 

%% add average CSM crossing time for ADAPT

% Guard against the case where AD_NoCSM is missing or empty
if isstruct(AD_NoCSM) && ~isempty(fieldnames(AD_NoCSM))

    reach_ind = fieldnames(AD_NoCSM);
    thresh    = 0.1;                % (head-fixed reaches, velocity-based)

    % Pre-allocate and loop
    Events_ADNoCSM.csmT2 = nan(1, numel(reach_ind));

    for reachnum = 1:numel(reach_ind)
        reach = AD_NoCSM.(reach_ind{reachnum}).kin;  % [time  vel  …]

        stim_inx = find(reach(:,2) >= thresh, 1, 'first');
        if isempty(stim_inx)                       % no crossing → skip
            continue
        end

        windowstart = reach(stim_inx,1);
        %windowend   = windowstart + 0.087;         % CSM-T2 (= +87 ms)
        windowend   = windowstart + 0.0815;

        [~, idx] = min( abs(reach(:,1) - windowend) );
        Events_ADNoCSM.csmT2(reachnum) = reach(idx,1);
    end

else
    % Structure absent or empty → return safe placeholder
    Events_ADNoCSM.csmT2 = [];
end

%% add average CSM crossing time aligned to T0: ADAPT
% % CSM REACHES
reach_ind = fieldnames(AD_CSM);

thresh=0.1; % For headfixed reaches 

for reachnum = 1:length(reach_ind)
    reach = AD_CSM.(reach_ind{reachnum}).kin;
    stim_inx = find(reach(:,2)>=thresh,1); 

    windowstart=reach(stim_inx,1); windowend=windowstart+(-0.055); % CSMT0
    
    differences = abs(reach(:,1) -windowend); % find value that is closest to 100 ms after stim 

    % Find the index of the minimum absolute difference
    [minDifference, index] = min(differences);
    
    Events_ADCSM.csmT0(reachnum) = (reach(index,1));
   
end 
%% add average CSM crossing time aligned to T0: ADAPT
% % NO CSM

% Guard against the case where AD_NoCSM is missing or empty
if isstruct(AD_NoCSM) && ~isempty(fieldnames(AD_NoCSM))

    reach_ind = fieldnames(AD_NoCSM);
    thresh    = 0.1;                % (head-fixed reaches, velocity-based)

    % Pre-allocate and loop
    Events_ADNoCSM.csmT0 = nan(1, numel(reach_ind));

    for reachnum = 1:numel(reach_ind)
        reach = AD_NoCSM.(reach_ind{reachnum}).kin;  % [time  vel  …]

        stim_inx = find(reach(:,2) >= thresh, 1, 'first');
        if isempty(stim_inx)                       % no crossing → skip
            continue
        end

        windowstart = reach(stim_inx,1);
        windowend   = windowstart +(-0.055);        

        [~, idx] = min( abs(reach(:,1) - windowend) );
        Events_ADNoCSM.csmT0(reachnum) = reach(idx,1);
    end

else
    % Structure absent or empty → return safe placeholder
    Events_ADNoCSM.csmT0 = [];
     Events_ADNoCSM.pos_thresh= [];
    Events_ADNoCSM.max_vx= [];
end

%% add average CSM crossing time for WASHOUT -T2
% % CSM REACHES
reach_ind = fieldnames(WO_CSM);

thresh=0.1; % For headfixed reaches 

for reachnum = 1:length(reach_ind)
    reach = WO_CSM.(reach_ind{reachnum}).kin;
    stim_inx = find(reach(:,2)>=thresh,1); 

    %windowstart=reach(stim_inx,1); windowend=windowstart+0.087; % CSMT2 
    windowstart=reach(stim_inx,1); windowend=windowstart+0.0815; % CSMT2 
    
    differences = abs(reach(:,1) -windowend); % find value that is closest to 100 ms after stim 

    % Find the index of the minimum absolute difference
    [minDifference, index] = min(differences);
    
    Events_WOCSM.csmT2(reachnum) = (reach(index,1));
   
end 

%% add average CSM crossing time for  WASHOUT EARLY

% Guard against the case where AD_NoCSM is missing or empty
if isstruct(WO_NoCSM) && ~isempty(fieldnames(WO_NoCSM))

    reach_ind = fieldnames(WO_NoCSM);
    thresh    = 0.1;                % (head-fixed reaches, velocity-based)

    % Pre-allocate and loop
    Events_WONoCSM.csmT2 = nan(1, numel(reach_ind));

    for reachnum = 1:numel(reach_ind)
        reach = WO_NoCSM.(reach_ind{reachnum}).kin;  % [time  vel  …]

        stim_inx = find(reach(:,2) >= thresh, 1, 'first');
        if isempty(stim_inx)                       % no crossing → skip
            continue
        end

        windowstart = reach(stim_inx,1);
        %windowend   = windowstart + 0.087;         % CSM-T2 (= +87 ms)
        windowend   = windowstart + 0.0815; 

        [~, idx] = min( abs(reach(:,1) - windowend) );
        Events_WONoCSM.csmT2(reachnum) = reach(idx,1);
    end

else
    % Structure absent or empty → return safe placeholder
    Events_WONoCSM.csmT2 = [];
end
%% add average CSM crossing time  to T0: WASHOUT
% % CSM REACHES
reach_ind = fieldnames(WO_CSM);

thresh=0.1; % For headfixed reaches 

for reachnum = 1:length(reach_ind)
    reach = WO_CSM.(reach_ind{reachnum}).kin;
    stim_inx = find(reach(:,2)>=thresh,1); 

    windowstart=reach(stim_inx,1); windowend=windowstart+(-0.055); % CSMT0
    
    differences = abs(reach(:,1) -windowend); % find value that is closest to 100 ms after stim 

    % Find the index of the minimum absolute difference
    [minDifference, index] = min(differences);
    
    Events_WOCSM.csmT0(reachnum) = (reach(index,1));
   
end 
%% add average CSM crossing time  to T0: WASHOUT
% % CSM REACHES
reach_ind = fieldnames(WO_NoCSM);

thresh=0.1; % For headfixed reaches 

for reachnum = 1:length(reach_ind)
    reach = WO_NoCSM.(reach_ind{reachnum}).kin;
    stim_inx = find(reach(:,2)>=thresh,1); 

    windowstart=reach(stim_inx,1); windowend=windowstart+(-0.055); % CSMT0
    
    differences = abs(reach(:,1) -windowend); % find value that is closest to 100 ms after stim 

    % Find the index of the minimum absolute difference
    [minDifference, index] = min(differences);
    
    Events_WONoCSM.csmT0(reachnum) = (reach(index,1));
   
end 
%% add average CSM crossing time for  WASHOUT EARLY

% Guard against the case where AD_NoCSM is missing or empty
if isstruct(WO_NoCSM) && ~isempty(fieldnames(WO_NoCSM))

    reach_ind = fieldnames(WO_NoCSM);
    thresh    = 0.1;                % (head-fixed reaches, velocity-based)

    % Pre-allocate and loop
    Events_WONoCSM.csmT0 = nan(1, numel(reach_ind));

    for reachnum = 1:numel(reach_ind)
        reach = WO_NoCSM.(reach_ind{reachnum}).kin;  % [time  vel  …]

        stim_inx = find(reach(:,2) >= thresh, 1, 'first');
        if isempty(stim_inx)                       % no crossing → skip
            continue
        end

        windowstart = reach(stim_inx,1);
        windowend   = windowstart +(-0.055);         % CSM-T2 (= +87 ms)

        [~, idx] = min( abs(reach(:,1) - windowend) );
        Events_WONoCSM.csmT0(reachnum) = reach(idx,1);
    end

else
    % Structure absent or empty → return safe placeholder
    Events_WONoCSM.csmT0 = [];
    Events_WONoCSM.max_vx = [];
    Events_WONoCSM.pos_thresh = [];
end
%% add average CSM crossing time for WASHOUT EARLY -T2
% % CSM REACHES
reach_ind = fieldnames(WOE_CSM);

thresh=0.1; % For headfixed reaches 

for reachnum = 1:length(reach_ind)
    reach = WOE_CSM.(reach_ind{reachnum}).kin;
    stim_inx = find(reach(:,2)>=thresh,1); 

   % windowstart=reach(stim_inx,1); windowend=windowstart+0.087; % CSMT2 
    windowstart=reach(stim_inx,1); windowend=windowstart+0.0815; % CSMT2 
    
    differences = abs(reach(:,1) -windowend); % find value that is closest to 100 ms after stim 

    % Find the index of the minimum absolute difference
    [minDifference, index] = min(differences);
    
    Events_WOECSM.csmT2(reachnum) = (reach(index,1));
   
end 


%% add average CSM crossing time for  WASHOUT EARLY

% Guard against the case where AD_NoCSM is missing or empty
if isstruct(WOE_NoCSM) && ~isempty(fieldnames(WOE_NoCSM))

    reach_ind = fieldnames(WOE_NoCSM);
    thresh    = 0.1;                % (head-fixed reaches, velocity-based)

    % Pre-allocate and loop
    Events_WOENoCSM.csmT2 = nan(1, numel(reach_ind));

    for reachnum = 1:numel(reach_ind)
        reach = WOE_NoCSM.(reach_ind{reachnum}).kin;  % [time  vel  …]

        stim_inx = find(reach(:,2) >= thresh, 1, 'first');
        if isempty(stim_inx)                       % no crossing → skip
            continue
        end

        windowstart = reach(stim_inx,1);
        %windowend   = windowstart + 0.087;         % CSM-T2 (= +87 ms)
        windowend   = windowstart + 0.0815;  

        [~, idx] = min( abs(reach(:,1) - windowend) );
        Events_WOENoCSM.csmT2(reachnum) = reach(idx,1);
    end

else
    % Structure absent or empty → return safe placeholder
    Events_WOENoCSM.csmT2 = [];
end


%% add average CSM crossing time for NoCSM aligned to T0: WASHOUT EARLY
% % CSM REACHES 
reach_ind = fieldnames(WOE_CSM);

thresh=0.1; % For headfixed reaches 

for reachnum = 1:length(reach_ind)
    reach = WOE_CSM.(reach_ind{reachnum}).kin;
    stim_inx = find(reach(:,2)>=thresh,1); 

    windowstart=reach(stim_inx,1); windowend=windowstart+(-0.055); % CSMT0
    
    differences = abs(reach(:,1) -windowend); % find value that is closest to 100 ms after stim 

    % Find the index of the minimum absolute difference
    [minDifference, index] = min(differences);
    
    Events_WOECSM.csmT0(reachnum) = (reach(index,1));
   
end 

%% add average CSM crossing time for NoCSM aligned to T0: WASHOUT EARLY


% Guard against the case where AD_NoCSM is missing or empty
if isstruct(WOE_NoCSM) && ~isempty(fieldnames(WOE_NoCSM))

    reach_ind = fieldnames(WOE_NoCSM);
    thresh    = 0.1;                % (head-fixed reaches, velocity-based)

    % Pre-allocate and loop
    Events_WOENoCSM.csmT0 = nan(1, numel(reach_ind));

    for reachnum = 1:numel(reach_ind)
        reach = WOE_NoCSM.(reach_ind{reachnum}).kin;  % [time  vel  …]

        stim_inx = find(reach(:,2) >= thresh, 1, 'first');
        if isempty(stim_inx)                       % no crossing → skip
            continue
        end

        windowstart = reach(stim_inx,1);
        windowend   = windowstart + (-0.055);         % CSM-T2 (= +87 ms)

        [~, idx] = min( abs(reach(:,1) - windowend) );
        Events_WOENoCSM.csmT0(reachnum) = reach(idx,1);
    end

else
    % Structure absent or empty → return safe placeholder
    Events_WOENoCSM.csmT0 = [];
    Events_WOENoCSM.pos_thresh= [];
    Events_WOENoCSM.max_vx= [];
end
%% ALL STIM OFF 
% Get Events for Baseline Trials
% % ALL CSM REACHES


% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["pos_thresh"]; 
%event_name = ["end_point","pos_thresh","reach_start","max_vx"]; 
enames = fieldnames(Events_StimoffCSM);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_StimoffCSM.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    % 
    data.units(c).StimoffCSM_pos_thresh(enum).Referenced = Ref;
    data.units(c).StimoffCSM_pos_thresh(enum).EventTrial = EvTr;
    data.units(c).StimoffCSM_pos_thresh(enum).Raster = BinMast;

    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).StimoffCSM_pos_thresh(enum).PSTH_raw = psth_raw;
    data.units(c).StimoffCSM_pos_thresh(enum).PSTH_smooth = psth_smooth;  

    % Compute PSTH using new computePSTH method (10ms bins + smoothing)
    bin_size = 10; % 10ms bin size
    smoothing_bins = 3; % Moving average over 3 bins

    [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);

    % Store results for the new method
    data.units(c).StimoffCSM_pos_thresh(enum).PSTH_rawbin = psth_rawbin;
    data.units(c).StimoffCSM_pos_thresh(enum).PSTH_smoothbin = psth_smoothbin;
    data.units(c).StimoffCSM_pos_thresh(enum).time_bins = time_bins;

end
end


%% Get Events for Baseline Trials- max v 
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["max_vx"]; 
%event_name = ["end_point","pos_thresh","reach_start","max_vx"]; 
enames = fieldnames(Events_StimoffCSM);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_StimoffCSM.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    % 
    data.units(c).StimoffCSM_max_vx(enum).Referenced = Ref;
    data.units(c).StimoffCSM_max_vx(enum).EventTrial = EvTr;
    data.units(c).StimoffCSM_max_vx(enum).Raster = BinMast;

    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).StimoffCSM_max_vx(enum).PSTH_raw = psth_raw;
    data.units(c).StimoffCSM_max_vx(enum).PSTH_smooth = psth_smooth;  

    % Compute PSTH using new computePSTH method (10ms bins + smoothing)
    bin_size = 10; % 10ms bin size
    smoothing_bins = 3; % Moving average over 3 bins

    [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);

    % Store results for the new method
    data.units(c).StimoffCSM_max_vx(enum).PSTH_rawbin = psth_rawbin;
    data.units(c).StimoffCSM_max_vx(enum).PSTH_smoothbin = psth_smoothbin;
    data.units(c).StimoffCSM_max_vx(enum).time_bins = time_bins;

end
end


%% BASELINE- T2 
% % CSM REACHES
% Testing out new code to plot a different type of PSTH 
baseline_epoch = 1000;
event_epoch = 1000;

% Choose event to align to
event_name = ["csmT2"]; 
enames = fieldnames(Events_StimoffCSM);
evnt = find(ismember(enames, event_name));

% Loop through all units
for c = 1:length(data.units)
    spike_times = data.units(c).spike_times * 1000; % Convert to ms
    
    % Loop through each event type
    for enum = 1:length(evnt)
        event_times = Events_StimoffCSM.(enames{evnt(enum)}) * 1000;
        window = [baseline_epoch; event_epoch];

        % Compute referenced spikes and raster binning
        [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
        data.units(c).StimoffCSM_csmT2(enum).Referenced = Ref;
        data.units(c).StimoffCSM_csmT2(enum).EventTrial = EvTr;
        data.units(c).StimoffCSM_csmT2(enum).Raster = BinMast;

        % Compute PSTH using original method
        [psth_raw, psth_smooth] = Psth1(BinMast, 5); % Was 5 before
        data.units(c).StimoffCSM_csmT2(enum).PSTH_raw = psth_raw;
        data.units(c).StimoffCSM_csmT2(enum).PSTH_smooth = psth_smooth;

        % Compute PSTH using new computePSTH method (10ms bins + smoothing)
        bin_size = 10; % 10ms bin size
        smoothing_bins = 3; % Moving average over 3 bins

       [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);
        
        % Store results for the new method
        data.units(c).StimoffCSM_csmT2(enum).PSTH_rawbin = psth_rawbin;
        data.units(c).StimoffCSM_csmT2(enum).PSTH_smoothbin = psth_smoothbin;
        data.units(c).StimoffCSM_csmT2(enum).time_bins = time_bins;

    end
end

%% BASELINE- T0
% Testing out new code to plot a different type of PSTH 
baseline_epoch = 1000;
event_epoch = 1000;

% Choose event to align to
event_name = ["csmT0"]; 
enames = fieldnames(Events_StimoffCSM);
evnt = find(ismember(enames, event_name));

% Loop through all units
for c = 1:length(data.units)
    spike_times = data.units(c).spike_times * 1000; % Convert to ms
    
    % Loop through each event type
    for enum = 1:length(evnt)
        event_times = Events_StimoffCSM.(enames{evnt(enum)}) * 1000;
        window = [baseline_epoch; event_epoch];

        % Compute referenced spikes and raster binning
        [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
        data.units(c).StimoffCSM_csmT0(enum).Referenced = Ref;
        data.units(c).StimoffCSM_csmT0(enum).EventTrial = EvTr;
        data.units(c).StimoffCSM_csmT0(enum).Raster = BinMast;

        % Compute PSTH using original method
        [psth_raw, psth_smooth] = Psth1(BinMast, 5); % Was 5 before
        data.units(c).StimoffCSM_csmT0(enum).PSTH_raw = psth_raw;
        data.units(c).StimoffCSM_csmT0(enum).PSTH_smooth = psth_smooth;

        % Compute PSTH using new computePSTH method (10ms bins + smoothing)
        bin_size = 10; % 10ms bin size
        smoothing_bins = 3; % Moving average over 3 bins

       [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);
        
        % Store results for the new method
        data.units(c).StimoffCSM_csmT0(enum).PSTH_rawbin = psth_rawbin;
        data.units(c).StimoffCSM_csmT0(enum).PSTH_smoothbin = psth_smoothbin;
        data.units(c).StimoffCSM_csmT0(enum).time_bins = time_bins;

    end
end
%% NO CSM- for all STIM OFF REACHES

% Get Events for Baseline Trials
% % ALL CSM REACHES


% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["pos_thresh"]; 
%event_name = ["end_point","pos_thresh","reach_start","max_vx"]; 
enames = fieldnames(Events_StimoffNoCSM);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_StimoffNoCSM.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    % 
    data.units(c).StimoffNoCSM_pos_thresh(enum).Referenced = Ref;
    data.units(c).StimoffNoCSM_pos_thresh(enum).EventTrial = EvTr;
    data.units(c).StimoffNoCSM_pos_thresh(enum).Raster = BinMast;

    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).StimoffNoCSM_pos_thresh(enum).PSTH_raw = psth_raw;
    data.units(c).StimoffNoCSM_pos_thresh(enum).PSTH_smooth = psth_smooth;  

    % Compute PSTH using new computePSTH method (10ms bins + smoothing)
    bin_size = 10; % 10ms bin size
    smoothing_bins = 3; % Moving average over 3 bins

    [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);

    % Store results for the new method
    data.units(c).StimoffNoCSM_pos_thresh(enum).PSTH_rawbin = psth_rawbin;
    data.units(c).StimoffNoCSM_pos_thresh(enum).PSTH_smoothbin = psth_smoothbin;
    data.units(c).StimoffNoCSM_pos_thresh(enum).time_bins = time_bins;

end
end


%% Get Events for Baseline Trials- max v 
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["max_vx"]; 
%event_name = ["end_point","pos_thresh","reach_start","max_vx"]; 
enames = fieldnames(Events_StimoffNoCSM);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_StimoffNoCSM.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    % 
    data.units(c).StimoffNoCSM_max_vx(enum).Referenced = Ref;
    data.units(c).StimoffNoCSM_max_vx(enum).EventTrial = EvTr;
    data.units(c).StimoffNoCSM_max_vx(enum).Raster = BinMast;

    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).StimoffNoCSM_max_vx(enum).PSTH_raw = psth_raw;
    data.units(c).StimoffNoCSM_max_vx(enum).PSTH_smooth = psth_smooth;  

    % Compute PSTH using new computePSTH method (10ms bins + smoothing)
    bin_size = 10; % 10ms bin size
    smoothing_bins = 3; % Moving average over 3 bins

    [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);

    % Store results for the new method
    data.units(c).StimoffNoCSM_max_vx(enum).PSTH_rawbin = psth_rawbin;
    data.units(c).StimoffNoCSM_max_vx(enum).PSTH_smoothbin = psth_smoothbin;
    data.units(c).StimoffNoCSM_max_vx(enum).time_bins = time_bins;

end
end


%% BASELINE- T2 
% % CSM REACHES
% Testing out new code to plot a different type of PSTH 
baseline_epoch = 1000;
event_epoch = 1000;

% Choose event to align to
event_name = ["csmT2"]; 
enames = fieldnames(Events_StimoffNoCSM);
evnt = find(ismember(enames, event_name));

% Loop through all units
for c = 1:length(data.units)
    spike_times = data.units(c).spike_times * 1000; % Convert to ms
    
    % Loop through each event type
    for enum = 1:length(evnt)
        event_times = Events_StimoffNoCSM.(enames{evnt(enum)}) * 1000;
        window = [baseline_epoch; event_epoch];

        % Compute referenced spikes and raster binning
        [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
        data.units(c).StimoffNoCSM_csmT2(enum).Referenced = Ref;
        data.units(c).StimoffNoCSM_csmT2(enum).EventTrial = EvTr;
        data.units(c).StimoffNoCSM_csmT2(enum).Raster = BinMast;

        % Compute PSTH using original method
        [psth_raw, psth_smooth] = Psth1(BinMast, 5); % Was 5 before
        data.units(c).StimoffNoCSM_csmT2(enum).PSTH_raw = psth_raw;
        data.units(c).StimoffNoCSM_csmT2(enum).PSTH_smooth = psth_smooth;

        % Compute PSTH using new computePSTH method (10ms bins + smoothing)
        bin_size = 10; % 10ms bin size
        smoothing_bins = 3; % Moving average over 3 bins

       [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);
        
        % Store results for the new method
        data.units(c).StimoffNoCSM_csmT2(enum).PSTH_rawbin = psth_rawbin;
        data.units(c).StimoffNoCSM_csmT2(enum).PSTH_smoothbin = psth_smoothbin;
        data.units(c).StimoffNoCSM_csmT2(enum).time_bins = time_bins;

    end
end

%% BASELINE- T0
% Testing out new code to plot a different type of PSTH 
baseline_epoch = 1000;
event_epoch = 1000;

% Choose event to align to
event_name = ["csmT0"]; 
enames = fieldnames(Events_StimoffNoCSM);
evnt = find(ismember(enames, event_name));

% Loop through all units
for c = 1:length(data.units)
    spike_times = data.units(c).spike_times * 1000; % Convert to ms
    
    % Loop through each event type
    for enum = 1:length(evnt)
        event_times = Events_StimoffNoCSM.(enames{evnt(enum)}) * 1000;
        window = [baseline_epoch; event_epoch];

        % Compute referenced spikes and raster binning
        [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
        data.units(c).StimoffNoCSM_csmT0(enum).Referenced = Ref;
        data.units(c).StimoffNoCSM_csmT0(enum).EventTrial = EvTr;
        data.units(c).StimoffNoCSM_csmT0(enum).Raster = BinMast;

        % Compute PSTH using original method
        [psth_raw, psth_smooth] = Psth1(BinMast, 5); % Was 5 before
        data.units(c).StimoffNoCSM_csmT0(enum).PSTH_raw = psth_raw;
        data.units(c).StimoffNoCSM_csmT0(enum).PSTH_smooth = psth_smooth;

        % Compute PSTH using new computePSTH method (10ms bins + smoothing)
        bin_size = 10; % 10ms bin size
        smoothing_bins = 3; % Moving average over 3 bins

       [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);
        
        % Store results for the new method
        data.units(c).StimoffNoCSM_csmT0(enum).PSTH_rawbin = psth_rawbin;
        data.units(c).StimoffNoCSM_csmT0(enum).PSTH_smoothbin = psth_smoothbin;
        data.units(c).StimoffNoCSM_csmT0(enum).time_bins = time_bins;

    end
end
%% Get Events for Baseline Trials
% % ALL CSM REACHES


% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["pos_thresh"]; 
%event_name = ["end_point","pos_thresh","reach_start","max_vx"]; 
enames = fieldnames(Events_BLCSM);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_BLCSM.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    % 
    data.units(c).BLCSM_pos_thresh(enum).Referenced = Ref;
    data.units(c).BLCSM_pos_thresh(enum).EventTrial = EvTr;
    data.units(c).BLCSM_pos_thresh(enum).Raster = BinMast;

    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).BLCSM_pos_thresh(enum).PSTH_raw = psth_raw;
    data.units(c).BLCSM_pos_thresh(enum).PSTH_smooth = psth_smooth;  

    % Compute PSTH using new computePSTH method (10ms bins + smoothing)
    bin_size = 10; % 10ms bin size
    smoothing_bins = 3; % Moving average over 3 bins

    [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);

    % Store results for the new method
    data.units(c).BLCSM_pos_thresh(enum).PSTH_rawbin = psth_rawbin;
    data.units(c).BLCSM_pos_thresh(enum).PSTH_smoothbin = psth_smoothbin;
    data.units(c).BLCSM_pos_thresh(enum).time_bins = time_bins;

end
end


%% Get Events for Baseline Trials- max v 
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["max_vx"]; 
%event_name = ["end_point","pos_thresh","reach_start","max_vx"]; 
enames = fieldnames(Events_BLCSM);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_BLCSM.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    % 
    data.units(c).BLCSM_max_vx(enum).Referenced = Ref;
    data.units(c).BLCSM_max_vx(enum).EventTrial = EvTr;
    data.units(c).BLCSM_max_vx(enum).Raster = BinMast;

    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).BLCSM_max_vx(enum).PSTH_raw = psth_raw;
    data.units(c).BLCSM_max_vx(enum).PSTH_smooth = psth_smooth;  

    % Compute PSTH using new computePSTH method (10ms bins + smoothing)
    bin_size = 10; % 10ms bin size
    smoothing_bins = 3; % Moving average over 3 bins

    [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);

    % Store results for the new method
    data.units(c).BLCSM_max_vx(enum).PSTH_rawbin = psth_rawbin;
    data.units(c).BLCSM_max_vx(enum).PSTH_smoothbin = psth_smoothbin;
    data.units(c).BLCSM_max_vx(enum).time_bins = time_bins;

end
end


%% BASELINE- T2 
% % CSM REACHES
% Testing out new code to plot a different type of PSTH 
baseline_epoch = 1000;
event_epoch = 1000;

% Choose event to align to
event_name = ["csmT2"]; 
enames = fieldnames(Events_BLCSM);
evnt = find(ismember(enames, event_name));

% Loop through all units
for c = 1:length(data.units)
    spike_times = data.units(c).spike_times * 1000; % Convert to ms
    
    % Loop through each event type
    for enum = 1:length(evnt)
        event_times = Events_BLCSM.(enames{evnt(enum)}) * 1000;
        window = [baseline_epoch; event_epoch];

        % Compute referenced spikes and raster binning
        [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
        data.units(c).BLCSM_csmT2(enum).Referenced = Ref;
        data.units(c).BLCSM_csmT2(enum).EventTrial = EvTr;
        data.units(c).BLCSM_csmT2(enum).Raster = BinMast;

        % Compute PSTH using original method
        [psth_raw, psth_smooth] = Psth1(BinMast, 5); % Was 5 before
        data.units(c).BLCSM_csmT2(enum).PSTH_raw = psth_raw;
        data.units(c).BLCSM_csmT2(enum).PSTH_smooth = psth_smooth;

        % Compute PSTH using new computePSTH method (10ms bins + smoothing)
        bin_size = 10; % 10ms bin size
        smoothing_bins = 3; % Moving average over 3 bins

       [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);
        
        % Store results for the new method
        data.units(c).BLCSM_csmT2(enum).PSTH_rawbin = psth_rawbin;
        data.units(c).BLCSM_csmT2(enum).PSTH_smoothbin = psth_smoothbin;
        data.units(c).BLCSM_csmT2(enum).time_bins = time_bins;

    end
end

%% BASELINE- T0
% Testing out new code to plot a different type of PSTH 
baseline_epoch = 1000;
event_epoch = 1000;

% Choose event to align to
event_name = ["csmT0"]; 
enames = fieldnames(Events_BLCSM);
evnt = find(ismember(enames, event_name));

% Loop through all units
for c = 1:length(data.units)
    spike_times = data.units(c).spike_times * 1000; % Convert to ms
    
    % Loop through each event type
    for enum = 1:length(evnt)
        event_times = Events_BLCSM.(enames{evnt(enum)}) * 1000;
        window = [baseline_epoch; event_epoch];

        % Compute referenced spikes and raster binning
        [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
        data.units(c).BLCSM_csmT0(enum).Referenced = Ref;
        data.units(c).BLCSM_csmT0(enum).EventTrial = EvTr;
        data.units(c).BLCSM_csmT0(enum).Raster = BinMast;

        % Compute PSTH using original method
        [psth_raw, psth_smooth] = Psth1(BinMast, 5); % Was 5 before
        data.units(c).BLCSM_csmT0(enum).PSTH_raw = psth_raw;
        data.units(c).BLCSM_csmT0(enum).PSTH_smooth = psth_smooth;

        % Compute PSTH using new computePSTH method (10ms bins + smoothing)
        bin_size = 10; % 10ms bin size
        smoothing_bins = 3; % Moving average over 3 bins

       [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);
        
        % Store results for the new method
        data.units(c).BLCSM_csmT0(enum).PSTH_rawbin = psth_rawbin;
        data.units(c).BLCSM_csmT0(enum).PSTH_smoothbin = psth_smoothbin;
        data.units(c).BLCSM_csmT0(enum).time_bins = time_bins;

    end
end

%% Get Events for Baseline Trials
% % ALL NO CSM REACHES


% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["pos_thresh"]; 
%event_name = ["end_point","pos_thresh","reach_start","max_vx"]; 
enames = fieldnames(Events_BLNoCSM);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_BLNoCSM.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    % 
    data.units(c).BLNoCSM_pos_thresh(enum).Referenced = Ref;
    data.units(c).BLNoCSM_pos_thresh(enum).EventTrial = EvTr;
    data.units(c).BLNoCSM_pos_thresh(enum).Raster = BinMast;

    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).BLNoCSM_pos_thresh(enum).PSTH_raw = psth_raw;
    data.units(c).BLNoCSM_pos_thresh(enum).PSTH_smooth = psth_smooth;  

    % Compute PSTH using new computePSTH method (10ms bins + smoothing)
    bin_size = 10; % 10ms bin size
    smoothing_bins = 3; % Moving average over 3 bins

    [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);

    % Store results for the new method
    data.units(c).BLNoCSM_pos_thresh(enum).PSTH_rawbin = psth_rawbin;
    data.units(c).BLNoCSM_pos_thresh(enum).PSTH_smoothbin = psth_smoothbin;
    data.units(c).BLNoCSM_pos_thresh(enum).time_bins = time_bins;

end
end


%% Get Events for Baseline Trials- max v 
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["max_vx"]; 
%event_name = ["end_point","pos_thresh","reach_start","max_vx"]; 
enames = fieldnames(Events_BLNoCSM);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_BLNoCSM.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    % 
    data.units(c).BLNoCSM_max_vx(enum).Referenced = Ref;
    data.units(c).BLNoCSM_max_vx(enum).EventTrial = EvTr;
    data.units(c).BLNoCSM_max_vx(enum).Raster = BinMast;

    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).BLNoCSM_max_vx(enum).PSTH_raw = psth_raw;
    data.units(c).BLNoCSM_max_vx(enum).PSTH_smooth = psth_smooth;  

    % Compute PSTH using new computePSTH method (10ms bins + smoothing)
    bin_size = 10; % 10ms bin size
    smoothing_bins = 3; % Moving average over 3 bins

    [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);

    % Store results for the new method
    data.units(c).BLNoCSM_max_vx(enum).PSTH_rawbin = psth_rawbin;
    data.units(c).BLNoCSM_max_vx(enum).PSTH_smoothbin = psth_smoothbin;
    data.units(c).BLNoCSM_max_vx(enum).time_bins = time_bins;

end
end


%% BASELINE- T2 
% % CSM REACHES
% Testing out new code to plot a different type of PSTH 
baseline_epoch = 1000;
event_epoch = 1000;

% Choose event to align to
event_name = ["csmT2"]; 
enames = fieldnames(Events_BLNoCSM);
evnt = find(ismember(enames, event_name));

% Loop through all units
for c = 1:length(data.units)
    spike_times = data.units(c).spike_times * 1000; % Convert to ms
    
    % Loop through each event type
    for enum = 1:length(evnt)
        event_times = Events_BLNoCSM.(enames{evnt(enum)}) * 1000;
        window = [baseline_epoch; event_epoch];

        % Compute referenced spikes and raster binning
        [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
        data.units(c).BLNoCSM_csmT2(enum).Referenced = Ref;
        data.units(c).BLNoCSM_csmT2(enum).EventTrial = EvTr;
        data.units(c).BLNoCSM_csmT2(enum).Raster = BinMast;

        % Compute PSTH using original method
        [psth_raw, psth_smooth] = Psth1(BinMast, 5); % Was 5 before
        data.units(c).BLNoCSM_csmT2(enum).PSTH_raw = psth_raw;
        data.units(c).BLNoCSM_csmT2(enum).PSTH_smooth = psth_smooth;

        % Compute PSTH using new computePSTH method (10ms bins + smoothing)
        bin_size = 10; % 10ms bin size
        smoothing_bins = 3; % Moving average over 3 bins

       [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);
        
        % Store results for the new method
        data.units(c).BLNoCSM_csmT2(enum).PSTH_rawbin = psth_rawbin;
        data.units(c).BLNoCSM_csmT2(enum).PSTH_smoothbin = psth_smoothbin;
        data.units(c).BLNoCSM_csmT2(enum).time_bins = time_bins;

    end
end

%% BASELINE- T0
% Testing out new code to plot a different type of PSTH 
baseline_epoch = 1000;
event_epoch = 1000;

% Choose event to align to
event_name = ["csmT0"]; 
enames = fieldnames(Events_BLNoCSM);
evnt = find(ismember(enames, event_name));

% Loop through all units
for c = 1:length(data.units)
    spike_times = data.units(c).spike_times * 1000; % Convert to ms
    
    % Loop through each event type
    for enum = 1:length(evnt)
        event_times = Events_BLNoCSM.(enames{evnt(enum)}) * 1000;
        window = [baseline_epoch; event_epoch];

        % Compute referenced spikes and raster binning
        [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
        data.units(c).BLNoCSM_csmT0(enum).Referenced = Ref;
        data.units(c).BLNoCSM_csmT0(enum).EventTrial = EvTr;
        data.units(c).BLNoCSM_csmT0(enum).Raster = BinMast;

        % Compute PSTH using original method
        [psth_raw, psth_smooth] = Psth1(BinMast, 5); % Was 5 before
        data.units(c).BLNoCSM_csmT0(enum).PSTH_raw = psth_raw;
        data.units(c).BLNoCSM_csmT0(enum).PSTH_smooth = psth_smooth;

        % Compute PSTH using new computePSTH method (10ms bins + smoothing)
        bin_size = 10; % 10ms bin size
        smoothing_bins = 3; % Moving average over 3 bins

       [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);
        
        % Store results for the new method
        data.units(c).BLNoCSM_csmT0(enum).PSTH_rawbin = psth_rawbin;
        data.units(c).BLNoCSM_csmT0(enum).PSTH_smoothbin = psth_smoothbin;
        data.units(c).BLNoCSM_csmT0(enum).time_bins = time_bins;

    end
end
%%
%% Get Events for Stim trials
% % ALL CSM ADAPT TRIALS 
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["pos_thresh"]; 
%event_name = ["end_point","pos_thresh","reach_start","max_vx"]; 
enames = fieldnames(Events_ADCSM);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_ADCSM.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
  

    data.units(c).ADCSM_pos_thresh(enum).Referenced = Ref;
    data.units(c).ADCSM_pos_thresh(enum).EventTrial = EvTr;
    data.units(c).ADCSM_pos_thresh(enum).Raster = BinMast;
    
    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).ADCSM_pos_thresh(enum).PSTH_raw = psth_raw;
    data.units(c).ADCSM_pos_thresh(enum).PSTH_smooth = psth_smooth;

    % Compute PSTH using new computePSTH method (10ms bins + smoothing)
    bin_size = 10; % 10ms bin size
    smoothing_bins = 3; % Moving average over 3 bins

    [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);

    % Store results for the new method
    data.units(c).ADCSM_pos_thresh(enum).PSTH_rawbin = psth_rawbin;
    data.units(c).ADCSM_pos_thresh(enum).PSTH_smoothbin = psth_smoothbin;
    data.units(c).ADCSM_pos_thresh(enum).time_bins = time_bins;
    
    
end
end

%% Get Events for ADAPT trials- max v 
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["max_vx"]; 
%event_name = ["end_point","pos_thresh","reach_start","max_vx"]; 
enames = fieldnames(Events_ADCSM);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_ADCSM.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    % 
    data.units(c).ADCSM_max_vx(enum).Referenced = Ref;
    data.units(c).ADCSM_max_vx(enum).EventTrial = EvTr;
    data.units(c).ADCSM_max_vx(enum).Raster = BinMast;

    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).ADCSM_max_vx(enum).PSTH_raw = psth_raw;
    data.units(c).ADCSM_max_vx(enum).PSTH_smooth = psth_smooth;  

    % Compute PSTH using new computePSTH method (10ms bins + smoothing)
    bin_size = 10; % 10ms bin size
    smoothing_bins = 3; % Moving average over 3 bins

    [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);

    % Store results for the new method
    data.units(c).ADCSM_max_vx(enum).PSTH_rawbin = psth_rawbin;
    data.units(c).ADCSM_max_vx(enum).PSTH_smoothbin = psth_smoothbin;
    data.units(c).ADCSM_max_vx(enum).time_bins = time_bins;

end
end


%% ADAPT- T2
% Testing out new code to plot a different type of PSTH 
baseline_epoch = 1000;
event_epoch = 1000;

% Choose event to align to
event_name = ["csmT2"]; 
enames = fieldnames(Events_ADCSM);
evnt = find(ismember(enames, event_name));

% Loop through all units
for c = 1:length(data.units)
    spike_times = data.units(c).spike_times * 1000; % Convert to ms
    
    % Loop through each event type
    for enum = 1:length(evnt)
        event_times = Events_ADCSM.(enames{evnt(enum)}) * 1000;
        window = [baseline_epoch; event_epoch];

        % Compute referenced spikes and raster binning
        [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
        data.units(c).ADCSM_csmT2(enum).Referenced = Ref;
        data.units(c).ADCSM_csmT2(enum).EventTrial = EvTr;
        data.units(c).ADCSM_csmT2(enum).Raster = BinMast;

        % Compute PSTH using original method
        [psth_raw, psth_smooth] = Psth1(BinMast, 5); % Was 5 before
        data.units(c).ADCSM_csmT2(enum).PSTH_raw = psth_raw;
        data.units(c).ADCSM_csmT2(enum).PSTH_smooth = psth_smooth;

        % Compute PSTH using new computePSTH method (10ms bins + smoothing)
        bin_size = 10; % 10ms bin size
        smoothing_bins = 3; % Moving average over 3 bins

       [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);
        
        % Store results for the new method
        data.units(c).ADCSM_csmT2(enum).PSTH_rawbin = psth_rawbin;
        data.units(c).ADCSM_csmT2(enum).PSTH_smoothbin = psth_smoothbin;
        data.units(c).ADCSM_csmT2(enum).time_bins = time_bins;

    end
end

%% BASELINE- T0
% Testing out new code to plot a different type of PSTH 
baseline_epoch = 1000;
event_epoch = 1000;

% Choose event to align to
event_name = ["csmT0"]; 
enames = fieldnames(Events_ADCSM);
evnt = find(ismember(enames, event_name));

% Loop through all units
for c = 1:length(data.units)
    spike_times = data.units(c).spike_times * 1000; % Convert to ms
    
    % Loop through each event type
    for enum = 1:length(evnt)
        event_times = Events_ADCSM.(enames{evnt(enum)}) * 1000;
        window = [baseline_epoch; event_epoch];

        % Compute referenced spikes and raster binning
        [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
        data.units(c).ADCSM_csmT0(enum).Referenced = Ref;
        data.units(c).ADCSM_csmT0(enum).EventTrial = EvTr;
        data.units(c).ADCSM_csmT0(enum).Raster = BinMast;

        % Compute PSTH using original method
        [psth_raw, psth_smooth] = Psth1(BinMast, 5); % Was 5 before
        data.units(c).ADCSM_csmT0(enum).PSTH_raw = psth_raw;
        data.units(c).ADCSM_csmT0(enum).PSTH_smooth = psth_smooth;

        % Compute PSTH using new computePSTH method (10ms bins + smoothing)
        bin_size = 10; % 10ms bin size
        smoothing_bins = 3; % Moving average over 3 bins

       [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);
        
        % Store results for the new method
        data.units(c).ADCSM_csmT0(enum).PSTH_rawbin = psth_rawbin;
        data.units(c).ADCSM_csmT0(enum).PSTH_smoothbin = psth_smoothbin;
        data.units(c).ADCSM_csmT0(enum).time_bins = time_bins;

    end
end
%% NOCSM ADAPT TRIALS

%% Get Events for Stim trials
% % ALL NO CSM ADAPT TRIALS 
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["pos_thresh"]; 
%event_name = ["end_point","pos_thresh","reach_start","max_vx"]; 
enames = fieldnames(Events_ADNoCSM);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_ADNoCSM.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
  

    data.units(c).ADNoCSM_pos_thresh(enum).Referenced = Ref;
    data.units(c).ADNoCSM_pos_thresh(enum).EventTrial = EvTr;
    data.units(c).ADNoCSM_pos_thresh(enum).Raster = BinMast;
    
    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).ADNoCSM_pos_thresh(enum).PSTH_raw = psth_raw;
    data.units(c).ADNoCSM_pos_thresh(enum).PSTH_smooth = psth_smooth;

    % Compute PSTH using new computePSTH method (10ms bins + smoothing)
    bin_size = 10; % 10ms bin size
    smoothing_bins = 3; % Moving average over 3 bins

    [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);

    % Store results for the new method
    data.units(c).ADNoCSM_pos_thresh(enum).PSTH_rawbin = psth_rawbin;
    data.units(c).ADNoCSM_pos_thresh(enum).PSTH_smoothbin = psth_smoothbin;
    data.units(c).ADNoCSM_pos_thresh(enum).time_bins = time_bins;
    
    
end
end

%% Get Events for ADAPT trials- max v 
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["max_vx"]; 
%event_name = ["end_point","pos_thresh","reach_start","max_vx"]; 
enames = fieldnames(Events_ADNoCSM);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_ADNoCSM.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    % 
    data.units(c).ADNoCSM_max_vx(enum).Referenced = Ref;
    data.units(c).ADNoCSM_max_vx(enum).EventTrial = EvTr;
    data.units(c).ADNoCSM_max_vx(enum).Raster = BinMast;

    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).ADNoCSM_max_vx(enum).PSTH_raw = psth_raw;
    data.units(c).ADNoCSM_max_vx(enum).PSTH_smooth = psth_smooth;  

    % Compute PSTH using new computePSTH method (10ms bins + smoothing)
    bin_size = 10; % 10ms bin size
    smoothing_bins = 3; % Moving average over 3 bins

    [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);

    % Store results for the new method
    data.units(c).ADNoCSM_max_vx(enum).PSTH_rawbin = psth_rawbin;
    data.units(c).ADNoCSM_max_vx(enum).PSTH_smoothbin = psth_smoothbin;
    data.units(c).ADNoCSM_max_vx(enum).time_bins = time_bins;

end
end


%% ADAPT- T2
% Testing out new code to plot a different type of PSTH 
baseline_epoch = 1000;
event_epoch = 1000;

% Choose event to align to
event_name = ["csmT2"]; 
enames = fieldnames(Events_ADNoCSM);
evnt = find(ismember(enames, event_name));

% Loop through all units
for c = 1:length(data.units)
    spike_times = data.units(c).spike_times * 1000; % Convert to ms
    
    % Loop through each event type
    for enum = 1:length(evnt)
        event_times = Events_ADNoCSM.(enames{evnt(enum)}) * 1000;
        window = [baseline_epoch; event_epoch];

        % Compute referenced spikes and raster binning
        [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
        data.units(c).ADNoCSM_csmT2(enum).Referenced = Ref;
        data.units(c).ADNoCSM_csmT2(enum).EventTrial = EvTr;
        data.units(c).ADNoCSM_csmT2(enum).Raster = BinMast;

        % Compute PSTH using original method
        [psth_raw, psth_smooth] = Psth1(BinMast, 5); % Was 5 before
        data.units(c).ADNoCSM_csmT2(enum).PSTH_raw = psth_raw;
        data.units(c).ADNoCSM_csmT2(enum).PSTH_smooth = psth_smooth;

        % Compute PSTH using new computePSTH method (10ms bins + smoothing)
        bin_size = 10; % 10ms bin size
        smoothing_bins = 3; % Moving average over 3 bins

       [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);
        
        % Store results for the new method
        data.units(c).ADNoCSM_csmT2(enum).PSTH_rawbin = psth_rawbin;
        data.units(c).ADNoCSM_csmT2(enum).PSTH_smoothbin = psth_smoothbin;
        data.units(c).ADNoCSM_csmT2(enum).time_bins = time_bins;

    end
end

%% BASELINE- T0
% Testing out new code to plot a different type of PSTH 
baseline_epoch = 1000;
event_epoch = 1000;

% Choose event to align to
event_name = ["csmT0"]; 
enames = fieldnames(Events_ADNoCSM);
evnt = find(ismember(enames, event_name));

% Loop through all units
for c = 1:length(data.units)
    spike_times = data.units(c).spike_times * 1000; % Convert to ms
    
    % Loop through each event type
    for enum = 1:length(evnt)
        event_times = Events_ADNoCSM.(enames{evnt(enum)}) * 1000;
        window = [baseline_epoch; event_epoch];

        % Compute referenced spikes and raster binning
        [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
        data.units(c).ADNoCSM_csmT0(enum).Referenced = Ref;
        data.units(c).ADNoCSM_csmT0(enum).EventTrial = EvTr;
        data.units(c).ADNoCSM_csmT0(enum).Raster = BinMast;

        % Compute PSTH using original method
        [psth_raw, psth_smooth] = Psth1(BinMast, 5); % Was 5 before
        data.units(c).ADNoCSM_csmT0(enum).PSTH_raw = psth_raw;
        data.units(c).ADNoCSM_csmT0(enum).PSTH_smooth = psth_smooth;

        % Compute PSTH using new computePSTH method (10ms bins + smoothing)
        bin_size = 10; % 10ms bin size
        smoothing_bins = 3; % Moving average over 3 bins

       [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);
        
        % Store results for the new method
        data.units(c).ADNoCSM_csmT0(enum).PSTH_rawbin = psth_rawbin;
        data.units(c).ADNoCSM_csmT0(enum).PSTH_smoothbin = psth_smoothbin;
        data.units(c).ADNoCSM_csmT0(enum).time_bins = time_bins;

    end
end
%% Get Events for Washout Trials- all trials 
% % ALL CSM REACHES
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["pos_thresh"];  
%event_name = ["end_point","pos_thresh","reach_start","max_vx"]; 
enames = fieldnames(Events_WOCSM);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_WOCSM.(enames{evnt(enum)})*1000;
    % event_times = event_times(1:5);
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    % 
    data.units(c).WOCSM_pos_thresh(enum).Referenced = Ref;
    data.units(c).WOCSM_pos_thresh(enum).EventTrial = EvTr;
    data.units(c).WOCSM_pos_thresh(enum).Raster = BinMast;

    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).WOCSM_pos_thresh(enum).PSTH_raw = psth_raw;
    data.units(c).WOCSM_pos_thresh(enum).PSTH_smooth = psth_smooth; 

    % Compute PSTH using new computePSTH method (10ms bins + smoothing)
    bin_size = 10; % 10ms bin size
    smoothing_bins = 3; % Moving average over 3 bins

    [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);

    % Store results for the new method
    data.units(c).WOCSM_pos_thresh(enum).PSTH_rawbin = psth_rawbin;
    data.units(c).WOCSM_pos_thresh(enum).PSTH_smoothbin = psth_smoothbin;
    data.units(c).WOCSM_pos_thresh(enum).time_bins = time_bins;
end
end

%% Get Events for WASHOUT trials- max v 
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["max_vx"]; 
%event_name = ["end_point","pos_thresh","reach_start","max_vx"]; 
enames = fieldnames(Events_WOCSM);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_WOCSM.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    % 
    data.units(c).WOCSM_max_vx(enum).Referenced = Ref;
    data.units(c).WOCSM_max_vx(enum).EventTrial = EvTr;
    data.units(c).WOCSM_max_vx(enum).Raster = BinMast;

    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).WOCSM_max_vx(enum).PSTH_raw = psth_raw;
    data.units(c).WOCSM_max_vx(enum).PSTH_smooth = psth_smooth;  

    % Compute PSTH using new computePSTH method (10ms bins + smoothing)
    bin_size = 10; % 10ms bin size
    smoothing_bins = 3; % Moving average over 3 bins

    [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);

    % Store results for the new method
    data.units(c).WOCSM_max_vx(enum).PSTH_rawbin = psth_rawbin;
    data.units(c).WOCSM_max_vx(enum).PSTH_smoothbin = psth_smoothbin;
    data.units(c).WOCSM_max_vx(enum).time_bins = time_bins;

end
end


%% WASHOUT- T2
% Testing out new code to plot a different type of PSTH 
baseline_epoch = 1000;
event_epoch = 1000;

% Choose event to align to
event_name = ["csmT2"]; 
enames = fieldnames(Events_WOCSM);
evnt = find(ismember(enames, event_name));

% Loop through all units
for c = 1:length(data.units)
    spike_times = data.units(c).spike_times * 1000; % Convert to ms
    
    % Loop through each event type
    for enum = 1:length(evnt)
        event_times = Events_WOCSM.(enames{evnt(enum)}) * 1000;
        window = [baseline_epoch; event_epoch];

        % Compute referenced spikes and raster binning
        [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
        data.units(c).WOCSM_csmT2(enum).Referenced = Ref;
        data.units(c).WOCSM_csmT2(enum).EventTrial = EvTr;
        data.units(c).WOCSM_csmT2(enum).Raster = BinMast;

        % Compute PSTH using original method
        [psth_raw, psth_smooth] = Psth1(BinMast, 5); % Was 5 before
        data.units(c).WOCSM_csmT2(enum).PSTH_raw = psth_raw;
        data.units(c).WOCSM_csmT2(enum).PSTH_smooth = psth_smooth;

        % Compute PSTH using new computePSTH method (10ms bins + smoothing)
        bin_size = 10; % 10ms bin size
        smoothing_bins = 3; % Moving average over 3 bins

       [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);
        
        % Store results for the new method
        data.units(c).WOCSM_csmT2(enum).PSTH_rawbin = psth_rawbin;
        data.units(c).WOCSM_csmT2(enum).PSTH_smoothbin = psth_smoothbin;
        data.units(c).WOCSM_csmT2(enum).time_bins = time_bins;

    end
end

%% WASHOUT - T0
% Testing out new code to plot a different type of PSTH 
baseline_epoch = 1000;
event_epoch = 1000;

% Choose event to align to
event_name = ["csmT0"]; 
enames = fieldnames(Events_WOCSM);
evnt = find(ismember(enames, event_name));

% Loop through all units
for c = 1:length(data.units)
    spike_times = data.units(c).spike_times * 1000; % Convert to ms
    
    % Loop through each event type
    for enum = 1:length(evnt)
        event_times = Events_WOCSM.(enames{evnt(enum)}) * 1000;
        window = [baseline_epoch; event_epoch];

        % Compute referenced spikes and raster binning
        [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
        data.units(c).WOCSM_csmT0(enum).Referenced = Ref;
        data.units(c).WOCSM_csmT0(enum).EventTrial = EvTr;
        data.units(c).WOCSM_csmT0(enum).Raster = BinMast;

        % Compute PSTH using original method
        [psth_raw, psth_smooth] = Psth1(BinMast, 5); % Was 5 before
        data.units(c).WOCSM_csmT0(enum).PSTH_raw = psth_raw;
        data.units(c).WOCSM_csmT0(enum).PSTH_smooth = psth_smooth;

        % Compute PSTH using new computePSTH method (10ms bins + smoothing)
        bin_size = 10; % 10ms bin size
        smoothing_bins = 3; % Moving average over 3 bins

       [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);
        
        % Store results for the new method
        data.units(c).WOCSM_csmT0(enum).PSTH_rawbin = psth_rawbin;
        data.units(c).WOCSM_csmT0(enum).PSTH_smoothbin = psth_smoothbin;
        data.units(c).WOCSM_csmT0(enum).time_bins = time_bins;

    end
end

%% WO NO CSM TRIALS
%% Get Events for Washout Trials- all NO CSM trials 
% % ALL CSM REACHES
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["pos_thresh"];  
%event_name = ["end_point","pos_thresh","reach_start","max_vx"]; 
enames = fieldnames(Events_WONoCSM);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_WONoCSM.(enames{evnt(enum)})*1000;
    % event_times = event_times(1:5);
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    % 
    data.units(c).WONoCSM_pos_thresh(enum).Referenced = Ref;
    data.units(c).WONoCSM_pos_thresh(enum).EventTrial = EvTr;
    data.units(c).WONoCSM_pos_thresh(enum).Raster = BinMast;

    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).WONoCSM_pos_thresh(enum).PSTH_raw = psth_raw;
    data.units(c).WONoCSM_pos_thresh(enum).PSTH_smooth = psth_smooth; 

    % Compute PSTH using new computePSTH method (10ms bins + smoothing)
    bin_size = 10; % 10ms bin size
    smoothing_bins = 3; % Moving average over 3 bins

    [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);

    % Store results for the new method
    data.units(c).WONoCSM_pos_thresh(enum).PSTH_rawbin = psth_rawbin;
    data.units(c).WONoCSM_pos_thresh(enum).PSTH_smoothbin = psth_smoothbin;
    data.units(c).WONoCSM_pos_thresh(enum).time_bins = time_bins;
end
end

%% Get Events for WASHOUT trials- max v 
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["max_vx"]; 
%event_name = ["end_point","pos_thresh","reach_start","max_vx"]; 
enames = fieldnames(Events_WONoCSM);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_WONoCSM.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    % 
    data.units(c).WONoCSM_max_vx(enum).Referenced = Ref;
    data.units(c).WONoCSM_max_vx(enum).EventTrial = EvTr;
    data.units(c).WONoCSM_max_vx(enum).Raster = BinMast;

    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).WONoCSM_max_vx(enum).PSTH_raw = psth_raw;
    data.units(c).WONoCSM_max_vx(enum).PSTH_smooth = psth_smooth;  

    % Compute PSTH using new computePSTH method (10ms bins + smoothing)
    bin_size = 10; % 10ms bin size
    smoothing_bins = 3; % Moving average over 3 bins

    [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);

    % Store results for the new method
    data.units(c).WONoCSM_max_vx(enum).PSTH_rawbin = psth_rawbin;
    data.units(c).WONoCSM_max_vx(enum).PSTH_smoothbin = psth_smoothbin;
    data.units(c).WONoCSM_max_vx(enum).time_bins = time_bins;

end
end


%% WASHOUT- T2
% Testing out new code to plot a different type of PSTH 
baseline_epoch = 1000;
event_epoch = 1000;

% Choose event to align to
event_name = ["csmT2"]; 
enames = fieldnames(Events_WONoCSM);
evnt = find(ismember(enames, event_name));

% Loop through all units
for c = 1:length(data.units)
    spike_times = data.units(c).spike_times * 1000; % Convert to ms
    
    % Loop through each event type
    for enum = 1:length(evnt)
        event_times = Events_WONoCSM.(enames{evnt(enum)}) * 1000;
        window = [baseline_epoch; event_epoch];

        % Compute referenced spikes and raster binning
        [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
        data.units(c).WONoCSM_csmT2(enum).Referenced = Ref;
        data.units(c).WONoCSM_csmT2(enum).EventTrial = EvTr;
        data.units(c).WONoCSM_csmT2(enum).Raster = BinMast;

        % Compute PSTH using original method
        [psth_raw, psth_smooth] = Psth1(BinMast, 5); % Was 5 before
        data.units(c).WONoCSM_csmT2(enum).PSTH_raw = psth_raw;
        data.units(c).WONoCSM_csmT2(enum).PSTH_smooth = psth_smooth;

        % Compute PSTH using new computePSTH method (10ms bins + smoothing)
        bin_size = 10; % 10ms bin size
        smoothing_bins = 3; % Moving average over 3 bins

       [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);
        
        % Store results for the new method
        data.units(c).WONoCSM_csmT2(enum).PSTH_rawbin = psth_rawbin;
        data.units(c).WONoCSM_csmT2(enum).PSTH_smoothbin = psth_smoothbin;
        data.units(c).WONoCSM_csmT2(enum).time_bins = time_bins;

    end
end

%% WASHOUT - T0
% Testing out new code to plot a different type of PSTH 
baseline_epoch = 1000;
event_epoch = 1000;

% Choose event to align to
event_name = ["csmT0"]; 
enames = fieldnames(Events_WONoCSM);
evnt = find(ismember(enames, event_name));

% Loop through all units
for c = 1:length(data.units)
    spike_times = data.units(c).spike_times * 1000; % Convert to ms
    
    % Loop through each event type
    for enum = 1:length(evnt)
        event_times = Events_WONoCSM.(enames{evnt(enum)}) * 1000;
        window = [baseline_epoch; event_epoch];

        % Compute referenced spikes and raster binning
        [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
        data.units(c).WONoCSM_csmT0(enum).Referenced = Ref;
        data.units(c).WONoCSM_csmT0(enum).EventTrial = EvTr;
        data.units(c).WONoCSM_csmT0(enum).Raster = BinMast;

        % Compute PSTH using original method
        [psth_raw, psth_smooth] = Psth1(BinMast, 5); % Was 5 before
        data.units(c).WONoCSM_csmT0(enum).PSTH_raw = psth_raw;
        data.units(c).WONoCSM_csmT0(enum).PSTH_smooth = psth_smooth;

        % Compute PSTH using new computePSTH method (10ms bins + smoothing)
        bin_size = 10; % 10ms bin size
        smoothing_bins = 3; % Moving average over 3 bins

       [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);
        
        % Store results for the new method
        data.units(c).WONoCSM_csmT0(enum).PSTH_rawbin = psth_rawbin;
        data.units(c).WONoCSM_csmT0(enum).PSTH_smoothbin = psth_smoothbin;
        data.units(c).WONoCSM_csmT0(enum).time_bins = time_bins;

    end
end
%% Get Events for Washout EARLY Trials- all trials - ALL CSM TRIALS 
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["pos_thresh"];  
%event_name = ["end_point","pos_thresh","reach_start","max_vx"]; 
enames = fieldnames(Events_WOECSM);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_WOECSM.(enames{evnt(enum)})*1000;
    % event_times = event_times(1:5);
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    % 
    data.units(c).WOECSM_pos_thresh(enum).Referenced = Ref;
    data.units(c).WOECSM_pos_thresh(enum).EventTrial = EvTr;
    data.units(c).WOECSM_pos_thresh(enum).Raster = BinMast;

    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).WOECSM_pos_thresh(enum).PSTH_raw = psth_raw;
    data.units(c).WOECSM_pos_thresh(enum).PSTH_smooth = psth_smooth; 

    % Compute PSTH using new computePSTH method (10ms bins + smoothing)
    bin_size = 10; % 10ms bin size
    smoothing_bins = 3; % Moving average over 3 bins

    [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);

    % Store results for the new method
    data.units(c).WOECSM_pos_thresh(enum).PSTH_rawbin = psth_rawbin;
    data.units(c).WOECSM_pos_thresh(enum).PSTH_smoothbin = psth_smoothbin;
    data.units(c).WOECSM_pos_thresh(enum).time_bins = time_bins;
end
end

%% Get Events for WASHOUT EARLY trials- max v 
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["max_vx"]; 
%event_name = ["end_point","pos_thresh","reach_start","max_vx"]; 
enames = fieldnames(Events_WOECSM);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_WOECSM.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    % 
    data.units(c).WOECSM_max_vx(enum).Referenced = Ref;
    data.units(c).WOECSM_max_vx(enum).EventTrial = EvTr;
    data.units(c).WOECSM_max_vx(enum).Raster = BinMast;

    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).WOECSM_max_vx(enum).PSTH_raw = psth_raw;
    data.units(c).WOECSM_max_vx(enum).PSTH_smooth = psth_smooth;  

    % Compute PSTH using new computePSTH method (10ms bins + smoothing)
    bin_size = 10; % 10ms bin size
    smoothing_bins = 3; % Moving average over 3 bins

    [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);

    % Store results for the new method
    data.units(c).WOECSM_max_vx(enum).PSTH_rawbin = psth_rawbin;
    data.units(c).WOECSM_max_vx(enum).PSTH_smoothbin = psth_smoothbin;
    data.units(c).WOECSM_max_vx(enum).time_bins = time_bins;

end
end


%% WASHOUT EARLY- T2
% Testing out new code to plot a different type of PSTH 
baseline_epoch = 1000;
event_epoch = 1000;

% Choose event to align to
event_name = ["csmT2"]; 
enames = fieldnames(Events_WOECSM);
evnt = find(ismember(enames, event_name));

% Loop through all units
for c = 1:length(data.units)
    spike_times = data.units(c).spike_times * 1000; % Convert to ms
    
    % Loop through each event type
    for enum = 1:length(evnt)
        event_times = Events_WOECSM.(enames{evnt(enum)}) * 1000;
        window = [baseline_epoch; event_epoch];

        % Compute referenced spikes and raster binning
        [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
        data.units(c).WOECSM_csmT2(enum).Referenced = Ref;
        data.units(c).WOECSM_csmT2(enum).EventTrial = EvTr;
        data.units(c).WOECSM_csmT2(enum).Raster = BinMast;

        % Compute PSTH using original method
        [psth_raw, psth_smooth] = Psth1(BinMast, 5); % Was 5 before
        data.units(c).WOECSM_csmT2(enum).PSTH_raw = psth_raw;
        data.units(c).WOECSM_csmT2(enum).PSTH_smooth = psth_smooth;

        % Compute PSTH using new computePSTH method (10ms bins + smoothing)
        bin_size = 10; % 10ms bin size
        smoothing_bins = 3; % Moving average over 3 bins

       [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);
        
        % Store results for the new method
        data.units(c).WOECSM_csmT2(enum).PSTH_rawbin = psth_rawbin;
        data.units(c).WOECSM_csmT2(enum).PSTH_smoothbin = psth_smoothbin;
        data.units(c).WOECSM_csmT2(enum).time_bins = time_bins;

    end
end

%% WASHOUT EARLY - T0
% Testing out new code to plot a different type of PSTH 
baseline_epoch = 1000;
event_epoch = 1000;

% Choose event to align to
event_name = ["csmT0"]; 
enames = fieldnames(Events_WOECSM);
evnt = find(ismember(enames, event_name));

% Loop through all units
for c = 1:length(data.units)
    spike_times = data.units(c).spike_times * 1000; % Convert to ms
    
    % Loop through each event type
    for enum = 1:length(evnt)
        event_times = Events_WOECSM.(enames{evnt(enum)}) * 1000;
        window = [baseline_epoch; event_epoch];

        % Compute referenced spikes and raster binning
        [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
        data.units(c).WOECSM_csmT0(enum).Referenced = Ref;
        data.units(c).WOECSM_csmT0(enum).EventTrial = EvTr;
        data.units(c).WOECSM_csmT0(enum).Raster = BinMast;

        % Compute PSTH using original method
        [psth_raw, psth_smooth] = Psth1(BinMast, 5); % Was 5 before
        data.units(c).WOECSM_csmT0(enum).PSTH_raw = psth_raw;
        data.units(c).WOECSM_csmT0(enum).PSTH_smooth = psth_smooth;

        % Compute PSTH using new computePSTH method (10ms bins + smoothing)
        bin_size = 10; % 10ms bin size
        smoothing_bins = 3; % Moving average over 3 bins

       [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);
        
        % Store results for the new method
        data.units(c).WOECSM_csmT0(enum).PSTH_rawbin = psth_rawbin;
        data.units(c).WOECSM_csmT0(enum).PSTH_smoothbin = psth_smoothbin;
        data.units(c).WOECSM_csmT0(enum).time_bins = time_bins;

    end
end
%% 
% all no CSM TRIALS 

%% Get Events for Washout EARLY Trials- all trials - ALL CSM TRIALS 
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["pos_thresh"];  
%event_name = ["end_point","pos_thresh","reach_start","max_vx"]; 
enames = fieldnames(Events_WOENoCSM);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_WOENoCSM.(enames{evnt(enum)})*1000;
    % event_times = event_times(1:5);
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    % 
    data.units(c).WOENoCSM_pos_thresh(enum).Referenced = Ref;
    data.units(c).WOENoCSM_pos_thresh(enum).EventTrial = EvTr;
    data.units(c).WOENoCSM_pos_thresh(enum).Raster = BinMast;

    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).WOENoCSM_pos_thresh(enum).PSTH_raw = psth_raw;
    data.units(c).WOENoCSM_pos_thresh(enum).PSTH_smooth = psth_smooth; 

    % Compute PSTH using new computePSTH method (10ms bins + smoothing)
    bin_size = 10; % 10ms bin size
    smoothing_bins = 3; % Moving average over 3 bins

    [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);

    % Store results for the new method
    data.units(c).WOENoCSM_pos_thresh(enum).PSTH_rawbin = psth_rawbin;
    data.units(c).WOENoCSM_pos_thresh(enum).PSTH_smoothbin = psth_smoothbin;
    data.units(c).WOENoCSM_pos_thresh(enum).time_bins = time_bins;
end
end

%% Get Events for WASHOUT EARLY trials- max v 
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["max_vx"]; 
%event_name = ["end_point","pos_thresh","reach_start","max_vx"]; 
enames = fieldnames(Events_WOENoCSM);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_WOENoCSM.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    % 
    data.units(c).WOENoCSM_max_vx(enum).Referenced = Ref;
    data.units(c).WOENoCSM_max_vx(enum).EventTrial = EvTr;
    data.units(c).WOENoCSM_max_vx(enum).Raster = BinMast;

    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).WOENoCSM_max_vx(enum).PSTH_raw = psth_raw;
    data.units(c).WOENoCSM_max_vx(enum).PSTH_smooth = psth_smooth;  

    % Compute PSTH using new computePSTH method (10ms bins + smoothing)
    bin_size = 10; % 10ms bin size
    smoothing_bins = 3; % Moving average over 3 bins

    [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);

    % Store results for the new method
    data.units(c).WOENoCSM_max_vx(enum).PSTH_rawbin = psth_rawbin;
    data.units(c).WOENoCSM_max_vx(enum).PSTH_smoothbin = psth_smoothbin;
    data.units(c).WOENoCSM_max_vx(enum).time_bins = time_bins;

end
end


%% WASHOUT EARLY- T2
% Testing out new code to plot a different type of PSTH 
baseline_epoch = 1000;
event_epoch = 1000;

% Choose event to align to
event_name = ["csmT2"]; 
enames = fieldnames(Events_WOENoCSM);
evnt = find(ismember(enames, event_name));

% Loop through all units
for c = 1:length(data.units)
    spike_times = data.units(c).spike_times * 1000; % Convert to ms
    
    % Loop through each event type
    for enum = 1:length(evnt)
        event_times = Events_WOENoCSM.(enames{evnt(enum)}) * 1000;
        window = [baseline_epoch; event_epoch];

        % Compute referenced spikes and raster binning
        [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
        data.units(c).WOENoCSM_csmT2(enum).Referenced = Ref;
        data.units(c).WOENoCSM_csmT2(enum).EventTrial = EvTr;
        data.units(c).WOENoCSM_csmT2(enum).Raster = BinMast;

        % Compute PSTH using original method
        [psth_raw, psth_smooth] = Psth1(BinMast, 5); % Was 5 before
        data.units(c).WOENoCSM_csmT2(enum).PSTH_raw = psth_raw;
        data.units(c).WOENoCSM_csmT2(enum).PSTH_smooth = psth_smooth;

        % Compute PSTH using new computePSTH method (10ms bins + smoothing)
        bin_size = 10; % 10ms bin size
        smoothing_bins = 3; % Moving average over 3 bins

       [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);
        
        % Store results for the new method
        data.units(c).WOENoCSM_csmT2(enum).PSTH_rawbin = psth_rawbin;
        data.units(c).WOENoCSM_csmT2(enum).PSTH_smoothbin = psth_smoothbin;
        data.units(c).WOENoCSM_csmT2(enum).time_bins = time_bins;

    end
end

%% WASHOUT EARLY - T0
% Testing out new code to plot a different type of PSTH 
baseline_epoch = 1000;
event_epoch = 1000;

% Choose event to align to
event_name = ["csmT0"]; 
enames = fieldnames(Events_WOENoCSM);
evnt = find(ismember(enames, event_name));

% Loop through all units
for c = 1:length(data.units)
    spike_times = data.units(c).spike_times * 1000; % Convert to ms
    
    % Loop through each event type
    for enum = 1:length(evnt)
        event_times = Events_WOENoCSM.(enames{evnt(enum)}) * 1000;
        window = [baseline_epoch; event_epoch];

        % Compute referenced spikes and raster binning
        [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
        data.units(c).WOENoCSM_csmT0(enum).Referenced = Ref;
        data.units(c).WOENoCSM_csmT0(enum).EventTrial = EvTr;
        data.units(c).WOENoCSM_csmT0(enum).Raster = BinMast;

        % Compute PSTH using original method
        [psth_raw, psth_smooth] = Psth1(BinMast, 5); % Was 5 before
        data.units(c).WOENoCSM_csmT0(enum).PSTH_raw = psth_raw;
        data.units(c).WOENoCSM_csmT0(enum).PSTH_smooth = psth_smooth;

        % Compute PSTH using new computePSTH method (10ms bins + smoothing)
        bin_size = 10; % 10ms bin size
        smoothing_bins = 3; % Moving average over 3 bins

       [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);
        
        % Store results for the new method
        data.units(c).WOENoCSM_csmT0(enum).PSTH_rawbin = psth_rawbin;
        data.units(c).WOENoCSM_csmT0(enum).PSTH_smoothbin = psth_smoothbin;
        data.units(c).WOENoCSM_csmT0(enum).time_bins = time_bins;

    end
end
%% Get Events for all CONTROL TRIALS
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["pos_thresh"];  
%event_name = ["end_point","pos_thresh","reach_start","max_vx"]; 
enames = fieldnames(Events_Ctrl);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_Ctrl.(enames{evnt(enum)})*1000;
    % event_times = event_times(1:5);
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    % 
    data.units(c).Ctrl_pos_thresh(enum).Referenced = Ref;
    data.units(c).Ctrl_pos_thresh(enum).EventTrial = EvTr;
    data.units(c).Ctrl_pos_thresh(enum).Raster = BinMast;

    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).Ctrl_pos_thresh(enum).PSTH_raw = psth_raw;
    data.units(c).Ctrl_pos_thresh(enum).PSTH_smooth = psth_smooth; 

    % Compute PSTH using new computePSTH method (10ms bins + smoothing)
    bin_size = 10; % 10ms bin size
    smoothing_bins = 3; % Moving average over 3 bins

    [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);

    % Store results for the new method
    data.units(c).Ctrl_pos_thresh(enum).PSTH_rawbin = psth_rawbin;
    data.units(c).Ctrl_pos_thresh(enum).PSTH_smoothbin = psth_smoothbin;
    data.units(c).Ctrl_pos_thresh(enum).time_bins = time_bins;
end
end

%% Get Events for all CONTROL TRIALS
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["end_point"];  
%event_name = ["end_point","pos_thresh","reach_start","max_vx"]; 
enames = fieldnames(Events_Ctrl);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_Ctrl.(enames{evnt(enum)})*1000;
    % event_times = event_times(1:5);
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    % 
    data.units(c).Ctrl_endpt(enum).Referenced = Ref;
    data.units(c).Ctrl_endpt(enum).EventTrial = EvTr;
    data.units(c).Ctrl_endpt(enum).Raster = BinMast;

    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).Ctrl_endpt(enum).PSTH_raw = psth_raw;
    data.units(c).Ctrl_endpt(enum).PSTH_smooth = psth_smooth; 

    % Compute PSTH using new computePSTH method (10ms bins + smoothing)
    bin_size = 10; % 10ms bin size
    smoothing_bins = 3; % Moving average over 3 bins

    [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);

    % Store results for the new method
    data.units(c).Ctrl_endpt(enum).PSTH_rawbin = psth_rawbin;
    data.units(c).Ctrl_endpt(enum).PSTH_smoothbin = psth_smoothbin;
    data.units(c).Ctrl_endpt(enum).time_bins = time_bins;
end
end

%% Get Events for all CONTROL TRIALS- maxvx 
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["max_vx"];  
%event_name = ["end_point","pos_thresh","reach_start","max_vx"]; 
enames = fieldnames(Events_Ctrl);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000; % converting to ms 
for enum = 1:length(evnt)
    event_times = Events_Ctrl.(enames{evnt(enum)})*1000;
    % event_times = event_times(1:5);
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    % 
    data.units(c).Ctrl_max_vx(enum).Referenced = Ref;
    data.units(c).Ctrl_max_vx(enum).EventTrial = EvTr;
    data.units(c).Ctrl_max_vx(enum).Raster = BinMast;

    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).Ctrl_max_vx(enum).PSTH_raw = psth_raw;
    data.units(c).Ctrl_max_vx(enum).PSTH_smooth = psth_smooth; 

    % Compute PSTH using new computePSTH method (10ms bins + smoothing)
    bin_size = 10; % 10ms bin size
    smoothing_bins = 3; % Moving average over 3 bins

    [psth_rawbin, psth_smoothbin, time_bins] = computePSTH(BinMast, bin_size, smoothing_bins);

    % Store results for the new method
    data.units(c).Ctrl_max_vx(enum).PSTH_rawbin = psth_rawbin;
    data.units(c).Ctrl_max_vx(enum).PSTH_smoothbin = psth_smoothbin;
    data.units(c).Ctrl_max_vx(enum).time_bins = time_bins;
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
    
    [psth_raw,psth_smooth] = Psth1opto(BinMast,5); % was 5 before 
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
    %close;
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

%% Get modulated intA cells 
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

%% Get Modulated IntA Neurons that are also optotagged 
% Find overlap between optotagged cells and reach-modulated IntA cells
overlap_indices = intersect(Salt_index, filtered_modulated_indices);

overlap_cluster_IDs = [data.units(overlap_indices).cluster_ID];
disp('Cluster IDs of overlapping cells:');
disp(overlap_cluster_IDs);

% Display the results
disp('Indices of cells that are both optotagged and modulated:');
disp(overlap_indices);
length(overlap_indices)

%% other way to potentially look at this 
% Find and store indices of cells that are both optotagged and modulated
optotagged_modulated_indices = intersect(Salt_index, filtered_modulated_indices);

% Display results
disp('Indices of optotagged and modulated cells:');
disp(optotagged_modulated_indices);
length(optotagged_modulated_indices)

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
%%
% % Plot Kinematics, Rasters and PSTHs 
% 
% epo = 500;
% 
% 
% folder_path="C:\Users\PersonLab\OneDrive - The University of Colorado Denver\Documents\NeuropixelsResults\StimAnimals\20220520_SAM010\POSTHRESH_BLOCKALL";
% 
% for c = 1:length(data.units)
% %for c = filtered_modulated_indices
% 
% 
% 
% 
%     figure('Position', [100, 100, 800, 900]); hold on
% 
%     % % 
%     subplot(5,1,1); hold on;
%     dim=8;
%     StimvelPlot(BL,dim,BL_color,0);
%     StimvelPlot(S,dim,Adaptcolor,0);
%     StimvelPlot(WO,dim,WO_color,0);
%     StimvelPlot(WOE,dim,WOE_color,0);
%     ylim('auto')
%     xlim([-.5,.5])
% 
% 
%     subplot(5,1,2); hold on
%     h = title(['cluster' ' ' num2str(data.units(c).cluster_ID) ' ' 'depth' ' ' num2str(data.units(c).depth)]) ;
%     set(h,'interpreter','none')
% 
%     % Plots Raster from binary matrix
%     tick_length = 1; % length of tick marks for spikes
%     spike_tick_width = 0.5; % width of lines marking spikes
%     event_tick_width = 2; % width of lines marking secondary events
% 
%     %baseline NoRR\
% 
%     rastr = data.units(c).BL_pos_thresh.Raster(:,(1000-epo):(2000-(1000-epo)));
%     %rastr = data.units(c).BL_csm.Raster(:,(1000-epo):(2000-(1000-epo)));
%     [y_co,x_co] = find(rastr);
%     x_co = x_co - epo;
%     line([x_co';x_co'],[y_co'-(tick_length/2);y_co'+(tick_length/2)],'Color',BL_color, 'LineWidth', spike_tick_width);
%     line([zeros(1,size(rastr,1));zeros(1,size(rastr,1))],[(1:1:size(rastr,1))-.5;(1:1:size(rastr,1))+.5],'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);
%     axis ij;
%     ylim([0,size(rastr,1)])
%     xlim([-epo,epo])
%     ylabel('Baseline ')
% 
%     subplot(5,1,3); hold on
%     % adapt all
%     rastr_stim = data.units(c).AD_pos_thresh.Raster(:,(1000-epo):(2000-(1000-epo)));
% 
%     [y_co_stim,x_co_stim] = find(rastr_stim);
%     x_co_stim = x_co_stim - epo;
%     line([x_co_stim';x_co_stim'],[y_co_stim'-(tick_length/2);y_co_stim'+(tick_length/2)],'Color',Adaptcolor, 'LineWidth', spike_tick_width);
%     line([zeros(1,size(rastr_stim,1));zeros(1,size(rastr_stim,1))],[(1:1:size(rastr_stim,1))-.5;(1:1:size(rastr_stim,1))+.5],'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);
%     axis ij;
%     ylim([0,size(rastr_stim,1)])
%     xlim([-epo,epo])
%     ylabel('Stim Trial')
% 
%     subplot(5,1,4); hold on
%     % washout 
%     rastr_wo = data.units(c).WOE_pos_thresh.Raster(:,(1000-epo):(2000-(1000-epo)));
%     %rastr_wo = data.units(c).WO_NoRR_pos_thresh.Raster(1:5,(1000-epo):(2000-(1000-epo)));
%     [y_co_wo,x_co_wo] = find(rastr_wo);
%     x_co_wo = x_co_wo - epo;
%     line([x_co_wo';x_co_wo'],[y_co_wo'-(tick_length/2);y_co_wo'+(tick_length/2)],'Color',WOE_color, 'LineWidth', spike_tick_width);
%     line([zeros(1,size(rastr_wo,1));zeros(1,size(rastr_wo,1))],[(1:1:size(rastr_wo,1))-.5;(1:1:size(rastr_wo,1))+.5],'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);
%     axis ij;
%     ylim([0,size(rastr_wo,1)])
%     xlim([-epo,epo])
%     ylabel('Washout Early ')
% 
% 
%     subplot(5,1,5); hold on
%     title('peth')
% 
%     psth_smooth = data.units(c).BL_pos_thresh.PSTH_smoothbin;
%     psth_raw = data.units(c).BL_pos_thresh.PSTH_rawbin;
% 
% 
%     psth_smooth_stim = data.units(c).AD_pos_thresh.PSTH_smoothbin;
%     psth_raw_stim = data.units(c).AD_pos_thresh.PSTH_rawbin;
% 
%     psth_smooth_wo = data.units(c).WOE_pos_thresh.PSTH_smoothbin;
%     psth_raw_wo = data.units(c).WOE_pos_thresh.PSTH_rawbin;
% 
%     % Generate time vector (assuming same binning and length)
%     num_bins = size(psth_smooth, 2);
%     time_vector = linspace(-1000, 1000, num_bins);  % Update if your window or binning is different
% 
%     plot(time_vector,psth_smooth,'Color', BL_color, 'LineWidth', 2);
%     hold on
%     plot(time_vector,psth_smooth_stim,'Color', Adaptcolor, 'LineWidth', 2);
%     hold on
%     plot(time_vector,psth_smooth_wo,'Color', WOE_color, 'LineWidth', 2);
% 
% 
%     %ylim([0 1]);
%     yline(0,'--')
%     xlim([-500 500])
%     xlabel('Time from positional threshold');
%     ylabel('Firing rate (Hz)');
%     % 
%     filename=sprintf('pos_thresh_Blockr%d.svg',c);
% 
%     full_path=fullfile(folder_path,filename);
%     saveas(gcf,full_path);
% 
%     close;
% % % % 
% % if mod(c,10) == 0
% %     drawnow
% %     input("press enter when done looking at figures")
% %     close all
% % end
% 
%  end

% %% Now do heatmaps 
% 
% % calc time of max firing rate 
% maxfr_inx = NaN(length(filtered_modulated_indices),1);
% 
% % Preallocate for speed using the number of filtered units
% num_filtered = numel(filtered_modulated_indices);
% maxfr_inx = zeros(1, num_filtered);
% 
% for i = 1:num_filtered
%     c = filtered_modulated_indices(i);
%     psth = data.units(c).BL_pos_thresh.PSTH_smoothbin;
%     [~, maxfr_inx(i)] = max(psth);
% end
% 
% % for i = 1:num_filtered
% %     c = filtered_modulated_indices(i);
% %     psth = data.units(c).BL_pos_thresh.PSTH_smooth;
% %     [~, maxfr_inx(i)] = max(psth);
% % end
% 
% % Now sort the max firing rate indices
% [maxfr_t, inx] = sort(maxfr_inx);
% 
% 
% %organizing into master matrix
% psth_masterBL = NaN(length(data.units(filtered_modulated_indices)),length(data.units(filtered_modulated_indices(inx(1))).BL_pos_thresh.PSTH_smoothbin));
% psth_masterAD = NaN(length(data.units(filtered_modulated_indices)),length(data.units(filtered_modulated_indices(inx(1))).AD_pos_thresh.PSTH_smoothbin));
% psth_masterWO = NaN(length(data.units(filtered_modulated_indices)),length(data.units(filtered_modulated_indices(inx(1))).WO_pos_thresh.PSTH_smoothbin));
% psth_masterWOE = NaN(length(data.units(filtered_modulated_indices)),length(data.units(filtered_modulated_indices(inx(1))).WOE_pos_thresh.PSTH_smoothbin));
% % 
% % psth_masterBL = NaN(length(data.units(filtered_modulated_indices)),length(data.units(filtered_modulated_indices(inx(1))).BL_pos_thresh.PSTH_smooth));
% % psth_masterAD = NaN(length(data.units(filtered_modulated_indices)),length(data.units(filtered_modulated_indices(inx(1))).AD_pos_thresh.PSTH_smooth));
% % psth_masterWO = NaN(length(data.units(filtered_modulated_indices)),length(data.units(filtered_modulated_indices(inx(1))).WO_pos_thresh.PSTH_smooth));
% % psth_masterWOE = NaN(length(data.units(filtered_modulated_indices)),length(data.units(filtered_modulated_indices(inx(1))).WOE_pos_thresh.PSTH_smooth));
% 
% for c =1:length(filtered_modulated_indices)
% 
%     psth_masterBL(c,:) = data.units(inx(c)).BL_pos_thresh.PSTH_smoothbin;
%     psth_masterAD(c,:) = data.units(inx(c)).AD_pos_thresh.PSTH_smoothbin;
%     psth_masterWO(c,:) = data.units(inx(c)).WO_pos_thresh.PSTH_smoothbin;
%     psth_masterWOE(c,:) = data.units(inx(c)).WOE_pos_thresh.PSTH_smoothbin;
% 
%     % psth_masterBL(c,:) = data.units(inx(c)).BL_pos_thresh.PSTH_smooth;
%     % psth_masterAD(c,:) = data.units(inx(c)).AD_pos_thresh.PSTH_smooth;
%     % psth_masterWO(c,:) = data.units(inx(c)).WO_pos_thresh.PSTH_smooth;
%     % psth_masterWOE(c,:) = data.units(inx(c)).WOE_pos_thresh.PSTH_smooth;
% 
% end
% 
% 
% % Now sort by peak firing rate 
% sorted_psth_BL = psth_masterBL(inx, :);
% sorted_psth_AD = psth_masterAD(inx, :);
% sorted_psth_WO = psth_masterWO(inx, :);
% sorted_psth_WOE = psth_masterWOE(inx, :);
%% Now plot heatmaps 
% 
% % Generate time vector
% num_bins = size(sorted_psth_BL, 2);
% time_vector = linspace(-1000, 1000, num_bins);  % assuming 10 ms bins over 200 bins
% 
% % Set tick marks for x-axis every 250 ms
% tick_times = -1000:250:1000;  
% [~, tick_indices] = min(abs(time_vector' - tick_times), [], 1);  % find closest indices
% 
% % Get global min and max for color scaling
% all_data = [sorted_psth_BL(:); sorted_psth_AD(:); sorted_psth_WOE(:); sorted_psth_WO(:)];
% vmin = min(all_data);
% vmax = max(all_data);
% 
% % Plot the sorted NoRRT2 PSTH heatmap
% figure('Position', [100, 100, 800, 900]); hold on
% 
% subplot(4,1,1); hold on;
% imagesc(sorted_psth_BL);
% colormap(slanCM('GnBu'));
% caxis([vmin vmax]);
% colorbar;
% axis tight;
% axis ij;
% set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times);
% xlabel('Time (ms)');
% ylabel('Units (sorted by peak firing time)');
% title('Baseline');
% 
% subplot(4,1,2); hold on;
% imagesc(sorted_psth_AD);
% colormap(slanCM('GnBu'));
% caxis([vmin vmax]);
% colorbar;
% axis tight;
% axis ij;
% set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times);
% xlabel('Time (ms)');
% ylabel('Units (sorted by peak firing time)');
% title('Adapt-Stim Block');
% 
% subplot(4,1,3); hold on;
% imagesc(sorted_psth_WOE);
% colormap(slanCM('GnBu'));
% caxis([vmin vmax]);
% colorbar;
% axis tight;
% axis ij;
% set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times);
% xlabel('Time (ms)');
% ylabel('Units (sorted by peak firing time)');
% title('Washout Early');
% 
% subplot(4,1,4); hold on;
% imagesc(sorted_psth_WO);
% colormap(slanCM('GnBu'));
% caxis([vmin vmax]);
% colorbar;
% axis tight;
% axis ij;
% set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times);
% xlabel('Time (ms)');
% ylabel('Units (sorted by peak firing time)');
% title('Washout All');
% %% Now plot heatmaps 
% % Generate time vector
% num_bins = size(sorted_psth_BL, 2);
% time_vector = linspace(-1000, 1000, num_bins);  % assuming 10 ms bins over 200 bins
% 
% % Set tick marks for x-axis every 250 ms
% tick_times = -1000:250:1000;  
% [~, tick_indices] = min(abs(time_vector' - tick_times), [], 1);  % find closest indices
% 
% % Get global min and max for color scaling
% all_data = [sorted_psth_BL(:); sorted_psth_AD(:); sorted_psth_WOE(:); sorted_psth_WO(:)];
% vmin = min(all_data);
% vmax = max(all_data);
% 
% % Create figure with tiled layout
% figure('Position', [100, 100, 800, 900]);
% tiledlayout(4,1, 'Padding', 'compact', 'TileSpacing', 'compact');
% 
% % Plot each condition
% nexttile; hold on;
% imagesc(sorted_psth_BL);
% colormap(slanCM('GnBu'));
% caxis([vmin vmax]);
% colorbar;
% axis tight;
% axis ij;
% set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times, 'FontSize', 10);
% xlabel('Time (ms)');
% ylabel('Units (sorted by peak firing time)');
% title('Baseline');
% 
% nexttile; hold on;
% imagesc(sorted_psth_AD);
% colormap(slanCM('GnBu'));
% caxis([vmin vmax]);
% colorbar;
% axis tight;
% axis ij;
% set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times, 'FontSize', 10);
% xlabel('Time (ms)');
% ylabel('Units (sorted by peak firing time)');
% title('Adapt-Stim Block');
% 
% nexttile; hold on;
% imagesc(sorted_psth_WOE);
% colormap(slanCM('GnBu'));
% caxis([vmin vmax]);
% colorbar;
% axis tight;
% axis ij;
% set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times, 'FontSize', 10);
% xlabel('Time (ms)');
% ylabel('Units (sorted by peak firing time)');
% title('Washout Early');
% 
% nexttile; hold on;
% imagesc(sorted_psth_WO);
% colormap(slanCM('GnBu'));
% caxis([vmin vmax]);
% colorbar;
% axis tight;
% axis ij;
% set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times, 'FontSize', 10);
% xlabel('Time (ms)');
% ylabel('Units (sorted by peak firing time)');
% title('Washout All');
% 
% 
% %% Take Difference 
% psth_diff = sorted_psth_WOE - sorted_psth_BL;
% %psth_diff = sorted_psth_AD - sorted_psth_BL;
% 
% % Compute symmetric color limits centered on 0
% max_val = max(abs(psth_diff(:)));
% 
% figure;
% imagesc(psth_diff);
% colormap(flipud(slanCM('RdBu')));  % Red=positive, Blue=negative
% clim([-max_val, max_val]);         % Set zero to be white in the center
% colorbar;
% 
% set(gca, 'XTick', tick_indices, 'XTickLabel', tick_times, 'FontSize', 10);
% xlabel('Time (ms)');
% ylabel('Units (sorted by peak firing time)');
% %title('WOE - BL Difference');
% title('AD - BL Difference');
% 
% %% Avg psth across units 
% 
% % Compute average PSTHs (mean across units)
% avg_psth_BL = mean(sorted_psth_BL, 1);  % 1 = average across rows (units)
% avg_psth_AD = mean(sorted_psth_AD, 1);
% avg_psth_WO = mean(sorted_psth_WO, 1);
% avg_psth_WOE = mean(sorted_psth_WOE, 1);
% 
% % Optional: Compute standard error of the mean (SEM)
% sem_psth_BL = std(sorted_psth_BL, 0, 1) / sqrt(size(sorted_psth_BL, 1));
% sem_psth_AD = std(sorted_psth_AD, 0, 1) / sqrt(size(sorted_psth_AD, 1));
% sem_psth_WO = std(sorted_psth_WO, 0, 1) / sqrt(size(sorted_psth_WO, 1));
% sem_psth_WOE = std(sorted_psth_WOE, 0, 1) / sqrt(size(sorted_psth_WOE, 1));
% 
% 
% % Generate time vector (assuming same binning and length)
% num_bins = size(sorted_psth_BL, 2);
% time_vector = linspace(-1000, 1000, num_bins);  % Update if your window or binning is different
% 
% % Plot average PSTHs
% figure;
% hold on;
% plot(time_vector, avg_psth_BL,'Color', BL_color, 'LineWidth', 2);
% plot(time_vector, avg_psth_AD, 'Color',Adaptcolor, 'LineWidth', 2);
% plot(time_vector, avg_psth_WO,'Color', WO_color, 'LineWidth', 2);
% plot(time_vector, avg_psth_WOE,'Color', WOE_color, 'LineWidth', 2);
% 
% % Optional: Add shaded error bars for SEM
% % fill([time_vector, fliplr(time_vector)], ...
% %      [avg_psth_BL + sem_psth_BL, fliplr(avg_psth_BL - sem_psth_BL)], ...
% %      BL_color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% % fill([time_vector, fliplr(time_vector)], ...
% %      [avg_psth_AD + sem_psth_AD, fliplr(avg_psth_AD - sem_psth_AD)], ...
% %      Adaptcolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% % fill([time_vector, fliplr(time_vector)], ...
% %      [avg_psth_WO + sem_psth_WO, fliplr(avg_psth_WO - sem_psth_WO)], ...
% %      WO_color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% % fill([time_vector, fliplr(time_vector)], ...
% %      [avg_psth_WOE + sem_psth_WOE, fliplr(avg_psth_WOE - sem_psth_WOE)], ...
% %      WOE_color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% %ylim([0 1]);
% xlim([-500 500])
% xlabel('Time (ms)');
% ylabel('Firing rate (Hz)');
% 
% title('Average PSTH across Units across Blocks');
% 
% %% Diff
% diffcolor=rgb('LightSeaGreen');
% % Compute average PSTHs (mean across units)
% avg_psth_diff = mean(psth_diff, 1);  % 1 = average across rows (units)
% 
% % Optional: Compute standard error of the mean (SEM)
% sem_psth_diff = std(psth_diff, 0, 1) / sqrt(size(psth_diff, 1));
% 
% % Generate time vector (assuming same binning and length)
% num_bins = size(psth_diff, 2);
% time_vector = linspace(-1000, 1000, num_bins);  % Update if your window or binning is different
% 
% % Plot average PSTHs
% figure;
% hold on;
% plot(time_vector, avg_psth_diff,'Color',diffcolor, 'LineWidth', 2);
% 
% % Optional: Add shaded error bars for SEM
% fill([time_vector, fliplr(time_vector)], ...
%      [avg_psth_diff + sem_psth_diff, fliplr(avg_psth_diff - sem_psth_diff)], ...
%      diffcolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% 
% %ylim([0 1]);
% yline(0,'--')
% %xlim([-250 250])
% xlim([-500 500])
% xlabel('Time (ms)');
% ylabel('Firing rate (Hz)');
% title('Av. Diff in FR CSM-NoCSM');

%%

cd("D:\Neural Data Alignment\HFBlockStim_CSMvsNoCSM_induced")
dataname = '220515_SAM012_Block_CSMvsNoCSM';
save(dataname, 'data','Reaches','ExcitRR_BL','ExcitRR_AD','ExcitRR_WO','ExcitRR_WOE', 'blk_track','blk_len','dta','suc','Salt_index','opto_index','Q','Q_frs','Q_reach','Kins','FRBins','FRCorrStats_edit','filtered_modulated_indices','optotagged_modulated_indices','pos_corr_cells','neg_corr_cells','-v7.3');