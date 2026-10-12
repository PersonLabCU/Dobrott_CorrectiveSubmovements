%% ChR2 structure with my data and Matt's data

% Need to first look at Medial reaches - separate out those for CID001 and
% CID002
% first need to rename structure 
names=fieldnames(ChR2);
an=names(7:8);

for i = 1:length(an)
    Reaches=ChR2.(an{i}).Reaches;
    [Reaches,bad]=rfilts_standard(Reaches,0.16,5,0.8);
      

  % Separate out into Med and Lat reaches 
    [Med,Lat]=TargetFilterML(Reaches);

    [Stim,Control]=StimFilter(Med);

    [S_C, F_C] = SuccessFilter(Control);
    [S_S, F_S] = SuccessFilter(Stim);
    
    ChR2.(an{i}).successrate.control = length(fieldnames(S_C)) / length(fieldnames(Control));
    ChR2.(an{i}).successrate.stim = length(fieldnames(S_S)) / length(fieldnames(Stim));

    ChR2.(an{i}).Reaches=Med;
    ChR2.(an{i}).LateralR=Lat;

end
%% Success Rate 
an=fieldnames(ChR2);
for a = 1:length(an)
    SR_control(a,1) = ChR2.(an{a}).successrate.control;
    SR_stim(a,1) = ChR2.(an{a}).successrate.stim;
    
end
%% 
an=fieldnames(ChR2);

for a = 1:length(an)
 
   Reaches= ChR2.(an{a}).Reaches;

    
    %setup variables 
    ag = 6;
    filt_state = 1;
    pct = .8;
    leng = 5;
    thresh = .16;
     dirdir=1;
   
    rnum = fieldnames(Reaches);
    [Stim, Control]=StimFilter(Reaches);
    
  % Get kinematic features for ctrl reaches 
   [ctrlRfilt,bad]=rfilts_standard(Control,0.16,5,0.8);
    CtrlR=accel_toreaches2(ctrlRfilt); ctrlKin=BasicKinematics(CtrlR);
    Ctrlacc=mean(ctrlKin.aavg_accel(:,2));
    Ctrldec=mean(ctrlKin.aavg_decel_70(:,2));
    C = AlignThresh(ctrlRfilt,thresh,dirdir,500);
    Ca=stimvel_accel(C);
    Cend=endpoint(ctrlRfilt);
    Cendmean=mean(Cend);
    
    C_CSM=AlignThreshCSM(ctrlRfilt,thresh,dirdir,400);
    C_CSM=stimvel_accel(C_CSM);

    % Get kinematic featrures for stim reaches 
    [stimRfilt,bad]=rfilts_standard(Stim,0.16,5,0.8);
    StimR=accel_toreaches2(stimRfilt);stimKin=BasicKinematics(StimR);
    Stimacc=mean(stimKin.aavg_accel(:,2));
    Stimdec=mean(stimKin.aavg_decel_70(:,2));
    S = AlignStim(stimRfilt,500);
    Sa=stimvel_accel(S);
    Send=endpoint(Stim);
    Sendmean=mean(Send);

    S_CSM=AlignCSM(stimRfilt,400);
    S_CSM=stimvel_accel(S_CSM);


    % Velocity percentiles 
     % now sort reaches based off of velocity percentiles 
         % Sort based off of change in distance during stim 
    % if using 250, stim start=26, stim end =31
% if using 500, star=51, end =56
    Sreaches=S.stimvel;
    fields=fieldnames(StimR);
        
    StimEndPos=[];
    StimStartPos=[];
    for i=1:length(fields)
        endpos=Sreaches(56,2,i);
        startpos=Sreaches(51,2,i); % stim start
        StimEndPos=[StimEndPos;endpos];
        StimStartPos=[StimStartPos;startpos];
    end
    %  Change in position
    StimChange=(StimEndPos-StimStartPos)*10;
    VmaxStim=StimChange;
    
   % VmaxStim=stimKin.vmax(:,2);
    
    p25=prctile(VmaxStim,25);
    p50=prctile(VmaxStim,50);
    p75=prctile(VmaxStim,75);
    ChR2.(an{a}).VmaxPrct.p25=p25;
    ChR2.(an{a}).VmaxPrct.p50=p50;
    ChR2.(an{a}).VmaxPrct.p75=p75;
    %find indices for elements
    below_p25_indices = find(VmaxStim <= p25);
    between_p25_p50_indices = find(VmaxStim > p25 & VmaxStim <= p50);
    between_p50_p75_indices = find(VmaxStim > p50 & VmaxStim <= p75);
    above_p75_indices = find(VmaxStim > p75);
    

    % Get field names from the structure
    Stimfields = fieldnames(StimR);
    stim25=Stimfields(below_p25_indices);
    stim50=Stimfields( between_p25_p50_indices);
    stim75=Stimfields( between_p50_p75_indices);
    stim100=Stimfields(above_p75_indices);

    % Find the indices of desired field names
    [~, below_p25_indices] = ismember(stim25, Stimfields);
    [~,between_p25_p50_indices ] = ismember(stim50, Stimfields);
    [~,between_p50_p75_indices ] = ismember(stim75, Stimfields);
     [~, above_p75_indices] = ismember(stim100, Stimfields);
    % Initialize a structure to hold extracted data
    Stim25 = struct(); Stim50=struct();Stim75=struct();Stim100=struct();

    % Extract data based on desired field names
    for j = 1:length(below_p25_indices)
        index = below_p25_indices(j);

        if index > 0
            fieldName = stim25{j};
            Stim25.(fieldName) = Extract(StimR.(fieldName));
        end
    end
      
     for p = 1:length(between_p25_p50_indices)
        index50 = between_p25_p50_indices(p);

        if index50 > 0
            fieldName50 = stim50{p};
            Stim50.(fieldName50) = Extract(StimR.(fieldName50));
        end
      end
    
       for w = 1:length(between_p50_p75_indices)
        index75 = between_p50_p75_indices(w);

        if index75 > 0
            fieldName75 = stim75{w};
            Stim75.(fieldName75) = Extract(StimR.(fieldName75));
        end
       end
       
     for c = 1:length(above_p75_indices)
        index100 = above_p75_indices(c);

        if index100 > 0
            fieldName100 = stim100{c};
            Stim100.(fieldName100) = Extract(StimR.(fieldName100));
        end
     end
   
    % Now save each structure 
    ChR2.(an{a}).Stim25Reaches=Stim25;
    ChR2.(an{a}).Stim50Reaches=Stim50;
    ChR2.(an{a}).Stim75Reaches=Stim75;
    ChR2.(an{a}).Stim100Reaches=Stim100;
   
    %make stimvel for each subset of reaches
    S25 = AlignStim(Stim25,500);S25=stimvel_accel(S25);
    S50 = AlignStim(Stim50,500);S50=stimvel_accel(S50);
    S75 = AlignStim(Stim75,500);S75=stimvel_accel(S75);
    S100 = AlignStim(Stim100,500);S100=stimvel_accel(S100);
    %save stimvel
    ChR2.(an{a}).S25 = S25;
    ChR2.(an{a}).S50 = S50;
    ChR2.(an{a}).S75 = S75;
    ChR2.(an{a}).S100 = S100;
  
      % Now add to Chr2 structure 
    ChR2.(an{a}).Control=CtrlR;
    ChR2.(an{a}).CtrlKin=ctrlKin; 
    ChR2.(an{a}).CtrlAcc=Ctrlacc;
    ChR2.(an{a}).CtrlDec=Ctrldec;
    ChR2.(an{a}).C=Ca;
    ChR2.(an{a}).C_CSM=C_CSM;
    ChR2.(an{a}).Cend=Cendmean;
    
    % add stim valu\es 
    ChR2.(an{a}).Stim=StimR;
    ChR2.(an{a}).StimKin=stimKin; 
    ChR2.(an{a}).StimAcc=Stimacc;
    ChR2.(an{a}).StimDec=Stimdec;
    ChR2.(an{a}).S=Sa;
    ChR2.(an{a}).S_CSM=S_CSM;
    ChR2.(an{a}).Send=Sendmean;
    

end   
%% Use ReReach Filter - Control Reaches 
ExcitRR_Ctrl=[]; % Savethis structure to use for graphing rereach timing in sep. script
for a=1:length(an)

    Creaches=ChR2.(an{a}).Control;

 
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
    
    ExcitRR_Ctrl.(an{a}).successrate.ReReach = length(fieldnames(S_C)) / length(fieldnames(Creaches));
    ExcitRR_Ctrl.(an{a}).successrate.NoReReach = length(fieldnames(S_S)) / length(fieldnames(Creaches));
    ExcitRR_Ctrl.(an{a}).CtrlReaches=Creaches;
    ExcitRR_Ctrl.(an{a}).C=C;
    ExcitRR_Ctrl.(an{a}).ReReaches=ReReaches;
    ExcitRR_Ctrl.(an{a}).RRKin=RRKin;
    ExcitRR_Ctrl.(an{a}).NoReReach=NoReach;
    ExcitRR_Ctrl.(an{a}).NRKin=NRKin;
    ExcitRR_Ctrl.(an{a}).CN=CN;
    ExcitRR_Ctrl.(an{a}).StimRR=StimRR;
    ExcitRR_Ctrl.(an{a}).MaxVel=MaxVel;
    % get % of reaches with a rereach
    cnum=length(fields1);
    rrnum=length(fields2);
    percentrr=rrnum/cnum;
    ExcitRR_Ctrl.(an{a}).Cnum=cnum;
    ExcitRR_Ctrl.(an{a}).RRnum=rrnum;
    ExcitRR_Ctrl.(an{a}).PercentRR=percentrr;
  
end 

%% Use ReReach Filter- STIM reaches 
ExcitRR=[]; % Savethis structure to use for graphing rereach timing in sep. script
for a=1:length(an)
    Sreaches=ChR2.(an{a}).Stim;
  

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
    
    ExcitRR.(an{a}).successrate.ReReach = length(fieldnames(S_C)) / length(fieldnames(Sreaches));
    ExcitRR.(an{a}).successrate.NoReReach = length(fieldnames(S_S)) / length(fieldnames(Sreaches));
    ExcitRR.(an{a}).ReReaches=ReReaches;
    ExcitRR.(an{a}).RRKin=RRKin;
    ExcitRR.(an{a}).S=S;
    ExcitRR.(an{a}).NoReReach=NoReach;
    ExcitRR.(an{a}).NRKin=NRKin;
    ExcitRR.(an{a}).SN=SN;
    ExcitRR.(an{a}).StimRR=StimRR;
    ExcitRR.(an{a}).MaxVel=MaxVel;
    snum=length(fields1);
    rrnum=length(fields2);
    percentrr=rrnum/snum;
    ExcitRR.(an{a}).Snum=snum;
    ExcitRR.(an{a}).RRnum=rrnum;
    ExcitRR.(an{a}).PercentRR=percentrr;
end 



%% Add accel column to SN data - not all animals have NoReReach - this is to add accel and to figure out what animals have none
NoReReachStimvel=[];
indices_to_delete = [];

for a=1:length(an)
    NRR=ExcitRR.(an{a}).SN;
    if isempty(NRR.stimvel)
        NoReReachStimvel.(an{a}).SN=[];
        % Add the index of the empty NRR.stimvel to the list
        indices_to_delete(end+1) = a;
    else
    NRR=stimvel_accel(NRR);
    NoReReachStimvel.(an{a}).SN=NRR;
    NoReReachStimvel.(an{a}).C= ExcitRR.(an{a}).C;
    end
end 

% Delete structures where NRR.stimvel is empty
for i = length(indices_to_delete):-1:1
    NoReReachStimvel = rmfield(NoReReachStimvel, an{indices_to_delete(i)});
end

%%
stimveldiffX=[];stimveldiffY=[];
for a=1:length(an)
   %  CavgX=ChR2.(an{a}).C.Average(26,5);
   % SavgX=ChR2.(an{a}).S.Average(26,5);

   CavgX=ChR2.(an{a}).C.Average(56,5); % from 500
   SavgX=ChR2.(an{a}).S.Average(56,5);
  

    CavgY=ChR2.(an{a}).C.Average(56,6);%from 500
    SavgY=ChR2.(an{a}).S.Average(56,6);

    DiffX=SavgX-CavgX;
    DiffY=SavgY-CavgY;
    stimveldiffX(a,1)=DiffX;
    stimveldiffY(a,1)=DiffY;
end 

AllAnimalsDiff=[stimveldiffX,stimveldiffY];
%% Plot stimvel plots 
%get average combined stimvel plot for 'out significant' animals
an=fieldnames(ChR2);

C_big = struct;
S_big = struct; 


for a = 1:length(an)
   
    C_big.stimvel(:,:,a) = ChR2.(an{a}).C.Average;
    S_big.stimvel(:,:,a) = ChR2.(an{a}).S.Average;
    
    
end


C_big.Average = nanmean(C_big.stimvel,3);
C_big.Std = nanstd(C_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(C_big.Average),1); 
for nnum = 1:length(C_big.Average)
    samp = squeeze(C_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
C_big.N = N;


S_big.Average = nanmean(S_big.stimvel,3);
S_big.Std = nanstd(S_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(C_big.Average),1); 
for nnum = 1:length(C_big.Average)
    samp = squeeze(S_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
S_big.N = N;

%% Stimvel for CSM

an=fieldnames(ChR2);

C_CSM_big = struct;
S_CSM_big = struct; 


for a = 1:length(an)
   
    C_CSM_big.stimvel(:,:,a) = ChR2.(an{a}).C_CSM.Average;
    S_CSM_big.stimvel(:,:,a) = ChR2.(an{a}).S_CSM.Average;
    
    
end


C_CSM_big.Average = nanmean(C_CSM_big.stimvel,3);
C_CSM_big.Std = nanstd(C_CSM_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(C_CSM_big.Average),1); 
for nnum = 1:length(C_CSM_big.Average)
    samp = squeeze(C_CSM_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
C_CSM_big.N = N;


S_CSM_big.Average = nanmean(S_CSM_big.stimvel,3);
S_CSM_big.Std = nanstd(S_CSM_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(C_CSM_big.Average),1); 
for nnum = 1:length(C_CSM_big.Average)
    samp = squeeze(S_CSM_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
S_CSM_big.N = N;


%% Get Stimvel plots for rereaches
SRR_big = struct;
SRR_CSM_big=struct;
ar=fieldnames(ExcitRR);
for a = 1:length(ar) 
    SRR_big.stimvel(:,:,a) = ExcitRR.(ar{a}).S.Average;
    SRR_CSM_big.stimvel(:,:,a) = ExcitRR.(ar{a}).S_CSM.Average
end

%ReReach
SRR_big.Average = nanmean(SRR_big.stimvel,3);
SRR_big.Std = nanstd(SRR_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(C_big.Average),1); 
for nnum = 1:length(C_big.Average)
    samp = squeeze(SRR_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
SRR_big.N = N;

%ReReach CSM algined
SRR_CSM_big.Average = nanmean(SRR_CSM_big.stimvel,3);
SRR_CSM_big.Std = nanstd(SRR_CSM_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(C_big.Average),1); 
for nnum = 1:length(C_big.Average)
    samp = squeeze(SRR_CSM_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
SRR_CSM_big.N = N;

%% No ReReach Structure
NRR_big=struct;

ar=fieldnames(NoReReachStimvel);
for a = 1:length(ar) 
    NRR_big.stimvel(:,:,a) = NoReReachStimvel.(ar{a}).SN.Average;
end
%NoReReach
NRR_big.Average = nanmean(NRR_big.stimvel,3);
NRR_big.Std = nanstd(NRR_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(NRR_big.Average),1); 
for nnum = 1:length(NRR_big.Average)
    samp = squeeze(NRR_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
NRR_big.N = N;
%% Plot stimvel
x1=0 ;x2=0.05;
y_upper=0.25; y_lower=0;
%Add acceleration column to stimvel
C_big=stimvel_accel(C_big);
S_big=stimvel_accel(S_big);


stimcolor= rgb('Teal');
linecolor=rgb('DarkGray');

dim =2;

StimvelPlot(C_big,dim,'k',1)
StimvelPlot(S_big,dim,stimcolor,0)

%StimvelPlot(P_big,dim,pluscolor,0)
hold on
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
xlim([-0.2 0.2])
%ylim([-0.5 1.5])
%ylim([-30 20])
ylim([0 0.25])
pbaspect([1 2 1])
%ylim([0 0.3])
yline(0,'--')
xlabel('Time from Stim (s)')
ylabel('Outward Velocity (cm/s)')
%ylabel('Outward acceleration')



% figure; hold on
% StimvelPlot2(C_big,dim,[.5 .5 .5],0)
% shadedErrorBars(C_big.Average(:,1),C_big.Average(:,dim),C_big.Std(:,dim)./sqrt(C_big.N),[.5 .5 .5])
% StimvelPlot2(S_big,dim,[.2 .6 .9],0)
% shadedErrorBars(S_big.Average(:,1),S_big.Average(:,dim),S_big.Std(:,dim)./sqrt(S_big.N),[.2 .6 .9])
% StimvelPlot2(P_big,dim,[0 .5 0],0)
% shadedErrorBars(P_big.Average(:,1),P_big.Average(:,dim),P_big.Std(:,dim)./sqrt(P_big.N),[.2 .6 .9])
%% stimvel aligned to CSM
C_CSM_big=stimvel_accel(C_CSM_big);
S_CSM_big=stimvel_accel(S_CSM_big);


stimcolor= rgb('Teal');
linecolor=rgb('DarkGray');

dim =2;

StimvelPlot(C_CSM_big,dim,'k',1)
StimvelPlot(SRR_CSM_big,dim,stimcolor,0)

xlim([-0.2 0.2])
%ylim([-1 1.5])
%ylim([-30 20])
ylim('auto')
pbaspect([1 2 1])
%ylim([0 0.3])
yline(0,'--')
xlabel('Time from CSM (s)')
ylabel('Outward Velocity (cm/s)')
%% PLot stimvel of rereaches 
norereach=rgb('MediumAquamarine');
linecolor=rgb('DarkGray');

dim = 2;
x1=0 ;x2=0.05;
y_upper=0.25; y_lower=0;
StimvelPlot(C_big,dim,'k',1)
StimvelPlot(SRR_big,dim,stimcolor,0)
StimvelPlot(NRR_big,dim,norereach,0)
hold on 
% Fill the area between the vertical lines
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');

yline(0,'--k')
xlim([-0.2 0.2])
ylim([0 0.25])
%ylim([-15 10])
%ylim('auto')
%yline(0,'--')
xlabel('Time from Stim (s)')
ylabel('Outward Velocity (cm/s)')
pbaspect([ 1 2 1])

%% PLot poistion 
x1=0 ;x2=0.05;
y_upper=3; y_lower=0;
dim = 2;

figure; hold on
StimvelPlot(C_big,dim,'k',1)
StimvelPlot(S_big,dim,stimcolor,0)
hold on
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
xlim([-0.2 0.2])
ylim([0 0.25])
pbaspect([ 1 2 1]) 
xlabel('Time from Stim (s)')
ylabel('Outward Position(cm)')
%% Look at reaches based off of velocity 
S25_big = struct; S50_big=struct;S75_big=struct; S100_big=struct;

for a = 1:length(an)
   
    S25_big.stimvel(:,:,a) = ChR2.(an{a}).S25.Average;
    S50_big.stimvel(:,:,a) = ChR2.(an{a}).S50.Average;
    S75_big.stimvel(:,:,a) = ChR2.(an{a}).S75.Average;
    S100_big.stimvel(:,:,a) = ChR2.(an{a}).S100.Average;
    
end

S25_big.Average = nanmean(S25_big.stimvel,3);
S25_big.Std = nanstd(S25_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(C_big.Average),1); 
for nnum = 1:length(C_big.Average)
    samp = squeeze(S25_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
S25_big.N = N;

S50_big.Average = nanmean(S50_big.stimvel,3);
S50_big.Std = nanstd(S50_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(C_big.Average),1); 
for nnum = 1:length(C_big.Average)
    samp = squeeze(S50_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
S50_big.N = N;

S75_big.Average = nanmean(S75_big.stimvel,3);
S75_big.Std = nanstd(S75_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(C_big.Average),1); 
for nnum = 1:length(C_big.Average)
    samp = squeeze(S75_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
S75_big.N = N;

S100_big.Average = nanmean(S100_big.stimvel,3);
S100_big.Std = nanstd(S100_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(C_big.Average),1); 
for nnum = 1:length(C_big.Average)
    samp = squeeze(S100_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
S100_big.N = N;

%% Now plot bbased off of velocity
stimcolor25= rgb('LightBlue');
stimcolor50=rgb('Dodgerblue');
stimcolor75=rgb('MediumBlue');
stimcolor100=rgb('MidnightBlue');
linecolor=rgb('DarkGray');


C_big=stimvel_accel(C_big);
S25_big=stimvel_accel(S25_big);
S50_big=stimvel_accel(S50_big);
S75_big=stimvel_accel(S75_big);
S100_big=stimvel_accel(S100_big);

dim = 8;
x1=0 ;x2=0.05;
y_upper=40; y_lower=-40;
StimvelPlot(C_big,dim,'k',1)
StimvelPlot(S25_big,dim,stimcolor25,0)
StimvelPlot(S50_big,dim,stimcolor50,0)
StimvelPlot(S75_big,dim,stimcolor75,0)
StimvelPlot(S100_big,dim,stimcolor100,0)
hold on 
% Fill the area between the vertical lines
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');

yline(0,'--k')
xlim([-0.2 0.2])
%ylim([-0.5 1.5])
%ylim([0 .25])
ylim('auto')
xlabel('Time from Stim (s)')
%ylabel('Outward Velocity (cm/s)')
%ylabel('Outward Position (cm)')
ylabel('Outward Acceleration (cm/s^2)')
%pbaspect([ 1 2 1])

figure; hold on
StimvelPlot2(S25_big,dim,stimcolor25,0)
shadedErrorBars(S25_big.Average(:,1),S25_big.Average(:,dim),S25_big.Std(:,dim)./sqrt(S25_big.N),stimcolor25)
StimvelPlot2(S50_big,dim,stimcolor50,0)
shadedErrorBars(S50_big.Average(:,1),S50_big.Average(:,dim),S50_big.Std(:,dim)./sqrt(S50_big.N),stimcolor50)
StimvelPlot2(S75_big,dim,stimcolor75,0)
shadedErrorBars(S75_big.Average(:,1),S75_big.Average(:,dim),S75_big.Std(:,dim)./sqrt(S75_big.N),stimcolor75)
StimvelPlot2(S100_big,dim,stimcolor100,0)
shadedErrorBars(S100_big.Average(:,1),S100_big.Average(:,dim),S100_big.Std(:,dim)./sqrt(S100_big.N),stimcolor100)
xlim([-0.2 0.2])
%ylim([-0.5 1.5])
%ylim([0 .25])
ylim('auto')
yline(0,'--')
xlabel('Time from Stim (s)')
%ylabel('Outward Velocity (cm/s)')
%ylabel('Outward Position (cm)')
ylabel('Outward Acceleration (cm/s^2)')
%pbaspect([ 1 2 1

%% Plot position only 
figure;

for j=1:11
    Cpos=C_big.stimvel(:,2,j);
    Spos=S_big.stimvel(:,2,j);
    x=[0.05 0.25];
    % plot(Cpos,(1:numel(Cpos)),'k')
    % hold on
    % plot(Spos,(1:numel(Cpos)),'b')
    % hold on
    plot(Cpos,'k')
    hold on
    plot(Spos,'b')
    hold on
    
end
%xlim([5 40])
ylim([0 0.3])
pbaspect([ 1 2 1])


%% Get mangnitude of effect for RR 
an=fieldnames(ExcitRR);
RRstimveldiffX=[];RRstimveldiffY=[];
for a=1:length(an)
    CavgX=ExcitRR.(an{a}).C.Average(26,5);
    SavgX=ExcitRR.(an{a}).S.Average(26,5);
    
    CavgY=ExcitRR.(an{a}).C.Average(26,6);
    SavgY=ExcitRR.(an{a}).S.Average(26,6);

    DiffX=SavgX-CavgX;
    DiffY=SavgY-CavgY;
    RRstimveldiffX(a,1)=DiffX;
    RRstimveldiffY(a,1)=DiffY;
end 

AllAnimalsDiffRR=[RRstimveldiffX,RRstimveldiffY];
%% get magnitude for NRR
nr=fieldnames(NoReReachStimvel);
RRstimveldiffX=[];RRstimveldiffY=[];
for a=1:length(nr)
    CavgX=NoReReachStimvel.(nr{a}).C.Average(26,5);
    SavgX=NoReReachStimvel.(nr{a}).SN.Average(26,5);

    CavgY=NoReReachStimvel.(nr{a}).C.Average(26,6);
    SavgY=NoReReachStimvel.(nr{a}).SN.Average(26,6);

    DiffX=SavgX-CavgX;
    DiffY=SavgY-CavgY;
    RRstimveldiffX(a,1)=DiffX;
    RRstimveldiffY(a,1)=DiffY;
end 

AllAnimalsDiffNRR=[RRstimveldiffX,RRstimveldiffY];
%% Get endpoint info 
EndptCtrl=[];
EndptStim=[];
for a=1:length(an)
    endptCtrl=ChR2.(an{a}).endpoint.control(:,1);
    endptStim=ChR2.(an{a}).endpoint.stim(:,1);
    Cmean=mean(endptCtrl);
    Smean=mean(endptStim);

    EndptCtrl(a)=Cmean;
    EndptStim(a)=Smean;
end 

%% Look at deceleration and peak velocity
figure; hold on
for a=1:length(ar)
    vmaxRR=ExcitRR.(ar{a}).RRKin.vmax(:,2);vmaxRR=vmaxRR*10;
    vmaxNR=ExcitRR.(ar{a}).NRKin.vmax(:,2);vmaxNR=vmaxNR*10;

    % decelRR=ExcitRR.(ar{a}).RRKin.aavg_decel_stim50(:,2);decelRR=decelRR*10;
    % decelNR=ExcitRR.(ar{a}).NRKin.aavg_decel_stim50(:,2);decelNR=decelNR*10;

    decelRR=ExcitRR.(ar{a}).RRKin.aavg_accel(:,2);decelRR=decelRR*10;
    decelNR=ExcitRR.(ar{a}).NRKin.aavg_accel(:,2);decelNR=decelNR*10;
    vmaxRRm=mean(vmaxRR);vmaxNRm=mean(vmaxNR);
    decelRRm=mean(decelRR); decelNRm=mean(decelNR);

    %linear regression
    RR_reg=fitlm(vmaxRR,decelRR)
    NR_reg=fitlm(vmaxNR,decelNR)
    plot(RR_reg)
    hold on
    plot(NR_reg,'Color','k')
    hold on
    scatter(vmaxRR,decelRR,'MarkerFaceColo',[0 0.5 0.5])
    hold on
    scatter(vmaxNR,decelNR,'MarkerFaceColor',[0.3984    0.8008    0.6641])
    scatter(vmaxRRm,decelRRm,50,'MarkerFaceColor',[0 0.5 0.5],'MarkerEdgeColor','k')
    hold on
    scatter(vmaxNRm,decelNRm,50,'MarkerFaceColor',[0.3984    0.8008    0.6641],'MarkerEdgeColor','k')
end 
 
%%

%corr per animal

vd_corr = []; vd_corrNR=[];
vd_corr_p = []; vd_corr_p=[];

vd_corr_rand = [];vd_corr_randNR = [];
vd_corr_p_rand = [];vd_corr_p_randNR = [];

fits = [];fitsNR=[];



RRcolor=rgb('Teal');
NRRcolor=rgb('MediumAquamarine');
%for plotting all linear correlations
figure; hold on
ar=fieldnames(ExcitRR)
for a = 1:length(ar)
    
    vmaxRR=ExcitRR.(ar{a}).RRKin.vmax(:,2); rpeakv=vmaxRR*10;
    vmaxNR=ExcitRR.(ar{a}).NRKin.vmax(:,2);rpeakvNR=vmaxNR*10;

    decelRR=ExcitRR.(ar{a}).RRKin.aavg_accel(:,2);rdecel=decelRR*10;
    decelNR=ExcitRR.(ar{a}).NRKin.aavg_accel(:,2);rdecelNR=decelNR*10;

    %decelRR=ExcitRR.(ar{a}).RRKin.aavg_decel_stim50(:,2);rdecel=decelRR*10;
    %decelNR=ExcitRR.(ar{a}).NRKin.aavg_decel_stim50(:,2);rdecelNR=decelNR*10;
  
    %compute correlation coefficient and pvalues
    [vd_corr(a), vd_corr_p(a)] = corr(rpeakv,rdecel);
    
    %for no rereach
     [vd_corrNR(a), vd_corr_pNR(a)] = corr(rpeakvNR,rdecelNR);
    
    %shuffle control for ReREach
    corr_temp = NaN(1,1000); p_temp = NaN(1,1000);
    for j = 1:1000
        rand_inx = randsample(length(rpeakv),length(rpeakv));
        [corr_temp(j), p_temp(j)] = corr(rpeakv(rand_inx),rdecel);
    end
    vd_corr_rand(a) = mean(corr_temp);
    vd_corr_p_rand(a) = mean(p_temp);

    % Shuffle control for no rereach
      corr_tempNR = NaN(1,1000); p_tempNR = NaN(1,1000);
    for h = 1:1000
        rand_inxNR = randsample(length(rpeakvNR),length(rpeakvNR));
        [corr_tempNR(h), p_tempNR(h)] = corr(rpeakvNR(rand_inxNR),rdecelNR);
    end
    vd_corr_randNR(a) = mean(corr_tempNR);
    vd_corr_p_randNR(a) = mean(p_tempNR);
    
    
    %for RR
    pf = polyfit(rpeakv,rdecel,1);
    fits(a,:) = pf;
    % for no re reach
    pfNR = polyfit(rpeakvNR,rdecelNR,1);
    fitsNR(a,:) = pfNR;
    
    %for plotting all linear correlations
    xs = [min(rpeakv) max(rpeakv)];
    pv = polyval(pf,xs);
    if vd_corr_p(a) <= .05
        col = [0 0.5 0.5];
    else
        col = [0.5 .5 .5 ];
    end
    plot(xs,pv,'Color',col,'LineWidth',1.5)

    % for no rereach
    %for plotting all linear correlations
    xsNR = [min(rpeakvNR) max(rpeakvNR)];
    pvNR = polyval(pfNR,xsNR);
    if vd_corr_pNR(a) <= .05
        colNR = [ 0.3984    0.8008    0.6641 ];
    else
        colNR = [0 0 0 ];
    end
    plot(xsNR,pvNR,'Color',colNR,'LineWidth',1.5)
    
end

%for plotting all linear correlations and mean
xs = [0 40];
pv = polyval(mean(fits),xs);
plot(xs,pv,'Color',[0 0.5 0.5],'LineWidth',4)
xlim(xs)
%ylim([-500 0])


% for plotting no re reach
xsNR = [0 40];
pvNR = polyval(mean(fitsNR),xsNR);
plot(xsNR,pvNR,'Color',[ 0.3984    0.8008    0.6641],'LineWidth',4)
xlim(xs)

figure; hold on
scatter(ones(length(ar),1),vd_corr,30,bmap(1,:),'filled','jitter','on','jitterAmount',.15,'MarkerFaceAlpha',.5)
scatter(ones(length(ar),1)*1.5,vd_corr_rand,30,[.8 .8 .8],'filled','jitter','on','jitterAmount',.15)
scatter(1,mean(vd_corr),170,bmap(1,:),'filled')
scatter(1.5,mean(vd_corr_rand),170,[.6 .6 .6],'filled')
% errorbar([1 1.5],[mean(dd_corr) mean(dd_corr_rand)],[std(dd_corr) std(dd_corr_rand)],...
%     'LineStyle','none','Color','k','CapSize',14)
% xlim([.7 1.8])
% ylim([-1 .1])
set(gcf,'Position',[100 100 300 400])


numcorr = sum(vd_corr_p <= .05)
meancorr = mean(vd_corr)
stecorr = std(vd_corr)/sqrt(26)

meancorr_rand = mean(vd_corr_rand)

%% Find min and max of accel during CSM time window ( 100ms after stim ) 
for a=1:length(an)
    stimvel=ExcitRR.(an{a}).S.stimvel;
    Accelmin=[];
    Accelmax=[];
   stimvelsize=size(stimvel);
   num=(stimvelsize(3));

    for b=1:num
        reach=stimvel(:,:,b);
        %window=reach(21:41,:);% row numbers for time of stim on to 100ms after - only for reaches with 200ms 
        window=reach(51:61,:);%
        accelmin=min(window(:,8));
        accelmax=max(window(:,8));
        Accelmin=[Accelmin;accelmin];
        Accelmax=[Accelmax;accelmax];
    end

    ExcitRR.(an{a}).AccelMin=Accelmin;
    ExcitRR.(an{a}).AccelMax=Accelmax;
    clear Accelmin Accelmax 
end



%% Plot Peaks of decel and CSM accel
AccelMax=[];
AccelMin=[];
DecelAll=[];
VMax=[];
RMSE=[];
RSquared=[];
figure;hold on
for a=1:length(an)
    Accelmin=ExcitRR.(an{a}).AccelMin;
    Accelmax=ExcitRR.(an{a}).AccelMax;
    Vmax=ExcitRR.(an{a}).RRKin.vmax(:,2);
    Decel=ExcitRR.(an{a}).RRKin.aavg_decel(:,2);
    AccelMax=[AccelMax;Accelmax];
    AccelMin=[AccelMin;Accelmin];
    DecelAll=[DecelAll;Decel];
    VMax=[VMax;Vmax];
    mdl=fitlm(Accelmin,Accelmax,'linear');
    rmse=mdl.RMSE;
    RMSE=[RMSE;rmse];
    rsquared=mdl.Rsquared.Adjusted;
    RSquared=[RSquared;rsquared];
    disp(mdl)
    %figure; hold on
    scatter(Accelmin, Accelmax,[],Vmax,'filled')
    colormap(brewermap([],('YlGnBu')))
    colorbar
    hold on
    h= plot(mdl)
    delete(h(1))
    set(h(2),'Color','k')
    set(h(3),'Color','k')
    set(h(4),'Color','k')
    
    xlim([-80 0])
    ylim([0 80])
    %set(gca,'XDir','reverse')
end
ylabel('AccelMax')
xlabel('AccelMin')
title('IntA- All')
%% Plot Peaks of decel and CSM accel
figure;hold on
for a=1:length(an)
    Accelmin=ExcitRR.(an{a}).AccelMin;
    Accelmax=ExcitRR.(an{a}).AccelMax;

    scatter(Accelmin, Accelmax)
end

%% 
KinVars_All=[VMax,AccelMin,AccelMax];

scatter(VMax,DecelAll)

%% Plot control and CSM -ACCEL
% Load KinVarsVariable for IntA-RN behav
% positive correlation
stimcolor= rgb('LightSteelBlue');
ctrlreach=rgb('Black');

poscolor=rgb('Tomato');negcolor=rgb('SteelBlue');

%neg variables
x=FR_Ctrlaccel_neg(:,2);
y=FR_Stimaccel_neg(:,2);
%pos vairables 
x1=FR_Ctrlaccel_pos(:,2);
y1=FR_Stimaccel_pos(:,2);

%CorrCoef
[r, p] = corr(x,y); [r1, p1] = corr(x1,y1);

% Fit a linear model
mdl = fitlm(x, y);
mdl_1=fitlm(x1,y1);

% Sort the x data for plotting purposes
[x_sorted, idx] = sort(x); [x_sorted1, idx1] = sort(x1);
y_sorted = y(idx); y_sorted1 = y1(idx1); 


% Predict the response and confidence intervals for sorted x data
[y_pred, y_ci] = predict(mdl, x_sorted, 'Alpha', 0.05); % 95% confidence intervals
[y_pred1, y_ci1] = predict(mdl_1, x_sorted1, 'Alpha', 0.05);

figure;
scatter(x, y,30,'filled','MarkerEdgeColor','k','MarkerFaceColor',negcolor);
hold on;
scatter(x1, y1,30,'filled','MarkerEdgeColor','k','MarkerFaceColor',poscolor);
hold on;
% Plot the fitted regression line
plot(x_sorted, y_pred,  'Color',negcolor, 'LineWidth', 1.5); % Red regression line
hold on
plot(x_sorted1, y_pred1,  'Color',poscolor, 'LineWidth', 1.5);
hold on
% Plot the confidence intervals
fill([x_sorted; flipud(x_sorted)], [y_ci(:,1); flipud(y_ci(:,2))], ...
    [ 0.6875    0.7656    0.8672], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on
fill([x_sorted1; flipud(x_sorted1)], [y_ci1(:,1); flipud(y_ci1(:,2))], ...
    [ 1.0000    0.3867    0.2773], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 

text(min(x), max(y), sprintf('r-neg = %.2f\n95%% p-value: [%.2f]', r, p), 'FontSize', 12, 'VerticalAlignment', 'top');
text(min(x1), max(y1), sprintf('r-pos = %.2f\n95%% p-value: [%.2f]', r1, p1), 'FontSize', 12, 'VerticalAlignment', 'bottom');
unity_line = plot(xlim, xlim, 'k--'); % Black dashed unity line

xline(0,'--')
yline(0,'--')
% Adding labels and title
xlabel('Control (Hz)');
ylabel('Stim (Hz)');
xlim([0 300])
ylim([0 400])
title('Peak FR during CSM-accel');
grid off;
hold off;

%%
mdl=fitlm(KinVars(:,2),KinVars(:,3))
figure; hold on
scatter(KinVars(:,2),KinVars(:,3),[],KinVars(:,1),'filled')
colormap(brewermap([],('YlGnBu')))
colorbar
hold on
h= plot(mdl)
delete(h(1))
set(h(2),'Color','k')
set(h(3),'Color','k')
set(h(4),'Color','k')

xlim([-80 0])
ylim([0 80])
%%
%corr per animal

vd_corr = []; 
vd_corr_p = []; 

vd_corr_rand = [];
vd_corr_p_rand = [];

fits = [];



RRcolor=rgb('Teal');

%for plotting all linear correlations
figure; hold on

for a = 1:length(an)
    
    Accelmin=ExcitRR.(an{a}).AccelMin; Accelmin=Accelmin*10;
    Accelmax=ExcitRR.(an{a}).AccelMax;Accelmax=Accelmax*10;

 
    %compute correlation coefficient and pvalues
    [vd_corr(a), vd_corr_p(a)] = corr(Accelmin,Accelmax);
    
    
    %shuffle control for ReREach
    corr_temp = NaN(1,1000); p_temp = NaN(1,1000);
    for j = 1:1000
        rand_inx = randsample(length(Accelmin),length(Accelmin));
        [corr_temp(j), p_temp(j)] = corr(Accelmin(rand_inx),Accelmax);
    end
    vd_corr_rand(a) = mean(corr_temp);
    vd_corr_p_rand(a) = mean(p_temp);

   
    
    %for RR
    pf = polyfit(Accelmin,Accelmax,1);
    fits(a,:) = pf;
  
    
    %for plotting all linear correlations
    xs = [min(Accelmin) max(Accelmax)];
    pv = polyval(pf,xs);
    if vd_corr_p(a) <= .05
        col = [0 0.5 0.5 0.5];
    else
        col = [0.5 .5 .5 ];
    end
    plot(xs,pv,'Color',col,'LineWidth',1.5)

   
end

%for plotting all linear correlations and mean
xs = [-400 0];
pv = polyval(mean(fits),xs);
plot(xs,pv,'Color',[0 0.5 0.5],'LineWidth',4)
%xlim(xs)
ylim([-500 0])
ylim([0 500])
ylabel('CSM peak accel')
xlabel('Peak Decel')

figure; hold on
scatter(ones(length(an),1),vd_corr,15,RRcolor,'filled','jitter','on','jitterAmount',.15,'MarkerFaceAlpha',.5)
scatter(ones(length(an),1)*1.5,vd_corr_rand,30,[.8 .8 .8],'filled','jitter','on','jitterAmount',.15)
scatter(1,mean(vd_corr),170,RRcolor,'filled')
scatter(1.5,mean(vd_corr_rand),170,[.6 .6 .6],'filled')
% errorbar([1 1.5],[mean(dd_corr) mean(dd_corr_rand)],[std(dd_corr) std(dd_corr_rand)],...
%     'LineStyle','none','Color','k','CapSize',14)
% xlim([.7 1.8])
ylim([-0.7 .1])
set(gcf,'Position',[100 100 300 400])


numcorr = sum(vd_corr_p <= .05)
meancorr = mean(vd_corr)
stecorr = std(vd_corr)/sqrt(26)

meancorr_rand = mean(vd_corr_rand)

%% Plot average accelertion and peak velocity : linear regression for each animal
%plot avg. accel and peak velocity linear regression per animal 
CtrlRValues=[]; CtrlRSquared=[]; CtrlSlope=[]; CtrlIntercept=[];
StimRValues=[]; StimRSquared=[];StimSlope=[];StimIntercept=[];
PlusRValues=[]; PlusRSquared=[];PlusSlope=[];PlusIntercept=[];

pcolor=rgb('DarkSeaGreen');
Pcolor=rgb('DarkGreen');
scolor=rgb('LightSteelBlue');
Scolor=rgb('RoyalBlue');

for a = 1:length(an)
    % get aavg accel for each trial 
    Ctrlacc= ChR2.(an{a}).CtrlKin.aavg_decel(:,2); Ctrlacc=Ctrlacc*10;
    Stimacc=ChR2.(an{a}).StimKin.aavg_accel(:,2);
    Plus1acc=ChR2.(an{a}).Plus1Kin.aavg_decel(:,2);Plus1acc=Plus1acc*10;
    
    % Calc mean of acc
    Ctrlacc_mean=mean(Ctrlacc);
    Stimacc_mean=mean(Stimacc);
    Plus1acc_mean=mean(Plus1acc);
    
    %get vmax 
    Ctrlvmax= ChR2.(an{a}).CtrlKin.vmax(:,2);Ctrlvmax=Ctrlvmax*10;
    Stimvmax=ChR2.(an{a}).StimKin.vmax(:,2);
    Plus1vmax=ChR2.(an{a}).Plus1Kin.vmax(:,2);Plus1vmax=Plus1vmax*10;
    
    % calculate R linear regression for each one and plot on separate graph
    RR_reg=polyfit(Ctrlvmax,Ctrlacc,1);
    CtrlyFit=polyval(Ctrlcoefficients,Ctrlvmax);
    Ctrlrvalue = corrcoef(Ctrlvmax,Ctrlacc);
    CtrlrValue = Ctrlrvalue(1, 2);
    CtrlR_squared=Ctrlrvalue(1,2)^2;
    CtrlRsqaured(a)=CtrlR_squared;
    CtrlRValues(a)=CtrlrValue;
    
    % get linear regression model and save 
    Ctrlslope=Ctrlcoefficients(1);
    Ctrlintercept=Ctrlcoefficients(2);
    CtrlSlope(a)=Ctrlslope;
    CtrlIntercept(a)=Ctrlintercept;
    
    % calculate R linear regression for each one and plot on separate graph
    Stimcoefficients=polyfit(Stimvmax,Stimacc,1);
    StimyFit=polyval(Stimcoefficients,Stimvmax);
    Stimrvalue = corrcoef(Stimvmax,Stimacc);
    StimrValue = Stimrvalue(1, 2);
    StimR_squared=Stimrvalue(1,2)^2;
    StimRsqaured(a)=StimR_squared;
    StimRValues(a)=StimrValue;
    
    % get linear regression model and save 
    Stimslope=Stimcoefficients(1);
    Stimintercept=Stimcoefficients(2);
    StimSlope(a)=Stimslope;
    StimIntercept(a)=Stimintercept;
    
    % calculate R for plus 1 trials 
    Pluscoefficients=polyfit(Plus1vmax,Plus1acc,1);
    PlusyFit=polyval(Pluscoefficients,Plus1vmax);
    Plusrvalue = corrcoef(Plus1vmax,Plus1acc);
    PlusrValue = Plusrvalue(1, 2);
    PlusR_squared=Plusrvalue(1,2)^2;
    PlusRsqaured(a)=PlusR_squared;
    PlusRValues(a)=PlusrValue;
    
    
    Plusslope=Pluscoefficients(1);
    Plusintercept=Pluscoefficients(2);
    PlusSlope(a)=Plusslope;
    PlusIntercept(a)=Plusintercept;
    
    
    std_X=std(Ctrlacc);
    std_Y=std(Stimacc);
    std_Y=std(Plus1acc);
    sem_X = std(Ctrlacc) / sqrt(length(Ctrlacc)); % SEM for X
    sem_Y = std(Plus1acc) / sqrt(length(Plus1acc)); % SEM for Y
    %sem_Y = std(Stimacc) / sqrt(length(Stimacc)); % SEM for Y
    
    %scatter(Ctrlacc_mean,Plus1acc_mean,'filled','MarkerFaceColor',[0.5 0.5 0.5])
    %scatter(Ctrlacc_mean,Stimacc_mean,'filled')
    %hold on
   % %errorbar(Ctrlaccmean, Plusaccmean, std_X, std_X, 'r.', 'MarkerSize', 10, 'LineWidth', 2); % Red error bars
   % errorbar(Ctrlacc_mean, Plus1acc_mean, sem_X, sem_X, sem_Y, sem_Y,'o','Color','k','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize', 10, 'LineWidth', 1); % Red error bars
   % errorbar(Ctrlacc_mean, Stimacc_mean, sem_X, sem_X, sem_Y, sem_Y, 'r.', 'MarkerSize', 10, 'LineWidth', 2);
    %labels=string(an);


   % [h,ext]=labelpoints(Ctrlacc_mean,Plus1acc_mean,labels(a))
%     %[h,ext]=labelpoints(Ctrlacc_mean,Stimacc_mean,labels(a))
%       xlim([0 10])
%       ylim([0 10])
% %     ylabel('Plus1acc')
% %     xlabel('Ctrlacc')
%     ylabel('Plus1decel')
%     xlabel('Ctrldecel')
%     slope=1;
%     intercept=0;
%     refline(slope,intercept)
%     hold on 
    
    
     scatter(Ctrlvmax,Ctrlacc,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','k')
%      hold on
 %     plot(Ctrlvmax, CtrlyFit,'Color',[0.5 0.5 0.5], 'LineWidth', 2); % Fitted line
%     hold on;
%     scatter(Stimvmax,Stimacc,'MarkerFaceColor',[0.6875    0.7656    0.8672],'MarkerEdgeColor','k')
%     hold on
%     plot(Stimvmax, StimyFit,'Color',[0.6875    0.7656    0.8672], 'LineWidth', 2); %
     hold on
     scatter(Plus1vmax,Plus1acc,'MarkerFaceColor',[ 0.5586    0.7344    0.5586],'MarkerEdgeColor','k')
%     hold on
%     plot(Plus1vmax, PlusyFit,'Color',[ 0.5586    0.7344    0.5586], 'LineWidth', 2); %
%     hold on
    
%     [Ctrlh,Ctrlp]=kstest2(Ctrlacc,Plus1acc)
% % Make Cd plot 
%     figure()
%     cdfplot(Ctrlacc)
% 
%     hold on
%     cdfplot(Plus1acc)
%     
end 


%calculate grand mean for Ctrl trials 
Ctrlgrand_mean_slope = mean(CtrlSlope);
Ctrlgrand_mean_intercept = mean(CtrlIntercept);
%CtrlX_range = [min([Ctrlvmax]):0.1:max([Ctrlvmax])];
CtrlX_range=[0 40];
% Calculate corresponding Y values for the grand mean line
CtrlY_grand_mean = Ctrlgrand_mean_slope * CtrlX_range +Ctrlgrand_mean_intercept;

%Calculate grand mean for Stim 
Stimgrand_mean_slope = mean(StimSlope);
Stimgrand_mean_intercept = mean(StimIntercept);
StimX_range = [min([Stimvmax]):0.1:max([Stimvmax])];
% Calculate corresponding Y values for the grand mean line
StimY_grand_mean = Stimgrand_mean_slope * StimX_range + Stimgrand_mean_intercept;

%Calculate grand mean for PLus 1 
Plusgrand_mean_slope = mean(PlusSlope);
Plusgrand_mean_intercept = mean(PlusIntercept);
%PlusX_range = [min([Plus1vmax]):0.1:max([Plus1vmax])];
PlusX_range=[0 40];
% Calculate corresponding Y values for the grand mean line
PlusY_grand_mean = Plusgrand_mean_slope * PlusX_range + Plusgrand_mean_intercept;

% Overlay the grand mean line
%figure()
hold on
plot(CtrlX_range, CtrlY_grand_mean,'Color',[0 0 0], 'LineWidth', 3, 'DisplayName', 'Grand Mean Line');
hold on
%plot(StimX_range, StimY_grand_mean,'Color',[ 0.2539    0.4102    0.8789], 'LineWidth', 3, 'DisplayName', 'Grand Mean Line');
%hold on
plot(PlusX_range, PlusY_grand_mean,'Color',[ 0.5586    0.7344    0.5586], 'LineWidth', 3, 'DisplayName', 'Grand Mean Line');
xlim([0 40])
%ylim([0 300])
ylim([-120 0])


xlabel('Peak Outward Velocity')
ylabel('Avg. Decel.')
%ylabel('Avg. Accel')
legend({['Ctrl slope = ' num2str(Ctrlgrand_mean_slope, '%.4f')], ['Plus1 slope = ' num2str(Plusgrand_mean_slope, '%.4f')]}, 'Location', 'SouthEast');

hold off;

%% get rid of zeros 
excitam = fieldnames(ExcitRR);
for b=1:length(excitam)
    RR_Excit=ExcitRR.(excitam{b}).StimRR;
 
    nonZeroRows_ex = [];

    % Iterate over each row of the array
    for i = 1:size(RR_Excit, 1)
        % Check if the current row contains any zeros
        if any(RR_Excit(i,:) == 0)
            % If zeros are present, skip this row
            continue;
        end
        
        % Append the current row to nonZeroRows if no zeros are present
        nonZeroRows_ex = [nonZeroRows_ex; RR_Excit(i,:)];
    end
    
        % Display the resulting array
        disp(nonZeroRows_ex);
        ExcitRR.(excitam{b}).StimRREdit=nonZeroRows_ex;
end

%% PLot histogram of RR timing for each animal
allcolor=rgb('Teal');

folder_path="C:\Users\Person Lab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Random stim\AllNeuron_MattandCourtney Data_11\RRTiming_IndividualAnimals";

for a=1:length(an)
    rrt=ExcitRR.(an{a}).StimRREdit;
    rrt=rrt(:,1); %timing of RR
    % num_binsall = ceil(sqrt(length(rrt)));
    num_binsall=20;
    figure; hold on
    histogram(rrt,num_binsall,'Normalization','Probability','EdgeColor','none','FaceColor',allcolor)
    xlim([0 0.25])
    xlabel('Time from Stim Onset (s)')
    title(an(a))
    
    filename=sprintf('ReReachTiming%d.tif',a);
    full_path=fullfile(folder_path,filename);
    saveas(gcf,full_path);
    
    close;
end

%% Get position after stim off for all stim reaches 

for a=1:length(an)
    Sreaches=ExcitRR.(an{a}).S.stimvel;
    fields=fieldnames(ExcitRR.(an{a}).ReReaches);
    
    StimEndPos=[];
    for i=1:length(fields)
        %endpos=Sreaches(46,2,i)
         endpos=Sreaches(41,2,i)
        StimEndPos=[StimEndPos;endpos];
    end
    ExcitRR.(an{a}).StimEndPos=StimEndPos;
end

%% Plot CSM timing and EndPos
%colors=brewermap(11,'Set3')
colors=colormap(slanCM('Greens',11))
%folder_path="C:\Users\Person Lab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Random stim\AllNeuron_MattandCourtney Data_11\RRTiming_IndividualAnimals\CSMvsEndStimPos";

figure; hold on
for i=1:length(an)
    csmtime=ExcitRR.(an{i}).StimRR(:,1);
    endpos=ExcitRR.(an{i}).StimEndPos;

    %figure; hold on
    scatter(csmtime,endpos,20,'filled','MarkerFaceColor',colors(i,:),'MarkerEdgeColor','k')
    % hold on

    % xlim([0.03 0.25])
    % ylim([0.1 0.26])
    % xline(0.054,'k','LineWidth',1.5)
    % xline(0.0864,'k','LineWidth',1.5)
    % xline(0.0972,'b','LineWidth',1.5)
    % xline(0.1296,'b--','LineWidth',2)
    % xline(0.1296,'r','LineWidth',1.5)
    % xline(0.162,'r','LineWidth',1.5)
    % yline(0.16,'k','LineWidth',1.5)
    % xlabel('CSM time ')
    % ylabel('Stim End Position')
    % 
    % filename=sprintf('CSMvsStimPos%d.tif',i);
    % full_path=fullfile(folder_path,filename);
    % saveas(gcf,full_path);
    % 
    % close;

   
end
xlim([0.03 0.25])
%ylim([0.1 0.26])
xline(0.054,'k','LineWidth',1.5)
xline(0.0864,'k','LineWidth',1.5)
xline(0.0972,'b','LineWidth',1.5)
xline(0.1296,'b--','LineWidth',2)
xline(0.1296,'r','LineWidth',1.5)
xline(0.162,'r','LineWidth',1.5)
yline(0.16,'k','LineWidth',1.5)
xlabel('CSM time ')
ylabel('Stim Start Position')
%% Separate out CSMs based off of timing of the CSM 

for a=1:length(an)
    Rreaches=ExcitRR.(an{a}).ReReaches;
    rrt=ExcitRR.(an{a}).StimRR; rrt=rrt(:,1);
 
    % Timings taken from histo from all animals 
    T1Reaches=[];
    T2Reaches=[];
    T3Reaches=[];
    T4Reaches=[];
    T5Reaches=[];

    %T1=find(rrt<=0.032);
    % T2=find(rrt>0.032 & rrt<=0.065);
    % T3=find(rrt>0.065 & rrt<=0.108);
    % T4=find(rrt>0.108 & rrt<=0.140);
    % T5=find(rrt>0.140);

    % new times based off of bin edges
    T1=find(rrt<=0.054);
    T2=find(rrt>0.054 & rrt<=0.0864);
    T3=find(rrt>=0.0972 & rrt<=0.1296);
    T4=find(rrt>0.1296 & rrt<=0.162);
    T5=find(rrt>0.162);
    rnum=fieldnames(Rreaches);

    for f=1:length(T1)
        idx1=T1(f);
        T1Reaches.(rnum{idx1})= ExcitRR.(an{a}).ReReaches.(rnum{idx1});  
    end 

  
    if isempty(T1Reaches)
        T1Kin=[];
        ST1=[];
    else
         T1Kin=BasicKinematics(T1Reaches);
        ST1=AlignStim(T1Reaches,400);
        %ST1=stimvel_accel(ST1);
    end 

    for h=1:length(T2)
        idx2=T2(h);
        T2Reaches.(rnum{idx2})= ExcitRR.(an{a}).ReReaches.(rnum{idx2});  
    end 
    T2Kin=BasicKinematics(T2Reaches);
    ST2=AlignStim(T2Reaches,400);
    ST2=stimvel_accel(ST2);

    for k=1:length(T3)
        idx3=T3(k);
        T3Reaches.(rnum{idx3})= ExcitRR.(an{a}).ReReaches.(rnum{idx3});  
    end 
    T3Kin=BasicKinematics(T3Reaches);
    ST3=AlignStim(T3Reaches,400);
    ST3=stimvel_accel(ST3);

     for j=1:length(T4)
        idx4=T4(j);
        T4Reaches.(rnum{idx4})= ExcitRR.(an{a}).ReReaches.(rnum{idx4});  
     end 
    if isempty(T4Reaches)
        T4Kin=[];
        ST4=[];
    else
         T4Kin=BasicKinematics(T4Reaches);
        ST4=AlignStim(T4Reaches,400);
        %ST1=stimvel_accel(ST4);
    end 

 


     for p=1:length(T5)
        idx5=T5(p);
        T5Reaches.(rnum{idx5})= ExcitRR.(an{a}).ReReaches.(rnum{idx5});  
     end 

     if isempty(T5Reaches)
        T5Kin=[];
        ST5=[];
    else
         T5Kin=BasicKinematics(T5Reaches);
        ST5=AlignStim(T5Reaches,400);
        %ST5=stimvel_accel(ST5);
    end 
   


    

    
    % Reaches 
     ExcitRR.(an{a}).RRT1=T1Reaches;
     ExcitRR.(an{a}).RRT2=T2Reaches;
     ExcitRR.(an{a}).RRT3=T3Reaches;
     ExcitRR.(an{a}).RRT4=T4Reaches;
     ExcitRR.(an{a}).RRT5=T5Reaches;

     % Kinematics
     ExcitRR.(an{a}).T1Kin=T1Kin;
     ExcitRR.(an{a}).T2Kin=T2Kin;
     ExcitRR.(an{a}).T3Kin=T3Kin;
     ExcitRR.(an{a}).T4Kin=T4Kin;
     ExcitRR.(an{a}).T5Kin=T5Kin;

     %Stimvel
     ExcitRR.(an{a}).ST1=ST1;
     ExcitRR.(an{a}).ST2=ST2;
     ExcitRR.(an{a}).ST3=ST3;
     ExcitRR.(an{a}).ST4=ST4;
     ExcitRR.(an{a}).ST5=ST5;

end 

% ST1, ST4, ST5 - not all animals have these 
%%
% Take out T1 from animal 3 

ST1_Stimvel=[];
indices_to_delete = [];

for a=1:length(an)
    ST1=ExcitRR.(an{a}).ST1;
    if isempty(ST1)
        ST1_Stimvel.(an{a}).ST1=[];
        % Add the index of the empty NRR.stimvel to the list
        indices_to_delete(end+1) = a;
    else
    ST1=stimvel_accel(ST1);
    ST1_Stimvel.(an{a}).ST1=ST1;
    ST1_Stimvel.(an{a}).C= ExcitRR.(an{a}).C;
    end
end 

% Delete structures where NRR.stimvel is empty
for i = length(indices_to_delete):-1:1
    ST1_Stimvel = rmfield(ST1_Stimvel, an{indices_to_delete(i)});
end

%% 
ST4_Stimvel=[];
indices_to_delete = [];

for a=1:length(an)
    ST4=ExcitRR.(an{a}).ST4;
    if isempty(ST4)
        ST4_Stimvel.(an{a}).ST4=[];
        % Add the index of the empty NRR.stimvel to the list
        indices_to_delete(end+1) = a;
    else
    ST4=stimvel_accel(ST4);
    ST4_Stimvel.(an{a}).ST4=ST4;
    ST4_Stimvel.(an{a}).C= ExcitRR.(an{a}).C;
    end
end 

% Delete structures where NRR.stimvel is empty
for i = length(indices_to_delete):-1:1
    ST4_Stimvel = rmfield(ST4_Stimvel, an{indices_to_delete(i)});
end

%% 
ST5_Stimvel=[];
indices_to_delete = [];

for a=1:length(an)
    ST5=ExcitRR.(an{a}).ST5;
    if isempty(ST5)
        ST5_Stimvel.(an{a}).ST5=[];
        % Add the index of the empty NRR.stimvel to the list
        indices_to_delete(end+1) = a;
    else
    ST5=stimvel_accel(ST5);
    ST5_Stimvel.(an{a}).ST5=ST5;
    ST5_Stimvel.(an{a}).C= ExcitRR.(an{a}).C;
    end
end 

% Delete structures where NRR.stimvel is empty
for i = length(indices_to_delete):-1:1
    ST5_Stimvel = rmfield(ST5_Stimvel, an{indices_to_delete(i)});
end
  %% Combine Stimvel and look at plots 

T2_big=struct;
T3_big = struct;
% T4_big=struct;
% T5_big=struct;


ar=fieldnames(ExcitRR);
for a = 1:length(ar) 
  
   T2_big.stimvel(:,:,a) = ExcitRR.(ar{a}).ST2.Average;
   T3_big.stimvel(:,:,a) = ExcitRR.(ar{a}).ST3.Average;
   % T4_big.stimvel(:,:,a) = ExcitRR.(ar{a}).ST4.Average;
   % T5_big.stimvel(:,:,a) = ExcitRR.(ar{a}).ST5.Average;
    
end


T2_big.Average = nanmean(T2_big.stimvel,3);
T2_big.Std = nanstd(T2_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T2_big.Average),1); 
for nnum = 1:length(T2_big.Average)
    samp = squeeze(T2_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
T2_big.N = N;

%
T3_big.Average = nanmean(T3_big.stimvel,3);
T3_big.Std = nanstd(T3_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T3_big.Average),1); 
for nnum = 1:length(T3_big.Average)
    samp = squeeze(T3_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
T3_big.N = N;

%
% T4_big.Average = nanmean(T4_big.stimvel,3);
% T4_big.Std = nanstd(T4_big.stimvel,0,3);
% 
% %loop through to count the number of reaches contributing to each time pt
% N = zeros(length(T4_big.Average),1); 
% for nnum = 1:length(T4_big.Average)
%     samp = squeeze(T4_big.stimvel(nnum,2,:));
%     N(nnum) = length(samp(~isnan(samp)));
% 
% end
% T4_big.N = N;
% 
% %
% T5_big.Average = nanmean(T5_big.stimvel,3);
% T5_big.Std = nanstd(T5_big.stimvel,0,3);
% 
% %loop through to count the number of reaches contributing to each time pt
% N = zeros(length(T5_big.Average),1); 
% for nnum = 1:length(T5_big.Average)
%     samp = squeeze(T5_big.stimvel(nnum,2,:));
%     N(nnum) = length(samp(~isnan(samp)));
% 
% end
% T5_big.N = N;



%% No tiing for certain animals 
T1_big=struct;

ar=fieldnames(ST1_Stimvel);
for a = 1:length(ar) 
    T1_big.stimvel(:,:,a) = ST1_Stimvel.(ar{a}).ST1.Average;
end
%NoReReach
T1_big.Average = nanmean(T1_big.stimvel,3);
T1_big.Std = nanstd(T1_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T1_big.Average),1); 
for nnum = 1:length(T1_big.Average)
    samp = squeeze(T1_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
T1_big.N = N;

%% 
T4_big=struct;

ar=fieldnames(ST4_Stimvel);
for a = 1:length(ar) 
    T4_big.stimvel(:,:,a) = ST4_Stimvel.(ar{a}).ST4.Average;
end
%NoReReach
T4_big.Average = nanmean(T4_big.stimvel,3);
T4_big.Std = nanstd(T4_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T4_big.Average),1); 
for nnum = 1:length(T4_big.Average)
    samp = squeeze(T4_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
T4_big.N = N;
%%
T5_big=struct;

ar=fieldnames(ST5_Stimvel);
for a = 1:length(ar) 
    T5_big.stimvel(:,:,a) = ST5_Stimvel.(ar{a}).ST5.Average;
end
%NoReReach
T5_big.Average = nanmean(T5_big.stimvel,3);
T5_big.Std = nanstd(T5_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T5_big.Average),1); 
for nnum = 1:length(T5_big.Average)
    samp = squeeze(T5_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
T5_big.N = N;


%% Plot stimvel
x1=0 ;x2=0.05;
y_upper= 20; y_lower=-40;
%Add acceleration column to stimvel
T1color=rgb('OrangeRed');
T2color=rgb('LightBlue');
T3color=rgb('DodgerBlue');
T4color=rgb('MediumBlue');
T5color=rgb('MidnightBlue');




stimcolor= rgb('Teal');
linecolor=rgb('DarkGray');

dim =8;

StimvelPlot(C_big,dim,'k',1)
StimvelPlot(T1_big,dim,T1color,0)
StimvelPlot(T2_big,dim,T2color,0)
StimvelPlot(T3_big,dim,T3color,0)
StimvelPlot(T4_big,dim,T4color,0)
StimvelPlot(T5_big,dim,T5color,0)


hold on
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
%xlim([-0.2 0.2])
ylim('auto')
%ylim([-0.5 1.5])
ylim([-40 20])
%pbaspect([1 2 1])
%ylim([0 0.3])
yline(0,'--')
xlabel('Time from Stim (s)')
ylabel('Outward Velocity (cm/s)')
%ylabel('Outward acceleration')
%% Find position at time=0 and time=50ms aligned to stimn
T1_start=[]; T1_end=[];T1_endpts=[];
% T1 - one animal does not have this 
ar=fieldnames(ST1_Stimvel);
for a= 1:length(ar)
    t1_stimvel=ExcitRR.(ar{a}).ST1.Average;
    t1_start=t1_stimvel(41,2);
    t1_end=t1_stimvel(46,2);
    T1_start=[T1_start;t1_start];
    T1_end=[T1_end;t1_end];

    % endpoint
    t1_end=ExcitRR.(ar{a}).T1Kin.endpt;
    t1_endpts=mean(t1_end(:,1));
    T1_endpts=[T1_endpts; t1_endpts];
end

T1_dif=T1_end-T1_start;
%% all other 
T2_start=[]; T2_end=[];T2_endpts=[];
T3_start=[]; T3_end=[];T3_endpts=[];
% T4_start=[]; T4_end=[];T4_endpts=[];
% T5_start=[]; T5_end=[];T5_endpts=[];

for a=1:length(an)
    %T2
    t2_stimvel=ExcitRR.(an{a}).ST2.Average;
    t2_start=t2_stimvel(41,2);
    t2_end=t2_stimvel(46,2);
    T2_start=[T2_start;t2_start];
    T2_end=[T2_end;t2_end];

    % endpoint
    t2_end=ExcitRR.(an{a}).T2Kin.endpt;
    t2_endpts=mean(t2_end(:,1));
    T2_endpts=[T2_endpts; t2_endpts];

    % T3
    t3_stimvel=ExcitRR.(an{a}).ST3.Average;
    t3_start=t3_stimvel(41,2);
    t3_end=t3_stimvel(46,2);
    T3_start=[T3_start;t3_start];
    T3_end=[T3_end;t3_end];

    % endpoint
    t3_end=ExcitRR.(an{a}).T3Kin.endpt;
    t3_endpts=mean(t3_end(:,1));
    T3_endpts=[T3_endpts; t3_endpts];


    % % T4
    % t4_stimvel=ExcitRR.(an{a}).ST4.Average;
    % t4_start=t4_stimvel(41,2);
    % t4_end=t4_stimvel(46,2);
    % T4_start=[T4_start;t4_start];
    % T4_end=[T4_end;t4_end];
    % 
    % 
    % % endpoint
    % t4_end=ExcitRR.(an{a}).T4Kin.endpt;
    % t4_endpts=mean(t4_end(:,1));
    % T4_endpts=[T4_endpts; t4_endpts];
    % 
    % 
    % %T5
    % t5_stimvel=ExcitRR.(an{a}).ST5.Average;
    % t5_start=t5_stimvel(41,2);
    % t5_end=t5_stimvel(46,2);
    % T5_start=[T5_start;t5_start];
    % T5_end=[T5_end;t5_end];
    % 
    % 
    % % endpoint
    % t5_end=ExcitRR.(an{a}).T5Kin.endpt;
    % t5_endpts=mean(t5_end(:,1));
    % T5_endpts=[T5_endpts; t5_endpts];



end 

T2_dif=T2_end-T2_start;
T3_dif=T3_end-T3_start;
%% 
T4_start=[]; T4_end=[];T4_endpts=[];
% T1 - one animal does not have this 
ar=fieldnames(ST4_Stimvel);
for a= 1:length(ar)
    t4_stimvel=ExcitRR.(ar{a}).ST4.Average;
    t4_start=t4_stimvel(41,2);
    t4_end=t4_stimvel(46,2);
    T4_start=[T4_start;t4_start];
    T4_end=[T4_end;t4_end];

    % endpoint
    t4_end=ExcitRR.(ar{a}).T4Kin.endpt;
    t4_endpts=mean(t4_end(:,1));
    T4_endpts=[T4_endpts; t4_endpts];
end

T4_dif=T4_end-T4_start;
%% 
T5_start=[]; T5_end=[];T5_endpts=[];
% T1 - one animal does not have this 
ar=fieldnames(ST5_Stimvel);
for a= 1:length(ar)
    t5_stimvel=ExcitRR.(ar{a}).ST5.Average;
    t5_start=t5_stimvel(41,2);
    t5_end=t5_stimvel(46,2);
    T5_start=[T5_start;t5_start];
    T5_end=[T5_end;t5_end];

    % endpoint
    t5_end=ExcitRR.(ar{a}).T5Kin.endpt;
    t5_endpts=mean(t5_end(:,1));
    T5_endpts=[T5_endpts; t5_endpts];
end

T5_dif=T5_end-T5_start;
%% Control reaches 
% Find position at time=0 and time=50ms aligned to stimn
C_startAll=[]; C_endAll=[];

ar=fieldnames(ExcitRR);
for a= 1:length(ar)
    c_stimvel=ExcitRR.(ar{a}).C.Average;
    c_start=c_stimvel(41,2);
    c_end=c_stimvel(46,2);
    C_startAll=[C_startAll;c_start];
    C_endAll=[C_endAll;c_end];

end

C_dif=C_endAll-C_startAll;
%%  Plot position 

T1color=rgb('OrangeRed');
T2color=rgb('LightBlue');
T3color=rgb('DodgerBlue');
T4color=rgb('MediumBlue');
T5color=rgb('MidnightBlue');

figure; hold on
scatter(T1_start,T1_end,30,T1color,'filled')
hold on
scatter(T2_start,T2_end,30,T2color,'filled')
hold on
scatter(T3_start,T3_end,30,T3color,'filled')
hold on
scatter(T4_start,T4_end,30,T4color,'filled')
hold on
scatter(T5_start,T5_end,30,T5color,'filled')
hold on
scatter(C_startAll,C_endAll,30,'k','filled')

xlabel('Start stim position')
ylabel('End stim position')

%% Mean and SEM  
T1_start_m=mean(T1_start); T1_end_m=mean(T1_end);
T1_sem_X = std(T1_start) / sqrt(length(T1_start)); % SEM for X
T1_sem_Y = std(T1_end) / sqrt(length(T1_end)); % SEM for Y

T2_start_m=mean(T2_start); T2_end_m=mean(T2_end);
T2_sem_X = std(T2_start) / sqrt(length(T2_start)); % SEM for X
T2_sem_Y = std(T2_end) / sqrt(length(T2_end)); % SEM for Y

T3_start_m=mean(T3_start); T3_end_m=mean(T3_end);
T3_sem_X = std(T3_start) / sqrt(length(T3_start)); % SEM for X
T3_sem_Y = std(T3_end) / sqrt(length(T3_end)); % SEM for Y

T4_start_m=mean(T4_start); T4_end_m=mean(T4_end);
T4_sem_X = std(T4_start) / sqrt(length(T4_start)); % SEM for X
T4_sem_Y = std(T4_end) / sqrt(length(T4_end)); % SEM for Y

T5_start_m=mean(T5_start); T5_end_m=mean(T5_end);
T5_sem_X = std(T5_start) / sqrt(length(T5_start)); % SEM for X
T5_sem_Y = std(T5_end) / sqrt(length(T5_end)); % SEM for Y

%ctrl
C_start_m=mean(C_start); C_end_m=mean(C_end);
C_sem_X = std(C_start) / sqrt(length(C_start)); % SEM for X
C_sem_Y = std(C_end) / sqrt(length(C_end)); % SEM for Y


%% Plot mean with error bars 
figure; hold on
scatter(T1_start_m,T1_end_m,80,T1color,'filled','MarkerEdgeColor','k')
hold on
errorbar(T1_start_m, T1_end_m, T1_sem_X, T1_sem_X, T1_sem_Y, T1_sem_Y,'k', 'LineWidth', 1); 
hold on
scatter(T2_start_m,T2_end_m,80,T2color,'filled','MarkerEdgeColor','k')
hold on
errorbar(T2_start_m, T2_end_m, T2_sem_X, T2_sem_X, T2_sem_Y, T2_sem_Y,'k', 'LineWidth', 1); 
hold on
scatter(T3_start_m,T3_end_m,80,T3color,'filled','MarkerEdgeColor','k')
hold on
errorbar(T3_start_m, T3_end_m, T3_sem_X, T3_sem_X, T3_sem_Y, T3_sem_Y,'k','LineWidth', 1); 
hold on
scatter(T4_start_m,T4_end_m,80,T4color,'filled','MarkerEdgeColor','k')
hold on
errorbar(T4_start_m, T4_end_m, T4_sem_X, T4_sem_X, T4_sem_Y, T4_sem_Y,'k', 'LineWidth', 1); 
hold on
scatter(T5_start_m,T5_end_m,80,T5color,'filled','MarkerEdgeColor','k')
hold on
errorbar(T5_start_m, T5_end_m, T5_sem_X, T5_sem_X, T5_sem_Y, T5_sem_Y,'k', 'LineWidth', 1); 
hold on
scatter(C_start_m,C_end_m,80,'k','filled','MarkerEdgeColor','k')
hold on
errorbar(C_start_m, C_end_m, C_sem_X, C_sem_X, C_sem_Y, C_sem_Y,'k','LineWidth', 1); 

ylabel('End Stim Position')
xlabel('Start Stim Position')
ylim([0.185 0.22])
 