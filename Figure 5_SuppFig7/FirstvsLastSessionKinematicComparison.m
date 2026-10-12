% First session vs. Last session kinematic comparison 

% RUN through sesssions and align to average time of CSM 

% ExcitRR_Ctrl1
an=fieldnames(ExcitRR_Ctrl1);
thresh=0.16;
dir=1;
for a=1:length(an)
    % First Session
    CtrlR=ExcitRR_Ctrl1.(an{a}).CtrlReaches;
    %C_CSM=AlignThreshCSM(CtrlR,thresh, dir, 500);
    C_CSM=AlignThresh(CtrlR,thresh, dir, 500);
    C_CSM=stimvel_accel(C_CSM); C_CSM=stimvel_accel_jerk(C_CSM);

    StimR=ExcitRR.(an{a}).ReReaches;
    %S_CSM=AlignThreshCSM(StimR,thresh, dir, 500);
    S_CSM=AlignStim(StimR,500);
    S_CSM=stimvel_accel(S_CSM); S_CSM=stimvel_accel_jerk(S_CSM);

    %Last Session
    CtrlR_last=ExcitRR_Ctrlend.(an{a}).CtrlReaches;
    %Cend_CSM=AlignThreshCSM(CtrlR_last,thresh, dir, 500);
    Cend_CSM=AlignThresh(CtrlR_last,thresh, dir, 500);
    Cend_CSM=stimvel_accel(Cend_CSM); Cend_CSM=stimvel_accel_jerk(Cend_CSM);

    StimR_last=ExcitRR_end.(an{a}).ReReaches;
    %Send_CSM=AlignThreshCSM(StimR_last,thresh, dir, 500);
    Send_CSM=AlignStim(StimR_last,500);
    Send_CSM=stimvel_accel(Send_CSM); Send_CSM=stimvel_accel_jerk(Send_CSM);

    ExcitRR_Ctrl1.(an{a}).C_CSM=C_CSM;
    ExcitRR_Ctrlend.(an{a}).Cend_CSM=Cend_CSM;

    ExcitRR.(an{a}).S_CSM=S_CSM;
    ExcitRR_end.(an{a}).Send_CSM=Send_CSM;

end

%%
an = fieldnames(ExcitRR_Ctrl1);
thresh = 0.16;
dir = 1;
K = 5;

for a = 1:length(an)
    name = an{a};

    % ---- FIRST SESSION ----
    CtrlR = ExcitRR_Ctrl1.(name).CtrlReaches;
    StimR = ExcitRR.(name).ReReaches;

    subsets = {'first','middle','last'};

    for s = 1:length(subsets)
        subName = subsets{s};

        C_sub = subsetReachStruct(CtrlR, subName, K);
        S_sub = subsetReachStruct(StimR, subName, K);

        %ExcitRR_Ctrl1.(name).(['C_CSM_' subName]) = stimvel_accel_jerk( ...
            %stimvel_accel(AlignThreshCSM(C_sub, thresh, dir, 500)) );

        ExcitRR_Ctrl1.(name).(['C_CSM_' subName]) = stimvel_accel_jerk( ...
            stimvel_accel(AlignThresh(C_sub, thresh, dir, 500)) );

        % ExcitRR.(name).(['S_CSM_' subName]) = stimvel_accel_jerk( ...
        %     stimvel_accel(AlignThreshCSM(S_sub, thresh, dir, 500)) );

        ExcitRR.(name).(['S_CSM_' subName]) = stimvel_accel_jerk( ...
            stimvel_accel(AlignStim(S_sub,500)) );
    end

    % ---- LAST SESSION ----
    CtrlR_last = ExcitRR_Ctrlend.(name).CtrlReaches;
    StimR_last = ExcitRR_end.(name).ReReaches;

    for s = 1:length(subsets)
        subName = subsets{s};

        Cend_sub = subsetReachStruct(CtrlR_last, subName, K);
        Send_sub = subsetReachStruct(StimR_last, subName, K);

        % ExcitRR_Ctrlend.(name).(['Cend_CSM_' subName]) = stimvel_accel_jerk( ...
        %     stimvel_accel(AlignThreshCSM(Cend_sub, thresh, dir, 500)) );

        ExcitRR_Ctrlend.(name).(['Cend_CSM_' subName]) = stimvel_accel_jerk( ...
            stimvel_accel(AlignThresh(Cend_sub, thresh, dir, 500)) );

        ExcitRR_end.(name).(['Send_CSM_' subName]) = stimvel_accel_jerk( ...
            stimvel_accel(AlignThreshCSM(Send_sub, thresh, dir, 500)) );

        ExcitRR_end.(name).(['Send_CSM_' subName]) = stimvel_accel_jerk( ...
            stimvel_accel(AlignStim(Send_sub,500)) );
    end
end

%% Now align Sham reaches 

am=fieldnames(ExcitRR_Sham);

for l=1:length(am)
    ShamR=ExcitRR_Sham.(am{l}).CtrlReaches;
    Sham_CSM=AlignThreshCSM(ShamR,thresh, dir, 500);
    Sham_CSM=stimvel_accel(Sham_CSM); Sham_CSM=stimvel_accel_jerk(Sham_CSM);
    
    ExcitRR_Sham.(am{l}).Sham_CSM=Sham_CSM;

end

 %% Combine Stimvel and look at plots 
C_big=struct;
T0_big=struct;
T1_big=struct;
T2_big = struct;
%T3_big=struct;

ar=fieldnames(ExcitRR_Ctrl1);
for a = 1:length(ar) 
   %C_big.stimvel(:,:,a) = ExcitRR_Ctrl1.(ar{a}).C.Average;
    C_big.stimvel(:,:,a) = ExcitRR_Ctrl1.(ar{a}).C_CSM.Average;
   T0_big.stimvel(:,:,a) = ExcitRR_Ctrl1.(ar{a}).ST0.Average;
   T1_big.stimvel(:,:,a) = ExcitRR_Ctrl1.(ar{a}).ST1.Average;
   T2_big.stimvel(:,:,a) = ExcitRR_Ctrl1.(ar{a}).ST2.Average;
  % T3_big.stimvel(:,:,a) = ExcitRR_Ctrl1.(ar{a}).ST3.Average;
 
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



T0_big.Average = nanmean(T0_big.stimvel,3);
T0_big.Std = nanstd(T0_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T0_big.Average),1); 
for nnum = 1:length(T0_big.Average)
    samp = squeeze(T0_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
T0_big.N = N;

T1_big.Average = nanmean(T1_big.stimvel,3);
T1_big.Std = nanstd(T1_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T1_big.Average),1); 
for nnum = 1:length(T1_big.Average)
    samp = squeeze(T1_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
T1_big.N = N;


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
% T3_big.Average = nanmean(T3_big.stimvel,3);
% T3_big.Std = nanstd(T3_big.stimvel,0,3);
% 
% %loop through to count the number of reaches contributing to each time pt
% N = zeros(length(T3_big.Average),1); 
% for nnum = 1:length(T3_big.Average)
%     samp = squeeze(T3_big.stimvel(nnum,2,:));
%     N(nnum) = length(samp(~isnan(samp)));
% 
% end
% T3_big.N = N;


%% Combine Stimvel and look at plots - for contorl first session
C_big=struct;
T0_big=struct;
T1_big=struct;
T2_big=struct;

% Add first/mid/last structures for C only
C_big_first=struct;
C_big_mid=struct;
C_big_last=struct;

ar=fieldnames(ExcitRR_Ctrl1);
for a = 1:length(ar) 
   % === All trials ===
   C_big.stimvel(:,:,a) = ExcitRR_Ctrl1.(ar{a}).C_CSM.Average;
   T0_big.stimvel(:,:,a) = ExcitRR_Ctrl1.(ar{a}).ST0.Average;
   T1_big.stimvel(:,:,a) = ExcitRR_Ctrl1.(ar{a}).ST1.Average;
   T2_big.stimvel(:,:,a) = ExcitRR_Ctrl1.(ar{a}).ST2.Average;

   % === First 5 ===
   C_big_first.stimvel(:,:,a) = ExcitRR_Ctrl1.(ar{a}).C_CSM_first.Average;
   % === Middle 5 ===
   C_big_mid.stimvel(:,:,a)   = ExcitRR_Ctrl1.(ar{a}).C_CSM_middle.Average;
   % === Last 5 ===
   C_big_last.stimvel(:,:,a)  = ExcitRR_Ctrl1.(ar{a}).C_CSM_last.Average;
end

% ----- Helper for average/std/N -----
calc_big = @(x) struct( ...
    'Average', nanmean(x.stimvel,3), ...
    'Std', nanstd(x.stimvel,0,3), ...
    'N', arrayfun(@(i) sum(~isnan(squeeze(x.stimvel(i,2,:)))), 1:size(x.stimvel,1))' ...
);

% ----- Compute group-level summaries -----
C_big = calc_big(C_big);
C_big_first = calc_big(C_big_first);
C_big_mid = calc_big(C_big_mid);
C_big_last = calc_big(C_big_last);

T0_big = calc_big(T0_big);
T1_big = calc_big(T1_big);
T2_big = calc_big(T2_big);

%% For control reaches last session

Cend_big=struct;
T0end_big=struct;
T1end_big=struct;
T2end_big=struct;

% Add first/mid/last structures for C only
Cend_big_first=struct;
Cend_big_mid=struct;
Cend_big_last=struct;

ar=fieldnames(ExcitRR_Ctrlend);
for a = 1:length(ar) 
   % === All trials ===
   Cend_big.stimvel(:,:,a) = ExcitRR_Ctrlend.(ar{a}).C_CSM.Average;
   T0end_big.stimvel(:,:,a) = ExcitRR_Ctrlend.(ar{a}).ST0.Average;
   T1end_big.stimvel(:,:,a) = ExcitRR_Ctrlend.(ar{a}).ST1.Average;
   T2end_big.stimvel(:,:,a) = ExcitRR_Ctrlend.(ar{a}).ST2.Average;

   % === First 5 ===
   Cend_big_first.stimvel(:,:,a) = ExcitRR_Ctrlend.(ar{a}).C_CSM_first.Average;
   % === Middle 5 ===
   Cend_big_mid.stimvel(:,:,a)   = ExcitRR_Ctrlend.(ar{a}).C_CSM_middle.Average;
   % === Last 5 ===
   Cend_big_last.stimvel(:,:,a)  = ExcitRR_Ctrl.(ar{a}).C_CSM_last.Average;
end

% ----- Helper for average/std/N -----
calc_big = @(x) struct( ...
    'Average', nanmean(x.stimvel,3), ...
    'Std', nanstd(x.stimvel,0,3), ...
    'N', arrayfun(@(i) sum(~isnan(squeeze(x.stimvel(i,2,:)))), 1:size(x.stimvel,1))' ...
);

% ----- Compute group-level summaries -----
C_big = calc_big(C_big);
C_big_first = calc_big(C_big_first);
C_big_mid = calc_big(C_big_mid);
C_big_last = calc_big(C_big_last);

T0_big = calc_big(T0_big);
T1_big = calc_big(T1_big);
T2_big = calc_big(T2_big);




%% For T0 reahces where 1 animal does not have any T0 reaces 
% Add accel column to SN data - not all animals have NoReReach - this is to add accel and to figure out what animals have none
NoT0Stimvel=[];
indices_to_delete = [];
an=fieldnames(ExcitRR_Ctrlend);
for a=1:length(an)
    NRR=ExcitRR_Ctrlend.(an{a}).ST0;
    if isempty(NRR)
        NoT0Stimvel.(an{a}).ST0=[];
        % Add the index of the empty NRR.stimvel to the list
        indices_to_delete(end+1) = a;
    else
   
    NoT0Stimvel.(an{a}).ST0=NRR;
    NoT0Stimvel.(an{a}).Cnum=ExcitRR_Ctrlend.(an{a}).Cnum;
    NoT0Stimvel.(an{a}).RRnum=ExcitRR_Ctrlend.(an{a}).RRnum;
    NoT0Stimvel.(an{a}).RRT0=ExcitRR_Ctrlend.(an{a}).RRT0;
    NoT0Stimvel.(an{a}).rrnum0=ExcitRR_Ctrlend.(an{a}).rrnum0;
    NoT0Stimvel.(an{a}).T0Kin=ExcitRR_Ctrlend.(an{a}).T0Kin;
    NoT0Stimvel.(an{a}).ST0V=ExcitRR_Ctrlend.(an{a}).ST0V;
     NoT0Stimvel.(an{a}).ST0_CSM=ExcitRR_Ctrlend.(an{a}).ST0_CSM;
   
    end

end 

% Delete structures where NRR.stimvel is empty
for i = length(indices_to_delete):-1:1
    NoT0Stimvel = rmfield(NoT0Stimvel, an{indices_to_delete(i)});
end

at=fieldnames(NoT0Stimvel);
for a = 1:length(at) 
 
   T0_bigend.stimvel(:,:,a) = NoT0Stimvel.(at{a}).ST0.Average;
  
end

%% For T3 reahces where 1 animal does not have any T3 reaches in the first session
% Add accel column to SN data - not all animals have NoReReach - this is to add accel and to figure out what animals have none
NoT3Stimvel=[];
indices_to_delete = [];
an=fieldnames(ExcitRR_Ctrl1);
for a=1:length(an)
    NRR=ExcitRR_Ctrl1.(an{a}).ST3;
    if isempty(NRR)
        NoT3Stimvel.(an{a}).ST3=[];
        % Add the index of the empty NRR.stimvel to the list
        indices_to_delete(end+1) = a;
    else
   
    NoT3Stimvel.(an{a}).ST3=NRR;
    NoT3Stimvel.(an{a}).Cnum=ExcitRR_Ctrl1.(an{a}).Cnum;
    NoT3Stimvel.(an{a}).RRnum=ExcitRR_Ctrl1.(an{a}).RRnum;
    NoT3Stimvel.(an{a}).RRT3=ExcitRR_Ctrl1.(an{a}).RRT3;
    NoT3Stimvel.(an{a}).rrnum3=ExcitRR_Ctrl1.(an{a}).rrnum3;
    NoT3Stimvel.(an{a}).T3Kin=ExcitRR_Ctrl1.(an{a}).T3Kin;
    NoT3Stimvel.(an{a}).ST3V=ExcitRR_Ctrl1.(an{a}).ST3V;
    NoT3Stimvel.(an{a}).ST3_CSM=ExcitRR_Ctrl1.(an{a}).ST3_CSM;
   
    end

end 

% Delete structures where NRR.stimvel is empty
for i = length(indices_to_delete):-1:1
    NoT3Stimvel = rmfield(NoT3Stimvel, an{indices_to_delete(i)});
end

at=fieldnames(NoT3Stimvel);
for a = 1:length(at) 
 
   T3_big.stimvel(:,:,a) = NoT3Stimvel.(at{a}).ST3.Average;
  
end
%% use only if an animal doesn't have any reaches in a category 
T0_bigend.Average = nanmean(T0_bigend.stimvel,3);
T0_bigend.Std = nanstd(T0_bigend.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T0_bigend.Average),1); 
for nnum = 1:length(T0_bigend.Average)
    samp = squeeze(T0_bigend.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
T0_bigend.N = N;
%% use only if an animal doesn't have any reaches in a category 
T3_big.Average = nanmean(T3_big.stimvel,3);
T3_big.Std = nanstd(T3_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T3_big.Average),1); 
for nnum = 1:length(T3_big.Average)
    samp = squeeze(T3_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
T3_big.N = N;
  %% Combine Stimvel and look at plots 

T1_bigstim=struct;
T2_bigstim = struct;
S_big=struct;


ar=fieldnames(ExcitRR);
for a = 1:length(ar) 
   T1_bigstim.stimvel(:,:,a) = ExcitRR.(ar{a}).ST1.Average;
   T2_bigstim.stimvel(:,:,a) = ExcitRR.(ar{a}).ST2.Average;
   %S_big.stimvel(:,:,a) = ExcitRR.(ar{a}).S.Average;
   S_big.stimvel(:,:,a) = ExcitRR.(ar{a}).S_CSM.Average;
  
end



S_big.Average = nanmean(S_big.stimvel,3);
S_big.Std = nanstd(S_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(S_big.Average),1); 
for nnum = 1:length(S_big.Average)
    samp = squeeze(S_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
S_big.N = N;



T1_bigstim.Average = nanmean(T1_bigstim.stimvel,3);
T1_bigstim.Std = nanstd(T1_bigstim.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T1_bigstim.Average),1); 
for nnum = 1:length(T1_bigstim.Average)
    samp = squeeze(T1_bigstim.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
T1_bigstim.N = N;


T2_bigstim.Average = nanmean(T2_bigstim.stimvel,3);
T2_bigstim.Std = nanstd(T2_bigstim.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T2_bigstim.Average),1); 
for nnum = 1:length(T2_bigstim.Average)
    samp = squeeze(T2_bigstim.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
T2_bigstim.N = N;
%%

T1_bigstimend=struct;
T2_bigstimend = struct;
S_bigend=struct;


ar=fieldnames(ExcitRR_end);
for a = 1:length(ar) 
   T1_bigstimend.stimvel(:,:,a) = ExcitRR_end.(ar{a}).ST1.Average;
   T2_bigstimend.stimvel(:,:,a) = ExcitRR_end.(ar{a}).ST2.Average;
  % S_bigend.stimvel(:,:,a) = ExcitRR_end.(ar{a}).S.Average;
   S_bigend.stimvel(:,:,a) = ExcitRR_end.(ar{a}).Send_CSM.Average;
  
end



S_bigend.Average = nanmean(S_bigend.stimvel,3);
S_bigend.Std = nanstd(S_bigend.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(S_bigend.Average),1); 
for nnum = 1:length(S_bigend.Average)
    samp = squeeze(S_bigend.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
S_bigend.N = N;



T1_bigstimend.Average = nanmean(T1_bigstimend.stimvel,3);
T1_bigstimend.Std = nanstd(T1_bigstimend.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T1_bigstimend.Average),1); 
for nnum = 1:length(T1_bigstimend.Average)
    samp = squeeze(T1_bigstimend.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
T1_bigstimend.N = N;


T2_bigstimend.Average = nanmean(T2_bigstimend.stimvel,3);
T2_bigstimend.Std = nanstd(T2_bigstimend.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T2_bigstimend.Average),1); 
for nnum = 1:length(T2_bigstimend.Average)
    samp = squeeze(T2_bigstimend.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
T2_bigstimend.N = N;

%% Sham 
Sham_big=struct;

ar=fieldnames(ExcitRR_Sham);
for a = 1:length(ar) 
   %Sham_big.stimvel(:,:,a) = ExcitRR_Sham.(ar{a}).C.Average;
   Sham_big.stimvel(:,:,a) = ExcitRR_Sham.(ar{a}).Sham_CSM.Average;
  
end

Sham_big.Average = nanmean(Sham_big.stimvel,3);
Sham_big.Std = nanstd(Sham_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(Sham_big.Average),1); 
for nnum = 1:length(Sham_big.Average)
    samp = squeeze(Sham_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
Sham_big.N = N;

%% PLot simple stimvel - first, last and sham alongside it 

Ctrlfirst=rgb('Orange');
Ctrllast=rgb('OrangeRed');
Shamcolor=rgb('Gray');

dim =2;

StimvelPlot(C_big,dim,Ctrlfirst,1)
StimvelPlot(C_bigend,dim,Ctrllast,0)
StimvelPlot(Sham_big,dim,Shamcolor,0)


xlim([-0.2 0.2])
%xticks(-0.2:0.01:0.2)
ylim('auto')
%ylim([-0.5 1.5])
%ylim([0 0.25])
%ylim([-20 20])
%pbaspect([1 2 1])
%ylim([0 0.3])
xline(0,'--')
%xline(0.05,'--')
yline(0,'--')
% the different CSM alignment times 
% xline(-0.055,'y')
% xline(0.02,'r')%T1
% xline(0.095,'b')%T2
% xline(0.165,'c')%T3
% xline(0.23,'m')
%xlabel('Time from Pos. Thresh.(s)')
xlabel('Time from CSM(s)')

%ylabel('Outward Position (cm)')
ylabel('Outward Velocity (cm/s)')
%ylabel('Outward acceleration')

%% PLot simple stimvel - first, last and sham alongside it 

Stimfirst=rgb('LightSteelBlue');
Stimlast=rgb('MediumBlue');

dim =8;

StimvelPlot(S_big,dim,Stimfirst,1)
StimvelPlot(S_bigend,dim,Stimlast,0)
StimvelPlot(Sham_big,dim,Shamcolor,0)


xlim([-0.2 0.2])
%xticks(-0.2:0.01:0.2)
ylim('auto')
%ylim([-0.5 1.5])
%ylim([0 0.25])
%ylim([-20 20])
%pbaspect([1 2 1])
%ylim([0 0.3])
xline(0,'--')
xline(0.05,'--')
yline(0,'--')
% the different CSM alignment times 
% xline(-0.055,'y')
% xline(0.02,'r')%T1
% xline(0.095,'b')%T2
% xline(0.165,'c')%T3
% xline(0.23,'m')
xlabel('Time from threshold crossing (s)')

%ylabel('Outward Position (cm)')
ylabel('Outward Velocity (cm/s)')
%ylabel('Outward acceleration')

%%
% PLot all together


Stimfirst=rgb('LightSteelBlue');
Stimlast=rgb('MediumBlue');

dim =2;

StimvelPlot(S_big,dim,Stimfirst,1)
StimvelPlot(S_bigend,dim,Stimlast,0)
StimvelPlot(Sham_big,dim,Shamcolor,0)
StimvelPlot(C_big,dim,Ctrlfirst,0)
StimvelPlot(C_bigend,dim,Ctrllast,0)


xlim([-0.2 0.2])
%xticks(-0.2:0.01:0.2)
ylim('auto')
%ylim([-0.5 1.5])
%ylim([0 0.25])
%ylim([-20 20])
%pbaspect([1 2 1])
%ylim([0 0.3])
xline(0,'--')
% xline(0.05,'--')
yline(0,'--')
% the different CSM alignment times 
% xline(-0.055,'y')
% xline(0.02,'r')%T1
% xline(0.095,'b')%T2
% xline(0.165,'c')%T3
% xline(0.23,'m')
%xlabel('Time from threshold crossing (s)')

xlabel('Time from CSM(s)')

%ylabel('Outward Position (cm)')
ylabel('Outward Velocity (cm/s)')
%ylabel('Outward acceleration')
%% Plot stimvel

%Add acceleration column to stimvel
T0color=rgb('Orange');
T1color=rgb('OrangeRed');
T2color=rgb('Crimson');
T3color=rgb('DarkRed');


dim =2;

StimvelPlot(T0_big,dim,T0color,1)
StimvelPlot(T1_big,dim,T1color,0)
StimvelPlot(T2_big,dim,T2color,0)
StimvelPlot(T3_big,dim,T3color,0)
StimvelPlot(C_big,dim,'k',0)

xlim([-0.2 0.2])
%xticks(-0.2:0.01:0.2)
ylim('auto')
%ylim([-0.5 1.5])
%ylim([0 0.27])
%ylim([-20 20])
pbaspect([1 2 1])
%ylim([0 0.3])
xline(0,'--')
xline(0.05,'--')
yline(0,'--')
% the different CSM alignment times 
% xline(-0.055,'y')
% xline(0.02,'r')%T1
% xline(0.095,'b')%T2
% xline(0.165,'c')%T3
% xline(0.23,'m')
xlabel('Time from Pos. Thresh.(s)')

%ylabel('Outward Position (cm)')
ylabel('Outward Velocity (cm/s)')
%ylabel('Outward acceleration')

%% Plot stimvel

T0colorstim=rgb('LightBlue');
T1colorstim=rgb('SteelBlue');
T2colorstim=rgb('MediumBlue');

stimcolor=rgb('DodgerBlue');
ctrlcolor=rgb('Black');
x1=0 ;x2=0.05;
y_upper=0.3; y_lower=0;

dim =2;

StimvelPlot(T1_bigstim,dim,T1colorstim,1)
StimvelPlot(T2_bigstim,dim,T2colorstim,0)
StimvelPlot(S_big,dim,stimcolor,0)

hold on
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
xlim([-0.2 0.2])
%xticks(-0.2:0.01:0.2)
ylim('auto')
%ylim([-0.5 1.5])
%ylim([0 0.27])
%ylim([-20 20])
pbaspect([1 2 1])
%ylim([0 0.3])
xline(0,'--')
xline(0.05,'--')
yline(0,'--')
% the different CSM alignment times 
% xline(-0.055,'y')
% xline(0.02,'r')%T1
% xline(0.095,'b')%T2
% xline(0.165,'c')%T3
% xline(0.23,'m')
xlabel('Time from Pos. Thresh.(s)')
%ylabel('Outward Position (cm)')
ylabel('Outward Velocity (cm/s)')
%ylabel('Outward acceleration')
%%
T0color=rgb('Thistle');
T1color=rgb('MediumOrchid');
T2color=rgb('SlateBlue');
T3color=rgb('DarkSlateBlue');


dim =2;

StimvelPlot(T0_bigend,dim,T0color,1)
StimvelPlot(T1_bigend,dim,T1color,0)
StimvelPlot(T2_bigend,dim,T2color,0)
StimvelPlot(T3_bigend,dim,T3color,0)
StimvelPlot(C_bigend,dim,'k',0)

xlim([-0.2 0.2])
%xticks(-0.2:0.01:0.2)
ylim('auto')
%ylim([-0.6 1.4])
%ylim([-0.05 0.3])
%ylim([-20 30])
pbaspect([1 2 1])
%ylim([0 0.3])
xline(0,'--')
xline(0.05,'--')
yline(0,'--')
% the different CSM alignment times 
% % xline(-0.055,'y')
% xline(0.02,'r')%T1
% xline(0.095,'b')%T2
% xline(0.165,'c')%T3
% xline(0.23,'m')
xlabel('Time from Pos. Thresh.(s)')

%ylabel('Outward Position (cm)')
%ylabel('Outward Velocity (cm/s)')

ylabel('Outward acceleration')

%% Plot stimvel

T1colorstimend=rgb('LightSkyBlue');
T2colorstimend=rgb('Navy');

stimcolor=rgb('DodgerBlue');
ctrlcolor=rgb('Black');
x1=0 ;x2=0.05;
y_upper=0.3; y_lower=0;

dim =2;

StimvelPlot(T1_bigstimend,dim,T1colorstimend,1)
StimvelPlot(T2_bigstimend,dim,T2colorstimend,0)
StimvelPlot(S_bigend,dim,stimcolor,0)

hold on
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
xlim([-0.2 0.2])
%xticks(-0.2:0.01:0.2)
ylim('auto')
%ylim([-0.5 1.5])
%ylim([0 0.27])
%ylim([-20 20])
pbaspect([1 2 1])
%ylim([0 0.3])
xline(0,'--')
xline(0.05,'--')
yline(0,'--')
% the different CSM alignment times 
% xline(-0.055,'y')
% xline(0.02,'r')%T1
% xline(0.095,'b')%T2
% xline(0.165,'c')%T3
% xline(0.23,'m')
xlabel('Time from Pos. Thresh.(s)')
%ylabel('Outward Position (cm)')
ylabel('Outward Velocity (cm/s)')
%ylabel('Outward acceleration')