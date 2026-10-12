%name your master struct
%navigate to folder with combined reaches in a single structure
% REaches Combined-strong stim 
ChR2 = struct;


fname = uigetdir;
direc = dir(strcat(fname,'\*.mat'));
%%
% filter reaches based off of position 

for i = 1:length(direc)
    %load reaches
    load(fullfile(fname,direc(i).name));
    
    %get name
    fullname = direc(i).name;
    mouse_name = regexprep(fullname,'.mat','');
    
    %save reaches
    thresh=0.1;dirdir=1;
    Reaches=anim_data.Reaches;
    %[Reaches,bad]=rfilts_standard(Reaches,0.16,5,0.8);
    %[Reaches,bad]=rfilts_standard(Reaches,0.1,5,0.8); % for headfixed reaching 
    [Reaches,bad]=rfilts_standard(Reaches,thresh,5,0.8); % new threshold 
    ChR2.(mouse_name).Reaches = Reaches;
    [SucAll,FailAll]=SuccessFilter(Reaches);

    AllRSuc= AlignThresh(SucAll,thresh,dirdir,200);
    AllRSuc=stimvel_accel(AllRSuc);

    AllRFail= AlignThresh(FailAll,thresh,dirdir,200);
    AllRFail=stimvel_accel(AllRFail);

    ChR2.(mouse_name).SucAll = SucAll;
    ChR2.(mouse_name).FailAll = FailAll;
    ChR2.(mouse_name).AllRSuc = AllRSuc;
    ChR2.(mouse_name).AllRFail = AllRFail;
    
    % % % do stimvel
    % get stim good for averaging directly off of stim point instead of thresh point
   % thresh = .16; dirdir = 1;
    %thresh = .1; 
    dirdir = 1;
    [Stim, Control] = StimFilter(Reaches);
    [Stim2] = stimfilt_thresh(Stim,thresh);
    [Stim3] = stimfilt_firstreach(Stim2, thresh, dirdir);
    [SucStim,FailStim]=SuccessFilter(Stim3);

    S = AlignStim(Stim3,200);
    S=stimvel_accel(S);

    S_CSM=AlignCSM(Stim3,500);
    S_CSM=stimvel_accel(S_CSM);


    Sv=AlignPeakV(Stim3,7,500);
    %S=stimvel_accel_jerk(S); %adding jerk column
    Ssuc=AlignStim(SucStim,200);
    Ssuc=stimvel_accel(Ssuc);
    Sfail=AlignStim(FailStim,200);
    Sfail=stimvel_accel(Sfail);

    [CtrlSucc,CtrlFaill]=SuccessFilter(Control);
    C = AlignThresh(Control,thresh,dirdir,200);
    C=stimvel_accel(C);
    Csuc = AlignThresh(CtrlSucc,thresh,dirdir,200);
    Csuc=stimvel_accel(Csuc);
    Cfail = AlignThresh(CtrlFaill,thresh,dirdir,200);
    Cfail=stimvel_accel(Cfail);

    C_CSM=AlignThreshCSM(Control,thresh,dirdir,500);
    C_CSM=stimvel_accel(C_CSM);

    Cv=AlignPeakV(Control,7,500);
    %C=stimvel_accel_jerk(C); %adding jerk column
    
    ChR2.(mouse_name).C = C;
    ChR2.(mouse_name).Csuc = Csuc;
    ChR2.(mouse_name).Cfail = Cfail;
    ChR2.(mouse_name).S = S;
    ChR2.(mouse_name).Ssuc = Ssuc;
    ChR2.(mouse_name).Sfail = Sfail;
    ChR2.(mouse_name).Cv = Cv;
    ChR2.(mouse_name).Sv = Sv;
    ChR2.(mouse_name).C_CSM = C_CSM;
    ChR2.(mouse_name).S_CSM = S_CSM;
    
    %get latency of max euclid distance in velocity (dim 0)
    [~,~,latency] = StimvelMag(C,S,0);
    %get magnitude of difference in each dimension, and sd
    [magnitude,deviation,p,z] = StimvelLatencyMag(C,S,latency);
    
    ChR2.(mouse_name).stimvel.latency = latency;
    ChR2.(mouse_name).stimvel.magnitude = magnitude;
    ChR2.(mouse_name).stimvel.sem = deviation;
    ChR2.(mouse_name).stimvel.p = p;
    ChR2.(mouse_name).stimvel.z = z;
    
    %get difference between control and stim, and integrate
    %this does it over the 100 ms post stim (including stim)
    SC_Difb = S.Average(:,2:end) - C.Average(:,2:end);
    SC_Dif = [C.Average(:,1) SC_Difb];
    areadiff = trapz(SC_Dif(21:31,1),abs(SC_Dif(21:31,5:7)));
    
    ChR2.(mouse_name).stimvel.areadiff = areadiff;
    
    % % % get basic kinematic features
    %Controlclipped = localmaxclip(Control); % not using for analysis looking at reacceleration
    %Controlclipped = shortfilter(Controlclipped,6); % using
    %rfilts_standard above
    %Controlclipped=shortfilter(Control,5);
    
    %ControlClipped=accel_toreaches(Control);% add acceleration column and jerk column
    ControlClipped=Control; 

    [CtrlSucc,CtrlFaill]=SuccessFilter(ControlClipped);
%      Stimclipped = localmaxclip(Stim3);
%      Stimclipped = shortfilter(Stimclipped,6);
    
    %Stimclipped=shortfilter(Stim3,5);

    %StimClipped=accel_toreaches(Stim3); % headfixed already has the accel category
    StimClipped=Stim3;
    [StimSucc,StimFaill]=SuccessFilter(StimClipped);
    
    % Save CTRL and STIM 
    ChR2.(mouse_name).CtrlReaches =ControlClipped;
    ChR2.(mouse_name).CtrlSuc =CtrlSucc;
    ChR2.(mouse_name).CtrlFail =CtrlFaill;
    ChR2.(mouse_name).StimReaches =StimClipped;
    ChR2.(mouse_name).StimSuc=StimSucc;
    ChR2.(mouse_name).StimFail =StimFaill;
    
    
    
    %Save kinematic variables 
    ChR2.(mouse_name).endpoint.control = endpoint(ControlClipped);
    ChR2.(mouse_name).endpoint.stim = endpoint(StimClipped);
    ChR2.(mouse_name).txmax.control = txmax(ControlClipped,thresh);
    ChR2.(mouse_name).txmax.stim = txmax(StimClipped,thresh);
    %Get Basic Kinematics
    Ctrlkin= BasicKinematics(ControlClipped);
    Ctrlendpt=mean(Ctrlkin.endpt);  
    Ctrlvmaxp=mean(Ctrlkin.vmax_p);
    ChR2.(mouse_name).CtrlKin = Ctrlkin;
    ChR2.(mouse_name).CtrlEndpt = Ctrlendpt;
    ChR2.(mouse_name).CtrlVmaxp = Ctrlvmaxp;
    %Success Control
    CtrlSkin= BasicKinematics(CtrlSucc);
    CtrlSendpt=mean(CtrlSkin.endpt);  
    ChR2.(mouse_name).CtrlSKin = CtrlSkin;
    ChR2.(mouse_name).CtrlSEndpt = CtrlSendpt;
    %Fail Control
    CtrlFkin= BasicKinematics(CtrlFaill);
    CtrlFendpt=mean(CtrlFkin.endpt);  
    ChR2.(mouse_name).CtrlFKin = CtrlFkin;
    ChR2.(mouse_name).CtrlFEndpt = CtrlFendpt;
    % All stim
    Stimkin=BasicKinematics(StimClipped);
    Stimendpt=mean(Stimkin.endpt);
    Stimvmaxp=mean(Stimkin.vmax_p);
    ChR2.(mouse_name).StimKin = Stimkin;
    ChR2.(mouse_name).StimEndpt = Stimendpt;
    ChR2.(mouse_name).StimVmaxp = Stimvmaxp;
    % Success Stim
    StimSkin=BasicKinematics(StimSucc);
    StimSendpt=mean(StimSkin.endpt);
    ChR2.(mouse_name).StimSKin = StimSkin;
    ChR2.(mouse_name).StimSEndpt = StimSendpt;

    % Fail Stim
    StimFkin=BasicKinematics(StimFaill);
    StimFendpt=mean(StimFkin.endpt);
    ChR2.(mouse_name).StimFKin = StimFkin;
    ChR2.(mouse_name).StimFEndpt = StimFendpt;

    % Sort based off of change in distance during stim 
    % if using 250, stim start=26, stim end =31
% if using 500, star=51, end =56
    Sreaches=S.stimvel;
    fields=fieldnames(StimClipped);
        
    % StimEndPos=[];
    % StimStartPos=[];
    % for i=1:length(fields)
    %     endpos=Sreaches(56,2,i);
    %     startpos=Sreaches(51,2,i); % stim start
    %     StimEndPos=[StimEndPos;endpos];
    %     StimStartPos=[StimStartPos;startpos];
    % end
    %  Change in position
    % StimChange=(StimEndPos-StimStartPos)*10;
    % VmaxStim=StimChange;

    % now sort reaches based off of velocity percentiles 
   % VmaxCtrl=Ctrlkin.vmax(:,2);
    VmaxCtrl=Ctrlkin.vmax(:,2);
    
    p25=prctile(VmaxCtrl,25);
    p50=prctile(VmaxCtrl,50);
    p75=prctile(VmaxCtrl,75);
    ChR2.(mouse_name).VmaxPrct.p25=p25;
    ChR2.(mouse_name).VmaxPrct.p50=p50;
    ChR2.(mouse_name).VmaxPrct.p75=p75;
    %find indices for elements
    below_p25_indices = find(VmaxCtrl <= p25);
    between_p25_p50_indices = find(VmaxCtrl > p25 & VmaxCtrl <= p50);
    between_p50_p75_indices = find(VmaxCtrl > p50 & VmaxCtrl <= p75);
    above_p75_indices = find(VmaxCtrl > p75);
    

    % Get field names from the structure
    Stimfields = fieldnames(ControlClipped);
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
            Stim25.(fieldName) = Extract(ControlClipped.(fieldName));
        end
    end
      
     for p = 1:length(between_p25_p50_indices)
        index50 = between_p25_p50_indices(p);

        if index50 > 0
            fieldName50 = stim50{p};
            Stim50.(fieldName50) = Extract(ControlClipped.(fieldName50));
        end
      end
    
       for w = 1:length(between_p50_p75_indices)
        index75 = between_p50_p75_indices(w);

        if index75 > 0
            fieldName75 = stim75{w};
            Stim75.(fieldName75) = Extract(ControlClipped.(fieldName75));
        end
       end
       
     for c = 1:length(above_p75_indices)
        index100 = above_p75_indices(c);

        if index100 > 0
            fieldName100 = stim100{c};
            Stim100.(fieldName100) = Extract(ControlClipped.(fieldName100));
        end
     end
   
    % Now save each structure 
    ChR2.(mouse_name).Stim25Reaches=Stim25;
    ChR2.(mouse_name).Stim50Reaches=Stim50;
    ChR2.(mouse_name).Stim75Reaches=Stim75;
    ChR2.(mouse_name).Stim100Reaches=Stim100;
    %make stimvel for each subset of reaches
     S25= AlignThresh(Stim25,thresh,dirdir,500);
    S50 = AlignThresh(Stim50,thresh,dirdir,500);
    S75 = AlignThresh(Stim75,thresh,dirdir,500);
    S100 = AlignThresh(Stim100,thresh,dirdir,500);
    %save stimvel
    ChR2.(mouse_name).S25 = S25;
    ChR2.(mouse_name).S50 = S50;
    ChR2.(mouse_name).S75 = S75;
    ChR2.(mouse_name).S100 = S100;
    
    % Separate out success and control 
    [S_C, F_C] = SuccessFilter(Control);
    [S_S, F_S] = SuccessFilter(Stim);
    
    ChR2.(mouse_name).successrate.control = length(fieldnames(S_C)) / length(fieldnames(Control));
    ChR2.(mouse_name).successrate.stim = length(fieldnames(S_S)) / length(fieldnames(Stim));
    
end
%% Get diff in stimvel 
an=fieldnames(ChR2);
stimveldiffX=[];stimveldiffY=[];
for a=1:length(an)
    CavgX=ChR2.(an{a}).C.Average(26,5);
    SavgX=ChR2.(an{a}).S.Average(26,5);

    CavgY=ChR2.(an{a}).C.Average(26,6);
    SavgY=ChR2.(an{a}).S.Average(26,6);

    DiffX=SavgX-CavgX;
    DiffY=SavgY-CavgY;
    stimveldiffX(a,1)=DiffX;
    stimveldiffY(a,1)=DiffY;
end 

ExcitAnimalsDiff=[stimveldiffX,stimveldiffY];





%%
%get basic kinematic parameters of control and stim reaches
an=fieldnames(ChR2);
%if you want the data at each animal's point of max euclid distance
for a = 1:length(an)
   
    chr2_endpt_c(a,:) = mean(ChR2.(an{a}).endpoint.control);
    chr2_endpt_s(a,:) = mean(ChR2.(an{a}).endpoint.stim);
    
    chr2_endpt_p(a) = ranksum(ChR2.(an{a}).endpoint.control(:,1),ChR2.(an{a}).endpoint.stim(:,1),'method','approximate');

    chr2_txmax_c(a,:) = mean(ChR2.(an{a}).txmax.control); 
    chr2_txmax_s(a,:) = mean(ChR2.(an{a}).txmax.stim);
    chr2_txmax_p(a) = ranksum(ChR2.(an{a}).txmax.control(:,1),ChR2.(an{a}).txmax.stim(:,1),'method','approximate');

    chr2_endpt_dist(a,1) = sqrt(sum(chr2_endpt_c(a,1:2).^2,2));
    chr2_endpt_dist(a,2) = sqrt(sum(chr2_endpt_s(a,1:2).^2,2));
    
    
    %get sds
    chr2_endpt_std_c(a,:) = std(ChR2.(an{a}).endpoint.control);
    chr2_endpt_std_s(a,:) = std(ChR2.(an{a}).endpoint.stim);
    
    endpts = [ChR2.(an{a}).endpoint.control(:,1) ; ChR2.(an{a}).endpoint.stim(:,1)];
    endpt_grps = ones(length(endpts),1);
    endpt_grps(length(ChR2.(an{a}).endpoint.control(:,1))+1:end) = 2;
    chr2_endpt_std_p(a) = vartestn(endpts,endpt_grps,'TestType','LeveneAbsolute','Display','off');
end

ranksum(chr2_endpt_std_c(:,1),chr2_endpt_std_s(:,1))

%%
%get stimvel values t ms after stim

t = .05;
dim = 2;

%replace first column of a with time
timediff = ChR2.(an{a}).C.Average(:,1) - t;
[~,mininx] = min(abs(timediff));

for a = 1:length(an)
   
    chr2_c(a,:) = nanmean(ChR2.(an{a}).C.stimvel(mininx,dim,:));
    chr2_s(a,:) = nanmean(ChR2.(an{a}).S.stimvel(mininx,dim,:));
    chr2_cs_p(a,:) = ranksum(squeeze(ChR2.(an{a}).C.stimvel(mininx,dim,:)),squeeze(ChR2.(an{a}).S.stimvel(mininx,dim,:)),'method','approximate');
end



%%
%get distance traveled from t to t2
t = .00;
t2 = .07;

dim = 2;

%replace first column of a with time
timediff = ChR2.(an{a}).C.Average(:,1) - t;
[~,mininx] = min(abs(timediff));

timediff2 = ChR2.(an{a}).C.Average(:,1) - t2;
[~,mininx2] = min(abs(timediff2));

for a = 1:length(an)
    stimvel_c = ChR2.(an{a}).C.stimvel;
    chr2_dist_c_all = squeeze(stimvel_c(mininx2,dim,:) - stimvel_c(mininx,dim,:));
    chr2_dist_c(a,1) = nanmean(chr2_dist_c_all);
    
    stimvel_s = ChR2.(an{a}).S.stimvel;
    chr2_dist_s_all = squeeze(stimvel_s(mininx2,dim,:) - stimvel_s(mininx,dim,:));
    chr2_dist_s(a,1) = nanmean(chr2_dist_s_all);
    
    chr2_dist_p(a) = ranksum(squeeze(chr2_dist_c_all),squeeze(chr2_dist_s_all),'method','approximate');
end
%%
%loop through distance traveled
%get distance traveled from t to t2
t = .00;
t2s = .01:.01:.2;
for aaa = 1:length(t2s)

dim = 2;
t2 = t2s(aaa);
%replace first column of a with time
timediff = ChR2.(an{a}).C.Average(:,1) - t;
[~,mininx] = min(abs(timediff));

timediff2 = ChR2.(an{a}).C.Average(:,1) - t2;
[~,mininx2] = min(abs(timediff2));

for a = 1:length(an)
    stimvel_c = ChR2.(an{a}).C.stimvel;
    chr2_dist_c_all = squeeze(stimvel_c(mininx2,dim,:) - stimvel_c(mininx,dim,:));
    chr2_dist_c(a,aaa) = nanmean(chr2_dist_c_all);
    
    stimvel_s = ChR2.(an{a}).S.stimvel;
    chr2_dist_s_all = squeeze(stimvel_s(mininx2,dim,:) - stimvel_s(mininx,dim,:));
    chr2_dist_s(a,aaa) = nanmean(chr2_dist_s_all);
    
    chr2_dist_p(a,aaa) = ranksum(squeeze(chr2_dist_c_all),squeeze(chr2_dist_s_all),'method','approximate');
end

end

chr2_dist_cmean = mean(chr2_dist_c);
chr2_dist_smean = mean(chr2_dist_s);
chr2_dist_cstd = std(chr2_dist_c);
chr2_dist_sstd= std(chr2_dist_s);
figure; hold on
plot(t2s,chr2_dist_cmean,'Color',[.5 .5 .5],'LineWidth',2)
plot(t2s,chr2_dist_smean,'Color',[.2 .5 .9],'LineWidth',2)
shadedErrorBars(t2s,chr2_dist_cmean,chr2_dist_cstd./sqrt(9),[.5 .5 .5])
shadedErrorBars(t2s,chr2_dist_smean,chr2_dist_sstd./sqrt(9),[.2 .5 .9])



%%
%get average combined stimvel plot for 'out significant' animals

C_big = struct;
Csuc_big = struct;
Cfail_big = struct;
S_big = struct;
Cv_big = struct;
Sv_big = struct;
C_CSM_big = struct;
S_CSM_big = struct;
Ssuc_big=struct;
Sfail_big=struct;
an=fieldnames(ChR2);
for a = 1:length(an)
   
    C_big.stimvel(:,:,a) = ChR2.(an{a}).C.Average;
    Csuc_big.stimvel(:,:,a) = ChR2.(an{a}).Csuc.Average;
    Cfail_big.stimvel(:,:,a) = ChR2.(an{a}).Cfail.Average;
    S_big.stimvel(:,:,a) = ChR2.(an{a}).S.Average;
    Cv_big.stimvel(:,:,a) = ChR2.(an{a}).Cv.Average;
    Sv_big.stimvel(:,:,a) = ChR2.(an{a}).Sv.Average;
    C_CSM_big.stimvel(:,:,a) = ChR2.(an{a}).C_CSM.Average;
    S_CSM_big.stimvel(:,:,a) = ChR2.(an{a}).S_CSM.Average;
    Ssuc_big.stimvel(:,:,a) = ChR2.(an{a}).Ssuc.Average;
    Sfail_big.stimvel(:,:,a) = ChR2.(an{a}).Sfail.Average;
    
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

Csuc_big.Average = nanmean(Csuc_big.stimvel,3);
Csuc_big.Std = nanstd(Csuc_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(Csuc_big.Average),1); 
for nnum = 1:length(Csuc_big.Average)
    samp = squeeze(Csuc_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
Csuc_big.N = N;

Cfail_big.Average = nanmean(Cfail_big.stimvel,3);
Cfail_big.Std = nanstd(Cfail_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(Cfail_big.Average),1); 
for nnum = 1:length(Cfail_big.Average)
    samp = squeeze(Cfail_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
Cfail_big.N = N;



S_big.Average = nanmean(S_big.stimvel,3);
S_big.Std = nanstd(S_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(C_big.Average),1); 
for nnum = 1:length(C_big.Average)
    samp = squeeze(S_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
S_big.N = N;

%For structures aligned to Velocity 
Cv_big.Average = nanmean(Cv_big.stimvel,3);
Cv_big.Std = nanstd(Cv_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(Cv_big.Average),1); 
for nnum = 1:length(Cv_big.Average)
    samp = squeeze(Cv_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
Cv_big.N = N;



Sv_big.Average = nanmean(Sv_big.stimvel,3);
Sv_big.Std = nanstd(Sv_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(Cv_big.Average),1); 
for nnum = 1:length(Cv_big.Average)
    samp = squeeze(Sv_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
Sv_big.N = N;

%For structures aligned to CSM
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

% %Succes Stim
Ssuc_big.Average = nanmean(Ssuc_big.stimvel,3);
Ssuc_big.Std = nanstd(Ssuc_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(C_big.Average),1); 
for nnum = 1:length(C_big.Average)
    samp = squeeze(Ssuc_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
Ssuc_big.N = N;

%fail
Sfail_big.Average = nanmean(Sfail_big.stimvel,3);
Sfail_big.Std = nanstd(Sfail_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(C_big.Average),1); 
for nnum = 1:length(C_big.Average)
    samp = squeeze(Sfail_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
Sfail_big.N = N;


%%
%stimcolor= rgb('DodgerBlue');
stimcolor= rgb('SteelBlue');
linecolor=rgb('DarkGray');
stimScolor=rgb('Teal');
stimFcolor=rgb('DarkRed');
ctrlScolor=rgb('Black');
ctrlFcolor=rgb('OrangeRed');

dim = 2;
x1=0 ;x2=0.05;
y_upper=0.2; y_lower=0;
StimvelPlot(C_big,dim,'k',1)
% StimvelPlot(Csuc_big,dim,'k',1)
% StimvelPlot(Cfail_big,dim, ctrlFcolor,0)
StimvelPlot(S_big,dim,stimcolor,0)
%StimvelPlot(S_big,11,stimcolor,0)
% StimvelPlot(Ssuc_big,dim,stimScolor,0)
% StimvelPlot(Sfail_big,dim,stimFcolor,0)
hold on 
% Fill the area between the vertical lines
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');

yline(0,'--k')
%xline(0,'--k')
xlim([-0.2 0.2])
ylim([0 0.2])
%ylim([-40 30])
%ylim('auto')
%yline(0,'--')
xlabel('Time from Stim (s)')
%ylabel('Outward Velocity (cm/s)')
ylabel('Outward Position (cm)')
%pbaspect([ 1 2 1])

% figure; hold on
% StimvelPlot2(C_big,dim,[.5 .5 .5],0)
% shadedErrorBars(C_big.Average(:,1),C_big.Average(:,dim),C_big.Std(:,dim)./sqrt(C_big.N),[.5 .5 .5])
% StimvelPlot2(S_big,dim,[.2 .6 .9],0)
% shadedErrorBars(S_big.Average(:,1),S_big.Average(:,dim),S_big.Std(:,dim)./sqrt(S_big.N),[.2 .6 .9])
% xlim([-0.2 0.2])
% %ylim([-0.5 0.5])
% %ylim([0 .1])
% %yline(0,'--')
% xlabel('Time from Stim (s)')
% ylabel('Outward Velocity (cm/s)')
% %pbaspect([ 1 2 1])
%%
dim = 5;

StimvelPlot(C,dim,'k',1)
StimvelPlot(S,dim,stimcolor,0)
%% PLot CSM
stimcolor= rgb('DodgerBlue');
linecolor=rgb('DarkGray');

dim = 2;

StimvelPlot(C_CSM_big,dim,'k',1)
StimvelPlot(S_CSM_big,dim,stimcolor,0)

yline(0,'--k')
xlim([-0.2 0.2])
%ylim([-0.5 1.5])
%ylim([-15 10])
ylim('auto')
%yline(0,'--')
xlabel('Time from CSM(s)')
ylabel('Outward Velocity (cm/s)')
 pbaspect([ 1 2 1])
%% Plot position only 
figure;

for j=1:5
    Cposx=C_big.stimvel(:,3,j);
    Sposx=S_big.stimvel(:,3,j);
    Cposy=C_big.stimvel(:,2,j);
    Sposy=S_big.stimvel(:,2,j);
    x=[0.05 0.25];
    % plot(Cpos,(1:numel(Cpos)),'k')
    % hold on
    % plot(Spos,(1:numel(Cpos)),'b')
    % hold on
    plot(Cposx,Cposy,'k')
    hold on
    plot(Sposx,Sposy,'b')
    hold on
    
end
%xlim([5 40])
pbaspect([ 1 2 1])

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

dim = 2;
x1=0 ;x2=0.05;
%y_upper=20; y_lower=-26;
StimvelPlot(C_big,dim,'k',1)
StimvelPlot(S25_big,dim,stimcolor25,0)
StimvelPlot(S50_big,dim,stimcolor50,0)
StimvelPlot(S75_big,dim,stimcolor75,0)
StimvelPlot(S100_big,dim,stimcolor100,0)
hold on 
% Fill the area between the vertical lines
%fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');

yline(0,'--k')
xlim([-0.2 0.5])
%ylim([-20 20])
ylim([0 .25])
%ylim('auto')
xlabel('Time from Stim (s)')
%ylabel('Outward Velocity (cm/s)')
ylabel('Outward Position (cm)')
%ylabel('Outward Acceleration (cm/s^2)')
pbaspect([ 1 2 1])

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
%pbaspect([ 1 2 1])

%% PLot single animal at different velocities 
a = 2;


C=ChR2.(an{a}).C;C=stimvel_accel(C);
S=ChR2.(an{a}).S;S=stimvel_accel(S);
S25 = ChR2.(an{a}).S25;S25=stimvel_accel(S25);
S50 = ChR2.(an{a}).S50;S50=stimvel_accel(S50);
S75 = ChR2.(an{a}).S75;S75=stimvel_accel(S75);
S100 = ChR2.(an{a}).S100;S100=stimvel_accel(S100);

dim = 8;
x1=0 ;x2=0.05;
%y_upper=0.25; y_lower=0;
y_upper=20; y_lower=-20;
StimvelPlot(C,dim,'k',1)
StimvelPlot(S,dim,stimcolor25,0)
% StimvelPlot(S25,dim,stimcolor25,0)
% StimvelPlot(S50,dim,stimcolor50,0)
% StimvelPlot(S75,dim,stimcolor75,0)
% StimvelPlot(S100,dim,stimcolor100,0)
hold on 
% Fill the area between the vertical lines
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');

yline(0,'--k')
%ylim([-1 1.8])
xlim([-0.25 0.25])
%ylim([0 .3])
ylim('auto')
xlabel('Time from Stim (s)')
%ylabel('Outward Velocity (cm/s)')
%ylabel('Outward Position (cm)')
ylabel('Outward Acceleration (cm/s^2)')
%pbaspect([ 1 2 1])


%% Use ReReach Filter -All Successful Reaches 
ExcitRR_AllSuc=[]; % Savethis structure to use for graphing rereach timing in sep. script
thresh=0.1;

an=fieldnames(ChR2);

for a=1:length(an)

    Creaches=ChR2.(an{a}).SucAll;
    %Creaches=ChR2.(an{a}).CtrlSuc;
    %Creaches=ChR2.(an{a}).CtrlFail;

 
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
    
    ExcitRR_AllSuc.(an{a}).successrate.ReReach = length(fieldnames(S_C)) / length(fieldnames(Creaches));
    ExcitRR_AllSuc.(an{a}).successrate.NoReReach = length(fieldnames(S_S)) / length(fieldnames(Creaches));
    ExcitRR_AllSuc.(an{a}).CtrlReaches=Creaches;
    ExcitRR_AllSuc.(an{a}).C=C;
    ExcitRR_AllSuc.(an{a}).ReReaches=ReReaches;
    %ExcitRR_AllSuc.(an{a}).PreThreshReReaches=PreThreshReReach;
    ExcitRR_AllSuc.(an{a}).RRKin=RRKin;
    ExcitRR_AllSuc.(an{a}).NoReReach=NoReach;
    ExcitRR_AllSuc.(an{a}).NRKin=NRKin;
    ExcitRR_AllSuc.(an{a}).CN=CN;
    ExcitRR_AllSuc.(an{a}).StimRR=StimRR;
    ExcitRR_AllSuc.(an{a}).MaxVel=MaxVel;
    % get % of reaches with a rereach
    cnum=length(fields1);
    rrnum=length(fields2);
    percentrr=rrnum/cnum;
    ExcitRR_AllSuc.(an{a}).Cnum=cnum;
    ExcitRR_AllSuc.(an{a}).RRnum=rrnum;
    ExcitRR_AllSuc.(an{a}).PercentRR=percentrr;
  
end 

%% Use ReReach Filter -All Successful Reaches 
ExcitRR_AllFail=[]; % Savethis structure to use for graphing rereach timing in sep. script
thresh=0.1;

an=fieldnames(ChR2);

for a=1:length(an)

    Creaches=ChR2.(an{a}).FailAll;
    %Creaches=ChR2.(an{a}).CtrlSuc;
    %Creaches=ChR2.(an{a}).CtrlFail;

 
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
    
    ExcitRR_AllFail.(an{a}).successrate.ReReach = length(fieldnames(S_C)) / length(fieldnames(Creaches));
    ExcitRR_AllFail.(an{a}).successrate.NoReReach = length(fieldnames(S_S)) / length(fieldnames(Creaches));
    ExcitRR_AllFail.(an{a}).CtrlReaches=Creaches;
    ExcitRR_AllFail.(an{a}).C=C;
    ExcitRR_AllFail.(an{a}).ReReaches=ReReaches;
    ExcitRR_AllFail.(an{a}).RRKin=RRKin;
    ExcitRR_AllFail.(an{a}).NoReReach=NoReach;
    ExcitRR_AllFail.(an{a}).NRKin=NRKin;
    ExcitRR_AllFail.(an{a}).CN=CN;
    ExcitRR_AllFail.(an{a}).StimRR=StimRR;
    ExcitRR_AllFail.(an{a}).MaxVel=MaxVel;
    % get % of reaches with a rereach
    cnum=length(fields1);
    rrnum=length(fields2);
    percentrr=rrnum/cnum;
    ExcitRR_AllFail.(an{a}).Cnum=cnum;
    ExcitRR_AllFail.(an{a}).RRnum=rrnum;
    ExcitRR_AllFail.(an{a}).PercentRR=percentrr;
  
end 
%% Use ReReach Filter - Control Reaches 
ExcitRR_Ctrl=[]; % Savethis structure to use for graphing rereach timing in sep. script
thresh=0.1;


for a=1:length(an)

    Creaches=ChR2.(an{a}).CtrlReaches;
    %Creaches=ChR2.(an{a}).CtrlSuc;
    %Creaches=ChR2.(an{a}).CtrlFail;

 
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
%% Use ReReach Filter - Control Reaches 
ExcitRR_CtrlFail=[]; % Savethis structure to use for graphing rereach timing in sep. script
%thresh=0.1;
for a=1:length(an)

   
    Creaches=ChR2.(an{a}).CtrlFail;

 
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
    
    ExcitRR_CtrlFail.(an{a}).successrate.ReReach = length(fieldnames(S_C)) / length(fieldnames(Creaches));
    ExcitRR_CtrlFail.(an{a}).successrate.NoReReach = length(fieldnames(S_S)) / length(fieldnames(Creaches));
    ExcitRR_CtrlFail.(an{a}).CtrlReaches=Creaches;
    ExcitRR_CtrlFail.(an{a}).C=C;
    ExcitRR_CtrlFail.(an{a}).ReReaches=ReReaches;
    ExcitRR_CtrlFail.(an{a}).RRKin=RRKin;
    ExcitRR_CtrlFail.(an{a}).NoReReach=NoReach;
    ExcitRR_CtrlFail.(an{a}).NRKin=NRKin;
    ExcitRR_CtrlFail.(an{a}).CN=CN;
    ExcitRR_CtrlFail.(an{a}).StimRR=StimRR;
    ExcitRR_CtrlFail.(an{a}).MaxVel=MaxVel;
    % get % of reaches with a rereach
    cnum=length(fields1);
    rrnum=length(fields2);

    percentrr=rrnum/cnum;
    ExcitRR_CtrlFail.(an{a}).Cnum=cnum;
    ExcitRR_CtrlFail.(an{a}).RRnum=rrnum;
    ExcitRR_CtrlFail.(an{a}).PercentRR=percentrr;
  
end 
%% Use ReReach Filter - Control Reaches -Success 
ExcitRR_CtrlSuc=[]; % Savethis structure to use for graphing rereach timing in sep. script
thresh=0.1;
for a=1:length(an)

  
    Creaches=ChR2.(an{a}).CtrlSuc;

    %Creaches=ChR2.(an{a}).CtrlFail;

 
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
    
    ExcitRR_CtrlSuc.(an{a}).successrate.ReReach = length(fieldnames(S_C)) / length(fieldnames(Creaches));
    ExcitRR_CtrlSuc.(an{a}).successrate.NoReReach = length(fieldnames(S_S)) / length(fieldnames(Creaches));
    ExcitRR_CtrlSuc.(an{a}).CtrlReaches=Creaches;
    ExcitRR_CtrlSuc.(an{a}).C=C;
    ExcitRR_CtrlSuc.(an{a}).ReReaches=ReReaches;
    ExcitRR_CtrlSuc.(an{a}).RRKin=RRKin;
    ExcitRR_CtrlSuc.(an{a}).NoReReach=NoReach;
    ExcitRR_CtrlSuc.(an{a}).NRKin=NRKin;
    ExcitRR_CtrlSuc.(an{a}).CN=CN;
    ExcitRR_CtrlSuc.(an{a}).StimRR=StimRR;
    ExcitRR_CtrlSuc.(an{a}).MaxVel=MaxVel;
    % get % of reaches with a rereach
    cnum=length(fields1);
    rrnum=length(fields2);
    percentrr=rrnum/cnum;
    ExcitRR_CtrlSuc.(an{a}).Cnum=cnum;
    ExcitRR_CtrlSuc.(an{a}).RRnum=rrnum;
    ExcitRR_CtrlSuc.(an{a}).PercentRR=percentrr;
  
end 
%% Use ReReach Filter- STIM reaches 
ExcitRR=[]; % Savethis structure to use for graphing rereach timing in sep. script
for a=1:length(an)
    Sreaches=ChR2.(an{a}).StimReaches;
  

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
       % SN=stimvel_accel(SN); %causing issues in the code 
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
%% Get position after stim off for all stim reaches 

for a=1:length(an)
    Sreaches=ExcitRR.(an{a}).S.stimvel;
    fields=fieldnames(ExcitRR.(an{a}).ReReaches);
    
    StimEndPos=[];
    for i=1:length(fields)
        %endpos=Sreaches(46,2,i);
        endpos=Sreaches(41,2,i); % stim start
        StimEndPos=[StimEndPos;endpos];
    end
    ExcitRR.(an{a}).StimEndPos=StimEndPos;
end

%% Plot CSM timing and EndPos
%colors=brewermap(11,'Set3')
colors=colormap(slanCM('Cool',2))
%folder_path="C:\Users\Person Lab\OneDrive - The University of Colorado Denver\Documents\Summary Data\Plots\Random stim\AllNeuron_MattandCourtney Data_11\RRTiming_IndividualAnimals\CSMvsEndStimPos";
an=fieldnames(ExcitRR);
figure; hold on
for i=1:length(an)
    csmtime=ExcitRR.(an{i}).StimRR(:,1);
    endpos=ExcitRR.(an{i}).StimEndPos;

    %figure; hold on
    scatter(csmtime,endpos,20,'filled','MarkerFaceColor',colors(i,:),'MarkerEdgeColor','k')
    % hold on
    % xline(0.054,'k')
    % xline(0.0864,'k')
    % xline(0.0972,'b')
    % xline(0.1296,'b')
    % xline(0.1296,'r')
    % xline(0.162,'r')
    % 
    % filename=sprintf('CSMvsStimPos%d.tif',i);
    % full_path=fullfile(folder_path,filename);
    % saveas(gcf,full_path);
    % 
    % close;

   
end
xlim([0.03 0.25])
ylim([0.1 0.26])
xline(0.054,'k','LineWidth',1.5)
xline(0.0864,'k','LineWidth',1.5)
xline(0.0972,'b','LineWidth',1.5)
xline(0.1296,'b--','LineWidth',2)
xline(0.1296,'r','LineWidth',1.5)
xline(0.162,'r','LineWidth',1.5)
yline(0.16,'k','LineWidth',1.5)
xlabel('CSM time ')
ylabel('Stim End Position')
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


%% Get Stimvel plots for rereaches
SRR_big = struct;


ar=fieldnames(ExcitRR);
for a = 1:length(ar) 
    SRR_big.stimvel(:,:,a) = ExcitRR.(ar{a}).S.Average;
  
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

%% Plot rereaches only 
stimcolor= rgb('DodgerBlue');
norereach=rgb('PowderBlue');
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
%ylim([-0.5 1.5])
%ylim([-15 10])
ylim('auto')
%yline(0,'--')
xlabel('Time from Stim (s)')
ylabel('Outward Velocity (cm/s)')
pbaspect([ 1 2 1])
%% Get mangnitude of effect for RR 
an=fieldnames(ChR2);
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

ExcitAnimalsDiffRR=[RRstimveldiffX,RRstimveldiffY];
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

ExcitAnimalsDiffNRR=[RRstimveldiffX,RRstimveldiffY];

%% Looks at timing of Stim for Ctrl and Stim reaches 
C_start=[]; C_end=[];
S_start=[]; S_end=[];


for a=1:length(an)
    %T2
    C_stimvel=ExcitRR.(an{a}).C.Average;
    c_start=C_stimvel(41,2);
    c_end=C_stimvel(46,2);
    C_start=[C_start;c_start];
    C_end=[C_end;c_end];

  
    % T3
    S_stimvel=ExcitRR.(an{a}).S.Average;
    s_start=S_stimvel(41,2);
    s_end=S_stimvel(46,2);
    S_start=[S_start;s_start];
    S_end=[S_end;s_end];


end 
%% Delta- difference in position 
C_dif_RN=C_end-C_start;
S_dif_RN=S_end-S_start;
%%  Plot position 

Scolor=rgb('Violet');


figure; hold on
scatter(C_start,C_end,30,'k','filled')
hold on
scatter(S_start,S_end,30,Scolor,'filled')


xlabel('Start stim position')
ylabel('End stim position')

%% Mean and SEM  
S_start_m=mean(S_start); S_end_m=mean(S_end);
S_sem_X = std(S_start) / sqrt(length(S_start)); % SEM for X
S_sem_Y = std(S_end) / sqrt(length(S_end)); % SEM for Y

C_start_m=mean(C_start); C_end_m=mean(C_end);
C_sem_X = std(C_start) / sqrt(length(C_start)); % SEM for X
C_sem_Y = std(C_end) / sqrt(length(C_end)); % SEM for Y


%% Plot mean with error bars 
figure; hold on
scatter(S_start_m,S_end_m,80,Scolor,'filled','MarkerEdgeColor','k')
hold on
errorbar(S_start_m, S_end_m, S_sem_X, S_sem_X, S_sem_Y, S_sem_Y,'k', 'LineWidth', 1); 
hold on
scatter(C_start_m,C_end_m,80,'k','filled','MarkerEdgeColor','k')
hold on
errorbar(C_start_m, C_end_m, C_sem_X, C_sem_X, C_sem_Y, C_sem_Y,'k','LineWidth', 1); 

ylabel('End Stim Position')
xlabel('Start Stim Position')
ylim([0.16 0.22])
xlim([0.154 0.174])
 
%% Find min and max of accel during CSM time window ( 100ms after stim ) 
for a=1:length(an)
    stimvel=ExcitRR.(an{a}).S.stimvel;
    Accelmin=[];
    Accelmax=[];
   stimvelsize=size(stimvel);
   num=(stimvelsize(3));

    for b=1:num
        reach=stimvel(:,:,b);
        window=reach(21:31,:);% row numbers for time of stim on to 100ms after - only for reaches with 200ms 
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
AccelMax_mean=[];
AccelMin_mean=[];
VMax=[];
RMSE=[];
RSquared=[];

figure;hold on
for a=1:length(an)
    Accelmin=ExcitRR.(an{a}).AccelMin; Accelmin_mean=mean(Accelmin);
    Accelmax=ExcitRR.(an{a}).AccelMax; Accelmax_mean=mean(Accelmax);
    Vmax=ExcitRR.(an{a}).RRKin.vmax(:,2);
    AccelMax=[AccelMax;Accelmax];
    AccelMin=[AccelMin;Accelmin];

    AccelMax_mean=[AccelMax_mean;Accelmax_mean];
    AccelMin_mean=[AccelMin_mean;Accelmin_mean];

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
%% 
KinVars=[VMax,AccelMin,AccelMax];
KinVars_mean=[AccelMin_mean,AccelMax_mean];
%% Plot control and CSM -ACCEL
% positive correlation
stimcolor= rgb('DodgerBlue');

%Accel vars 
x=KinVars(:,2);
y=KinVars(:,3);

% x=KinVars_mean(:,1);
% y=KinVars_mean(:,2);

%CorrCoef
[r, p] = corr(x,y); 

% Fit a linear model
mdl = fitlm(x, y);

% Sort the x data for plotting purposes
[x_sorted, idx] = sort(x); 
y_sorted = y(idx); 


% Predict the response and confidence intervals for sorted x data
[y_pred, y_ci] = predict(mdl, x_sorted, 'Alpha', 0.05); % 95% confidence intervals

figure;
% scatter(x, y,30,'filled','MarkerEdgeColor','k','MarkerFaceColor',stimcolor);
% hold on;
% Plot the fitted regression line
plot(x_sorted, y_pred,  'Color',stimcolor, 'LineWidth', 1.5); % Red regression line
hold on
% Plot the confidence intervals
fill([x_sorted; flipud(x_sorted)], [y_ci(:,1); flipud(y_ci(:,2))], ...
    [0.1172    0.5625    1.0000], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 


text(min(x), max(y), sprintf('r-neg = %.2f\n95%% p-value: [%.2f]', r, p), 'FontSize', 12, 'VerticalAlignment', 'top');

%unity_line = plot(xlim, xlim, 'k--'); % Black dashed unity line

% Adding labels and title
xlabel('AccelMin');
ylabel('AccelMax');
% xlim([0 300])
%ylim([0 80])
%set(gca,'XDir','reverse')
title('Accel Decel during CSM');
grid off;
hold off;

%%
mdl=fitlm(KinVars(:,2),KinVars(:,3))
figure; hold on
scatter(KinVars(:,2),KinVars(:,3),[],KinVars(:,1),'filled')
%scatter(KinVars(:,1),KinVars(:,2))
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
mdlall=fitlm(AccelMin,AccelMax)
figure; hold on
plot(mdlall)
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
xlim(xs)
%ylim([-500 0])
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

%% Plot velocity 
figure; hold on
for a=1:length(an)
    vels=ExcitRR.(an{a}).MaxVel; peak1=vels(:,1);peak2=vels(:,2);
    %figure; hold on
    scatter(peak1, peak2)
    hold on
end 
%% plot re-reach timing 
figure; hold on
colors=cool(5);
for a=1:length(an)
    RRtime=ExcitRR.(an{a}).StimRR(:,1);
    posRR=RRtime>=0;
    RRtimepos=RRtime(posRR);
    [cdf,x]=ecdf(RRtimepos);
    plot(x,cdf,'LineWidth',4,'Color',colors(a,:))
end
xlim([ 0 1])
%% Plot position per animal on individual reaches 
a = 1;
Sreaches=ChR2.(an{a}).StimReaches;
Creaches=ChR2.(an{a}).CtrlReaches;

x1=0 ;x2=0.05;y_upper=0.3; y_lower=0;

%scatcolor=rgb('CadetBlue');
scatcolor=rgb('DodgerBlue');

creachfields=fieldnames(Creaches);
sreachfields=fieldnames(Sreaches);


Sonly=AlignStim(Sreaches,1000); % need a long enough epoch to look at the re reach 
Conly=AlignThresh(Creaches,thresh,dirdir,1000);

Sonlyreach=Sonly.stimvel(:,2,34);% plot per reach
Conlyreach=Conly.stimvel(:,2,18); %5 is a good reach, also 19 %for double reach 38
stimes=Sonly.stimvel(:,1,2);
sz=20;
figure;
scatter(stimes,Sonlyreach,sz,'filled','MarkerFaceColor',[0.1172    0.5625    1.0000])
hold on
scatter(stimes,Conlyreach,sz,'filled','MarkerFaceColor',[0 0 0])
hold on
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
xlabel('Time from stim (s)')
%ylabel('Outward Position (cm)')
%ylabel('Upward Position (cm)')
ylabel('Lateral Position (cm)')
%ylim([0 0.1])
%ylim([0 0.25])
ylim([0 0.3])
xlim([-0.2 0.8])
%%  Plot all stim reaches for 1 animal
[SucStim,Failstim]=SuccessFilter(Sreaches);
Sucfields=fieldnames(SucStim);
Failfields=fieldnames(Failstim)

Sfailcolor=rgb('Crimson');
SucSonly=AlignStim(SucStim,1000)
FailSonly=AlignStim(Failstim,1000);

folder_path="C:\Users\Courtney\Documents\Summary Data\Plots\Random stim\Reacceleration Plots\Animal1_FailStim_positionPLots" ;

%for k=1:length(Sucfields)
for k=1:length(Failfields)
    %StimR=SucSonly.stimvel(:,2,k);
    FailR=FailSonly.stimvel(:,2,k);
    
    Conlyreach=Conly.stimvel(:,2,8); %8 is a Successful ctrl reach
    stimes=SucSonly.stimvel(:,1,2);
    
    sz=20;
    figure;
    %scatter(stimes,StimR,sz,'filled','MarkerFaceColor',[0.3711    0.6172    0.6250])
    scatter(stimes,FailR,sz,'filled','MarkerFaceColor',[0.8594    0.0781    0.2344])
    hold on
    scatter(stimes,Conlyreach,sz,'filled','MarkerFaceColor',[0 0 0])
    hold on
    fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
    xlabel('Time from stim (s)')
    ylabel('Outward Position (cm)')
    %title(Sucfields(k)) ;
    title(Failfields(k));
    %ylabel('Upward Position (cm)')
    %ylabel('Lateral Position (cm)')
    %ylim([0 0.1])
    %ylim([0 0.25])
    ylim([0 0.3])
    
  
    %filename=sprintf('SucStimReach%d.tif',k);
    filename=sprintf('FailStimReach%d.tif',k);
    full_path=fullfile(folder_path,filename);
    saveas(gcf,full_path);
    
    close;
end
 

%% Find velocity peaks 
for a=1:length(an)
    Sreaches=ChR2.(an{a}).Ssuc.stimvel; 
    Snums=size(Sreaches);Snum=Snums(3);

    Inx_peak=zeros(Snum,2);

    for reachnum=1:Snum  
        reach1 = Sreaches(1:21,5,reachnum); %velocity before stim for reaches aligned to stim
        reach1max=max(reach1);

        reach2 = Sreaches(26:end,5,reachnum); %velocity after stim
        reach2max=max(reach2);

       peaks=[reach1max,reach2max];
       Inx_peak(reachnum,:)=peaks;

    end

    ChR2.(an{a}).PeakV=Inx_peak;
end
%% Plot velocity peaks per animal 
%figure;hold on
for a=1:length(an)
    PeakVs=ChR2.(an{a}).PeakV;
  figure;hold on
    scatter(PeakVs(:,1),PeakVs(:,2))
    hold on 
end 

figure;hold on
for a=1:length(an)
    PeakVs=ChR2.(an{a}).PeakV;
    mdl=fitlm(PeakVs(:,1),PeakVs(:,2))
    scatter(PeakVs(:,1),PeakVs(:,2))
    hold on 
    plot(mdl)
    hold on
end 


figure;
for a=1:length(an)
    PeakVs=ChR2.(an{a}).PeakV;
    Peakmean=mean(PeakVs); Peaksd=std(PeakVs)
   
    plot(Peakmean(:,1),Peakmean(:,2))
    hold on 
    errorbar(Peakmean(:,1),Peakmean(:,2),Peaksd(:,2),'horizontal','LineStyle','none')
    hold on 
    errorbar(Peakmean(:,1),Peakmean(:,2),Peaksd(:,1),'horizontal','LineStyle','none')
    hold on
end

%% Get basic kinematics for these reaches 
Slowkin=BasicKinematics(StimRSlow);
Fastkin=BasicKinematics(StimRFast);

SS = AlignStim(StimRSlow,500);
%SS=stimvel_accel(SS);
SF=AlignStim(StimRFast,500)
SF=stimvel_accel(SF);

SFcolor=rgb('DarkCyan');
SScolor=rgb('LightSkyBlue');

dim=6;
StimvelPlot(SS,dim,SScolor,1)
StimvelPlot(SF,dim,SFcolor,0)
%pbaspect([ 1 2 1])
%ylim([0 0.3])
ylim('auto')

figure()
scatter(Slowkin.vmax(:,2),Slowkin.aavg_decel(:,2))
hold on
scatter(Fastkin.vmax(:,2),Fastkin.aavg_decel(:,2))

decS=mean(Slowkin.aavg_decel(:,2));
decF=mean(Fastkin.aavg_decel(:,2));
Dec=[decF,decS];
accS=mean(Slowkin.aavg_accel(:,2));
accF=mean(Fastkin.aavg_accel(:,2));
Acc=[accF,accS];
vmaxS=mean(Slowkin.vmax(:,2));
vmaxF=mean(Fastkin.vmax(:,2));
Vmax=[vmaxF,vmaxS];
%%
% Compute the empirical cumulative distribution function (ECDF)
Faccel=Fastkin.aavg_accel(:,2);
Saccel=Slowkin.aavg_accel(:,2);
[f, x] = ecdf(Faccel);
[fs,fx]=ecdf(Saccel);
[h,p,ksstat]=kstest2(Faccel,Saccel);

% Create the cumulative distribution plot
figure;
plot(x, f, 'b-', 'LineWidth', 2);
hold on
plot(fx, fs, 'g-', 'LineWidth', 2);
xlabel('X-axis Label');
ylabel('Cumulative Probability');
title('Cumulative Distribution Plot');

% Optionally, add a grid
grid on;

% Display the KS test results
fprintf('Kolmogorov-Smirnov Test:\n');
fprintf('KS Statistic = %.4f\n', ksstat);
fprintf('p-value = %.4f\n', p);

% Determine statistical significance
alpha = 0.05; % Set the significance level
if p < alpha
    fprintf('The two datasets are statistically different (reject null hypothesis).\n');
else
    fprintf('There is no significant difference between the two datasets (fail to reject null hypothesis).\n');
end
    
%% Plot value of Vmax percentile for each animal 
orange=rgb('Orange');
orangered=rgb('orangered');
red=rgb('Crimson');
figure;hold on
for a=1:length(an)
   Pctvalues=ChR2.(an{a}).VmaxPrct
   sz=60;
   scatter(25,Pctvalues.p25,sz,'filled','MarkerFaceColor',[1 0.6445 0])  
   hold on
  scatter(50,Pctvalues.p50,sz,'filled','MarkerFaceColor',[ 1.0000  0.2695  0])
  hold on
  scatter(75,Pctvalues.p75,sz,'filled','MarkerFaceColor',[  0.8594 0.0781 0.2344])
  hold on
  plot([25,50],[Pctvalues.p25,Pctvalues.p50],'k')
  hold on
  plot([50,75],[Pctvalues.p50,Pctvalues.p75],'k')
  %text(25, Pctvalues.p25, an{a}, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
  text(75, Pctvalues.p75, an{a}, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end 
  
customticks=[25,50,75];
xticks(customticks)
xticklabels({25, 50 75})
xlabel('Percentile')
ylabel('Outward Velocity')
title('Velocity Percentiles per animal')

%% Look at acceleration 
C_big_accel=stimvel_accel(C_big);
S_big_accel=stimvel_accel(S_big);

dim =8 ; %outward acceleration
x1=0 ;x2=0.05;
%y_upper=1.5; y_lower=-0.5;
y_upper=30; y_lower=-40;
StimvelPlot(C_big_accel,dim,'k',1)
StimvelPlot(S_big_accel,dim,stimcolor,0)
hold on 
% Fill the area between the vertical lines
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');

ylim('auto')
yline(0,'--k')
xlim([-0.2 0.2])
%ylim([-0.5 1.5])
%ylim([0 .25])
%yline(0,'--')
xlabel('Time from Stim (s)')
ylabel('Outward Velocity (cm/s)')
%pbaspect([ 1 2 1])


%%
%get average standard deviation stimvel plot for 'out significant' animals

C_std = struct;
S_std = struct;
Ssuc_std=struct;
Sfail_std=struct;
for a = 1:length(an)
   
    C_std.stimvel(:,:,a) = ChR2.(an{a}).C.Std;
    S_std.stimvel(:,:,a) = ChR2.(an{a}).S.Std;
    Ssuc_std.stimvel(:,:,a) = ChR2.(an{a}).Ssuc.Std;
    Sfail_std.stimvel(:,:,a) = ChR2.(an{a}).Sfail.Std;
    
end


C_std.Average = nanmean(C_std.stimvel,3);
C_std.Average(:,1) = C.Average(:,1);
C_std.Std = nanstd(C_std.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(C_std.Average),1); 
for nnum = 1:length(C_std.Average)
    samp = squeeze(C_std.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
C_std.N = N;



S_std.Average = nanmean(S_std.stimvel,3);
S_std.Average(:,1) = C.Average(:,1);
S_std.Std = nanstd(S_std.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(C_std.Average),1); 
for nnum = 1:length(C_std.Average)
    samp = squeeze(S_std.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end


S_std.N = N;


dim = 5;

StimvelPlot(C_std,dim,[.5 .5 .5],1)
StimvelPlot(S_std,dim,[.2 .6 .9],0)
ylim auto

figure; hold on
StimvelPlot2(C_std,dim,[.5 .5 .5],0)
shadedErrorBars(C_std.Average(:,1),C_std.Average(:,dim),C_std.Std(:,dim)./sqrt(C_std.N),[.5 .5 .5])
StimvelPlot2(S_std,dim,[.2 .6 .9],0)
shadedErrorBars(S_std.Average(:,1),S_std.Average(:,dim),S_std.Std(:,dim)./sqrt(S_std.N),[.2 .6 .9])
ylim auto

%%
%difference curves and pvalue heatmaps
Diff_stimvel = struct;
Pval_chr2 = []; Zval_chr2 = [];
for a = 1:length(an)
   
    C = ChR2.(an{a}).C;
    S = ChR2.(an{a}).S;
    
%     C = ChR2.(an{a}).Ssuc;
%     S = ChR2.(an{a}).Sfail;
    
    
    Diff_chr2.stimvel(:,:,a) = S.Average - C.Average;
    Diff_chr2.stimvel(:,1,a) = C.Average(:,1);
    
    
    %get pvalues at each time point for each animal
    
    
    parray = ones(size(C.stimvel,1),2);
    
    for pnum = 1: length(parray)
        a4b = squeeze(C.stimvel(pnum,10,:)); % change this when you want to look at positin or velocity 
        b4b = squeeze(S.stimvel(pnum,10,:));
        
        [parray(pnum,1),b,pstat] = ranksum(a4b,b4b,'method','approximate');
        parray(pnum,2) = pstat.zval;
    end

    Pval_chr2(a,:) = parray(:,1);
    Zval_chr2(a,:) = parray(:,2);
end

Diff_chr2.Average = nanmean(Diff_chr2.stimvel,3);
Diff_chr2.Std = nanstd(Diff_chr2.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(Diff_chr2.Average),1); 
for nnum = 1:length(Diff_chr2.Average)
    samp = squeeze(Diff_chr2.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
Diff_chr2.N = N;



% figure;
% hold on;

subplot(2,1,1); hold on
x1=0 ;x2=0.05;
%y_upper=1.5; y_lower=-0.5;
y_upper=40; y_lower=-30;

dmap = brewermap(13,'Spectral');
dmap = dmap([1:5 9:13],:); 
dim = 10; % make sure to change this 
for dnum = 1:size(Diff_chr2.stimvel,3)
    
    reach = Diff_chr2.stimvel(:,:,dnum);
    plot(reach(:,1),reach(:,dim),'Color',[dmap(dnum,:) .3],'LineWidth',2)
    
    [~,minx] = min(reach(:,dim));
    
    scatter(reach(minx,1),1.0,40,dmap(dnum,:),'filled')

end


StimvelPlot2(Diff_chr2,dim,[.2 .2 .2],0)
hold on
% Fill the area between the vertical lines
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
yline(0,'--k')
ylim([-20 20])
xlim([-.1 .2])
%ylim([-.05 0.05])

subplot(2,1,2); hold on

gmap = brewermap(96,'Greys');
rmap = brewermap(64,'Reds');
combomap = flipud([gmap(1:32,:) ; rmap(33:end,:)]);
%figure; hold on
imagesc(Pval_chr2)
colormap(combomap);
set(gca,'colorscale','log')
caxis([.005 .5])
%colorbar('Ticks',[.005 .01 .05 .1 .5]);
xlim([10.5 41.5])
%ylim([1 5])



%%
%3D endpoints

dmap = brewermap(13,'Spectral');
dmap = dmap([1:5 9:13],:);
figure; hold on
for a = 1:length(an)
    
    endpts_c = ChR2.(an{a}).endpoint.control;
    scatter3(endpts_c(:,1),endpts_c(:,2),endpts_c(:,3),15,dmap(2,:),'filled');
    
    
end
for a = 1:length(an)
    endpts_s = ChR2.(an{a}).endpoint.stim;
    scatter3(endpts_s(:,1),endpts_s(:,2),endpts_s(:,3),15,dmap(8,:),'filled');
    
end


figure; hold on
for a = 1:length(an)
    
    
    scatter3(chr2_endpt_c(:,1),chr2_endpt_c(:,2),chr2_endpt_c(:,3),70,dmap(2,:),'filled','MarkerFaceAlpha',.4);
    scatter3(chr2_endpt_s(:,1),chr2_endpt_s(:,2),chr2_endpt_s(:,3),70,dmap(8,:),'filled','MarkerFaceAlpha',.4);

    
end



%%

% Success Rate

for a = 1:length(an)
    SR_control(a,1) = ChR2.(an{a}).successrate.control;
    SR_stim(a,1) = ChR2.(an{a}).successrate.stim;
    
end

[p,h] = signrank(SR_control,SR_stim);

figure; hold on

scatter(ones(1,length(SR_control)),SR_control,45,'filled',...
    'jitter','on','jitterAmount',0.1);
scatter(ones(1,length(SR_stim))*2,SR_stim,45,'filled',...
    'jitter','on','jitterAmount',0.1);

title(['Success Rate: p = ',num2str(p)])
xlim([0,3])


%%

%plot single stimvel curve for animal

a = 1;
stimcolor= rgb('DodgerBlue');
linecolor=rgb('DarkGray');

C = ChR2.(an{a}).C;
S = ChR2.(an{a}).S;
dim = 2;

StimvelPlot(C,dim,'k',1)
StimvelPlot(S,dim,stimcolor,0)

xlim([-0.2 0.2])
ylim([0 0.3])
%ylim([-0.5 0.5])

%yline(0,'--')
xlabel('Time from Stim (s)')
%ylabel('Outward Velocity (cm/s)')
ylabel('Outward position(cm/s)')
%pbaspect([1 2 1])

%% 
% Plot variabiity of endpoitns 
StimScolor=rgb('DodgerBlue');
StimFcolor=rgb('OrangeRed');
CtrlScolor=rgb('DimGray');
CtrlFcolor=rgb('DarkRed');
CtrlSucend=[]; CtrlFailend=[];
CtrlSucSD=[];CtrlFailSD=[];

StimSucend=[]; StimFailend=[];
StimSucSD=[]; StimFailSD=[];


an=fieldnames(ChR2);
numbins=20;
for a=1:length(an) 
    endpts_cS = ChR2.(an{a}).CtrlSKin.endpt;
    endpts_cS=endpts_cS*10;
    CSucend=mean(endpts_cS);
    CtrlSucend=[CtrlSucend; CSucend];
    CSucSD=std(endpts_cS,0,1);
    CtrlSucSD=[CtrlSucSD; CSucSD];
  
    endpts_cF=ChR2.(an{a}).CtrlFKin.endpt;
    endpts_cF=endpts_cF*10;
    CFailend=mean(endpts_cF);
    CtrlFailend=[CtrlFailend; CFailend];
    CFailSD=std(endpts_cF,0,1);
    CtrlFailSD=[CtrlFailSD; CFailSD];

    endpts_sS = ChR2.(an{a}).StimSKin.endpt;
    endpts_sS=endpts_sS*10;
    SSucend=mean(endpts_sS);
    StimSucend=[StimSucend; SSucend];
    SSucSD=std(endpts_sS,0,1);
    StimSucSD=[StimSucSD; SSucSD];
  
    endpts_sF=ChR2.(an{a}).StimFKin.endpt;
    endpts_sF=endpts_sF*10;
    SFailend=mean(endpts_sF);
    StimFailend=[StimFailend; SFailend];
    SFailSD=std(endpts_sF,0,1);
    StimFailSD=[StimFailSD; SFailSD];

    figure;
    histogram(endpts_cS,numbins,'FaceColor',CtrlScolor)
    hold on 
    histogram(endpts_cF,numbins,'FaceColor',CtrlFcolor)
    hold on 
    histogram(endpts_sS,numbins,'FaceColor',StimScolor)
    hold on 
    histogram(endpts_sF,numbins,'FaceColor',StimFcolor)
end

%% plot endpoints
% Plot variabiity of endpoitns 
Stimcolor=rgb('DodgerBlue');

Ctrlcolor=rgb('DimGray');

Ctrlend=[]; 
Stimend=[]; 


an=fieldnames(ChR2);
%numbins=20;
figure; 
for a=1:length(an) 
    endpts_c = ChR2.(an{a}).CtrlKin.endpt;
    endpts_c=endpts_c*10;
    Cend=mean(endpts_c);
    Ctrlend=[Ctrlend; Cend];
  

    endpts_s= ChR2.(an{a}).StimKin.endpt;
    endpts_s=endpts_s*10;
    Send=mean(endpts_s);
    Stimend=[Stimend; Send];

    %figure;
    histogram(endpts_c,'FaceColor',Ctrlcolor)
    hold on 
    histogram(endpts_s,'FaceColor',Stimcolor)
 
end