%name your master struct

ChR2= struct;
%use strong stim folder - session 

fname = uigetdir;
direc = dir(strcat(fname,'\*.mat'));
%%

for i = 1:length(direc)
    %load reaches
    load(fullfile(fname,direc(i).name));
    
     %get name
    fullname = direc(i).name;
    mouse_name = regexprep(fullname,'.mat','');
    ChR2.(mouse_name).anim_data = anim_data;
    
    %setup variables 
    ag = 6;
    filt_state = 1;
    pct = .8;
    leng = 5;
    thresh = .16;
   
    dnum = fieldnames (anim_data);

 %First session
    Reaches=anim_data.(dnum{1}).Reaches;
    [Reaches,bad]=rfilts_standard(Reaches,thresh,5,0.8);
    Reaches=accel_toreaches(Reaches);
    ChR2.(mouse_name).Reaches_first = Reaches;
    
    % % % do stimvel
    % get stim good for averaging directly off of stim point instead of thresh point
    thresh = .16; dirdir = 1;
    [Stim, Control] = StimFilter(Reaches);
    [Stim2] = stimfilt_thresh(Stim,thresh);
    [Stim3] = stimfilt_firstreach(Stim2, thresh, dirdir);
    [SucStim,FailStim]=SuccessFilter(Stim3);
    StimKin=BasicKinematics(Stim3);
    CtrlKin=BasicKinematics(Control);

    S = AlignStim(Stim3,200);
    S=stimvel_accel(S);
    Ssuc=AlignStim(SucStim,200);
    Ssuc=stimvel_accel(Ssuc);
    Sfail=AlignStim(FailStim,200);
    Sfail=stimvel_accel(Sfail);

    C = AlignThresh(Control,thresh,dirdir,200);
    C=stimvel_accel(C);

    % basic kinematics 
    % Separate out success and control 
    [S_C, F_C] = SuccessFilter(Control);
    [S_S, F_S] = SuccessFilter(Stim3);
    
    ChR2.(mouse_name).successrate_first.control = length(fieldnames(S_C)) / length(fieldnames(Control));
    ChR2.(mouse_name).successrate_first.stim = length(fieldnames(S_S)) / length(fieldnames(Stim3));

     %first session
    ChR2.(mouse_name).Control=Control;
     ChR2.(mouse_name).CtrlKin=CtrlKin;
    ChR2.(mouse_name).C=C;
    ChR2.(mouse_name).Stim=Stim3;
    ChR2.(mouse_name).StimKin=StimKin;
    ChR2.(mouse_name).S=S;
    ChR2.(mouse_name).Ssuc=Ssuc;
    ChR2.(mouse_name).Sfail=Sfail;

    %Last session
    Reaches_end=anim_data.(dnum{end}).Reaches;
    [Reaches_end,bad_end]=rfilts_standard(Reaches_end,thresh,5,0.8);
    Reaches_end=accel_toreaches(Reaches_end);
    ChR2.(mouse_name).Reaches_end = Reaches_end;
    %stimvel
    [Stim_end, Control_end] = StimFilter(Reaches_end);
    [Stim2_end] = stimfilt_thresh(Stim_end,thresh);
    [Stim3_end] = stimfilt_firstreach(Stim2_end, thresh, dirdir);
    [SucStim_end,FailStim_end]=SuccessFilter(Stim3_end);
    StimEndKin=BasicKinematics(Stim3_end);
    CtrlEndKin=BasicKinematics(Control_end);

        %stimvel
    S_end = AlignStim(Stim3_end,200);
    S_end=stimvel_accel(S_end);
    Ssuc_end=AlignStim(SucStim_end,200);
    Ssuc_end=stimvel_accel(Ssuc_end);
    Sfail_end=AlignStim(FailStim_end,200);
    Sfail_end=stimvel_accel(Sfail_end);

    C_end = AlignThresh(Control_end,thresh,dirdir,200);
    C_end=stimvel_accel(C_end);
     
     % Separate out success and control 
    [S_Cend, F_Cend] = SuccessFilter(Control_end);
    [S_Send, F_Send] = SuccessFilter(Stim3_end);
    
    ChR2.(mouse_name).successrate_end.control = length(fieldnames(S_Cend)) / length(fieldnames(Control_end));
    ChR2.(mouse_name).successrate_end.stim = length(fieldnames(S_Send)) / length(fieldnames(Stim3_end));
     

        %last session
       ChR2.(mouse_name).Control_end=Control_end;
       ChR2.(mouse_name).CtrlEndKin=CtrlEndKin;
       ChR2.(mouse_name).C_end=C_end;
       ChR2.(mouse_name).Stim_end=Stim3_end;
       ChR2.(mouse_name).StimEndKin=StimEndKin;
       ChR2.(mouse_name).S_end=S_end;
       ChR2.(mouse_name).Ssuc=Ssuc;
       ChR2.(mouse_name).Sfail=Sfail;
   
         
    end
    
clearvars -except ChR2






%% Create structure that averages all stimvel structure
C1=struct;
S1=struct;
Cend=struct;
Send=struct;

an=fieldnames(ChR2);
for a=1:length(an)
        C1.stimvel(:,:,a) =ChR2.(an{a}).C.Average;
        S1.stimvel(:,:,a) =ChR2.(an{a}).S.Average;

        Cend.stimvel(:,:,a) =ChR2.(an{a}).C_end.Average;
        Send.stimvel(:,:,a) =ChR2.(an{a}).S_end.Average;
       
end

C1.Average = nanmean(C1.stimvel,3);
C1.Std = nanstd(C1.stimvel,0,3);

N = zeros(length(C1.Average),1); 
for nnum = 1:length(C1.Average)
    samp = squeeze(C1.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
C1.N = N;

S1.Average = nanmean(S1.stimvel,3);
S1.Std = nanstd(S1.stimvel,0,3);

N = zeros(length(S1.Average),1); 
for nnum = 1:length(S1.Average)
    samp = squeeze(S1.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
S1.N = N;


Cend.Average = nanmean(Cend.stimvel,3);
Cend.Std = nanstd(Cend.stimvel,0,3);

N = zeros(length(Cend.Average),1); 
for nnum = 1:length(Cend.Average)
    samp = squeeze(Cend.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
Cend.N = N;

Send.Average = nanmean(Send.stimvel,3);
Send.Std = nanstd(Send.stimvel,0,3);

N = zeros(length(Send.Average),1); 
for nnum = 1:length(Send.Average)
    samp = squeeze(Send.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));
end
Send.N = N;

%%
% stim1=rgb('DodgerBlue');
% ctrl1=rgb('DarkOrange');
% stimend=rgb('MediumBlue');
% ctrlend=rgb('MediumPurple');

stim1=rgb('LightSteelBlue');
ctrl1=rgb('DimGray');
stimend=rgb('Navy');
ctrlend=rgb('Black');

dim = 2;
x1=0 ;x2=0.05;
y_upper=0.3; y_lower=0;

StimvelPlot(C1,dim,ctrl1,1)
StimvelPlot(S1,dim,stim1,0)
StimvelPlot(Cend,dim,ctrlend,0)
StimvelPlot(Send,dim,stimend,0)
hold on 
% Fill the area between the vertical lines
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
yline(0,'--k')
xlim([-0.2 0.2])
%ylim([-0.5 1.5])
%ylim([-15 15])
ylim('auto')
%yline(0,'--')
xlabel('Time from Stim (s)')
ylabel('Outward Velocity (cm/s)')
pbaspect([ 1 2 1])   

%% Get RR info for first and last session
thresh=0.16;dirdir = 1;

% Use ReReach Filter - Control Reaches 
ExcitRR_Ctrl1=[]; % Savethis structure to use for graphing rereach timing in sep. script
for a=1:length(an)

    Creaches=ChR2.(an{a}).Control;


    % ReReach
    %[ReReaches,StimRR,MaxVel]=ReReachFilter_accel_ctrlreaches(Creaches); % currently looking within 200ms 
    [ReReaches,StimRR,MaxVel]=ReReachFilter_accel_ctrlreaches_edit(Creaches); 
    RRKin=BasicKinematics(ReReaches);

    C = AlignThresh(ReReaches,thresh,dirdir,200);
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
        CN = AlignThresh(NoReach,thresh,dirdir,200);
        %CN=stimvel_accel(CN); %causing issues in the code 
    end 

    [S_C, F_C] = SuccessFilter(ReReaches);
    [S_S, F_S] = SuccessFilter(NoReach);
    
    ExcitRR_Ctrl1.(an{a}).successrate.ReReach = length(fieldnames(S_C)) / length(fieldnames(Creaches));
    ExcitRR_Ctrl1.(an{a}).successrate.NoReReach = length(fieldnames(S_S)) / length(fieldnames(Creaches));
    ExcitRR_Ctrl1.(an{a}).CtrlReaches=Creaches;
    ExcitRR_Ctrl1.(an{a}).C=C;
    ExcitRR_Ctrl1.(an{a}).ReReaches=ReReaches;
    ExcitRR_Ctrl1.(an{a}).RRKin=RRKin;
    ExcitRR_Ctrl1.(an{a}).NoReReach=NoReach;
    ExcitRR_Ctrl1.(an{a}).NRKin=NRKin;
    ExcitRR_Ctrl1.(an{a}).CN=CN;
    ExcitRR_Ctrl1.(an{a}).StimRR=StimRR;
    ExcitRR_Ctrl1.(an{a}).MaxVel=MaxVel;
    % get % of reaches with a rereach
    cnum=length(fields1);
    rrnum=length(fields2);
    percentrr=rrnum/cnum;
    ExcitRR_Ctrl1.(an{a}).Cnum=cnum;
    ExcitRR_Ctrl1.(an{a}).RRnum=rrnum;
    ExcitRR_Ctrl1.(an{a}).PercentRR=percentrr;
  
end 

%%
% Use ReReach Filter - Control Reaches 
ExcitRR_Ctrlend=[]; % Savethis structure to use for graphing rereach timing in sep. script
for a=1:length(an)

    Creaches=ChR2.(an{a}).Control_end;

 
    % ReReach
   % [ReReaches,StimRR,MaxVel]=ReReachFilter_accel_ctrlreaches(Creaches); % currently looking within 200ms 
    [ReReaches,StimRR,MaxVel]=ReReachFilter_accel_ctrlreaches_edit(Creaches);
    RRKin=BasicKinematics(ReReaches);

    C = AlignThresh(ReReaches,thresh,dirdir,200);
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
        CN = AlignThresh(NoReach,thresh,dirdir,200);
        %CN=stimvel_accel(CN); %causing issues in the code 
    end 

    [S_C, F_C] = SuccessFilter(ReReaches);
    [S_S, F_S] = SuccessFilter(NoReach);
    
    ExcitRR_Ctrlend.(an{a}).successrate.ReReach = length(fieldnames(S_C)) / length(fieldnames(Creaches));
    ExcitRR_Ctrlend.(an{a}).successrate.NoReReach = length(fieldnames(S_S)) / length(fieldnames(Creaches));
    ExcitRR_Ctrlend.(an{a}).CtrlReaches=Creaches;
    ExcitRR_Ctrlend.(an{a}).C=C;
    ExcitRR_Ctrlend.(an{a}).ReReaches=ReReaches;
    ExcitRR_Ctrlend.(an{a}).RRKin=RRKin;
    ExcitRR_Ctrlend.(an{a}).NoReReach=NoReach;
    ExcitRR_Ctrlend.(an{a}).NRKin=NRKin;
    ExcitRR_Ctrlend.(an{a}).CN=CN;
    ExcitRR_Ctrlend.(an{a}).StimRR=StimRR;
    ExcitRR_Ctrlend.(an{a}).MaxVel=MaxVel;
    % get % of reaches with a rereach
    cnum=length(fields1);
    rrnum=length(fields2);
    percentrr=rrnum/cnum;
    ExcitRR_Ctrlend.(an{a}).Cnum=cnum;
    ExcitRR_Ctrlend.(an{a}).RRnum=rrnum;
    ExcitRR_Ctrlend.(an{a}).PercentRR=percentrr;
  
end 

%% Use ReReach Filter- STIM reaches 
ExcitRR=[]; % Savethis structure to use for graphing rereach timing in sep. script
for a=1:length(an)
    Sreaches=ChR2.(an{a}).Stim;
  

    % Stim Re
    [ReReaches,StimRR,MaxVel]=ReReachFilter_accel(Sreaches);
    RRKin=BasicKinematics(ReReaches);
    
    S = AlignStim(ReReaches,200);
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
        SN = AlignStim(NoReach,200);
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

%% Use ReReach Filter- STIM reaches 
ExcitRR_end=[]; % Savethis structure to use for graphing rereach timing in sep. script
for a=1:length(an)
    Sreaches=ChR2.(an{a}).Stim_end;
  

    % Stim Re
    [ReReaches,StimRR,MaxVel]=ReReachFilter_accel(Sreaches);
    RRKin=BasicKinematics(ReReaches);
    
    S = AlignStim(ReReaches,200);
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
        SN = AlignStim(NoReach,200);
        %SN=stimvel_accel(SN); %causing issues in the code 
    end 

    [S_C, F_C] = SuccessFilter(ReReaches);
    [S_S, F_S] = SuccessFilter(NoReach);
    
    ExcitRR_end.(an{a}).successrate.ReReach = length(fieldnames(S_C)) / length(fieldnames(Sreaches));
    ExcitRR_end.(an{a}).successrate.NoReReach = length(fieldnames(S_S)) / length(fieldnames(Sreaches));
    ExcitRR_end.(an{a}).ReReaches=ReReaches;
    ExcitRR_end.(an{a}).RRKin=RRKin;
    ExcitRR_end.(an{a}).S=S;
    ExcitRR_end.(an{a}).NoReReach=NoReach;
    ExcitRR_end.(an{a}).NRKin=NRKin;
    ExcitRR_end.(an{a}).SN=SN;
    ExcitRR_end.(an{a}).StimRR=StimRR;
    ExcitRR_end.(an{a}).MaxVel=MaxVel;
    snum=length(fields1);
    rrnum=length(fields2);
    percentrr=rrnum/snum;
    ExcitRR_end.(an{a}).Snum=snum;
    ExcitRR_end.(an{a}).RRnum=rrnum;
    ExcitRR_end.(an{a}).PercentRR=percentrr;
end 
%% Ctrl first
ctrl1 = fieldnames(ExcitRR_Ctrl1);
for b=1:length(ctrl1)
    RR_Excit=ExcitRR_Ctrl1.(ctrl1{b}).StimRR;
 
    nonZeroRows_ctrl1 = [];

    % Iterate over each row of the array
    for i = 1:size(RR_Excit, 1)
        % Check if the current row contains any zeros
        if any(RR_Excit(i,:) == 0)
            % If zeros are present, skip this row
            continue;
        end
        
        % Append the current row to nonZeroRows if no zeros are present
        nonZeroRows_ctrl1 = [nonZeroRows_ctrl1; RR_Excit(i,:)];
    end
    
        % Display the resulting array
        disp(nonZeroRows_ctrl1);
        ExcitRR_Ctrl1.(ctrl1{b}).StimRREdit=nonZeroRows_ctrl1;
end
%%
ctrlend = fieldnames(ExcitRR_Ctrlend);
for b=1:length(ctrlend)
    RR_Excit=ExcitRR_Ctrlend.(ctrlend{b}).StimRR;
 
    nonZeroRows_ctrlend = [];

    % Iterate over each row of the array
    for i = 1:size(RR_Excit, 1)
        % Check if the current row contains any zeros
        if any(RR_Excit(i,:) == 0)
            % If zeros are present, skip this row
            continue;
        end
        
        % Append the current row to nonZeroRows if no zeros are present
        nonZeroRows_ctrlend = [nonZeroRows_ctrlend; RR_Excit(i,:)];
    end
    
        % Display the resulting array
        disp(nonZeroRows_ctrlend);
        ExcitRR_Ctrlend.(ctrlend{b}).StimRREdit=nonZeroRows_ctrlend;
end
%%
stim1=fieldnames(ExcitRR)
for b=1:length(stim1)
    RR_Excit=ExcitRR.(stim1{b}).StimRR;
 
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
        ExcitRR.(stim1{b}).StimRREdit=nonZeroRows_ex;
end

%%
stimend=fieldnames(ExcitRR_end)
for b=1:length(stimend)
    RR_Excit=ExcitRR_end.(stimend{b}).StimRR;
 
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
        ExcitRR_end.(stimend{b}).StimRREdit=nonZeroRows_ex;
end

%% New Method for calculating CDF
% Define a common set of x-values (bin edges)
X_Min = 0;
X_Max = 0.2;

% Gather all data points to determine the bin edges
allData = [];
ctrl1=fieldnames(ExcitRR_Ctrl1);
for a = 1:length(ctrl1)
    allT = ExcitRR_Ctrl1.(ctrl1{a}).StimRREdit(:,1);
    posRRall = allT >= 0;
    RRTpos_all = allT(posRRall);
    allData = [allData; RRTpos_all];
end

ctrlend=fieldnames(ExcitRR_Ctrlend);
for m = 1:length(ctrlend)
    excitT = ExcitRR_Ctrlend.(ctrlend{m}).StimRREdit(:,1);
    posRRext = excitT > 0;
    RRTpos_ext = excitT(posRRext);
    allData = [allData; RRTpos_ext];
end

stim1=fieldnames(ExcitRR);
for t = 1:length(stim1)
    RNT = ExcitRR.(stim1{t}).StimRREdit(:,1);
    posRRextRN = RNT >= 0;
    RRTpos_extRN = RNT(posRRextRN);
    allData = [allData; RRTpos_extRN];
end

stimend=fieldnames(ExcitRR_end);
for p = 1:length(stimend)
    RNT1 = ExcitRR_end.(stimend{p}).StimRREdit(:,1);
    posRRextRN1 = RNT1 >= 0;
    RRTpos_extRN1 = RNT1(posRRextRN1);
    allData = [allData; RRTpos_extRN1];
end

%%
% Use Freedman-Diaconis rule to determine the bin width
iqrAll = iqr(allData); % Interquartile range
binWidth = 2 * iqrAll / (length(allData)^(1/3));
num_bins = ceil((X_Max - X_Min) / binWidth);

% Generate bin edges
common_x = linspace(X_Min, X_Max, num_bins);

% Initialize arrays to store individual CDFs
ctrl1CDFs = zeros(num_bins - 1, length(ctrl1));
ctrlendCDFs = zeros(num_bins - 1, length(ctrlend));
stim1CDFs = zeros(num_bins - 1, length(stim1));
stimendCDFs = zeros(num_bins - 1, length(stimend));

% Define colors
stim1color=rgb('LightSteelBlue');
ctrl1color=rgb('DimGray');
stimendcolor=rgb('Navy');
ctrlendcolor=rgb('Black');

figure;
hold on

%first ctrl reach
for a = 1:length(ctrl1)
    allT = ExcitRR_Ctrl1.(ctrl1{a}).StimRREdit(:,1);
    posRRall = allT >= 0;
    RRTpos_all = allT(posRRall);

    % Calculate the CDF using the common_x as bin edges
    cdf = histcounts(RRTpos_all, common_x, 'Normalization', 'cdf');

    % Store the individual CDFs
    ctrl1CDFs(:, a) = cdf;

    % Plot the individual subject's CDF
    plot(common_x(1:end-1), cdf, 'LineWidth', 1.5, 'Color', ctrl1color);
end

% end ctrl reach
for m = 1:length(ctrlend)
    excitT = ExcitRR_Ctrlend.(ctrlend{m}).StimRREdit(:,1);
    posRRext = excitT > 0;
    RRTpos_ext = excitT(posRRext);

    % Calculate the CDF using the common_x as bin edges
    cdfE = histcounts(RRTpos_ext, common_x, 'Normalization', 'cdf');

    % Store the individual CDFs
    ctrlendCDFs(:, m) = cdfE;

    % Plot the individual subject's CDF
    plot(common_x(1:end-1), cdfE, 'LineWidth', 1.5, 'Color', ctrlendcolor);
end

% stim start reaches 
for t = 1:length(stim1)
    RNT = ExcitRR.(stim1{t}).StimRREdit(:,1);
    posRRextRN = RNT >= 0;
    RRTpos_extRN = RNT(posRRextRN);

    % Calculate the CDF using the common_x as bin edges
    cdfRN = histcounts(RRTpos_extRN, common_x, 'Normalization', 'cdf');

    % Store the individual CDFs
    stim1CDFs(:, t) = cdfRN;

    % Plot the individual subject's CDF
    plot(common_x(1:end-1), cdfRN, 'LineWidth', 1.5, 'Color', stim1color);
end

% stim start reaches 
for p = 1:length(stimend)
    RNT1 = ExcitRR_end.(stimend{p}).StimRREdit(:,1);
    posRRextRN1 = RNT1 >= 0;
    RRTpos_extRN1 = RNT1(posRRextRN1);

    % Calculate the CDF using the common_x as bin edges
    cdfRN1 = histcounts(RRTpos_extRN1, common_x, 'Normalization', 'cdf');

    % Store the individual CDFs
    stimendCDFs(:, p) = cdfRN1;

    % Plot the individual subject's CDF
    plot(common_x(1:end-1), cdfRN1, 'LineWidth', 1.5, 'Color', stimendcolor);
end

% Calculate the average CDFs
averagectrl1CDF = mean(ctrl1CDFs, 2);
averagectrlendCDF = mean(ctrlendCDFs, 2);
averagestim1CDF = mean(stim1CDFs, 2);
averagestimendCDF = mean(stimendCDFs, 2);

% Plot the average CDFs
plot(common_x(1:end-1), averagectrl1CDF, 'LineWidth', 4, 'Color', ctrl1color);
plot(common_x(1:end-1), averagectrlendCDF, 'LineWidth', 4, 'Color', ctrlendcolor);
plot(common_x(1:end-1), averagestim1CDF, 'LineWidth', 4, 'Color', stim1color);
plot(common_x(1:end-1), averagestimendCDF, 'LineWidth', 4, 'Color', stimendcolor);


xlabel('Timing of corrective sub-movement(s)');
ylabel('Probability');
ylim([0 1]);
%xlim([0 0.2])

%%
BigRR_ctrl1=[];
BigRR_ctrl1mean=[];

an=fieldnames(ExcitRR_Ctrl1);
for a=1:length(an)
     %allRR =  ExcitRR_Ctrl1.(an{a}).StimRREdit(:,1);
     allRR =  ExcitRR_Ctrl1.(an{a}).StimRR(:,1);
     BigRR_ctrl1=[BigRR_ctrl1;allRR];
     allmean=mean(allRR);
     BigRR_ctrl1mean=[BigRR_ctrl1mean;allmean];
    
end
%%
BigRR_ctrlend=[];
BigRR_ctrlendmean=[];

an=fieldnames(ExcitRR_Ctrlend);
for a=1:length(an)
     %allRR =  ExcitRR_Ctrlend.(an{a}).StimRREdit(:,1);
     allRR =  ExcitRR_Ctrlend.(an{a}).StimRR(:,1);
     BigRR_ctrlend=[BigRR_ctrlend;allRR];
     allmean=mean(allRR);
     BigRR_ctrlendmean=[BigRR_ctrlendmean;allmean];
    
end
%%
BigRR_stim1=[];
BigRR_stim1mean=[];

an=fieldnames(ExcitRR);
for a=1:length(an)
     %allRR =  ExcitRR.(an{a}).StimRREdit(:,1);
     allRR =  ExcitRR.(an{a}).StimRR(:,1);
     BigRR_stim1=[BigRR_stim1;allRR];
     allmean=mean(allRR);
     BigRR_stim1mean=[BigRR_stim1mean;allmean];
    
end


%%
BigRR_stimend=[];
BigRR_stimendmean=[];

an=fieldnames(ExcitRR_end);
for a=1:length(an)
     %allRR =  ExcitRR_end.(an{a}).StimRREdit(:,1);
     allRR =  ExcitRR_end.(an{a}).StimRR(:,1);
     BigRR_stimend=[BigRR_stimend;allRR];
     allmean=mean(allRR);
     BigRR_stimendmean=[BigRR_stimendmean;allmean];
    
end
%% 
num_bins_ctrl1 = ceil(sqrt(length(BigRR_ctrl1)));
num_bins_ctrlend=ceil(sqrt(length(BigRR_ctrlend)));

num_bins_stim1 = ceil(sqrt(length(BigRR_stim1)));
num_bins_stimend=ceil(sqrt(length(BigRR_stimend)));


ctrl1_mean=mean(BigRR_ctrl1mean); 
ctrlend_mean=mean(BigRR_ctrlendmean); 
mean1=ctrl1_mean; 
mean2=ctrlend_mean;

stim1_mean=mean(BigRR_stim1mean); 
stimend_mean=mean(BigRR_stimendmean); 
mean3=stim1_mean; 
mean4=stimend_mean;


%%
% ctrl1color=rgb('DarkOrange')
% ctrlendcolor=rgb('MediumPurple');
figure; hold on
histogram(BigRR_ctrl1,'BinWidth',0.02,'Normalization','Probability','EdgeColor','none','FaceColor',ctrl1color)
hold on
histogram(BigRR_ctrlend,'BinWidth',0.02,'Normalization','Probability','EdgeColor','none','FaceColor',ctrlendcolor)
hold on
xline(mean1, 'Color', ctrl1color, 'LineWidth', 2);
hold on
xline(mean2, 'Color', ctrlendcolor, 'LineWidth', 2);
xlabel('Time from Stim Onset (s)')
ylabel('Probability')
title('Ctrl Reaches- First vs. Last)')

figure; hold on
histogram(BigRR_ctrl1,num_bins_ctrl1,'Normalization','Probability','EdgeColor','none','FaceColor',ctrl1color)
hold on
histogram(BigRR_ctrlend,num_bins_ctrlend,'Normalization','Probability','EdgeColor','none','FaceColor',ctrlendcolor)
hold on
xline(mean1, 'Color', ctrl1color, 'LineWidth', 2);
hold on
xline(mean2, 'Color', ctrlendcolor, 'LineWidth', 2);
xlabel('Time from Stim Onset (s)')
ylabel('Probability')


figure; hold on
histogram(BigRR_stim1,num_bins_stim1,'Normalization','Probability','EdgeColor','none','FaceColor',stim1color)
hold on
histogram(BigRR_stimend,num_bins_stimend,'Normalization','Probability','EdgeColor','none','FaceColor',stimendcolor)
hold on
xline(mean3, 'Color', stim1color, 'LineWidth', 2);
hold on
xline(mean4, 'Color', stimendcolor, 'LineWidth', 2);

xlabel('Time from Stim Onset (s)')
ylabel('Probability')

figure; hold on
histogram(BigRR_stim1,'BinWidth',0.02,'Normalization','Probability','EdgeColor','none','FaceColor',stim1color)
hold on
histogram(BigRR_stimend,'BinWidth',0.02,'Normalization','Probability','EdgeColor','none','FaceColor',stimendcolor)
hold on
xline(mean3, 'Color', stim1color, 'LineWidth', 2);
hold on
xline(mean4, 'Color', stimendcolor, 'LineWidth', 2);

xlabel('Time from Stim Onset (s)')
ylabel('Probability')
title('Stim Reaches- First vs. Last)')
%% 
%% Separate out CSMs based off of timing of the CSM 
thresh=0.16;

dirdir= 1;
an=fieldnames(ExcitRR_Ctrl1);
for a=1:length(an)
    Rreaches=ExcitRR_Ctrl1.(an{a}).ReReaches;
    rrt=ExcitRR_Ctrl1.(an{a}).StimRR; rrt=rrt(:,1);
    
    CtrlReaches=ExcitRR_Ctrl1.(an{a}).CtrlReaches;

 
    % Timings taken from histo from all animals 
    T0Reaches=[];
    T1Reaches=[];
    T2Reaches=[];
    T3Reaches=[];
  

    % new times based off of bin edges

    % T1=find(rrt>0 & rrt<=0.0588);
    % T2=find(rrt>=0.0662 & rrt<=0.1254);
    % T3=find(rrt>0.1328 & rrt<=0.2);
    
    % Bins for RR times 
    T0=find(rrt==0);
    T1=find(rrt>0 & rrt<=0.0588);
    T2=find(rrt>0.0588 & rrt<0.1328);
    T3=find(rrt>0.1328 & rrt<=0.2);

    rnum=fieldnames(Rreaches);
    
    % for all reaches with a 0 time 
   for s=1:length(T0)
        idx0=T0(s);
        T0Reaches.(rnum{idx0})= ExcitRR_Ctrl1.(an{a}).ReReaches.(rnum{idx0});  
    end 
    T0Kin=BasicKinematics(T0Reaches);
   % ST0 = AlignThresh(T0Reaches,thresh,dirdir,500);
    ST0 = AlignThresh(T0Reaches,thresh,dirdir,200) ; % for heatmap
    %ST0 = AlignThresh_speed(T0Reaches,thresh,dirdir,500) ; % for heatmap
    ST0=stimvel_accel(ST0);
    ST0V = AlignPeakV(T0Reaches,7,500);
    ST0V=stimvel_accel(ST0V);
    ST0_CSM=AlignThreshCSM_T0(T0Reaches,thresh,dirdir,200);
    ST0_CSM=stimvel_accel(ST0_CSM);

    % Ctrl to T0
    CST0_CSM=AlignThreshCSM_T0(CtrlReaches,thresh,dirdir,500);
    CST0_CSM=stimvel_accel(CST0_CSM);

    for f=1:length(T1)
        idx1=T1(f);
        T1Reaches.(rnum{idx1})= ExcitRR_Ctrl1.(an{a}).ReReaches.(rnum{idx1});  
    end 

  
    if isempty(T1Reaches)
        T1Kin=[];
        ST1=[];
    else
       
         T1Kin=BasicKinematics(T1Reaches);
         %ST1 = AlignThresh(T1Reaches,thresh,dirdir,500);
         ST1 = AlignThresh(T1Reaches,thresh,dirdir,200);
         %ST1 = AlignThresh_speed(T1Reaches,thresh,dirdir,500);
         ST1=stimvel_accel(ST1);
    end 
    ST1V = AlignPeakV(T1Reaches,7,500);
    ST1V=stimvel_accel(ST1V);
    ST1_CSM=AlignThreshCSM_T1(T1Reaches,thresh,dirdir,200);
    ST1_CSM=stimvel_accel(ST1_CSM);


    % Ctrl to T1
    CST1_CSM=AlignThreshCSM_T1(CtrlReaches,thresh,dirdir,500);
    CST1_CSM=stimvel_accel(CST1_CSM);

    for h=1:length(T2)
        idx2=T2(h);
        T2Reaches.(rnum{idx2})= ExcitRR_Ctrl1.(an{a}).ReReaches.(rnum{idx2});  
    end 
    
    T2Kin=BasicKinematics(T2Reaches);
    %ST2 = AlignThresh(T2Reaches,thresh,dirdir,500);
    ST2 = AlignThresh(T2Reaches,thresh,dirdir,200);
   % ST2 = AlignThresh_speed(T2Reaches,thresh,dirdir,500);
    ST2=stimvel_accel(ST2);
    ST2V = AlignPeakV(T2Reaches,7,500);
    ST2V=stimvel_accel(ST2V);
    ST2_CSM=AlignThreshCSM_T2(T2Reaches,thresh,dirdir,200);
    ST2_CSM=stimvel_accel(ST2_CSM);


    % Ctrl to T2
    CST2_CSM=AlignThreshCSM_T2(CtrlReaches,thresh,dirdir,500);
    CST2_CSM=stimvel_accel(CST2_CSM);

    for k=1:length(T3)
        idx3=T3(k);
        T3Reaches.(rnum{idx3})= ExcitRR_Ctrl1.(an{a}).ReReaches.(rnum{idx3});  
    end 
    
    if isempty(T3Reaches)
        T3Kin=[];
        ST3=[];
        ST3V=[];
        ST3_CSM=[];
    else
       
         T3Kin=BasicKinematics(T3Reaches);
         %ST1 = AlignThresh(T1Reaches,thresh,dirdir,500);
         ST3 = AlignThresh(T3Reaches,thresh,dirdir,200);
         %ST1 = AlignThresh_speed(T1Reaches,thresh,dirdir,500);
         ST3=stimvel_accel(ST3);
         ST3V = AlignPeakV(T3Reaches,7,500);
        ST3V=stimvel_accel(ST3V);
        ST3_CSM=AlignThreshCSM_T3(T3Reaches,thresh,dirdir,200);
        ST3_CSM=stimvel_accel(ST3_CSM);
    end 

    

    % Ctrl to T0
    CST3_CSM=AlignThreshCSM_T3(CtrlReaches,thresh,dirdir,500);
    CST3_CSM=stimvel_accel(CST3_CSM);

    % Reaches 
     ExcitRR_Ctrl1.(an{a}).RRT0=T0Reaches;
     ExcitRR_Ctrl1.(an{a}).RRT1=T1Reaches;
     ExcitRR_Ctrl1.(an{a}).RRT2=T2Reaches;
     ExcitRR_Ctrl1.(an{a}).RRT3=T3Reaches;
    
     NoReReach=ExcitRR_Ctrl1.(an{a}).NoReReach;
    nocsmfields=fieldnames(NoReReach);
    t0fields=fieldnames(T0Reaches);
    t1fields=fieldnames(T1Reaches);
    t2fields=fieldnames(T2Reaches);

    nocsmnum=length(nocsmfields);
    rrnum0=length(t0fields);
    rrnum1=length(t1fields);
    rrnum2=length(t2fields);
    

    if isempty(T3Reaches)
        t3fields=[];
        rrnum3=[];
    else
        t3fields=fieldnames(T3Reaches);
        rrnum3=length(t3fields);
    end

    ExcitRR_Ctrl1.(an{a}).nocsmnum=nocsmnum;
    ExcitRR_Ctrl1.(an{a}).rrnum0=rrnum0;
    ExcitRR_Ctrl1.(an{a}).rrnum1=rrnum1;
    ExcitRR_Ctrl1.(an{a}).rrnum2=rrnum2;
    ExcitRR_Ctrl1.(an{a}).rrnum3=rrnum3;

     % Kinematics
     ExcitRR_Ctrl1.(an{a}).T0Kin=T0Kin;
     ExcitRR_Ctrl1.(an{a}).T1Kin=T1Kin;
     ExcitRR_Ctrl1.(an{a}).T2Kin=T2Kin;
     ExcitRR_Ctrl1.(an{a}).T3Kin=T3Kin;


     %Stimvel
     ExcitRR_Ctrl1.(an{a}).ST0=ST0;
     ExcitRR_Ctrl1.(an{a}).ST1=ST1;
     ExcitRR_Ctrl1.(an{a}).ST2=ST2;
     ExcitRR_Ctrl1.(an{a}).ST3=ST3;

     ExcitRR_Ctrl1.(an{a}).ST0V=ST0V;
     ExcitRR_Ctrl1.(an{a}).ST1V=ST1V;
     ExcitRR_Ctrl1.(an{a}).ST2V=ST2V;
     ExcitRR_Ctrl1.(an{a}).ST3V=ST3V;

     ExcitRR_Ctrl1.(an{a}).ST0_CSM=ST0_CSM;
     ExcitRR_Ctrl1.(an{a}).ST1_CSM=ST1_CSM;
     ExcitRR_Ctrl1.(an{a}).ST2_CSM=ST2_CSM;
     ExcitRR_Ctrl1.(an{a}).ST3_CSM=ST3_CSM;

     ExcitRR_Ctrl1.(an{a}).CST0_CSM=CST0_CSM;
     ExcitRR_Ctrl1.(an{a}).CST1_CSM=CST1_CSM;
     ExcitRR_Ctrl1.(an{a}).CST2_CSM=CST2_CSM;
     ExcitRR_Ctrl1.(an{a}).CST3_CSM=CST3_CSM;
     % % Add NRR CSM alignment 
     % NRRreaches=ExcitRR_Ctrl.(an{a}).NoReReach;
     % NRR_CSM=AlignThreshCSM_T4(NRRreaches,thresh,dirdir,500);
     % NRR_CSM=stimvel_accel(NRR_CSM);
     % ExcitRR_Ctrl.(an{a}).NRR_CSM=NRR_CSM;
         
    
    %  % success rate 
    %  [S_T0, F_T0] = SuccessFilter(T0Reaches);
    %  [S_T1, F_T1] = SuccessFilter(T1Reaches);
    %  [S_T2, F_T2] = SuccessFilter(T2Reaches);
    %  %[S_T3, F_T3] = SuccessFilter(T3Reaches);
    % 
    %  ExcitRR_Ctrl1.(an{a}).successrate.T0 = length(fieldnames(S_T0)) / length(fieldnames(T0Reaches));
    %  ExcitRR_Ctrl1.(an{a}).successrate.T1 = length(fieldnames(S_T1)) / length(fieldnames(T1Reaches));
    %  ExcitRR_Ctrl1.(an{a}).successrate.T2 = length(fieldnames(S_T2)) / length(fieldnames(T2Reaches));
    %  %ExcitRR_Ctrl1.(an{a}).successrate.T3 = length(fieldnames(S_T3)) / length(fieldnames(T3Reaches));
    % 
    %    % Initialize variables
    % numReaches = zeros(1, 4);
    % successRates = zeros(1, 4);
    % 
    % % Check if each category exists before calculating success rates
    % if ~isempty(T0Reaches)
    %     numReaches(1) = length(fieldnames(T0Reaches));
    %     successRates(1) = length(fieldnames(S_T0)) / numReaches(1);
    % else
    %     numReaches(1) = 0; % No reaches in T0
    %     successRates(1) = NaN; % Set as NaN so it doesn't affect weighted avg
    % end
    % 
    % if ~isempty(T1Reaches)
    %     numReaches(2) = length(fieldnames(T1Reaches));
    %     successRates(2) = length(fieldnames(S_T1)) / numReaches(2);
    % else
    %     numReaches(2) = 0;
    %     successRates(2) = NaN;
    % end
    % 
    % if ~isempty(T2Reaches)
    %     numReaches(3) = length(fieldnames(T2Reaches));
    %     successRates(3) = length(fieldnames(S_T2)) / numReaches(3);
    % else
    %     numReaches(3) = 0;
    %     successRates(3) = NaN;
    % end
    % 
    % if ~isempty(T3Reaches)
    %     numReaches(4) = length(fieldnames(T3Reaches));
    %     successRates(4) = length(fieldnames(S_T3)) / numReaches(4);
    % else
    %     numReaches(4) = 0;
    %     successRates(4) = NaN;
    % end
    % 
    % 
    % % Store results
    % ExcitRR_Ctrl1.(an{a}).successrate.T0 = successRates(1);
    % ExcitRR_Ctrl1.(an{a}).successrate.T1 = successRates(2);
    % ExcitRR_Ctrl1.(an{a}).successrate.T2 = successRates(3);
    % ExcitRR_Ctrl1.(an{a}).successrate.T3 = successRates(4);
    % 
    % % Compute total reaches (excluding empty categories)
    % totalReaches = sum(numReaches);
    % 
    % % Compute weighted success rates per category
    % weightedSuccessRates = (successRates .* numReaches) / totalReaches;
    % 
    % % Handle NaN values (set to 0 for visualization)
    % weightedSuccessRates(isnan(weightedSuccessRates)) = 0;
    % 
    % % Store results
    % ExcitRR_Ctrl1.(an{a}).successrate.T0_weighted = weightedSuccessRates(1);
    % ExcitRR_Ctrl1.(an{a}).successrate.T1_weighted = weightedSuccessRates(2);
    % ExcitRR_Ctrl1.(an{a}).successrate.T2_weighted = weightedSuccessRates(3);
    % ExcitRR_Ctrl1.(an{a}).successrate.T3_weighted = weightedSuccessRates(4);

end 

% ST1, ST4, ST5 - not all animals have these 

%% Separate out CSMs based off of timing of the CSM - Last sessions
thresh=0.16;

dirdir= 1;
an=fieldnames(ExcitRR_Ctrlend);
for a=1:length(an)
    Rreaches=ExcitRR_Ctrlend.(an{a}).ReReaches;
    rrt=ExcitRR_Ctrlend.(an{a}).StimRR; rrt=rrt(:,1);
    
    CtrlReaches=ExcitRR_Ctrlend.(an{a}).CtrlReaches;

 
    % Timings taken from histo from all animals 
    T0Reaches=[];
    T1Reaches=[];
    T2Reaches=[];
    T3Reaches=[];
  

    % new times based off of bin edges

    % T1=find(rrt>0 & rrt<=0.0588);
    % T2=find(rrt>=0.0662 & rrt<=0.1254);
    % T3=find(rrt>0.1328 & rrt<=0.2);
    
    % Bins for RR times 
    T0=find(rrt==0);
    T1=find(rrt>0 & rrt<=0.0588);
    T2=find(rrt>0.0588 & rrt<0.1328);
    T3=find(rrt>0.1328 & rrt<=0.2);

    rnum=fieldnames(Rreaches);
    
    % for all reaches with a 0 time 
   for s=1:length(T0)
        idx0=T0(s);
        T0Reaches.(rnum{idx0})= ExcitRR_Ctrlend.(an{a}).ReReaches.(rnum{idx0});  
   end 

    if isempty(T0Reaches)
        T0Kin=[];
        ST0=[];
        ST0V=[];
        ST0_CSM=[];
    else
        T0Kin=BasicKinematics(T0Reaches);
       % ST0 = AlignThresh(T0Reaches,thresh,dirdir,500);
        ST0 = AlignThresh(T0Reaches,thresh,dirdir,200) ; % for heatmap
        %ST0 = AlignThresh_speed(T0Reaches,thresh,dirdir,500) ; % for heatmap
        ST0=stimvel_accel(ST0);
        ST0V = AlignPeakV(T0Reaches,7,500);
        ST0V=stimvel_accel(ST0V);
        ST0_CSM=AlignThreshCSM_T0(T0Reaches,thresh,dirdir,200);
        ST0_CSM=stimvel_accel(ST0_CSM);
    end


    % Ctrl to T0
    CST0_CSM=AlignThreshCSM_T0(CtrlReaches,thresh,dirdir,500);
    CST0_CSM=stimvel_accel(CST0_CSM);

    for f=1:length(T1)
        idx1=T1(f);
        T1Reaches.(rnum{idx1})= ExcitRR_Ctrlend.(an{a}).ReReaches.(rnum{idx1});  
    end 

  
    if isempty(T1Reaches)
        T1Kin=[];
        ST1=[];
    else
       
         T1Kin=BasicKinematics(T1Reaches);
         %ST1 = AlignThresh(T1Reaches,thresh,dirdir,500);
         ST1 = AlignThresh(T1Reaches,thresh,dirdir,200);
         %ST1 = AlignThresh_speed(T1Reaches,thresh,dirdir,500);
         ST1=stimvel_accel(ST1);
    end 
    ST1V = AlignPeakV(T1Reaches,7,500);
    ST1V=stimvel_accel(ST1V);
    ST1_CSM=AlignThreshCSM_T1(T1Reaches,thresh,dirdir,200);
    ST1_CSM=stimvel_accel(ST1_CSM);


    % Ctrl to T1
    CST1_CSM=AlignThreshCSM_T1(CtrlReaches,thresh,dirdir,500);
    CST1_CSM=stimvel_accel(CST1_CSM);

    for h=1:length(T2)
        idx2=T2(h);
        T2Reaches.(rnum{idx2})= ExcitRR_Ctrlend.(an{a}).ReReaches.(rnum{idx2});  
    end 
    T2Kin=BasicKinematics(T2Reaches);
    %ST2 = AlignThresh(T2Reaches,thresh,dirdir,500);
    ST2 = AlignThresh(T2Reaches,thresh,dirdir,200);
   % ST2 = AlignThresh_speed(T2Reaches,thresh,dirdir,500);
    ST2=stimvel_accel(ST2);
    ST2V = AlignPeakV(T2Reaches,7,500);
    ST2V=stimvel_accel(ST2V);
    ST2_CSM=AlignThreshCSM_T2(T2Reaches,thresh,dirdir,200);
    ST2_CSM=stimvel_accel(ST2_CSM);


    % Ctrl to T2
    CST2_CSM=AlignThreshCSM_T2(CtrlReaches,thresh,dirdir,500);
    CST2_CSM=stimvel_accel(CST2_CSM);

    for k=1:length(T3)
        idx3=T3(k);
        T3Reaches.(rnum{idx3})= ExcitRR_Ctrlend.(an{a}).ReReaches.(rnum{idx3});  
    end 
    
    T3Kin=BasicKinematics(T3Reaches);
    %ST3 = AlignThresh(T3Reaches,thresh,dirdir,500);
    ST3 = AlignThresh(T3Reaches,thresh,dirdir,200);
    %ST3 = AlignThresh_speed(T3Reaches,thresh,dirdir,500);
    ST3=stimvel_accel(ST3);
    ST3V = AlignPeakV(T3Reaches,7,500);
    ST3V=stimvel_accel(ST3V);
    ST3_CSM=AlignThreshCSM_T3(T3Reaches,thresh,dirdir,200);
    ST3_CSM=stimvel_accel(ST3_CSM);
    

    % Ctrl to T0
    CST3_CSM=AlignThreshCSM_T3(CtrlReaches,thresh,dirdir,500);
    CST3_CSM=stimvel_accel(CST3_CSM);

    % Reaches 
     ExcitRR_Ctrlend.(an{a}).RRT0=T0Reaches;
     ExcitRR_Ctrlend.(an{a}).RRT1=T1Reaches;
     ExcitRR_Ctrlend.(an{a}).RRT2=T2Reaches;
     ExcitRR_Ctrlend.(an{a}).RRT3=T3Reaches;
    

    if isempty(T0Reaches)
        t0fields=[];
        rrnum0=[];
    else
        t0fields=fieldnames(T0Reaches);
        rrnum0=length(t0fields);
    end
    
    NoReReach=ExcitRR_Ctrlend.(an{a}).NoReReach;
    nocsmfields=fieldnames(NoReReach);
    t1fields=fieldnames(T1Reaches);
    t2fields=fieldnames(T2Reaches);
    t3fields=fieldnames(T3Reaches);
    
    nocsmnum=length(nocsmfields);
    rrnum1=length(t1fields);
    rrnum2=length(t2fields);
    rrnum3=length(t3fields);
   
    ExcitRR_Ctrlend.(an{a}).nocsmnum=nocsmnum;
    ExcitRR_Ctrlend.(an{a}).rrnum0=rrnum0;
    ExcitRR_Ctrlend.(an{a}).rrnum1=rrnum1;
    ExcitRR_Ctrlend.(an{a}).rrnum2=rrnum2;
    ExcitRR_Ctrlend.(an{a}).rrnum3=rrnum3;

     % Kinematics
     ExcitRR_Ctrlend.(an{a}).T0Kin=T0Kin;
     ExcitRR_Ctrlend.(an{a}).T1Kin=T1Kin;
     ExcitRR_Ctrlend.(an{a}).T2Kin=T2Kin;
     ExcitRR_Ctrlend.(an{a}).T3Kin=T3Kin;


     %Stimvel
     ExcitRR_Ctrlend.(an{a}).ST0=ST0;
     ExcitRR_Ctrlend.(an{a}).ST1=ST1;
     ExcitRR_Ctrlend.(an{a}).ST2=ST2;
     ExcitRR_Ctrlend.(an{a}).ST3=ST3;

     ExcitRR_Ctrlend.(an{a}).ST0V=ST0V;
     ExcitRR_Ctrlend.(an{a}).ST1V=ST1V;
     ExcitRR_Ctrlend.(an{a}).ST2V=ST2V;
     ExcitRR_Ctrlend.(an{a}).ST3V=ST3V;

     ExcitRR_Ctrlend.(an{a}).ST0_CSM=ST0_CSM;
     ExcitRR_Ctrlend.(an{a}).ST1_CSM=ST1_CSM;
     ExcitRR_Ctrlend.(an{a}).ST2_CSM=ST2_CSM;
     ExcitRR_Ctrlend.(an{a}).ST3_CSM=ST3_CSM;

     ExcitRR_Ctrlend.(an{a}).CST0_CSM=CST0_CSM;
     ExcitRR_Ctrlend.(an{a}).CST1_CSM=CST1_CSM;
     ExcitRR_Ctrlend.(an{a}).CST2_CSM=CST2_CSM;
     ExcitRR_Ctrlend.(an{a}).CST3_CSM=CST3_CSM;
     % % Add NRR CSM alignment 
     % NRRreaches=ExcitRR_Ctrl.(an{a}).NoReReach;
     % NRR_CSM=AlignThreshCSM_T4(NRRreaches,thresh,dirdir,500);
     % NRR_CSM=stimvel_accel(NRR_CSM);
     % ExcitRR_Ctrl.(an{a}).NRR_CSM=NRR_CSM;
     
      if isempty(T0Reaches)
        S_T0=[];
        F_T0=[];
        ExcitRR_Ctrlend.(an{a}).successrate.T0 =[];
    else
        [S_T0, F_T0] = SuccessFilter(T0Reaches); 
        ExcitRR_Ctrlend.(an{a}).successrate.T0 = length(fieldnames(S_T0)) / length(fieldnames(T0Reaches));
    end

     % success rate 
     %[S_T0, F_T0] = SuccessFilter(T0Reaches);
     [S_T1, F_T1] = SuccessFilter(T1Reaches);
     [S_T2, F_T2] = SuccessFilter(T2Reaches);
     [S_T3, F_T3] = SuccessFilter(T3Reaches);
    
    % ExcitRR_Ctrlend.(an{a}).successrate.T0 = length(fieldnames(S_T0)) / length(fieldnames(T0Reaches));
     ExcitRR_Ctrlend.(an{a}).successrate.T1 = length(fieldnames(S_T1)) / length(fieldnames(T1Reaches));
     ExcitRR_Ctrlend.(an{a}).successrate.T2 = length(fieldnames(S_T2)) / length(fieldnames(T2Reaches));
     ExcitRR_Ctrlend.(an{a}).successrate.T3 = length(fieldnames(S_T3)) / length(fieldnames(T3Reaches));
    
        % Initialize variables
    numReaches = zeros(1, 4);
    successRates = zeros(1, 4);
    
    % Check if each category exists before calculating success rates
    if ~isempty(T0Reaches)
        numReaches(1) = length(fieldnames(T0Reaches));
        successRates(1) = length(fieldnames(S_T0)) / numReaches(1);
    else
        numReaches(1) = 0; % No reaches in T0
        successRates(1) = NaN; % Set as NaN so it doesn't affect weighted avg
    end
    
    if ~isempty(T1Reaches)
        numReaches(2) = length(fieldnames(T1Reaches));
        successRates(2) = length(fieldnames(S_T1)) / numReaches(2);
    else
        numReaches(2) = 0;
        successRates(2) = NaN;
    end
    
    if ~isempty(T2Reaches)
        numReaches(3) = length(fieldnames(T2Reaches));
        successRates(3) = length(fieldnames(S_T2)) / numReaches(3);
    else
        numReaches(3) = 0;
        successRates(3) = NaN;
    end
    
    if ~isempty(T3Reaches)
        numReaches(4) = length(fieldnames(T3Reaches));
        successRates(4) = length(fieldnames(S_T3)) / numReaches(4);
    else
        numReaches(4) = 0;
        successRates(4) = NaN;
    end
   
    
    % Store results
    ExcitRR_Ctrlend.(an{a}).successrate.T0 = successRates(1);
    ExcitRR_Ctrlend.(an{a}).successrate.T1 = successRates(2);
    ExcitRR_Ctrlend.(an{a}).successrate.T2 = successRates(3);
    ExcitRR_Ctrlend.(an{a}).successrate.T3 = successRates(4);
  
    % Compute total reaches (excluding empty categories)
    totalReaches = sum(numReaches);
    
    % Compute weighted success rates per category
    weightedSuccessRates = (successRates .* numReaches) / totalReaches;
    
    % Handle NaN values (set to 0 for visualization)
    weightedSuccessRates(isnan(weightedSuccessRates)) = 0;
    
    % Store results
    ExcitRR_Ctrlend.(an{a}).successrate.T0_weighted = weightedSuccessRates(1);
    ExcitRR_Ctrlend.(an{a}).successrate.T1_weighted = weightedSuccessRates(2);
    ExcitRR_Ctrlend.(an{a}).successrate.T2_weighted = weightedSuccessRates(3);
    ExcitRR_Ctrlend.(an{a}).successrate.T3_weighted = weightedSuccessRates(4);

end 

% ST1, ST4, ST5 - not all animals have these 

%% First session stim 
thresh=0.16;

dirdir= 1;
an=fieldnames(ExcitRR);
for a=1:length(an)
    Rreaches=ExcitRR.(an{a}).ReReaches;
    rrt=ExcitRR.(an{a}).StimRR; rrt=rrt(:,1);
    
   % CtrlReaches=ExcitRR.(an{a}).CtrlReaches;

 
    % Timings taken from histo from all animals 
    T0Reaches=[];
    T1Reaches=[];
    T2Reaches=[];
    T3Reaches=[];
  

    % new times based off of bin edges

    % T1=find(rrt>0 & rrt<=0.0588);
    % T2=find(rrt>=0.0662 & rrt<=0.1254);
    % T3=find(rrt>0.1328 & rrt<=0.2);
    
    % Bins for RR times 
    T0=find(rrt==0);
    T1=find(rrt>0 & rrt<=0.0588);
    T2=find(rrt>0.0588 & rrt<0.1328);
    T3=find(rrt>0.1328 & rrt<=0.2);

    rnum=fieldnames(Rreaches);
    
    % for all reaches with a 0 time 
   for s=1:length(T0)
        idx0=T0(s);
        T0Reaches.(rnum{idx0})= ExcitRR.(an{a}).ReReaches.(rnum{idx0});  
   end 

   if isempty(T0Reaches)
        T0Kin=[];
        ST0=[];
        ST0V=[];
        ST0_CSM=[];
    else
        T0Kin=BasicKinematics(T0Reaches);
        ST0 = AlignStim(T0Reaches,200);
        ST0=stimvel_accel(ST0);
        ST0V = AlignPeakV(T0Reaches,7,500);
        ST0V=stimvel_accel(ST0V);
        ST0_CSM=AlignCSM_T0(T0Reaches,200); % might have to change this to the Stim one ?? 
        ST0_CSM=stimvel_accel(ST0_CSM);
   end 
 

    % Ctrl to T0
    % CST0_CSM=AlignThreshCSM_T0(CtrlReaches,thresh,dirdir,500);
    % CST0_CSM=stimvel_accel(CST0_CSM);

    for f=1:length(T1)
        idx1=T1(f);
        T1Reaches.(rnum{idx1})= ExcitRR.(an{a}).ReReaches.(rnum{idx1});  
    end 

  
    if isempty(T1Reaches)
        T1Kin=[];
        ST1=[];
    else
       
         T1Kin=BasicKinematics(T1Reaches);
         ST1=AlignStim(T1Reaches,200);
         ST1=stimvel_accel(ST1);
    end 
    ST1V = AlignPeakV(T1Reaches,7,500);
    ST1V=stimvel_accel(ST1V);
    ST1_CSM=AlignCSM_T1(T1Reaches,200);
    ST1_CSM=stimvel_accel(ST1_CSM);


    % % Ctrl to T1
    % CST1_CSM=AlignThreshCSM_T1(CtrlReaches,thresh,dirdir,500);
    % CST1_CSM=stimvel_accel(CST1_CSM);

    for h=1:length(T2)
        idx2=T2(h);
        T2Reaches.(rnum{idx2})= ExcitRR.(an{a}).ReReaches.(rnum{idx2});  
    end 
    T2Kin=BasicKinematics(T2Reaches);
    
    ST2 = AlignStim(T2Reaches,200);
   
    ST2=stimvel_accel(ST2);
    ST2V = AlignPeakV(T2Reaches,7,500);
    ST2V=stimvel_accel(ST2V);
    ST2_CSM=AlignCSM_T2(T2Reaches,200);
    ST2_CSM=stimvel_accel(ST2_CSM);


    % % Ctrl to T2
    % CST2_CSM=AlignThreshCSM_T2(CtrlReaches,thresh,dirdir,500);
    % CST2_CSM=stimvel_accel(CST2_CSM);

    for k=1:length(T3)
        idx3=T3(k);
        T3Reaches.(rnum{idx3})= ExcitRR.(an{a}).ReReaches.(rnum{idx3});  
    end 
     
    if isempty(T3Reaches)
        T3Kin=[];
        ST3=[];
        ST3V=[];
        ST3_CSM=[];
    else
       
         T3Kin=BasicKinematics(T3Reaches);
         ST3=AlignStim(T3Reaches,200);
         ST3=stimvel_accel(ST3);

        ST3V = AlignPeakV(T3Reaches,7,500);
        ST3V=stimvel_accel(ST3V);
        ST3_CSM=AlignCSM_T3(T3Reaches,200);
        ST3_CSM=stimvel_accel(ST3_CSM);
    
    end 
    

   
    

    % % Ctrl to T0
    % CST3_CSM=AlignThreshCSM_T3(CtrlReaches,thresh,dirdir,500);
    % CST3_CSM=stimvel_accel(CST3_CSM);

    % Reaches 
     ExcitRR.(an{a}).RRT0=T0Reaches;
     ExcitRR.(an{a}).RRT1=T1Reaches;
     ExcitRR.(an{a}).RRT2=T2Reaches;
     ExcitRR.(an{a}).RRT3=T3Reaches;
    
   if isempty(T0Reaches)
        t0fields=[];
        rrnum0=[];
    else
        t0fields=fieldnames(T0Reaches);
        rrnum0=length(t0fields);
   end

   
    t1fields=fieldnames(T1Reaches);
    t2fields=fieldnames(T2Reaches);

    if isempty(T3Reaches)
        t3fields=[];
        rrnum3=[];
    else
        t3fields=fieldnames(T3Reaches);
        rrnum3=length(t3fields);
   end
    
    

    
    rrnum1=length(t1fields);
    rrnum2=length(t2fields);
    
   
    ExcitRR.(an{a}).rrnum0=rrnum0;
    ExcitRR.(an{a}).rrnum1=rrnum1;
    ExcitRR.(an{a}).rrnum2=rrnum2;
    ExcitRR.(an{a}).rrnum3=rrnum3;

     % Kinematics
     ExcitRR.(an{a}).T0Kin=T0Kin;
     ExcitRR.(an{a}).T1Kin=T1Kin;
     ExcitRR.(an{a}).T2Kin=T2Kin;
     ExcitRR.(an{a}).T3Kin=T3Kin;


     %Stimvel
     ExcitRR.(an{a}).ST0=ST0;
     ExcitRR.(an{a}).ST1=ST1;
     ExcitRR.(an{a}).ST2=ST2;
     ExcitRR.(an{a}).ST3=ST3;

     ExcitRR.(an{a}).ST0V=ST0V;
     ExcitRR.(an{a}).ST1V=ST1V;
     ExcitRR.(an{a}).ST2V=ST2V;
     ExcitRR.(an{a}).ST3V=ST3V;

     ExcitRR.(an{a}).ST0_CSM=ST0_CSM;
     ExcitRR.(an{a}).ST1_CSM=ST1_CSM;
     ExcitRR.(an{a}).ST2_CSM=ST2_CSM;
     ExcitRR.(an{a}).ST3_CSM=ST3_CSM;
     % 
     % ExcitRR.(an{a}).CST0_CSM=CST0_CSM;
     % ExcitRR.(an{a}).CST1_CSM=CST1_CSM;
     % ExcitRR.(an{a}).CST2_CSM=CST2_CSM;
     % ExcitRR.(an{a}).CST3_CSM=CST3_CSM;
     % % Add NRR CSM alignment 
     % NRRreaches=ExcitRR_Ctrl.(an{a}).NoReReach;
     % NRR_CSM=AlignThreshCSM_T4(NRRreaches,thresh,dirdir,500);
     % NRR_CSM=stimvel_accel(NRR_CSM);
     % ExcitRR_Ctrl.(an{a}).NRR_CSM=NRR_CSM;
         
    
     % success rate 
     % [S_T0, F_T0] = SuccessFilter(T0Reaches);
     % [S_T1, F_T1] = SuccessFilter(T1Reaches);
     % [S_T2, F_T2] = SuccessFilter(T2Reaches);
     % [S_T3, F_T3] = SuccessFilter(T3Reaches);
     % 
     % ExcitRR.(an{a}).successrate.T0 = length(fieldnames(S_T0)) / length(fieldnames(T0Reaches));
     % ExcitRR.(an{a}).successrate.T1 = length(fieldnames(S_T1)) / length(fieldnames(T1Reaches));
     % ExcitRR.(an{a}).successrate.T2 = length(fieldnames(S_T2)) / length(fieldnames(T2Reaches));
     % ExcitRR.(an{a}).successrate.T3 = length(fieldnames(S_T3)) / length(fieldnames(T3Reaches));

end 

% ST1, ST4, ST5 - not all animals have these 

%% LAST session stim 
thresh=0.16;

dirdir= 1;
an=fieldnames(ExcitRR_end);
for a=1:length(an)
    Rreaches=ExcitRR_end.(an{a}).ReReaches;
    rrt=ExcitRR_end.(an{a}).StimRR; rrt=rrt(:,1);
    
   % CtrlReaches=ExcitRR.(an{a}).CtrlReaches;

 
    % Timings taken from histo from all animals 
    T0Reaches=[];
    T1Reaches=[];
    T2Reaches=[];
    T3Reaches=[];
  

    % new times based off of bin edges

    % T1=find(rrt>0 & rrt<=0.0588);
    % T2=find(rrt>=0.0662 & rrt<=0.1254);
    % T3=find(rrt>0.1328 & rrt<=0.2);
    
    % Bins for RR times 
    T0=find(rrt==0);
    T1=find(rrt>0 & rrt<=0.0588);
    T2=find(rrt>0.0588 & rrt<0.1328);
    T3=find(rrt>0.1328 & rrt<=0.2);

    rnum=fieldnames(Rreaches);
    
    % for all reaches with a 0 time 
   for s=1:length(T0)
        idx0=T0(s);
        T0Reaches.(rnum{idx0})= ExcitRR_end.(an{a}).ReReaches.(rnum{idx0});  
   end 

   if isempty(T0Reaches)
        T0Kin=[];
        ST0=[];
        ST0V=[];
        ST0_CSM=[];
    else
        T0Kin=BasicKinematics(T0Reaches);
        ST0 = AlignStim(T0Reaches,200);
        ST0=stimvel_accel(ST0);
        ST0V = AlignPeakV(T0Reaches,7,500);
        ST0V=stimvel_accel(ST0V);
        ST0_CSM=AlignCSM_T0(T0Reaches,200); % might have to change this to the Stim one ?? 
        ST0_CSM=stimvel_accel(ST0_CSM);
   end 
 

    % Ctrl to T0
    % CST0_CSM=AlignThreshCSM_T0(CtrlReaches,thresh,dirdir,500);
    % CST0_CSM=stimvel_accel(CST0_CSM);

    for f=1:length(T1)
        idx1=T1(f);
        T1Reaches.(rnum{idx1})= ExcitRR_end.(an{a}).ReReaches.(rnum{idx1});  
    end 

  
    if isempty(T1Reaches)
        T1Kin=[];
        ST1=[];
    else
       
         T1Kin=BasicKinematics(T1Reaches);
         ST1=AlignStim(T1Reaches,200);
         ST1=stimvel_accel(ST1);
    end 
    ST1V = AlignPeakV(T1Reaches,7,500);
    ST1V=stimvel_accel(ST1V);
    ST1_CSM=AlignCSM_T1(T1Reaches,200);
    ST1_CSM=stimvel_accel(ST1_CSM);


    % % Ctrl to T1
    % CST1_CSM=AlignThreshCSM_T1(CtrlReaches,thresh,dirdir,500);
    % CST1_CSM=stimvel_accel(CST1_CSM);

    for h=1:length(T2)
        idx2=T2(h);
        T2Reaches.(rnum{idx2})= ExcitRR_end.(an{a}).ReReaches.(rnum{idx2});  
    end 
    T2Kin=BasicKinematics(T2Reaches);
    
    ST2 = AlignStim(T2Reaches,200);
   
    ST2=stimvel_accel(ST2);
    ST2V = AlignPeakV(T2Reaches,7,500);
    ST2V=stimvel_accel(ST2V);
    ST2_CSM=AlignCSM_T2(T2Reaches,200);
    ST2_CSM=stimvel_accel(ST2_CSM);


    % % Ctrl to T2
    % CST2_CSM=AlignThreshCSM_T2(CtrlReaches,thresh,dirdir,500);
    % CST2_CSM=stimvel_accel(CST2_CSM);

    for k=1:length(T3)
        idx3=T3(k);
        T3Reaches.(rnum{idx3})= ExcitRR_end.(an{a}).ReReaches.(rnum{idx3});  
    end 
     
    if isempty(T3Reaches)
        T3Kin=[];
        ST3=[];
        ST3V=[];
        ST3_CSM=[];
    else
       
         T3Kin=BasicKinematics(T3Reaches);
         ST3=AlignStim(T3Reaches,200);
         ST3=stimvel_accel(ST3);

        ST3V = AlignPeakV(T3Reaches,7,500);
        ST3V=stimvel_accel(ST3V);
        ST3_CSM=AlignCSM_T3(T3Reaches,200);
        ST3_CSM=stimvel_accel(ST3_CSM);
    
    end 
    

   
    

    % % Ctrl to T0
    % CST3_CSM=AlignThreshCSM_T3(CtrlReaches,thresh,dirdir,500);
    % CST3_CSM=stimvel_accel(CST3_CSM);

    % Reaches 
     ExcitRR_end.(an{a}).RRT0=T0Reaches;
     ExcitRR_end.(an{a}).RRT1=T1Reaches;
     ExcitRR_end.(an{a}).RRT2=T2Reaches;
     ExcitRR_end.(an{a}).RRT3=T3Reaches;
    
   if isempty(T0Reaches)
        t0fields=[];
        rrnum0=[];
    else
        t0fields=fieldnames(T0Reaches);
        rrnum0=length(t0fields);
   end

   
    t1fields=fieldnames(T1Reaches);
    t2fields=fieldnames(T2Reaches);

    if isempty(T3Reaches)
        t3fields=[];
        rrnum3=[];
    else
        t3fields=fieldnames(T3Reaches);
        rrnum3=length(t3fields);
   end
    
    

    
    rrnum1=length(t1fields);
    rrnum2=length(t2fields);
    
   
    ExcitRR_end.(an{a}).rrnum0=rrnum0;
    ExcitRR_end.(an{a}).rrnum1=rrnum1;
    ExcitRR_end.(an{a}).rrnum2=rrnum2;
    ExcitRR_end.(an{a}).rrnum3=rrnum3;

     % Kinematics
     ExcitRR_end.(an{a}).T0Kin=T0Kin;
     ExcitRR_end.(an{a}).T1Kin=T1Kin;
     ExcitRR_end.(an{a}).T2Kin=T2Kin;
     ExcitRR_end.(an{a}).T3Kin=T3Kin;


     %Stimvel
     ExcitRR_end.(an{a}).ST0=ST0;
     ExcitRR_end.(an{a}).ST1=ST1;
     ExcitRR_end.(an{a}).ST2=ST2;
     ExcitRR_end.(an{a}).ST3=ST3;

     ExcitRR_end.(an{a}).ST0V=ST0V;
     ExcitRR_end.(an{a}).ST1V=ST1V;
     ExcitRR_end.(an{a}).ST2V=ST2V;
     ExcitRR_end.(an{a}).ST3V=ST3V;

     ExcitRR_end.(an{a}).ST0_CSM=ST0_CSM;
     ExcitRR_end.(an{a}).ST1_CSM=ST1_CSM;
     ExcitRR_end.(an{a}).ST2_CSM=ST2_CSM;
     ExcitRR_end.(an{a}).ST3_CSM=ST3_CSM;
     % 
     % ExcitRR.(an{a}).CST0_CSM=CST0_CSM;
     % ExcitRR.(an{a}).CST1_CSM=CST1_CSM;
     % ExcitRR.(an{a}).CST2_CSM=CST2_CSM;
     % ExcitRR.(an{a}).CST3_CSM=CST3_CSM;
     % % Add NRR CSM alignment 
     % NRRreaches=ExcitRR_Ctrl.(an{a}).NoReReach;
     % NRR_CSM=AlignThreshCSM_T4(NRRreaches,thresh,dirdir,500);
     % NRR_CSM=stimvel_accel(NRR_CSM);
     % ExcitRR_Ctrl.(an{a}).NRR_CSM=NRR_CSM;
         
    
     % success rate 
     % [S_T0, F_T0] = SuccessFilter(T0Reaches);
     % [S_T1, F_T1] = SuccessFilter(T1Reaches);
     % [S_T2, F_T2] = SuccessFilter(T2Reaches);
     % [S_T3, F_T3] = SuccessFilter(T3Reaches);
     % 
     % ExcitRR.(an{a}).successrate.T0 = length(fieldnames(S_T0)) / length(fieldnames(T0Reaches));
     % ExcitRR.(an{a}).successrate.T1 = length(fieldnames(S_T1)) / length(fieldnames(T1Reaches));
     % ExcitRR.(an{a}).successrate.T2 = length(fieldnames(S_T2)) / length(fieldnames(T2Reaches));
     % ExcitRR.(an{a}).successrate.T3 = length(fieldnames(S_T3)) / length(fieldnames(T3Reaches));

end 

% ST1, ST4, ST5 - not all animals have these 
 %% Combine Stimvel and look at plots 
C_big=struct;
T0_big=struct;
T1_big=struct;
T2_big = struct;
T3_big=struct;

ar=fieldnames(ExcitRR_Ctrl1);
for a = 1:length(ar) 
   C_big.stimvel(:,:,a) = ExcitRR_Ctrl1.(ar{a}).C.Average;
   T0_big.stimvel(:,:,a) = ExcitRR_Ctrl1.(ar{a}).ST0.Average;
   T1_big.stimvel(:,:,a) = ExcitRR_Ctrl1.(ar{a}).ST1.Average;
   T2_big.stimvel(:,:,a) = ExcitRR_Ctrl1.(ar{a}).ST2.Average;
   T3_big.stimvel(:,:,a) = ExcitRR_Ctrl1.(ar{a}).ST3.Average;
 
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
C_bigend=struct;
T0_bigend=struct;
T1_bigend=struct;
T2_bigend = struct;
T3_bigend=struct;

ar=fieldnames(ExcitRR_Ctrlend);
for a = 1:length(ar) 
   C_bigend.stimvel(:,:,a) = ExcitRR_Ctrlend.(ar{a}).C.Average;
   %T0_bigend.stimvel(:,:,a) = ExcitRR_Ctrlend.(ar{a}).ST0.Average;
   T1_bigend.stimvel(:,:,a) = ExcitRR_Ctrlend.(ar{a}).ST1.Average;
   T2_bigend.stimvel(:,:,a) = ExcitRR_Ctrlend.(ar{a}).ST2.Average;
   T3_bigend.stimvel(:,:,a) = ExcitRR_Ctrlend.(ar{a}).ST3.Average;
 
end

C_bigend.Average = nanmean(C_bigend.stimvel,3);
C_bigend.Std = nanstd(C_bigend.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(C_bigend.Average),1); 
for nnum = 1:length(C_bigend.Average)
    samp = squeeze(C_bigend.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
C_bigend.N = N;



% T0_bigend.Average = nanmean(T0_bigend.stimvel,3);
% T0_bigend.Std = nanstd(T0_bigend.stimvel,0,3);
% 
% %loop through to count the number of reaches contributing to each time pt
% N = zeros(length(T0_bigend.Average),1); 
% for nnum = 1:length(T0_bigend.Average)
%     samp = squeeze(T0_bigend.stimvel(nnum,2,:));
%     N(nnum) = length(samp(~isnan(samp)));
% 
% end
% T0_bigend.N = N;

T1_bigend.Average = nanmean(T1_bigend.stimvel,3);
T1_bigend.Std = nanstd(T1_bigend.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T1_bigend.Average),1); 
for nnum = 1:length(T1_bigend.Average)
    samp = squeeze(T1_bigend.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
T1_bigend.N = N;


T2_bigend.Average = nanmean(T2_bigend.stimvel,3);
T2_bigend.Std = nanstd(T2_bigend.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T2_bigend.Average),1); 
for nnum = 1:length(T2_bigend.Average)
    samp = squeeze(T2_bigend.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
T2_bigend.N = N;

%
T3_bigend.Average = nanmean(T3_bigend.stimvel,3);
T3_bigend.Std = nanstd(T3_bigend.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T3_bigend.Average),1); 
for nnum = 1:length(T3_bigend.Average)
    samp = squeeze(T3_bigend.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
T3_bigend.N = N;

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
 
   T0_big.stimvel(:,:,a) = NoT0Stimvel.(at{a}).ST0.Average;
  
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
  %% Combine Stimvel and look at plots 

T1_bigstim=struct;
T2_bigstim = struct;
S_big=struct;


ar=fieldnames(ExcitRR);
for a = 1:length(ar) 
   T1_bigstim.stimvel(:,:,a) = ExcitRR.(ar{a}).ST1.Average;
   T2_bigstim.stimvel(:,:,a) = ExcitRR.(ar{a}).ST2.Average;
   S_big.stimvel(:,:,a) = ExcitRR.(ar{a}).S.Average;
  
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
   S_bigend.stimvel(:,:,a) = ExcitRR_end.(ar{a}).S.Average;
  
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