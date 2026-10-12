%% Separate out CSMs based off of timing of the CSM 
thresh=0.1;
dirdir= 1;
% ExcitRR_Ctrl=ExcitRR_AllSuc;

an=fieldnames(ExcitRR_Ctrl);
for a=1:length(an)
    Rreaches=ExcitRR_Ctrl.(an{a}).ReReaches;
    C= AlignThresh(Rreaches,thresh,dirdir,200);
    C=stimvel_accel(C);
    C=stimvel_accel_jerk(C);
    ExcitRR_Ctrl.(an{a}).C=C;
    rrt=ExcitRR_Ctrl.(an{a}).StimRR; rrt=rrt(:,1);

    CN=ExcitRR_Ctrl.(an{a}).CN;
    CN=stimvel_accel(CN);
    CN=stimvel_accel_jerk(CN);
    ExcitRR_Ctrl.(an{a}).CN=CN;

 
    % Timings taken from histo from all animals 
    T0Reaches=[];
    T1Reaches=[];
    T2Reaches=[];
    T3Reaches=[];
   
    
    % % % Bins covering all - from freely behaving
    % T0=find(rrt==0);
    % T1=find(rrt>0 & rrt<=0.0588);
    % T2=find(rrt>0.0588 & rrt<0.1328);
    % T3=find(rrt>0.1328 & rrt<=0.2);
    
    % using jerk as the binning windows 
    % T0=find(rrt>=0 & rrt<=0.02);
    % T1=find(rrt>0.02 & rrt<=0.04);
    % T2=find(rrt>0.04 & rrt<0.06);
    % T3=find(rrt>=0.06);



    % Bins covering all 
    % T0=find(rrt==0);
    % T1=find(rrt>0 & rrt<=0.0672);
    % T2=find(rrt>0.0672 & rrt<=0.108);
    % T3=find(rrt>0.108 & rrt<=0.2);

    % % Bins covering all - for headfixed reaches only 
    % T0=find(rrt==0);
    % T1=find(rrt>0 & rrt<=0.0522);
    % T2=find(rrt>0.0522 & rrt<=0.1218);
    % T3=find(rrt>0.1218 & rrt<=0.2);

     % % Bins covering all - for headfixed reaches only 
    T0=find(rrt==0);
    T1=find(rrt>0 & rrt<=0.06);
    T2=find(rrt>0.06 & rrt<=0.12);
    T3=find(rrt>0.12 & rrt<=0.2);




    rnum=fieldnames(Rreaches);
    
    % for all reaches with a 0 time 
   for s=1:length(T0)
        idx0=T0(s);
        T0Reaches.(rnum{idx0})= ExcitRR_Ctrl.(an{a}).ReReaches.(rnum{idx0});  
    end 
    T0Kin=BasicKinematics(T0Reaches);
    ST0 = AlignThresh(T0Reaches,thresh,dirdir,200);
    ST0=stimvel_accel(ST0);
    ST0=stimvel_accel_jerk(ST0);
    ST0V = AlignPeakV(T0Reaches,7,500);
    ST0V=stimvel_accel(ST0V);
    ST0_CSM=AlignThresh(T0Reaches,thresh,dirdir,500);
    ST0_CSM=stimvel_accel(ST0_CSM);

    for f=1:length(T1)
        idx1=T1(f);
        T1Reaches.(rnum{idx1})= ExcitRR_Ctrl.(an{a}).ReReaches.(rnum{idx1});  
    end 

  
    if isempty(T1Reaches)
        T1Kin=[];
        ST1=[];
    else
       
         T1Kin=BasicKinematics(T1Reaches);
         ST1 = AlignThresh(T1Reaches,thresh,dirdir,200);
         ST1=stimvel_accel(ST1);
         ST1=stimvel_accel_jerk(ST1);
    end 
    ST1V = AlignPeakV(T1Reaches,7,500);
    ST1V=stimvel_accel(ST1V);
    ST1_CSM=AlignThreshCSM_T1HF(T1Reaches,thresh,dirdir,500);
    ST1_CSM=stimvel_accel(ST1_CSM);

    for h=1:length(T2)
        idx2=T2(h);
        T2Reaches.(rnum{idx2})= ExcitRR_Ctrl.(an{a}).ReReaches.(rnum{idx2});  
    end 
    T2Kin=BasicKinematics(T2Reaches);
    ST2 = AlignThresh(T2Reaches,thresh,dirdir,200);
    ST2=stimvel_accel(ST2);
    ST2=stimvel_accel_jerk(ST2);
    ST2V = AlignPeakV(T2Reaches,7,500);
    ST2V=stimvel_accel(ST2V);
    ST2_CSM=AlignThreshCSM_T2HF(T2Reaches,thresh,dirdir,500);
    ST2_CSM=stimvel_accel(ST2_CSM);

    for k=1:length(T3)
        idx3=T3(k);
        T3Reaches.(rnum{idx3})= ExcitRR_Ctrl.(an{a}).ReReaches.(rnum{idx3});  
    end 

    if isempty(T3Reaches)
        T3Kin=[];
        ST3=[];
    else
       
         T3Kin=BasicKinematics(T3Reaches);
         ST3 = AlignThresh(T3Reaches,thresh,dirdir,200);
         ST3=stimvel_accel(ST3);
         ST3=stimvel_accel_jerk(ST3);
         ST3V = AlignPeakV(T3Reaches,7,500);
         ST3V=stimvel_accel(ST3V);
         ST3_CSM=AlignThreshCSM_T3HF(T3Reaches,thresh,dirdir,500);
         ST3_CSM=stimvel_accel(ST3_CSM);
    end 



    % Reaches 
     ExcitRR_Ctrl.(an{a}).RRT0=T0Reaches;
     ExcitRR_Ctrl.(an{a}).RRT1=T1Reaches;
     ExcitRR_Ctrl.(an{a}).RRT2=T2Reaches;
     ExcitRR_Ctrl.(an{a}).RRT3=T3Reaches;
   
     NoCSM=ExcitRR_Ctrl.(an{a}).NoReReach;
    nocsmfields=fieldnames(NoCSM); 
    t0fields=fieldnames(T0Reaches);
    t1fields=fieldnames(T1Reaches);
    t2fields=fieldnames(T2Reaches);

    if isempty(T3Reaches)
        t3fields=[];
        rrnum3=[];
    else
        t3fields=fieldnames(T3Reaches);
        rrnum3=length(t3fields);
    end
    
    nocsmnum=length(nocsmfields);
    rrnum0=length(t0fields);
    rrnum1=length(t1fields);
    rrnum2=length(t2fields);
   

    ExcitRR_Ctrl.(an{a}).nocsmnum=nocsmnum;
    ExcitRR_Ctrl.(an{a}).rrnum0=rrnum0;
    ExcitRR_Ctrl.(an{a}).rrnum1=rrnum1;
    ExcitRR_Ctrl.(an{a}).rrnum2=rrnum2;
    ExcitRR_Ctrl.(an{a}).rrnum3=rrnum3;
    % ExcitRR_Ctrl.(an{a}).rrnum4=rrnum4;
    % ExcitRR_Ctrl.(an{a}).rrnum5=rrnum5;

     % Kinematics
     ExcitRR_Ctrl.(an{a}).T0Kin=T0Kin;
     ExcitRR_Ctrl.(an{a}).T1Kin=T1Kin;
     ExcitRR_Ctrl.(an{a}).T2Kin=T2Kin;
     ExcitRR_Ctrl.(an{a}).T3Kin=T3Kin;
     % ExcitRR_Ctrl.(an{a}).T4Kin=T4Kin;
     % ExcitRR_Ctrl.(an{a}).T5Kin=T5Kin;


     %Stimvel
     ExcitRR_Ctrl.(an{a}).ST0=ST0;
     ExcitRR_Ctrl.(an{a}).ST1=ST1;
     ExcitRR_Ctrl.(an{a}).ST2=ST2;
     ExcitRR_Ctrl.(an{a}).ST3=ST3;
     % ExcitRR_Ctrl.(an{a}).ST4=ST4;
     % ExcitRR_Ctrl.(an{a}).ST5=ST5;

     ExcitRR_Ctrl.(an{a}).ST0V=ST0V;
     ExcitRR_Ctrl.(an{a}).ST1V=ST1V;
     ExcitRR_Ctrl.(an{a}).ST2V=ST2V;
     ExcitRR_Ctrl.(an{a}).ST3V=ST3V;
     % ExcitRR_Ctrl.(an{a}).ST4V=ST4V;
     % ExcitRR_Ctrl.(an{a}).ST5V=ST5V;

     ExcitRR_Ctrl.(an{a}).ST0_CSM=ST0_CSM;
     ExcitRR_Ctrl.(an{a}).ST1_CSM=ST1_CSM;
     ExcitRR_Ctrl.(an{a}).ST2_CSM=ST2_CSM;
     ExcitRR_Ctrl.(an{a}).ST3_CSM=ST3_CSM;

     % Add NRR CSM alignment 
     % NRRreaches=ExcitRR_Ctrl.(an{a}).NoReReach;
     % NRR_CSM=AlignThreshCSM_T4(NRRreaches,thresh,dirdir,500);
     % NRR_CSM=stimvel_accel(NRR_CSM);
     % ExcitRR_Ctrl.(an{a}).NRR_CSM=NRR_CSM;
     % 
    
     % success rate 
     % [S_T0, F_T0] = SuccessFilter(T0Reaches);
     % [S_T1, F_T1] = SuccessFilter(T1Reaches);
     % [S_T2, F_T2] = SuccessFilter(T2Reaches);
     % [S_T3, F_T3] = SuccessFilter(T3Reaches);
     % % [S_T4, F_T4] = SuccessFilter(T4Reaches);
     % % [S_T5, F_T5] = SuccessFilter(T5Reaches);
     % 
     % ExcitRR_Ctrl.(an{a}).successrate.T0 = length(fieldnames(S_T0)) / length(fieldnames(T0Reaches));
     % ExcitRR_Ctrl.(an{a}).successrate.T1 = length(fieldnames(S_T1)) / length(fieldnames(T1Reaches));
     % ExcitRR_Ctrl.(an{a}).successrate.T2 = length(fieldnames(S_T2)) / length(fieldnames(T2Reaches));
     % ExcitRR_Ctrl.(an{a}).successrate.T3 = length(fieldnames(S_T3)) / length(fieldnames(T3Reaches));
     %ExcitRR_Ctrl.(an{a}).successrate.T4 = length(fieldnames(S_T4)) / length(fieldnames(T4Reaches));
     %ExcitRR_Ctrl.(an{a}).successrate.T5 = length(fieldnames(S_T5)) / length(fieldnames(T5Reaches));

end 

% ST1, ST4, ST5 - not all animals have these 


%% Separate out CSMs based off of timing of the CSM 
thresh=0.1;

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
    T1=find(rrt>0 & rrt<=0.06);
    T2=find(rrt>0.06 & rrt<=0.12);
    T3=find(rrt>0.12 & rrt<=0.2);
   

    rnum=fieldnames(Rreaches);
    
    % for all reaches with a 0 time 
   for s=1:length(T0)
        idx0=T0(s);
        T0Reaches.(rnum{idx0})= ExcitRR.(an{a}).ReReaches.(rnum{idx0});  
   end 

   if isempty(T0Reaches)
        T0Kin=[];
        ST0=[];
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

  %% Combine Stimvel and look at plots 

T0_big=struct;
T1_big=struct;
T2_big = struct;
T3_big=struct;
C_big=struct;

NRR_big=struct;


ar=fieldnames(ExcitRR_Ctrl);
for a = 1:length(ar) 
   T0_big.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).ST0.Average;
   T1_big.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).ST1.Average;
   T2_big.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).ST2.Average;
   T3_big.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).ST3.Average;
   C_big.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).C.Average;
   NRR_big.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).CN.Average;
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
% % 

%
NRR_big.Average = nanmean(NRR_big.stimvel,3);
NRR_big.Std = nanstd(NRR_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(NRR_big.Average),1); 
for nnum = 1:length(NRR_big.Average)
    samp = squeeze(NRR_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
NRR_big.N = N;

%% For T0 reahces where 1 animal does not have any T0 reaces 
% Add accel column to SN data - not all animals have NoReReach - this is to add accel and to figure out what animals have none
NoT3Stimvel=[];
indices_to_delete = [];

for a=1:length(an)
    NRR=ExcitRR_Ctrl.(an{a}).ST3;
    if isempty(NRR)
        NoT3Stimvel.(an{a}).ST0=[];
        % Add the index of the empty NRR.stimvel to the list
        indices_to_delete(end+1) = a;
    else

    NoT3Stimvel.(an{a}).ST3=NRR;
    NoT3Stimvel.(an{a}).Cnum=ExcitRR_Ctrl.(an{a}).Cnum;
    NoT3Stimvel.(an{a}).RRnum=ExcitRR_Ctrl.(an{a}).RRnum;
    NoT3Stimvel.(an{a}).RRT3=ExcitRR_Ctrl.(an{a}).RRT3;
    NoT3Stimvel.(an{a}).rrnum3=ExcitRR_Ctrl.(an{a}).rrnum3;
    NoT3Stimvel.(an{a}).T3Kin=ExcitRR_Ctrl.(an{a}).T3Kin;
    NoT3Stimvel.(an{a}).ST3V=ExcitRR_Ctrl.(an{a}).ST3V;
     %NoT3Stimvel.(an{a}).ST3_CSM=ExcitRR_Ctrl.(an{a}).ST3_CSM;

    end

end 

% Delete structures where NRR.stimvel is empty
for i = length(indices_to_delete):-1:1
    NoT3Stimvel = rmfield(NoT3Stimvel, an{indices_to_delete(i)});
end
%

at=fieldnames(NoT3Stimvel);
for a = 1:length(at) 

   T3_big.stimvel(:,:,a) = NoT3Stimvel.(at{a}).ST3.Average;

end


T3_big.Average = nanmean(T3_big.stimvel,3);
T3_big.Std = nanstd(T3_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T3_big.Average),1); 
for nnum = 1:length(T3_big.Average)
    samp = squeeze(T3_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
T3_big.N = N;
%% PEak V
  % Combine Stimvel and look at plots 

T0_bigV=struct;
T1_bigV=struct;
T2_bigV = struct;
T3_bigV=struct;




ar=fieldnames(ExcitRR_Ctrl);
for a = 1:length(ar) 
   T0_bigV.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).ST0V.Average;
   T1_bigV.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).ST1V.Average;
   T2_bigV.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).ST2V.Average;
   T3_bigV.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).ST3V.Average;
    
end

T0_bigV.Average = nanmean(T0_bigV.stimvel,3);
T0_bigV.Std = nanstd(T0_bigV.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T0_bigV.Average),1); 
for nnum = 1:length(T0_bigV.Average)
    samp = squeeze(T0_bigV.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
T0_bigV.N = N;

T1_bigV.Average = nanmean(T1_bigV.stimvel,3);
T1_bigV.Std = nanstd(T1_bigV.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T1_bigV.Average),1); 
for nnum = 1:length(T1_bigV.Average)
    samp = squeeze(T1_bigV.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
T1_bigV.N = N;


T2_bigV.Average = nanmean(T2_bigV.stimvel,3);
T2_bigV.Std = nanstd(T2_bigV.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T2_bigV.Average),1); 
for nnum = 1:length(T2_bigV.Average)
    samp = squeeze(T2_bigV.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
T2_bigV.N = N;

%
T3_bigV.Average = nanmean(T3_bigV.stimvel,3);
T3_bigV.Std = nanstd(T3_bigV.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T3_bigV.Average),1); 
for nnum = 1:length(T3_bigV.Average)
    samp = squeeze(T3_bigV.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
T3_bigV.N = N;

 %% Combine Stimvel and look at plots 

T0_bigCSM=struct;
T1_bigCSM=struct;
T2_bigCSM= struct;
%T3_bigCSM=struct;
%NRR_bigCSM=struct;


ar=fieldnames(ExcitRR_Ctrl);
for a = 1:length(ar) 
   T0_bigCSM.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).ST0_CSM.Average;
   T1_bigCSM.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).ST1_CSM.Average;
   T2_bigCSM.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).ST2_CSM.Average;
   %T3_bigCSM.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).ST3_CSM.Average;
   %NRR_bigCSM.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).NRR_CSM.Average;
    
end

T0_bigCSM.Average = nanmean(T0_bigCSM.stimvel,3);
T0_bigCSM.Std = nanstd(T0_bigCSM.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T0_bigCSM.Average),1); 
for nnum = 1:length(T0_bigCSM.Average)
    samp = squeeze(T0_bigCSM.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
T0_bigCSM.N = N;

T1_bigCSM.Average = nanmean(T1_bigCSM.stimvel,3);
T1_bigCSM.Std = nanstd(T1_bigCSM.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T1_bigCSM.Average),1); 
for nnum = 1:length(T1_bigCSM.Average)
    samp = squeeze(T1_bigCSM.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
T1_bigCSM.N = N;


T2_bigCSM.Average = nanmean(T2_bigCSM.stimvel,3);
T2_bigCSM.Std = nanstd(T2_bigCSM.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T2_bigCSM.Average),1); 
for nnum = 1:length(T2_bigCSM.Average)
    samp = squeeze(T2_bigCSM.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
T2_bigCSM.N = N;

% %
% T3_bigCSM.Average = nanmean(T3_bigCSM.stimvel,3);
% T3_bigCSM.Std = nanstd(T3_bigCSM.stimvel,0,3);
% 
% %loop through to count the number of reaches contributing to each time pt
% N = zeros(length(T3_bigCSM.Average),1); 
% for nnum = 1:length(T3_bigCSM.Average)
%     samp = squeeze(T3_bigCSM.stimvel(nnum,2,:));
%     N(nnum) = length(samp(~isnan(samp)));
% 
% end
% T3_bigCSM.N = N;


%
% NRR_bigCSM.Average = nanmean(NRR_bigCSM.stimvel,3);
% NRR_bigCSM.Std = nanstd(NRR_bigCSM.stimvel,0,3);
% 
% %loop through to count the number of reaches contributing to each time pt
% N = zeros(length(NRR_bigCSM.Average),1); 
% for nnum = 1:length(NRR_bigCSM.Average)
%     samp = squeeze(NRR_bigCSM.stimvel(nnum,2,:));
%     N(nnum) = length(samp(~isnan(samp)));
% 

% end
% NRR_bigCSM.N = N;

%% Plot stimvel
x1=0 ;x2=0.05;
%y_upper= 1.5; y_lower=-0.5;
%Add acceleration column to stimvel
T0color=rgb('Thistle');
T1color=rgb('MediumOrchid');
T2color=rgb('SlateBlue');
T3color=rgb('DarkSlateBlue');
NRRcolor=rgb('OrangeRed');


dim =5;

StimvelPlot(T0_big,dim,T0color,1)
StimvelPlot(T1_big,dim,T1color,0)
StimvelPlot(T2_big,dim,T2color,0)
StimvelPlot(T3_big,dim,T3color,0)
StimvelPlot(C_big,dim,'k',0)
StimvelPlot(NRR_big,dim,NRRcolor,0)

%hold on
%fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
xlim([-0.2 0.2])
%xticks(-0.2:0.01:0.2)
ylim('auto')
%ylim([-0.6 1.4])
%ylim([0 0.25])
%ylim([-20 20])
%pbaspect([1 2 1])
%ylim([0 0.3])
xline(0,'--')
xline(0.05,'--')
yline(0,'--')

% xline(0.23,'m')
xlabel('Time from Pos. Thresh.(s)')
%ylabel('Outward Position (cm)')
ylabel('Outward Velocity (cm/s)')
%ylabel('Outward acceleration')

%% Plot stimvel aligned to peak V
x1=0 ;x2=0.05;
%y_upper= 1.5; y_lower=-0.5;
%Add acceleration column to stimvel
T0color=rgb('Thistle');
T1color=rgb('MediumOrchid');
T2color=rgb('SlateBlue');
T3color=rgb('DarkSlateBlue');


dim =5;

StimvelPlot(T0_bigV,dim,T0color,1)
StimvelPlot(T1_bigV,dim,T1color,0)
StimvelPlot(T2_bigV,dim,T2color,0)
StimvelPlot(T3_bigV,dim,T3color,0)


%hold on
%fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
xlim([-0.2 0.2])
ylim('auto')
%ylim([-0.6 1.4])
%ylim([0 0.25])
%ylim([-20 20])
%pbaspect([1 2 1])
%ylim([0 0.3])
xline(0,'--')
xline(0.05,'--')
yline(0,'--')
xlabel('Time from Peak Velocity(s)')
ylabel('Outward Position (cm)')
%ylabel('Outward Velocity (cm/s)')
%ylabel('Outward acceleration')

%% Plot stimvel for CSM
T0color=rgb('Thistle');
T1color=rgb('MediumOrchid');
T2color=rgb('SlateBlue');
T3color=rgb('DarkSlateBlue');
NRRcolor=rgb('OrangeRed');


dim =5;

StimvelPlot(T0_bigCSM,dim,T0color,1)
StimvelPlot(T1_bigCSM,dim,T1color,0)
StimvelPlot(T2_bigCSM,dim,T2color,0)
% StimvelPlot(T3_bigCSM,dim,T3color,0)
% StimvelPlot(NRR_bigCSM,dim,NRRcolor,0)

%xlim([-0.25 0.25])
xlim([-0.5 0.5])
%xticks(-0.2:0.01:0.2)
ylim('auto')
%ylim([-0.6 1.4])
%ylim([0 0.25])
%ylim([-20 20])
%pbaspect([1 2 1])
%ylim([0 0.3])
xline(0,'--')
yline(0,'--')

xlabel('Time from CSM(s)')
ylabel('Outward Position (cm)')
%ylabel('Outward Velocity (cm/s)')
%ylabel('Outward acceleration')
%% Plot position of max velocity for each time point 
 pcolor=rgb('Plum');
figure; hold on
for a=1:length(an)
    vmaxp=ExcitRR_Ctrl.(an{a}).RRKin.vmax_p(:,2);
    rrtime=ExcitRR_Ctrl.(an{a}).StimRR(:,1);
    sz=30;

    scatter(rrtime,vmaxp,sz,pcolor,'filled')
    xlim([0.05 0.2])
    ylim([0 0.25])

end
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
%% Get Endpoint 
T0Endpt=[];
T1Endpt=[];
T2Endpt=[];
T3Endpt=[];
T0Endpt_std=[];
T1Endpt_std=[];
T2Endpt_std=[];
T3Endpt_std=[];

for a=1:length(an)
    T0kin=ExcitRR_Ctrl.(an{a}).T0Kin.endpt(:,1);
    T0endpt=mean(T0kin);
    T0endpt_std=std(T0kin);

    T1kin=ExcitRR_Ctrl.(an{a}).T1Kin.endpt(:,1);
    T1endpt=mean(T1kin);
    T1endpt_std=std(T1kin);

    T2kin=ExcitRR_Ctrl.(an{a}).T2Kin.endpt(:,1);
    T2endpt=mean(T2kin);
    T2endpt_std=std(T2kin);

    T3kin=ExcitRR_Ctrl.(an{a}).T3Kin.endpt(:,1);
    T3endpt=mean(T3kin);
    T3endpt_std=std(T3kin);
    
    T0Endpt=[T0Endpt;T0endpt];
    T0Endpt_std=[T0Endpt_std;T0endpt_std];

    T1Endpt=[T1Endpt;T1endpt];
    T1Endpt_std=[T1Endpt_std;T1endpt_std];

    T2Endpt=[T2Endpt;T2endpt];
    T2Endpt_std=[T2Endpt_std;T2endpt_std];

    T3Endpt=[T3Endpt;T3endpt];
    T3Endpt_std=[T3Endpt_std;T3endpt_std];


end 
