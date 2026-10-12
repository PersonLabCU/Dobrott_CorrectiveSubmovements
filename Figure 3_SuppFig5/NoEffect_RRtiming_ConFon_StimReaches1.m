%% Separate out CSMs based off of timing of the CSM 
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
C_big=struct;
T0_big=struct;
T1_big=struct;
T2_big = struct;
T3_big=struct;
S_big=struct;


ar=fieldnames(ExcitRR);
for a = 1:length(ar) 
  %C_big.stimvel(:,:,a) = ExcitRR.(ar{a}).C.Average;

   %T0_big.stimvel(:,:,a) = ExcitRR.(ar{a}).ST0.Average;
   T1_big.stimvel(:,:,a) = ExcitRR.(ar{a}).ST1.Average;
   T2_big.stimvel(:,:,a) = ExcitRR.(ar{a}).ST2.Average;
   T3_big.stimvel(:,:,a) = ExcitRR.(ar{a}).ST3.Average;
   S_big.stimvel(:,:,a) = ExcitRR.(ar{a}).S.Average;
   %NRR_big.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).CN.Average;
end

% C_big.Average = nanmean(C_big.stimvel,3);
% C_big.Std = nanstd(C_big.stimvel,0,3);
% 
% %loop through to count the number of reaches contributing to each time pt
% N = zeros(length(C_big.Average),1); 
% for nnum = 1:length(C_big.Average)
%     samp = squeeze(C_big.stimvel(nnum,2,:));
%     N(nnum) = length(samp(~isnan(samp)));
% 
% end
% C_big.N = N;

S_big.Average = nanmean(S_big.stimvel,3);
S_big.Std = nanstd(S_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(S_big.Average),1); 
for nnum = 1:length(S_big.Average)
    samp = squeeze(S_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
S_big.N = N;



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
%%
ac=fieldnames(ExcitRR_Ctrl);
for a = 1:length(ac) 
 
   C_big.stimvel(:,:,a) = ExcitRR_Ctrl.(ac{a}).C.Average;
   
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
%% For T0 reahces where 1 animal does not have any T0 reaces 
% Add accel column to SN data - not all animals have NoReReach - this is to add accel and to figure out what animals have none
NoT0Stimvel=[];
indices_to_delete = [];

for a=1:length(an)
    NRR=ExcitRR.(an{a}).ST0;
    if isempty(NRR)
        NoT0Stimvel.(an{a}).ST0=[];
        % Add the index of the empty NRR.stimvel to the list
        indices_to_delete(end+1) = a;
    else
   
    NoT0Stimvel.(an{a}).ST0=NRR;
    NoT0Stimvel.(an{a}).Snum=ExcitRR.(an{a}).Snum;
    NoT0Stimvel.(an{a}).RRnum=ExcitRR.(an{a}).RRnum;
    NoT0Stimvel.(an{a}).RRT0=ExcitRR.(an{a}).RRT0;
    NoT0Stimvel.(an{a}).rrnum0=ExcitRR.(an{a}).rrnum0;
    NoT0Stimvel.(an{a}).T0Kin=ExcitRR.(an{a}).T0Kin;
    NoT0Stimvel.(an{a}).ST0V=ExcitRR.(an{a}).ST0V;
     NoT0Stimvel.(an{a}).ST0_CSM=ExcitRR.(an{a}).ST0_CSM;
   
    end

end 

% Delete structures where NRR.stimvel is empty
for i = length(indices_to_delete):-1:1
    NoT0Stimvel = rmfield(NoT0Stimvel, an{indices_to_delete(i)});
end
%%

at=fieldnames(NoT0Stimvel);
for a = 1:length(at) 
 
   T0_big.stimvel(:,:,a) = NoT0Stimvel.(at{a}).ST0.Average;
  
end


T0_big.Average = nanmean(T0_big.stimvel,3);
T0_big.Std = nanstd(T0_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(T0_big.Average),1); 
for nnum = 1:length(T0_big.Average)
    samp = squeeze(T0_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
T0_big.N = N;
%% PEak V
  % Combine Stimvel and look at plots 

% T0_bigV=struct;
% T1_bigV=struct;
% T2_bigV = struct;
% T3_bigV=struct;
% 
% 
% 
% ar=fieldnames(ExcitRR_Ctrl);
% for a = 1:length(ar) 
%    T0_bigV.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).ST0V.Average;
% 
%    T1_bigV.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).ST1V.Average;
%    T2_bigV.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).ST2V.Average;
%    T3_bigV.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).ST3V.Average;
% 
% end
% 
% T0_bigV.Average = nanmean(T0_bigV.stimvel,3);
% T0_bigV.Std = nanstd(T0_bigV.stimvel,0,3);
% 
% %loop through to count the number of reaches contributing to each time pt
% N = zeros(length(T0_bigV.Average),1); 
% for nnum = 1:length(T0_bigV.Average)
%     samp = squeeze(T0_bigV.stimvel(nnum,2,:));
%     N(nnum) = length(samp(~isnan(samp)));
% 
% end
% T0_bigV.N = N;
% 
% T1_bigV.Average = nanmean(T1_bigV.stimvel,3);
% T1_bigV.Std = nanstd(T1_bigV.stimvel,0,3);
% 
% %loop through to count the number of reaches contributing to each time pt
% N = zeros(length(T1_bigV.Average),1); 
% for nnum = 1:length(T1_bigV.Average)
%     samp = squeeze(T1_bigV.stimvel(nnum,2,:));
%     N(nnum) = length(samp(~isnan(samp)));
% 
% end
% T1_bigV.N = N;
% 
% 
% T2_bigV.Average = nanmean(T2_bigV.stimvel,3);
% T2_bigV.Std = nanstd(T2_bigV.stimvel,0,3);
% 
% %loop through to count the number of reaches contributing to each time pt
% N = zeros(length(T2_bigV.Average),1); 
% for nnum = 1:length(T2_bigV.Average)
%     samp = squeeze(T2_bigV.stimvel(nnum,2,:));
%     N(nnum) = length(samp(~isnan(samp)));
% 
% end
% T2_bigV.N = N;
% 
% %
% T3_bigV.Average = nanmean(T3_bigV.stimvel,3);
% T3_bigV.Std = nanstd(T3_bigV.stimvel,0,3);
% 
% %loop through to count the number of reaches contributing to each time pt
% N = zeros(length(T3_bigV.Average),1); 
% for nnum = 1:length(T3_bigV.Average)
%     samp = squeeze(T3_bigV.stimvel(nnum,2,:));
%     N(nnum) = length(samp(~isnan(samp)));
% 
% end
% T3_bigV.N = N;

 %% Combine Stimvel and look at plots 

% T0_bigCSM=struct;
% T1_bigCSM=struct;
% T2_bigCSM= struct;
% T3_bigCSM=struct;
% %NRR_bigCSM=struct;
% 
% 
% ar=fieldnames(ExcitRR_Ctrl);
% for a = 1:length(ar) 
%    T0_bigCSM.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).ST0_CSM.Average;
%    T1_bigCSM.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).ST1_CSM.Average;
%    T2_bigCSM.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).ST2_CSM.Average;
%    T3_bigCSM.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).ST3_CSM.Average;
%    %NRR_bigCSM.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).NRR_CSM.Average;
% 
% end
% 
% T0_bigCSM.Average = nanmean(T0_bigCSM.stimvel,3);
% T0_bigCSM.Std = nanstd(T0_bigCSM.stimvel,0,3);
% 
% %loop through to count the number of reaches contributing to each time pt
% N = zeros(length(T0_bigCSM.Average),1); 
% for nnum = 1:length(T0_bigCSM.Average)
%     samp = squeeze(T0_bigCSM.stimvel(nnum,2,:));
%     N(nnum) = length(samp(~isnan(samp)));
% 
% end
% T0_bigCSM.N = N;
% 
% T1_bigCSM.Average = nanmean(T1_bigCSM.stimvel,3);
% T1_bigCSM.Std = nanstd(T1_bigCSM.stimvel,0,3);
% 
% %loop through to count the number of reaches contributing to each time pt
% N = zeros(length(T1_bigCSM.Average),1); 
% for nnum = 1:length(T1_bigCSM.Average)
%     samp = squeeze(T1_bigCSM.stimvel(nnum,2,:));
%     N(nnum) = length(samp(~isnan(samp)));
% 
% end
% T1_bigCSM.N = N;
% 
% 
% T2_bigCSM.Average = nanmean(T2_bigCSM.stimvel,3);
% T2_bigCSM.Std = nanstd(T2_bigCSM.stimvel,0,3);
% 
% %loop through to count the number of reaches contributing to each time pt
% N = zeros(length(T2_bigCSM.Average),1); 
% for nnum = 1:length(T2_bigCSM.Average)
%     samp = squeeze(T2_bigCSM.stimvel(nnum,2,:));
%     N(nnum) = length(samp(~isnan(samp)));
% 
% end
% T2_bigCSM.N = N;
% 
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
% 
% % %
% % NRR_bigCSM.Average = nanmean(NRR_bigCSM.stimvel,3);
% % NRR_bigCSM.Std = nanstd(NRR_bigCSM.stimvel,0,3);
% % 
% % %loop through to count the number of reaches contributing to each time pt
% % N = zeros(length(NRR_bigCSM.Average),1); 
% % for nnum = 1:length(NRR_bigCSM.Average)
% %     samp = squeeze(NRR_bigCSM.stimvel(nnum,2,:));
% %     N(nnum) = length(samp(~isnan(samp)));
% % 
% % end
% % NRR_bigCSM.N = N;


%% Combine Stimvel and look at plots 

% CT0_bigCSM=struct;
% CT1_bigCSM=struct;
% CT2_bigCSM= struct;
% CT3_bigCSM=struct;
% %NRR_bigCSM=struct;
% 
% 
% ar=fieldnames(ExcitRR_Ctrl);
% for a = 1:length(ar) 
%    CT0_bigCSM.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).CST0_CSM.Average;
%    CT1_bigCSM.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).CST1_CSM.Average;
%    CT2_bigCSM.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).CST2_CSM.Average;
%    CT3_bigCSM.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).CST3_CSM.Average;
%    %NRR_bigCSM.stimvel(:,:,a) = ExcitRR_Ctrl.(ar{a}).NRR_CSM.Average;
% 
% end
% 
% CT0_bigCSM.Average = nanmean(CT0_bigCSM.stimvel,3);
% CT0_bigCSM.Std = nanstd(CT0_bigCSM.stimvel,0,3);
% 
% %loop through to count the number of reaches contributing to each time pt
% N = zeros(length(CT0_bigCSM.Average),1); 
% for nnum = 1:length(CT0_bigCSM.Average)
%     samp = squeeze(CT0_bigCSM.stimvel(nnum,2,:));
%     N(nnum) = length(samp(~isnan(samp)));
% 
% end
% CT0_bigCSM.N = N;
% 
% CT1_bigCSM.Average = nanmean(CT1_bigCSM.stimvel,3);
% CT1_bigCSM.Std = nanstd(CT1_bigCSM.stimvel,0,3);
% 
% %loop through to count the number of reaches contributing to each time pt
% N = zeros(length(CT1_bigCSM.Average),1); 
% for nnum = 1:length(CT1_bigCSM.Average)
%     samp = squeeze(CT1_bigCSM.stimvel(nnum,2,:));
%     N(nnum) = length(samp(~isnan(samp)));
% 
% end
% CT1_bigCSM.N = N;
% 
% 
% CT2_bigCSM.Average = nanmean(CT2_bigCSM.stimvel,3);
% CT2_bigCSM.Std = nanstd(CT2_bigCSM.stimvel,0,3);
% 
% %loop through to count the number of reaches contributing to each time pt
% N = zeros(length(CT2_bigCSM.Average),1); 
% for nnum = 1:length(CT2_bigCSM.Average)
%     samp = squeeze(CT2_bigCSM.stimvel(nnum,2,:));
%     N(nnum) = length(samp(~isnan(samp)));
% 
% end
% CT2_bigCSM.N = N;
% 
% %
% CT3_bigCSM.Average = nanmean(CT3_bigCSM.stimvel,3);
% CT3_bigCSM.Std = nanstd(CT3_bigCSM.stimvel,0,3);
% 
% %loop through to count the number of reaches contributing to each time pt
% N = zeros(length(CT3_bigCSM.Average),1); 
% for nnum = 1:length(CT3_bigCSM.Average)
%     samp = squeeze(CT3_bigCSM.stimvel(nnum,2,:));
%     N(nnum) = length(samp(~isnan(samp)));
% 
% end
% CT3_bigCSM.N = N;
%% Plot stimvel
x1=0 ;x2=0.05;
%y_upper= 1.5; y_lower=-0.5;
%Add acceleration column to stimvel
% T0color=rgb('Thistle');
% T1color=rgb('MediumOrchid');
% T2color=rgb('SlateBlue');
% T3color=rgb('DarkSlateBlue');

T0color=rgb('LightBlue');
T1color=rgb('SteelBlue');
T2color=rgb('MediumBlue');
T3color=rgb('MidnightBlue');

NRRcolor=rgb('OrangeRed');
stimcolor=rgb('DodgerBlue');
ctrlcolor=rgb('Black');

dim =2;

StimvelPlot(T0_big,dim,T0color,1)
StimvelPlot(T1_big,dim,T1color,0)
StimvelPlot(T2_big,dim,T2color,0)
StimvelPlot(T3_big,dim,T3color,0)
StimvelPlot(S_big,dim,stimcolor,0)
StimvelPlot(C_big,dim,'k',0)

%StimvelPlot(NRR_big,dim,NRRcolor,0)

%hold on
%fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
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
ylabel('Outward Position (cm)')
%ylabel('Outward Velocity (cm/s)')
%ylabel('Outward acceleration')

%% Plot stimvel aligned to peak V
% x1=0 ;x2=0.05;
% %y_upper= 1.5; y_lower=-0.5;
% %Add acceleration column to stimvel
% 
% 
% dim =2;
% 
% StimvelPlot(T0_bigV,dim,T0color,1)
% StimvelPlot(T1_bigV,dim,T1color,0)
% StimvelPlot(T2_bigV,dim,T2color,0)
% StimvelPlot(T3_bigV,dim,T3color,0)
% 
% 
% %hold on
% %fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
% xlim([-0.2 0.2])
% ylim('auto')
% %ylim([-0.6 1.4])
% %ylim([0 0.25])
% %ylim([-20 20])
% pbaspect([1 2 1])
% %ylim([0 0.3])
% xline(0,'--')
% xline(0.05,'--')
% yline(0,'--')
% xlabel('Time from Peak Velocity(s)')
% ylabel('Outward Position (cm)')
% %ylabel('Outward Velocity (cm/s)')
% %ylabel('Outward acceleration')

%% Plot stimvel for CSM


% dim =2;
% 
% StimvelPlot(T0_bigCSM,dim,T0color,1)
% StimvelPlot(T1_bigCSM,dim,T1color,0)
% StimvelPlot(T2_bigCSM,dim,T2color,0)
% StimvelPlot(T3_bigCSM,dim,T3color,0)
% StimvelPlot(NRR_bigCSM,dim,NRRcolor,0)
% 
% %xlim([-0.25 0.25])
% xlim([-0.5 0.5])
% %xticks(-0.2:0.01:0.2)
% ylim('auto')
% %ylim([-0.6 1.4])
% %ylim([0 0.25])
% %ylim([-20 20])
% %pbaspect([1 2 1])
% %ylim([0 0.3])
% xline(0,'--')
% yline(0,'--')
% 
% xlabel('Time from CSM(s)')
% ylabel('Outward Position (cm)')
%ylabel('Outward Velocity (cm/s)')
%ylabel('Outward acceleration')

%% Get Endpoint 
T0Endpt=[]; T0Endpt_std=[];
for a=1:length(at)
    T0kin=NoT0Stimvel.(at{a}).T0Kin.endpt(:,1);
    T0endpt=mean(T0kin);
    T0endpt_std=std(T0kin);

    T0Endpt=[T0Endpt;T0endpt];
    T0Endpt_std=[T0Endpt_std;T0endpt_std];
end
%%
an=fieldnames(ExcitRR);

T1Endpt=[];
T2Endpt=[];
T3Endpt=[];
T1Endpt_std=[];
T2Endpt_std=[];
T3Endpt_std=[];

for a=1:length(an)

    T1kin=ExcitRR.(an{a}).T1Kin.endpt(:,1);
    T1endpt=mean(T1kin);
    T1endpt_std=std(T1kin);

    T2kin=ExcitRR.(an{a}).T2Kin.endpt(:,1);
    T2endpt=mean(T2kin);
    T2endpt_std=std(T2kin);

    T3kin=ExcitRR.(an{a}).T3Kin.endpt(:,1);
    T3endpt=mean(T3kin);
    T3endpt_std=std(T3kin);

    T1Endpt=[T1Endpt;T1endpt];
    T1Endpt_std=[T1Endpt_std;T1endpt_std];

    T2Endpt=[T2Endpt;T2endpt];
    T2Endpt_std=[T2Endpt_std;T2endpt_std];

    T3Endpt=[T3Endpt;T3endpt];
    T3Endpt_std=[T3Endpt_std;T3endpt_std];


end 


%% Find min and max of accel during CSM time window ( 100ms after stim ) 
%T0- not all animals have a T0
at=fieldnames(NoT0Stimvel);

for a=1:length(at)
    stimvel=NoT0Stimvel.(at{a}).ST0_CSM.stimvel;
    Accelmin=[];
    Accelmax=[];
   stimvelsize=size(stimvel);
   num=(stimvelsize(3));

    for b=1:num
        reach=stimvel(:,:,b);
        decel_window=reach(1:21,:);
        accel_window=reach(21:41,:);
       
        accelmin=min(decel_window(:,8));
        accelmax=max(accel_window(:,8));
        Accelmin=[Accelmin;accelmin];
        Accelmax=[Accelmax;accelmax];
    end

    NoT0Stimvel.(at{a}).AccelMinT0=Accelmin;
    NoT0Stimvel.(at{a}).AccelMaxT0=Accelmax;
    clear Accelmin Accelmax 
end

%% Peak Velocity for nitial movement and CSM 
%T0
for a=1:length(at)
    stimvel=NoT0Stimvel.(at{a}).ST0_CSM.stimvel;
    Accelmin=[];
    Accelmax=[];
    Time1=[];
    Time2=[];
   stimvelsize=size(stimvel);
   num=(stimvelsize(3));

    for b=1:num
        reach=stimvel(:,:,b);
        decel_window=reach(1:21,:);
        accel_window=reach(21:41,:);
       
        % accelmin=max(decel_window(:,5));
        % accelmax=max(accel_window(:,5));

        [accelmin,idx_r]=max(decel_window(:,5));
        [accelmax,idx_csm]=max(accel_window(:,5));

        Accelmin=[Accelmin;accelmin];
        Accelmax=[Accelmax;accelmax];

        time1=decel_window(idx_r,1);
        time2=accel_window(idx_csm,1);
        Time1=[Time1;time1];
        Time2=[Time2;time2];

    end

    NoT0Stimvel.(at{a}).Vmax1T0=Accelmin;
    NoT0Stimvel.(at{a}).Vmax2T0=Accelmax;

    NoT0Stimvel.(at{a}).Time1T0=Time1;
    NoT0Stimvel.(at{a}).Time2T0=Time2;
    clear Accelmin Accelmax Time1 Time2
end

%%
for a = 1:length(an)
    stimvel = ExcitRR.(an{a}).ST1_CSM.stimvel;
    Accelmin = [];
    Accelmax = [];

    % Check the number of dimensions in stimvel
    if ndims(stimvel) < 3
        % Single instance case (2D array)
        num = 1;
    else
        % Multiple instances case (3D array)
        num = size(stimvel, 3);
    end

    % Process the data based on the value of num
    if num == 1
        % Single instance
        reach = stimvel;  % No need for indexing along the third dimension
        decel_window = reach(1:21, :);
        accel_window = reach(21:41, :);

        accelmin = min(decel_window(:, 8));
        accelmax = max(accel_window(:, 8));

        Accelmin = accelmin;
        Accelmax = accelmax;
    else
        % Multiple instances
        for b = 1:num
            reach = stimvel(:, :, b);
            decel_window = reach(1:21, :);
            accel_window = reach(21:41, :);

            accelmin = min(decel_window(:, 8));
            accelmax = max(accel_window(:, 8));
            Accelmin = [Accelmin; accelmin];
            Accelmax = [Accelmax; accelmax];
        end
    end

    % Save results back to the structure
    ExcitRR.(an{a}).AccelMinT1 = Accelmin;
    ExcitRR.(an{a}).AccelMaxT1 = Accelmax;

    % Clear temporary variables
    clear Accelmin Accelmax
end
%% PEak Velocity 
for a = 1:length(an)
    stimvel = ExcitRR.(an{a}).ST1_CSM.stimvel;
    Accelmin = [];
    Accelmax = [];
     Time1=[];
    Time2=[];

    % Check the number of dimensions in stimvel
    if ndims(stimvel) < 3
        % Single instance case (2D array)
        num = 1;
    else
        % Multiple instances case (3D array)
        num = size(stimvel, 3);
    end

    % Process the data based on the value of num
    if num == 1
        % Single instance
        reach = stimvel;  % No need for indexing along the third dimension
        decel_window = reach(1:21, :);
        accel_window = reach(21:41, :);

        [accelmin,idx_r]=max(decel_window(:,5));
        [accelmax,idx_csm]=max(accel_window(:,5));
        % get time of max 
        time1=decel_window(idx_r,1);
        time2=accel_window(idx_csm,1);

        Accelmin = accelmin;
        Accelmax = accelmax;

        Time1=[Time1;time1];
        Time2=[Time2;time2];
    else
        % Multiple instances
        for b = 1:num
            reach = stimvel(:, :, b);
            decel_window = reach(1:21, :);
            accel_window = reach(21:41, :);

         
             [accelmin,idx_r]=max(decel_window(:,5));
            [accelmax,idx_csm]=max(accel_window(:,5));

            % get time of max 
            time1=decel_window(idx_r,1);
            time2=accel_window(idx_csm,1);

            Accelmin = [Accelmin; accelmin];
            Accelmax = [Accelmax; accelmax];

            Time1=[Time1;time1];
            Time2=[Time2;time2];
        end
    end

    % Save results back to the structure
    ExcitRR.(an{a}).Vmax1T1 = Accelmin;
    ExcitRR.(an{a}).Vmax2T1 = Accelmax;

    ExcitRR.(an{a}).Time1T1=Time1;
    ExcitRR.(an{a}).Time2T1=Time2;

    % Clear temporary variables
    clear Accelmin Accelmax Time1 Time2
end


%% 
%T2
for a=1:length(an)
    stimvel=ExcitRR.(an{a}).ST2_CSM.stimvel;
    Accelmin=[];
    Accelmax=[];
   stimvelsize=size(stimvel);
   num=(stimvelsize(3));

    for b=1:num
        reach=stimvel(:,:,b);
        decel_window=reach(1:21,:);
        accel_window=reach(21:41,:);
       
        accelmin=min(decel_window(:,8));
        accelmax=max(accel_window(:,8));
        Accelmin=[Accelmin;accelmin];
        Accelmax=[Accelmax;accelmax];
    end

    ExcitRR.(an{a}).AccelMinT2=Accelmin;
    ExcitRR.(an{a}).AccelMaxT2=Accelmax;
    clear Accelmin Accelmax 
end
%% Peak Velocity for nitial movement and CSM 
%T2
for a=1:length(an)
    stimvel=ExcitRR.(an{a}).ST2_CSM.stimvel;
    Accelmin=[];
    Accelmax=[];
    Time1=[];
    Time2=[];
   stimvelsize=size(stimvel);
   num=(stimvelsize(3));

    for b=1:num
        reach=stimvel(:,:,b);
        decel_window=reach(1:21,:);
        accel_window=reach(21:41,:);
        % 
        % accelmin=max(decel_window(:,5));
        % accelmax=max(accel_window(:,5));
        [accelmin,idx_r]=max(decel_window(:,5));
        [accelmax,idx_csm]=max(accel_window(:,5));
        Accelmin=[Accelmin;accelmin];
        Accelmax=[Accelmax;accelmax];

        time1=decel_window(idx_r,1);
        time2=accel_window(idx_csm,1);
        Time1=[Time1;time1];
        Time2=[Time2;time2];
    end

    ExcitRR.(an{a}).Vmax1T2=Accelmin;
    ExcitRR.(an{a}).Vmax2T2=Accelmax;

    ExcitRR.(an{a}).Time1T2=Time1;
    ExcitRR.(an{a}).Time2T2=Time2;

    clear Accelmin Accelmax Time1 Time2
end
%% 
%T3
for a=1:length(an)
    stimvel=ExcitRR.(an{a}).ST3_CSM.stimvel;
    Accelmin=[];
    Accelmax=[];
   stimvelsize=size(stimvel);
   num=(stimvelsize(3));

    for b=1:num
        reach=stimvel(:,:,b);
        decel_window=reach(1:21,:);
        accel_window=reach(21:41,:);
       
        accelmin=min(decel_window(:,8));
        accelmax=max(accel_window(:,8));
        Accelmin=[Accelmin;accelmin];
        Accelmax=[Accelmax;accelmax];
    end

    ExcitRR.(an{a}).AccelMinT3=Accelmin;
    ExcitRR.(an{a}).AccelMaxT3=Accelmax;
    clear Accelmin Accelmax 
end

%% Peak Velocity for nitial movement and CSM 
%T3
for a=1:length(an)
    stimvel=ExcitRR.(an{a}).ST3_CSM.stimvel;
    Accelmin=[];
    Accelmax=[];
    Time1=[];
    Time2=[];
   stimvelsize=size(stimvel);
   num=(stimvelsize(3));

    for b=1:num
        reach=stimvel(:,:,b);
        decel_window=reach(1:21,:);
        accel_window=reach(21:41,:);
        % 
        % accelmin=max(decel_window(:,5));
        % accelmax=max(accel_window(:,5));
        [accelmin,idx_r]=max(decel_window(:,5));
        [accelmax,idx_csm]=max(accel_window(:,5));
        Accelmin=[Accelmin;accelmin];
        Accelmax=[Accelmax;accelmax];
        

        time1=decel_window(idx_r,1);
        time2=accel_window(idx_csm,1);
        Time1=[Time1;time1];
        Time2=[Time2;time2];

    end

    ExcitRR.(an{a}).Vmax1T3=Accelmin;
    ExcitRR.(an{a}).Vmax2T3=Accelmax;

    ExcitRR.(an{a}).Time1T3=Time1;
    ExcitRR.(an{a}).Time2T3=Time2;
    clear Accelmin Accelmax Time1 Time2
end
%% Combine time bins across animals 
AccelMaxT0=[];AccelMaxT1=[]; AccelMaxT2=[]; AccelMaxT3=[]; 
AccelMinT0=[];AccelMinT1=[];AccelMinT2=[];AccelMinT3=[];
%% T0 time bin
for a=1:length(at)
    AccelminT0=NoT0Stimvel.(at{a}).AccelMinT0;
    AccelmaxT0=NoT0Stimvel.(at{a}).AccelMaxT0; 
    AccelMaxT0=[AccelMaxT0;AccelmaxT0];
    AccelMinT0=[AccelMinT0;AccelminT0];
end
%%
for a=1:length(an)

    AccelminT1=ExcitRR.(an{a}).AccelMinT1;
    AccelmaxT1=ExcitRR.(an{a}).AccelMaxT1; 
    AccelMaxT1=[AccelMaxT1;AccelmaxT1];
    AccelMinT1=[AccelMinT1;AccelminT1];

    AccelminT2=ExcitRR.(an{a}).AccelMinT2;
    AccelmaxT2=ExcitRR.(an{a}).AccelMaxT2; 
    AccelMaxT2=[AccelMaxT2;AccelmaxT2];
    AccelMinT2=[AccelMinT2;AccelminT2];

    AccelminT3=ExcitRR.(an{a}).AccelMinT3;
    AccelmaxT3=ExcitRR.(an{a}).AccelMaxT3; 
    AccelMaxT3=[AccelMaxT3;AccelmaxT3];
    AccelMinT3=[AccelMinT3;AccelminT3];

end 
%% Combine time bins across animals 
Vmax1T0=[];Vmax1T1=[];Vmax1T2=[];Vmax1T3=[];
Vmax2T0=[];Vmax2T1=[];Vmax2T2=[];Vmax2T3=[];
%%
for a=1:length(at)
    vmax1T0=NoT0Stimvel.(at{a}).Vmax1T0;
    vmax2T0=NoT0Stimvel.(at{a}).Vmax2T0; 
    Vmax1T0=[Vmax1T0;vmax1T0];
    Vmax2T0=[Vmax2T0;vmax2T0];
end

%%
for a=1:length(an)

    vmax1T1=ExcitRR.(an{a}).Vmax1T1;
    vmax2T1=ExcitRR.(an{a}).Vmax2T1; 
    Vmax1T1=[Vmax1T1;vmax1T1];
    Vmax2T1=[Vmax2T1;vmax2T1];

    vmax1T2=ExcitRR.(an{a}).Vmax1T2;
    vmax2T2=ExcitRR.(an{a}).Vmax2T2; 
    Vmax1T2=[Vmax1T2;vmax1T2];
    Vmax2T2=[Vmax2T2;vmax2T2];

    vmax1T3=ExcitRR.(an{a}).Vmax1T3;
    vmax2T3=ExcitRR.(an{a}).Vmax2T3; 
    Vmax1T3=[Vmax1T3;vmax1T3];
    Vmax2T3=[Vmax2T3;vmax2T3];

end 
%% Combine time bins across animals 
Time1T0=[];Time1T1=[];Time1T2=[];Time1T3=[];
Time2T0=[];Time2T1=[];Time2T2=[];Time2T3=[];
%%
for a=1:length(at)
    t1T0=NoT0Stimvel.(at{a}).Time1T0;
    t2T0=NoT0Stimvel.(at{a}).Time2T0; 
    Time1T0=[Time1T0;t1T0];
    Time2T0=[Time2T0;t2T0];
end
%%
for a=1:length(an)

    t1T1=ExcitRR.(an{a}).Time1T1;
    t2T1=ExcitRR.(an{a}).Time2T1; 
    Time1T1=[Time1T1;t1T1];
    Time2T1=[Time2T1;t2T1];

    t1T2=ExcitRR.(an{a}).Time1T2;
    t2T2=ExcitRR.(an{a}).Time2T2; 
    Time1T2=[Time1T2;t1T2];
    Time2T2=[Time2T2;t2T2];

    t1T3=ExcitRR.(an{a}).Time1T3;
    t2T3=ExcitRR.(an{a}).Time2T3; 
    Time1T3=[Time1T3;t1T3];
    Time2T3=[Time2T3;t2T3];

end 
%% TIme 
TimeT0=[]; TimeT1=[];TimeT2=[]; TimeT3=[];
TimeT0=[Time1T0,Time2T0];
TimeT1=[Time1T1,Time2T1];
TimeT2=[Time1T2,Time2T2];
TimeT3=[Time1T3,Time2T3];
%% Combine into the same structure 
KinVarsT0=[];KinVarsT1=[]; KinVarsT2=[];KinVarsT3=[];

KinVarsT0=[AccelMinT0,AccelMaxT0];
KinVarsT1=[AccelMinT1,AccelMaxT1];
KinVarsT2=[AccelMinT2,AccelMaxT2];
KinVarsT3=[AccelMinT3,AccelMaxT3];

VMaxT0=[]; VMaxT1=[]; VMaxT2=[]; VMaxT3=[];
% Combine velocities - 2 peaks 
VMaxT0=[Vmax1T0,Vmax2T0];
VMaxT1=[Vmax1T1,Vmax2T1];
VMaxT2=[Vmax1T2,Vmax2T2];
VMaxT3=[Vmax1T3,Vmax2T3];

%% Remove outlier


row_vmax1=any(Vmax1T1>3,2);
Vmax1T1_edit=Vmax1T1;
Vmax1T1_edit(row_vmax1,:)=[];

row_vmax2=any(Vmax1T2>3,2);
Vmax1T2_edit=Vmax1T2;
Vmax1T2_edit(row_vmax2,:)=[];


row_vmax3=any(Vmax1T3>3,2);
Vmax1T3_edit=Vmax1T3;
Vmax1T3_edit(row_vmax3,:)=[];
%% 
row_vmax2_neg=any(Vmax2T2<-1,2);
Vmax2T2_edit=Vmax2T2;
Vmax2T2_edit(row_vmax2_neg,:)=[];


row_vmax3_neg=any(Vmax2T3<-1,2);
Vmax2T3_edit=Vmax2T3;
Vmax2T3_edit(row_vmax3_neg,:)=[];
%% Plot vmax
%initialc=rgb('Indigo');
%csmc=rgb('DarkCyan');
initialc=rgb('DarkMagenta');
csmc=rgb('MediumBlue');
numbin0=ceil(sqrt(length(Vmax1T0)));
numbin1=ceil(sqrt(length(Vmax1T1_edit)));
numbin2=ceil(sqrt(length(Vmax1T2_edit)));
numbin3=ceil(sqrt(length(Vmax1T3_edit)));

numbinsall=(numbin0+numbin1+numbin2+numbin3)/4

%% Calculate mean for each 
MeanT0=mean(VMaxT0);
MeanT1=mean(VMaxT1);
MeanT2=mean(VMaxT2);
MeanT3=mean(VMaxT3);

MeanT0t=mean(TimeT0);
MeanT1t=mean(TimeT1);
MeanT2t=mean(TimeT2);
MeanT3t=mean(TimeT3);
%% Time 
% Combine data to calculate bin edges based on the entire range
all_datat = [Time1T0; Time2T0];  % Use relevant datasets for each subplot
all_data1t = [Time1T1; Time2T1]; 
all_data2t = [Time1T2; Time2T2]; 
all_data3t = [Time1T3; Time2T3]; 

numbinsall=round(numbinsall)
% All have same axis and bin sizes 
all_data_combinedt=[all_datat;all_data1t;all_data2t;all_data3t];
bin_edges_allt=linspace(min(all_data_combinedt), max(all_data_combinedt), numbinsall + 1);

% Calculate bin edges uniformly for the combined data range
bin_edgest = linspace(min(all_datat), max(all_datat), numbin0 + 1);
bin_edges1t = linspace(min(all_data1t), max(all_data1t), numbin1 + 1);
bin_edges2t = linspace(min(all_data2t), max(all_data2t), numbin2 + 1);
bin_edges3t = linspace(min(all_data3t), max(all_data3t), numbin3 + 1);
%% T
% Combine data to calculate bin edges based on the entire range
all_data=[]; alldata1=[]; alldata2=[]; alldata3=[];
all_data = [Vmax1T0; Vmax2T0];  % Use relevant datasets for each subplot
all_data1 = [Vmax1T1_edit; Vmax2T1]; 
all_data2 = [Vmax1T2_edit; Vmax2T2_edit]; 
all_data3 = [Vmax1T3_edit; Vmax2T3_edit]; 

% All have same axis and bin sizes 
all_data_combined=[];
all_data_combined=[all_data;all_data1;all_data2;all_data3];
bin_edges_all=linspace(min(all_data_combined), max(all_data_combined), numbinsall + 1); % change num of bins 

% Calculate bin edges uniformly for the combined data range
bin_edges = linspace(min(all_data), max(all_data), numbin0 + 1);
bin_edges1 = linspace(min(all_data1), max(all_data1), numbin1 + 1);
bin_edges2 = linspace(min(all_data2), max(all_data2), numbin2 + 1);
bin_edges3 = linspace(min(all_data3), max(all_data3), numbin3 + 1);
%%
% Plot the histograms using the same bin edges
%figure;
subplot(2,2,1)
% Histogram 1
h1 = histogram(Vmax1T0, bin_edges, 'Normalization', 'count', ...
    'EdgeColor', 'none', 'FaceColor', initialc);
hold on;
xline(MeanT0(:,1),'Color',initialc)

% Histogram 2
h2 = histogram(Vmax2T0, bin_edges, 'Normalization', 'count', ...
    'EdgeColor', 'none', 'FaceColor', csmc);
xline(MeanT0(:,2),'Color',csmc)
title('T0 Bin');

subplot(2,2,2)
% Histogram 1
h3 = histogram(Vmax1T1_edit, bin_edges1, 'Normalization', 'count', ...
    'EdgeColor', 'none', 'FaceColor', initialc);
hold on;

% Histogram 2
h4 = histogram(Vmax2T1, bin_edges1, 'Normalization', 'count', ...
    'EdgeColor', 'none', 'FaceColor', csmc);
title('T1 Bin');


subplot(2,2,3)
% Histogram 1
h5 = histogram(Vmax1T2, bin_edges2, 'Normalization', 'count', ...
    'EdgeColor', 'none', 'FaceColor', initialc);
hold on;

% Histogram 2
h6 = histogram(Vmax2T2, bin_edges2, 'Normalization', 'count', ...
    'EdgeColor', 'none', 'FaceColor', csmc);
title('T2 Bin');

subplot(2,2,4)
% Histogram 1
h7 = histogram(Vmax1T3, bin_edges3, 'Normalization', 'count', ...
    'EdgeColor', 'none', 'FaceColor', initialc);
hold on;

% Histogram 2
h8 = histogram(Vmax2T3, bin_edges3, 'Normalization', 'count', ...
    'EdgeColor', 'none', 'FaceColor', csmc);
title('T3 Bin');

% Set axis labels
xlabel('Max Outward Velocity');
ylabel('Probability');

%% Using all the same X axis and bin number 
% Plot the histograms using the same bin edges
figure;
set(gcf,'Position',[100,100,400,800])
subplot(4,1,1)
% Histogram 1
h1 = histogram(Vmax1T0, bin_edges_all, 'Normalization', 'count', ...
    'EdgeColor', 'w', 'FaceColor', initialc);
hold on;
xline(MeanT0(:,1),'Color',initialc)

% Histogram 2
h2 = histogram(Vmax2T0, bin_edges_all, 'Normalization', 'count', ...
    'EdgeColor', 'w', 'FaceColor', csmc);
xline(MeanT0(:,2),'Color',csmc)
title('T0 Bin');


subplot(4,1,2)
% Histogram 1
h3 = histogram(Vmax1T1_edit, bin_edges_all, 'Normalization', 'count', ...
    'EdgeColor', 'w', 'FaceColor', initialc);
hold on;
xline(MeanT1(:,1),'Color',initialc)

% Histogram 2
h4 = histogram(Vmax2T1, bin_edges_all, 'Normalization', 'count', ...
    'EdgeColor', 'w', 'FaceColor', csmc);
title('T1 Bin');
xline(MeanT1(:,2),'Color',csmc)

subplot(4,1,3)
% Histogram 1
h5 = histogram(Vmax1T2_edit, bin_edges_all, 'Normalization', 'count', ...
    'EdgeColor', 'w', 'FaceColor', initialc);
hold on;
xline(MeanT2(:,1),'Color',initialc)

% Histogram 2
h6 = histogram(Vmax2T2_edit, bin_edges_all, 'Normalization', 'count', ...
    'EdgeColor', 'w', 'FaceColor', csmc);
title('T2 Bin');
xline(MeanT2(:,2),'Color',csmc)

subplot(4,1,4)
% Histogram 1
h7 = histogram(Vmax1T3_edit, bin_edges_all, 'Normalization', 'count', ...
    'EdgeColor', 'w', 'FaceColor', initialc);
hold on;
xline(MeanT3(:,1),'Color',initialc)

% Histogram 2
h8 = histogram(Vmax2T3_edit, bin_edges_all, 'Normalization', 'count', ...
    'EdgeColor', 'w', 'FaceColor', csmc);
title('T3 Bin');
xline(MeanT3(:,2),'Color',csmc)



% Set axis labels
xlabel('Max Outward Velocity');
ylabel('Probability');

%% TIME
% Plot the histograms using the same bin edges
figure;
set(gcf,'Position',[100,100,400,800])
subplot(4,1,1)
% Histogram 1
h1 = histogram(Time1T0, bin_edges_allt, 'Normalization', 'count', ...
    'EdgeColor', 'w', 'FaceColor', initialc);
hold on;
xline(MeanT0t(:,1),'Color',initialc)

% Histogram 2
h2 = histogram(Time2T0, bin_edges_allt, 'Normalization', 'count', ...
    'EdgeColor', 'w', 'FaceColor', csmc);
xline(MeanT0t(:,2),'Color',csmc)
title('T0 Bin');


subplot(4,1,2)
% Histogram 1
h3 = histogram(Time1T1, bin_edges_allt, 'Normalization', 'count', ...
    'EdgeColor', 'w', 'FaceColor', initialc);
hold on;
xline(MeanT1t(:,1),'Color',initialc)

% Histogram 2
h4 = histogram(Time2T1, bin_edges_allt, 'Normalization', 'count', ...
    'EdgeColor', 'w', 'FaceColor', csmc);
title('T1 Bin');
xline(MeanT1t(:,2),'Color',csmc)

subplot(4,1,3)
% Histogram 1
h5 = histogram(Time1T2, bin_edges_allt, 'Normalization', 'count', ...
    'EdgeColor', 'w', 'FaceColor', initialc);
hold on;
xline(MeanT2t(:,1),'Color',initialc)

% Histogram 2
h6 = histogram(Time2T2, bin_edges_allt, 'Normalization', 'count', ...
    'EdgeColor', 'w', 'FaceColor', csmc);
title('T2 Bin');
xline(MeanT2t(:,2),'Color',csmc)

subplot(4,1,4)
% Histogram 1
h7 = histogram(Time1T3, bin_edges_allt, 'Normalization', 'count', ...
    'EdgeColor', 'w', 'FaceColor', initialc);
hold on;
xline(MeanT3t(:,1),'Color',initialc)

% Histogram 2
h8 = histogram(Time2T3, bin_edges_allt, 'Normalization', 'count', ...
    'EdgeColor', 'w', 'FaceColor', csmc);
title('T3 Bin');
xline(MeanT3t(:,2),'Color',csmc)



% Set axis labels
xlabel('Time from CSM (s)');
ylabel('Count');

%% Find NaN values 
rows_with_nan0=any(isnan(KinVarsT0),2);
KinVarsT0_edit=KinVarsT0;
KinVarsT0_edit(rows_with_nan0,:)=[];

rows_with_nan1=any(isnan(KinVarsT1),2);
KinVarsT1_edit=KinVarsT1;
KinVarsT1_edit(rows_with_nan1,:)=[];

rows_with_nan2=any(isnan(KinVarsT2),2);
KinVarsT2_edit=KinVarsT2;
KinVarsT2_edit(rows_with_nan2,:)=[];

rows_with_nan3=any(isnan(KinVarsT3),2);
KinVarsT3_edit=KinVarsT3;
KinVarsT3_edit(rows_with_nan3,:)=[];

%% Find values out of range 
rows_to_delete0=any(KinVarsT0_edit>80,2);
KinVarsT0_edit(rows_to_delete0,:)=[];

rows_to_delete1=any(KinVarsT1_edit>80,2);
KinVarsT1_edit(rows_to_delete1,:)=[];

rows_to_delete2=any(KinVarsT2_edit>80,2);
KinVarsT2_edit(rows_to_delete2,:)=[];

rows_to_delete3=any(KinVarsT3_edit>80,2);
KinVarsT3_edit(rows_to_delete3,:)=[];
%%
rows_to_delete0=any(KinVarsT0_edit<-80,2);
KinVarsT0_edit(rows_to_delete0,:)=[];

rows_to_delete1=any(KinVarsT1_edit<-80,2);
KinVarsT1_edit(rows_to_delete1,:)=[];

rows_to_delete2=any(KinVarsT2_edit<-80,2);
KinVarsT2_edit(rows_to_delete2,:)=[];

rows_to_delete3=any(KinVarsT3_edit<-80,2);
KinVarsT3_edit(rows_to_delete3,:)=[];
%% For reaches with no decel before 
rows_to_delete0=any(KinVarsT0_edit(:,1)>0,2);
KinVarsT0_edit(rows_to_delete0,:)=[];
%% Get slope for each animal and each time bin 
an=fieldnames(ExcitRR)
for a=1:length(an)

    decelT1=ExcitRR.(an{a}).AccelMinT1;
    accelT1=ExcitRR.(an{a}).AccelMaxT1;

    decelT2=ExcitRR.(an{a}).AccelMinT2;
    accelT2=ExcitRR.(an{a}).AccelMaxT2;

    decelT3=ExcitRR.(an{a}).AccelMinT3;
    accelT3=ExcitRR.(an{a}).AccelMaxT3;

   
    mdl1 = fitlm(decelT1, accelT1);
    mdl2 = fitlm(decelT2, accelT2);
    mdl3 = fitlm(decelT3, accelT3);

    [r1, p1] = corr(decelT1,accelT1,'rows','Complete');
    [r2, p2] = corr(decelT2,accelT2,'rows','Complete'); 
    [r3, p3] = corr(decelT3,accelT3,'rows','Complete'); 

  
    slope1=mdl1.Coefficients.Estimate(2);
    slope2=mdl2.Coefficients.Estimate(2);
    slope3=mdl3.Coefficients.Estimate(2);

    ExcitRR.(an{a}).slope1=slope1;
    ExcitRR.(an{a}).slope2=slope2;
    ExcitRR.(an{a}).slope3=slope3;

    ExcitRR.(an{a}).r1=r1;
    ExcitRR.(an{a}).r2=r2;
    ExcitRR.(an{a}).r3=r3;

    ExcitRR.(an{a}).p1=p1;
    ExcitRR.(an{a}).p2=p2;
    ExcitRR.(an{a}).p3=p3;

end
%%
at= fieldnames(NoT0Stimvel);
for b=1:length(at)
    decelT0=NoT0Stimvel.(at{b}).AccelMinT0;
    accelT0=NoT0Stimvel.(at{b}).AccelMaxT0;

    mdl0 = fitlm(decelT0, accelT0);
    [r0, p0] = corr(decelT0,accelT0,'rows','Complete');
    slope0=mdl0.Coefficients.Estimate(2);

     NoT0Stimvel.(at{b}).slope0=slope0;
     NoT0Stimvel.(at{b}).r0= r0;
     NoT0Stimvel.(at{b}).p0= p0;

end
%% Plot control and CSM -ACCEL
% positive correlation

%Accel vars 
x0=KinVarsT0_edit(:,1);x1=KinVarsT1_edit(:,1);x2=KinVarsT2_edit(:,1);x3=KinVarsT3_edit(:,1);
y0=KinVarsT0_edit(:,2);y1=KinVarsT1_edit(:,2);y2=KinVarsT2_edit(:,2);y3=KinVarsT3_edit(:,2);


%CorrCoef
[r0, p0] = corr(x0,y0,'rows','Complete'); 
[r1, p1] = corr(x1,y1,'rows','Complete'); 
[r2, p2] = corr(x2,y2,'rows','Complete'); 
[r3, p3] = corr(x3,y3,'rows','Complete'); 

% Fit a linear model
mdl0 = fitlm(x0, y0);
mdl1 = fitlm(x1, y1);
mdl2 = fitlm(x2, y2);
mdl3 = fitlm(x3, y3);

% Sort the x data for plotting purposes
[x_sorted0, idx0] = sort(x0); 
y_sorted0 = y0(idx0); 

[x_sorted1, idx1] = sort(x1); 
y_sorted1 = y1(idx1); 

[x_sorted2, idx2] = sort(x2); 
y_sorted2 = y2(idx2); 

[x_sorted3, idx3] = sort(x3); 
y_sorted3 = y3(idx3); 

% Predict the response and confidence intervals for sorted x data
[y_pred0, y_ci0] = predict(mdl0, x_sorted0, 'Alpha', 0.05); % 95% confidence intervals
[y_pred1, y_ci1] = predict(mdl1, x_sorted1, 'Alpha', 0.05);
[y_pred2, y_ci2] = predict(mdl2, x_sorted2, 'Alpha', 0.05);
[y_pred3, y_ci3] = predict(mdl3, x_sorted3, 'Alpha', 0.05);

figure;
% scatter(x0, y0,20,'filled','MarkerEdgeColor','k','MarkerFaceColor',T0color);
% hold on;
% scatter(x1, y1,20,'filled','MarkerEdgeColor','k','MarkerFaceColor',T1color);
% hold on;
% scatter(x2, y2,20,'filled','MarkerEdgeColor','k','MarkerFaceColor',T2color);
% hold on;
% scatter(x3, y3,20,'filled','MarkerEdgeColor','k','MarkerFaceColor',T3color);
% hold on;
% Plot the fitted regression line
plot(x_sorted0, y_pred0,  'Color',T0color, 'LineWidth', 1.5); % Red regression line
hold on
plot(x_sorted1, y_pred1,  'Color',T1color, 'LineWidth', 1.5); % Red regression line
hold on
plot(x_sorted2, y_pred2,  'Color',T2color, 'LineWidth', 1.5); % Red regression line
hold on
plot(x_sorted3, y_pred3,  'Color',T3color, 'LineWidth', 1.5); % Red regression line
hold on
% Plot the confidence intervals
fill([x_sorted0; flipud(x_sorted0)], [y_ci0(:,1); flipud(y_ci0(:,2))], ...
    [0.6758 0.8438 0.8984], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 

hold on
% Plot the confidence intervals
fill([x_sorted1; flipud(x_sorted1)], [y_ci1(:,1); flipud(y_ci1(:,2))], ...
    [0.2734 0.5078 0.7031], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on
% Plot the confidence intervals
fill([x_sorted2; flipud(x_sorted2)], [y_ci2(:,1); flipud(y_ci2(:,2))], ...
    [0 0 0.8008], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on
% Plot the confidence intervals
fill([x_sorted3; flipud(x_sorted3)], [y_ci3(:,1); flipud(y_ci3(:,2))], ...
    [0.0977 0.0977 0.4375], 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 



text(-5,39, sprintf('T0 Bin r = %.2f  p-value: [%.2f]', r0, p0), 'FontSize', 8, 'VerticalAlignment', 'top');
text(-5,37, sprintf('T1 Bin r = %.2f  p-value: [%.2f]', r1, p1), 'FontSize', 8, 'VerticalAlignment', 'top');
text(-5,35, sprintf('T2 Bin r = %.2f  p-value: [%.2f]', r2, p2), 'FontSize', 8, 'VerticalAlignment', 'top');
text(-5,33, sprintf('T3 Bin r = %.2f  p-value: [%.2f]', r3, p3), 'FontSize', 8, 'VerticalAlignment', 'top');

%unity_line = plot(xlim, xlim, 'k--'); % Black dashed unity line

% Adding labels and title
xlabel('AccelMin');
ylabel('AccelMax');
% xlim([-70 0])
% ylim([0 40])
set(gca,'XDir','reverse')
title('Accel Decel during CSM');
grid off;
hold off;

%% CDF Kinematic PLots 
%% Different way to calculate CDF 


% % Initialize variables to store all CDFs
% all_cdfs_t0 = [];
% all_cdfs_t1= []; all_cdfs_t3 = [];
% all_cdfs_t2= []; 
% 
% global_min=inf; global_max=-inf; 
% % Loop through each animal to determine the global min and max
% for a = 1:length(an)
% 
%     t0 = ExcitRR.(an{a}).T0Kin.tort(:,1);
%     t1 = ExcitRR.(an{a}).T1Kin.tort(:,1);
%     t2 = ExcitRR.(an{a}).T2Kin.tort(:,1);
%     t3 = ExcitRR.(an{a}).T3Kin.tort(:,1);
% 
% 
%     % Update global min and max
%     global_min = min(global_min, min([t0; t1; t2; t3]));
%     global_max = max(global_max, max([t0; t1; t2; t3]));
% end
% 
% % Define common x_values based on global min and max
% x_values = linspace(global_min, global_max, 1000);
% %%
% % Plot CDFs for each animal and interpolate to common x_values
% figure;
% hold on;
% 
% % Loop through each animal
% for a = 1:length(an)
% 
%     t0t = ExcitRR.(an{a}).T0Kin.tort(:,1);
%     t1t = ExcitRR.(an{a}).T1Kin.tort(:,1);
%     t2t = ExcitRR.(an{a}).T2Kin.tort(:,1);
%     t3t = ExcitRR.(an{a}).T3Kin.tort(:,1);
% 
% 
%     % Sort data and calculate CDF
%     sorted_data_t0 = sort(t0t);
%     cdf_values_t0 = (1:length(sorted_data_t0)) / length(sorted_data_t0);
% 
%     sorted_data_t1 = sort(t1t);
%     cdf_values_t1 = (1:length(sorted_data_t1)) / length(sorted_data_t1);
% 
%     sorted_data_t2 = sort(t2t);
%     cdf_values_t2 = (1:length(sorted_data_t2)) / length(sorted_data_t2);
% 
%     sorted_data_t3 = sort(t3t);
%     cdf_values_t3 = (1:length(sorted_data_t3)) / length(sorted_data_t3);
% 
%     % Plot the CDF for the current animal
%    % Use the same figure for all CDFs
%     plot(sorted_data_t0, cdf_values_t0, 'LineWidth', 1.5,'Color',T0color);
%     hold on; % Hold the plot for the next animal's CDF
%     plot(sorted_data_t1, cdf_values_t1, 'LineWidth', 1.5,'Color',T1color);
%     hold on
%     plot(sorted_data_t2, cdf_values_t2, 'LineWidth', 1.5,'Color',T2color);
%     hold on
%     plot(sorted_data_t3, cdf_values_t3, 'LineWidth', 1.5,'Color',T3color);
% 
% 
%      % Interpolate CDFs to the common x_values
%     interpolated_cdf_t0 = interp1(sorted_data_t0, cdf_values_t0, x_values, 'linear', 'extrap');
%     interpolated_cdf_t1 = interp1(sorted_data_t1, cdf_values_t1, x_values, 'linear', 'extrap');
%     interpolated_cdf_t2 = interp1(sorted_data_t2, cdf_values_t2, x_values, 'linear', 'extrap');
%     interpolated_cdf_t3 = interp1(sorted_data_t3, cdf_values_t3, x_values, 'linear', 'extrap');
% 
% 
%     % Ensure CDF values are between 0 and 1
%     interpolated_cdf_t0 = min(max(interpolated_cdf_t0, 0), 1);
%     interpolated_cdf_t1 = min(max(interpolated_cdf_t1, 0), 1);
%     interpolated_cdf_t2 = min(max(interpolated_cdf_t2, 0), 1);
%     interpolated_cdf_t3 = min(max(interpolated_cdf_t3, 0), 1);
% 
%     % Store interpolated CDFs
%     all_cdfs_t0 = [all_cdfs_t0; interpolated_cdf_t0];
%     all_cdfs_t1 = [all_cdfs_t1; interpolated_cdf_t1];
%     all_cdfs_t2 = [all_cdfs_t2; interpolated_cdf_t2];
%     all_cdfs_t3 = [all_cdfs_t3; interpolated_cdf_t3];
% 
% end
% 
% 
% % Calculate the mean CDF across all animals
% mean_cdf_t0 = mean(all_cdfs_t0, 1);
% mean_cdf_t1 = mean(all_cdfs_t1, 1);
% mean_cdf_t2 = mean(all_cdfs_t2, 1);
% mean_cdf_t3 = mean(all_cdfs_t3, 1);
% hold on
% 
% % Plot the mean CDF on the same plot
% plot(x_values, mean_cdf_t0, 'LineWidth', 5,'Color',T0color); % Plot in black with a thicker line
% hold on
% plot(x_values, mean_cdf_t1, 'LineWidth', 5,'Color',T1color); 
% hold on
% plot(x_values, mean_cdf_t2, 'LineWidth', 5,'Color',T2color); 
% hold on
% plot(x_values, mean_cdf_t3, 'LineWidth', 5,'Color',T3color); 
% 
% 
% % Customize the plot
% %ylim([0 1])
% %xlim([-20 1])
% xlim([0 100])
% %set(gca, 'XDir', 'reverse'); % Reverse the x-axis direction
% xlabel('Tortuosity');
% ylabel('CDF');
% title('CDFs for All Animals and Mean CDF');
% 
% % [h, p, k] = kstest2(Ctrl_VmaxP,Plus1S_VmaxP);
% % [h1, p1, k1] = kstest2(Ctrl_VmaxP,Plus1F_VmaxP);
% % legend( ['p-value (ctrl,Plus1S) = ' num2str(p, '%.4f');'p-value (ctrl,Plus1F) = ' num2str(p1, '%.4f')], 'Location', 'SouthEast');
% %%
% 
% % Perform K-S tests between different conditions
% [data_t0, data_t1, data_t2, data_t3] = deal([]);
% for a = 1:length(an)
% 
%     data_t0 = [data_t0; ExcitRR.(an{a}).T0Kin.tort(:,1)];
%     data_t1 = [data_t1; ExcitRR.(an{a}).T1Kin.tort(:,1)];
%     data_t2 = [data_t2; ExcitRR.(an{a}).T2Kin.tort(:,1)];
%     data_t3 = [data_t3; ExcitRR.(an{a}).T3Kin.tort(:,1)];
% 
% end
% 
% % Run K-S tests between different conditions
% % [~, p_ctrl_t0] = kstest2(data_t0, data_t1);
% % [~, p_ctrl_sf] = kstest2(data_t0, data_t2);
% % [~, p_ctrl_p1] = kstest2(data_t0, data_t3);
% % [~, p_ctrl_p1f] = kstest2(data_t1, data_t0);
% % [~, p_ss_sf] = kstest2(data_t1, data_t2);
% % [~, p_ss_p1] = kstest2(data_t1, data_t2);
% % [~, p_ss_p1f] = kstest2(data_t2, data_t0);
% % [~, p_sf_p1] = kstest2(data_t2, data_t2);
% % [~, p_sf_p1f] = kstest2(data_t2, data_t3);
% % 
% 
% % Display results
% % fprintf('K-S Test Results:\n');
% % fprintf('Control vs StimSuc: p = %.4f\n', p_ctrl_ss);
% % fprintf('Control vs StimFail: p = %.4f\n', p_ctrl_sf);
% % fprintf('Control vs Plus1: p = %.4f\n', p_ctrl_p1);
% % fprintf('Control vs Plus1Fail: p = %.4f\n', p_ctrl_p1f);
% % fprintf('StimSuc vs StimFail: p = %.4f\n', p_ss_sf);
% % fprintf('StimSuc vs Plus1: p = %.4f\n', p_ss_p1);
% % fprintf('StimSuc vs Plus1Fail: p = %.4f\n', p_ss_p1f);
% % fprintf('StimFail vs Plus1: p = %.4f\n', p_sf_p1);
% % fprintf('StimFail vs Plus1Fail: p = %.4f\n', p_sf_p1f);
% % fprintf('Plus1 vs Plus1Fail: p = %.4f\n', p_p1_p1f);
% 
% % legend( ['p-value (ctrl,Plus1S) = ' num2str(p_ctrl_p1, '%.4f');'p-value (ctrl,Plus1F) = ' num2str( p_ctrl_p1f, '%.4f')], 'Location', 'NorthWest')
% % %legend( ['p-value (ctrl,StimS) = ' num2str(p_ctrl_ss, '%.4f');'p-value (ctrl,StimF) = ' num2str( p_ctrl_sf, '%.4f')], 'Location', 'NorthWest')
% hold off; % Release the plot
% 
% %% Ksdensity plots 
% T0_Kin=[];
% T1_Kin=[]; 
% T2_Kin=[]; 
% T3_Kin=[];
% 
% figure;
% hold on;
% 
% % Loop through each animal
% for a = 1:length(an)
%     % Extract data for the current animal
%     t0 = ExcitRR.(an{a}).T0Kin.tort(:,1);
%     t1 = ExcitRR.(an{a}).T1Kin.tort(:,1);
%     t2 = ExcitRR.(an{a}).T2Kin.tort(:,1);
%     t3 = ExcitRR.(an{a}).T3Kin.tort(:,1);
% 
%     T0_Kin=[T0_Kin;t0];
%     T1_Kin=[T1_Kin;t1];
%     T2_Kin=[T2_Kin;t2];
%     T3_Kin=[T3_Kin;t3];
% 
%     [f1,xi1]=ksdensity(t0);
%     [f2,xi2]=ksdensity(t1);
%     [f3,xi3]=ksdensity(t2);
%     [f4,xi4]=ksdensity(t3);
% 
%     % %plot individual animal ksdensirt 
%     % plot(xi1, f1, 'LineWidth', 1.5,'Color',T0color);
%     % hold on; 
%     % plot(xi2, f2, 'LineWidth', 1.5,'Color',T1color);
%     % hold on
%     % plot(xi3, f3, 'LineWidth', 1.5,'Color',T2color);
%     % hold on
%     % plot(xi4, f4, 'LineWidth', 1.5,'Color',T3color);
%     % 
% end
% 
% %average ksdensity across all anims 
% [f1_comb,xi1_comb]=ksdensity(T0_Kin);
% [f2_comb,xi2_comb]=ksdensity(T1_Kin);
% [f3_comb,xi3_comb]=ksdensity(T2_Kin);
% [f4_comb,xi4_comb]=ksdensity(T3_Kin);
% 
% t0mean=mean(T0_Kin);
% t1mean=mean(T1_Kin);
% t2mean=mean(T2_Kin);
% t3mean=mean(T3_Kin);
% 
% % Plot the mean ksdensity on the same plot
% plot(xi1_comb, f1_comb, 'LineWidth', 4,'Color',T0color);
% hold on; 
% plot(xi2_comb, f2_comb, 'LineWidth', 4,'Color',T1color);
% hold on
% plot(xi3_comb, f3_comb, 'LineWidth', 4,'Color',T2color);
% hold on
% plot(xi4_comb, f4_comb, 'LineWidth', 4,'Color',T3color);
% 
% 
% % Customize the plot
% 
% % Mark the mean on the plot
% % xline(cmean)
% % xline(ssmean,'Color',sscolor)
% % xline(sfmean,'Color',sfcolor)
% % xline(p1mean,'Color',p1color)
% % xline(p1fmean,'Color',p1fcolor)
% 
% xlabel('Tort');
% xlim([0 50])
% %xlim([0 6])
% %xlim([-30 5])
% 
% title('KSDensity for All Animals and Mean KSDensity');