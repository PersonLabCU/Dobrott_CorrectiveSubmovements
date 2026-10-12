%% Treat all ReReaches as T0 reaches - Baseline
% Load in Structure with TimeBins already saved 
T0_divide_BL=[];
thresh=0.1;
an=fieldnames(ExcitRR_BL);
for a=1:length(an)
   
    t0r=ExcitRR_BL.(an{a}).RRT0;
    if isempty(t0r)
        ReReach0=[];
        StimRR0=[];
        MaxVel0=[];
        ReachesWithoutSustainedJerk0=[];
    else
        t0r=accel_toreaches_jerk(t0r);
        [ReReach0, StimRR0, MaxVel0,ReachesWithoutSustainedJerk0] = CSMEarlyThresh_jerk_edit7(t0r,thresh);
    end
    
    t1r=ExcitRR_BL.(an{a}).RRT1;
    t1r=accel_toreaches_jerk(t1r);
    [ReReach1, StimRR1, MaxVel1,ReachesWithoutSustainedJerk1] = CSMEarlyThresh_jerk_edit7(t1r,thresh); % only looking for max and no threshold 
   
    t2r=ExcitRR_BL.(an{a}).RRT2;
    t2r=accel_toreaches_jerk(t2r);
    [ReReach2, StimRR2, MaxVel2,ReachesWithoutSustainedJerk2] = CSMEarlyThresh_jerk_edit7(t2r,thresh);

    t3r=ExcitRR_BL.(an{a}).RRT3;
    if isempty(t3r)
        ReReach3=[];
        StimRR3=[];
        MaxVel3=[];
        ReachesWithoutSustainedJerk3=[];
    else
        t3r=accel_toreaches_jerk(t3r);
        [ReReach3, StimRR3, MaxVel3,ReachesWithoutSustainedJerk3] = CSMEarlyThresh_jerk_edit7(t3r,thresh);
    end
   
    
    T0_divide_BL.(an{a}).ReReach0=ReReach0;
    T0_divide_BL.(an{a}).StimRR0=StimRR0;
    T0_divide_BL.(an{a}).NoJerk0=ReachesWithoutSustainedJerk0;

    T0_divide_BL.(an{a}).ReReach1=ReReach1;
    T0_divide_BL.(an{a}).StimRR1=StimRR1;
    T0_divide_BL.(an{a}).NoJerk1=ReachesWithoutSustainedJerk1;

    T0_divide_BL.(an{a}).ReReach2=ReReach2;
    T0_divide_BL.(an{a}).StimRR2=StimRR2;
    T0_divide_BL.(an{a}).NoJerk2=ReachesWithoutSustainedJerk2;
 
    T0_divide_BL.(an{a}).ReReach3=ReReach3;
    T0_divide_BL.(an{a}).StimRR3=StimRR3;
    T0_divide_BL.(an{a}).NoJerk3=ReachesWithoutSustainedJerk3;


    % Save to ExcitRR_Sham
    % Run through all ReReaches - easier than doing all reaches separated
    AllRR=ExcitRR_BL.(an{a}).ReReaches;
    AllRR=accel_toreaches_jerk(AllRR);
    [ReReachAll, StimRRAll, MaxVelAll,ReachesWithoutSustainedJerkAll] = CSMEarlyThresh_jerk_edit7(AllRR,thresh);
    ExcitRR_BL.(an{a}).StimRRAll_thresh=StimRRAll;
    ExcitRR_BL.(an{a}).ReReachAll_thresh=ReReachAll;

end

%% Treat all ReReaches as T0 reaches - Adapt stimulation block
% Load in Structure with TimeBins already saved 
T0_divide_AD=[];
thresh=0.1;
an=fieldnames(ExcitRR_AD);
for a=1:length(an)
    
    t0r=ExcitRR_AD.(an{a}).RRT0;
    if isempty(t0r)
        ReReach0=[];
        StimRR0=[];
        MaxVel0=[];
        ReachesWithoutSustainedJerk0=[];
    else
        t0r=accel_toreaches_jerk(t0r);
        [ReReach0, StimRR0, MaxVel0,ReachesWithoutSustainedJerk0] = CSMEarlyThresh_jerk_edit7(t0r,thresh);
    end


    t1r=ExcitRR_AD.(an{a}).RRT1;
    if isempty(t1r)
        ReReach1=[];
        StimRR1=[];
        MaxVel1=[];
        ReachesWithoutSustainedJerk1=[];
    else
        t1r=accel_toreaches_jerk(t1r);
        [ReReach1, StimRR1, MaxVel1,ReachesWithoutSustainedJerk1] = CSMEarlyThresh_jerk_edit7(t1r,thresh);
    end
   
    t2r=ExcitRR_AD.(an{a}).RRT2;
    if isempty(t2r)
        ReReach2=[];
        StimRR2=[];
        MaxVel2=[];
        ReachesWithoutSustainedJerk2=[];
    else
        t2r=accel_toreaches_jerk(t2r);
        [ReReach2, StimRR2, MaxVel2,ReachesWithoutSustainedJerk2] = CSMEarlyThresh_jerk_edit7(t2r,thresh);
    end

    t3r=ExcitRR_AD.(an{a}).RRT3;
    if isempty(t3r)
        ReReach3=[];
        StimRR3=[];
        MaxVel3=[];
        ReachesWithoutSustainedJerk3=[];
    else
        t3r=accel_toreaches_jerk(t3r);
        [ReReach3, StimRR3, MaxVel3,ReachesWithoutSustainedJerk3] = CSMEarlyThresh_jerk_edit7(t3r,thresh);
    end
    
    T0_divide_AD.(an{a}).ReReach0=ReReach0;
    T0_divide_AD.(an{a}).StimRR0=StimRR0;
    T0_divide_AD.(an{a}).NoJerk0=ReachesWithoutSustainedJerk0;

    T0_divide_AD.(an{a}).ReReach1=ReReach1;
    T0_divide_AD.(an{a}).StimRR1=StimRR1;
    T0_divide_AD.(an{a}).NoJerk1=ReachesWithoutSustainedJerk1;

    T0_divide_AD.(an{a}).ReReach2=ReReach2;
    T0_divide_AD.(an{a}).StimRR2=StimRR2;
    T0_divide_AD.(an{a}).NoJerk2=ReachesWithoutSustainedJerk2;
 
    T0_divide_AD.(an{a}).ReReach3=ReReach3;
    T0_divide_AD.(an{a}).StimRR3=StimRR3;
    T0_divide_AD.(an{a}).NoJerk3=ReachesWithoutSustainedJerk3;


    % Save to ExcitRR_Sham
    % Run through all ReReaches - easier than doing all reaches separated
    AllRR=ExcitRR_AD.(an{a}).ReReaches;
    AllRR=accel_toreaches_jerk(AllRR);
    [ReReachAll, StimRRAll, MaxVelAll,ReachesWithoutSustainedJerkAll] = CSMEarlyThresh_jerk_edit7(AllRR,thresh);
    ExcitRR_AD.(an{a}).StimRRAll_thresh=StimRRAll;
    ExcitRR_AD.(an{a}).ReReachAll_thresh=ReReachAll;

end


%% Treat all ReReaches as T0 reaches - Washout all
% Load in Structure with TimeBins already saved 
T0_divide_WO=[];
thresh=0.1;
an=fieldnames(ExcitRR_WO);
for a=1:length(an)
    
    t0r=ExcitRR_WO.(an{a}).RRT0;
    if isempty(t0r)
        ReReach0=[];
        StimRR0=[];
        MaxVel0=[];
        ReachesWithoutSustainedJerk0=[];
    else
        t0r=accel_toreaches_jerk(t0r);
        [ReReach0, StimRR0, MaxVel0,ReachesWithoutSustainedJerk0] = CSMEarlyThresh_jerk_edit7(t0r,thresh);
    end


    t1r=ExcitRR_WO.(an{a}).RRT1;
    if isempty(t1r)
        ReReach1=[];
        StimRR1=[];
        MaxVel1=[];
        ReachesWithoutSustainedJerk1=[];
    else
        t1r=accel_toreaches_jerk(t1r);
        [ReReach1, StimRR1, MaxVel1,ReachesWithoutSustainedJerk1] = CSMEarlyThresh_jerk_edit7(t1r,thresh);
    end
   
    t2r=ExcitRR_WO.(an{a}).RRT2;
    if isempty(t2r)
        ReReach2=[];
        StimRR2=[];
        MaxVel2=[];
        ReachesWithoutSustainedJerk2=[];
    else
        t2r=accel_toreaches_jerk(t2r);
        [ReReach2, StimRR2, MaxVel2,ReachesWithoutSustainedJerk2] = CSMEarlyThresh_jerk_edit7(t2r,thresh);
    end

    t3r=ExcitRR_WO.(an{a}).RRT3;
    if isempty(t3r)
        ReReach3=[];
        StimRR3=[];
        MaxVel3=[];
        ReachesWithoutSustainedJerk3=[];
    else
        t3r=accel_toreaches_jerk(t3r);
        [ReReach3, StimRR3, MaxVel3,ReachesWithoutSustainedJerk3] = CSMEarlyThresh_jerk_edit7(t3r,thresh);
    end
    
    T0_divide_WO.(an{a}).ReReach0=ReReach0;
    T0_divide_WO.(an{a}).StimRR0=StimRR0;
    T0_divide_WO.(an{a}).NoJerk0=ReachesWithoutSustainedJerk0;

    T0_divide_WO.(an{a}).ReReach1=ReReach1;
    T0_divide_WO.(an{a}).StimRR1=StimRR1;
    T0_divide_WO.(an{a}).NoJerk1=ReachesWithoutSustainedJerk1;

    T0_divide_WO.(an{a}).ReReach2=ReReach2;
    T0_divide_WO.(an{a}).StimRR2=StimRR2;
    T0_divide_WO.(an{a}).NoJerk2=ReachesWithoutSustainedJerk2;
 
    T0_divide_WO.(an{a}).ReReach3=ReReach3;
    T0_divide_WO.(an{a}).StimRR3=StimRR3;
    T0_divide_WO.(an{a}).NoJerk3=ReachesWithoutSustainedJerk3;


    % Save to ExcitRR_Sham
    % Run through all ReReaches - easier than doing all reaches separated
    AllRR=ExcitRR_WO.(an{a}).ReReaches;
    AllRR=accel_toreaches_jerk(AllRR);
    [ReReachAll, StimRRAll, MaxVelAll,ReachesWithoutSustainedJerkAll] = CSMEarlyThresh_jerk_edit7(AllRR,thresh);
    ExcitRR_WO.(an{a}).StimRRAll_thresh=StimRRAll;
    ExcitRR_WO.(an{a}).ReReachAll_thresh=ReReachAll;

end

%% Treat all ReReaches as T0 reaches - Washout early
% Load in Structure with TimeBins already saved 
T0_divide_WOE=[];
thresh=0.1;
an=fieldnames(ExcitRR_WO_E);
for a=1:length(an)
    
    t0r=ExcitRR_WO_E.(an{a}).RRT0;
    if isempty(t0r)
        ReReach0=[];
        StimRR0=[];
        MaxVel0=[];
        ReachesWithoutSustainedJerk0=[];
    else
        t0r=accel_toreaches_jerk(t0r);
        [ReReach0, StimRR0, MaxVel0,ReachesWithoutSustainedJerk0] = CSMEarlyThresh_jerk_edit7(t0r,thresh);
    end


    t1r=ExcitRR_WO_E.(an{a}).RRT1;
    if isempty(t1r)
        ReReach1=[];
        StimRR1=[];
        MaxVel1=[];
        ReachesWithoutSustainedJerk1=[];
    else
        t1r=accel_toreaches_jerk(t1r);
        [ReReach1, StimRR1, MaxVel1,ReachesWithoutSustainedJerk1] = CSMEarlyThresh_jerk_edit7(t1r,thresh);
    end
   
    t2r=ExcitRR_WO_E.(an{a}).RRT2;
    if isempty(t2r)
        ReReach2=[];
        StimRR2=[];
        MaxVel2=[];
        ReachesWithoutSustainedJerk2=[];
    else
        t2r=accel_toreaches_jerk(t2r);
        [ReReach2, StimRR2, MaxVel2,ReachesWithoutSustainedJerk2] = CSMEarlyThresh_jerk_edit7(t2r,thresh);
    end

    t3r=ExcitRR_WO_E.(an{a}).RRT3;
    if isempty(t3r)
        ReReach3=[];
        StimRR3=[];
        MaxVel3=[];
        ReachesWithoutSustainedJerk3=[];
    else
        t3r=accel_toreaches_jerk(t3r);
        [ReReach3, StimRR3, MaxVel3,ReachesWithoutSustainedJerk3] = CSMEarlyThresh_jerk_edit7(t3r,thresh);
    end
    
    T0_divide_WOE.(an{a}).ReReach0=ReReach0;
    T0_divide_WOE.(an{a}).StimRR0=StimRR0;
    T0_divide_WOE.(an{a}).NoJerk0=ReachesWithoutSustainedJerk0;

    T0_divide_WOE.(an{a}).ReReach1=ReReach1;
    T0_divide_WOE.(an{a}).StimRR1=StimRR1;
    T0_divide_WOE.(an{a}).NoJerk1=ReachesWithoutSustainedJerk1;

    T0_divide_WOE.(an{a}).ReReach2=ReReach2;
    T0_divide_WOE.(an{a}).StimRR2=StimRR2;
    T0_divide_WOE.(an{a}).NoJerk2=ReachesWithoutSustainedJerk2;
 
    T0_divide_WOE.(an{a}).ReReach3=ReReach3;
    T0_divide_WOE.(an{a}).StimRR3=StimRR3;
    T0_divide_WOE.(an{a}).NoJerk3=ReachesWithoutSustainedJerk3;


    % Save to ExcitRR_Sham
    % Run through all ReReaches - easier than doing all reaches separated
    AllRR=ExcitRR_WO_E.(an{a}).ReReaches;
    AllRR=accel_toreaches_jerk(AllRR);
    [ReReachAll, StimRRAll, MaxVelAll,ReachesWithoutSustainedJerkAll] = CSMEarlyThresh_jerk_edit7(AllRR,thresh);
    ExcitRR_WO_E.(an{a}).StimRRAll_thresh=StimRRAll;
    ExcitRR_WO_E.(an{a}).ReReachAll_thresh=ReReachAll;

end

%% Treat all ReReaches as T0 reaches - Adapt stimulation block
% Load in Structure with TimeBins already saved 
T0_divide_WOL=[];
thresh=0.1;
an=fieldnames(ExcitRR_WO_L);
for a=1:length(an)
    
    t0r=ExcitRR_WO_L.(an{a}).RRT0;
    if isempty(t0r)
        ReReach0=[];
        StimRR0=[];
        MaxVel0=[];
        ReachesWithoutSustainedJerk0=[];
    else
        t0r=accel_toreaches_jerk(t0r);
        [ReReach0, StimRR0, MaxVel0,ReachesWithoutSustainedJerk0] = CSMEarlyThresh_jerk_edit7(t0r,thresh);
    end


    t1r=ExcitRR_WO_L.(an{a}).RRT1;
    if isempty(t1r)
        ReReach1=[];
        StimRR1=[];
        MaxVel1=[];
        ReachesWithoutSustainedJerk1=[];
    else
        t1r=accel_toreaches_jerk(t1r);
        [ReReach1, StimRR1, MaxVel1,ReachesWithoutSustainedJerk1] = CSMEarlyThresh_jerk_edit7(t1r,thresh);
    end
   
    t2r=ExcitRR_WO_L.(an{a}).RRT2;
    if isempty(t2r)
        ReReach2=[];
        StimRR2=[];
        MaxVel2=[];
        ReachesWithoutSustainedJerk2=[];
    else
        t2r=accel_toreaches_jerk(t2r);
        [ReReach2, StimRR2, MaxVel2,ReachesWithoutSustainedJerk2] = CSMEarlyThresh_jerk_edit7(t2r,thresh);
    end

    t3r=ExcitRR_WO_L.(an{a}).RRT3;
    if isempty(t3r)
        ReReach3=[];
        StimRR3=[];
        MaxVel3=[];
        ReachesWithoutSustainedJerk3=[];
    else
        t3r=accel_toreaches_jerk(t3r);
        [ReReach3, StimRR3, MaxVel3,ReachesWithoutSustainedJerk3] = CSMEarlyThresh_jerk_edit7(t3r,thresh);
    end
    
    T0_divide_WOL.(an{a}).ReReach0=ReReach0;
    T0_divide_WOL.(an{a}).StimRR0=StimRR0;
    T0_divide_WOL.(an{a}).NoJerk0=ReachesWithoutSustainedJerk0;

    T0_divide_WOL.(an{a}).ReReach1=ReReach1;
    T0_divide_WOL.(an{a}).StimRR1=StimRR1;
    T0_divide_WOL.(an{a}).NoJerk1=ReachesWithoutSustainedJerk1;

    T0_divide_WOL.(an{a}).ReReach2=ReReach2;
    T0_divide_WOL.(an{a}).StimRR2=StimRR2;
    T0_divide_WOL.(an{a}).NoJerk2=ReachesWithoutSustainedJerk2;
 
    T0_divide_WOL.(an{a}).ReReach3=ReReach3;
    T0_divide_WOL.(an{a}).StimRR3=StimRR3;
    T0_divide_WOL.(an{a}).NoJerk3=ReachesWithoutSustainedJerk3;


    % Save to ExcitRR_Sham
    % Run through all ReReaches - easier than doing all reaches separated
    AllRR=ExcitRR_WO_L.(an{a}).ReReaches;
    AllRR=accel_toreaches_jerk(AllRR);
    [ReReachAll, StimRRAll, MaxVelAll,ReachesWithoutSustainedJerkAll] = CSMEarlyThresh_jerk_edit7(AllRR,thresh);
    ExcitRR_WO_L.(an{a}).StimRRAll_thresh=StimRRAll;
    ExcitRR_WO_L.(an{a}).ReReachAll_thresh=ReReachAll;

end

%% Now filter through , get total number within pre thresh window - FIRST SESSION- Ctrl
for a = 1:length(an)
    RRtimes1= T0_divide_BL.(an{a}).StimRR1;
    RRtimes2= T0_divide_BL.(an{a}).StimRR2;
    RRtimes3= T0_divide_BL.(an{a}).StimRR3;

    % Keep only values >= -0.08 and <= -0.02
    RRtimes1_filt = RRtimes1(RRtimes1 >= -0.08 & RRtimes1 <= -0.02);
    RRtimes2_filt = RRtimes2(RRtimes2 >= -0.08 & RRtimes2 <= -0.02);
    RRtimes3_filt = RRtimes3(RRtimes3 >= -0.08 & RRtimes3 <= -0.02);

    % (Optional) Save back the filtered values, if desired:
    ExcitRR_BL.(an{a}).StimRR1_filt=  RRtimes1_filt;
    ExcitRR_BL.(an{a}).StimRR2_filt=  RRtimes2_filt;
    ExcitRR_BL.(an{a}).StimRR3_filt=  RRtimes3_filt;

    AllbinsRR=ExcitRR_BL.(an{a}).AllBins;

    combinedData=[];
    combinedData = [AllbinsRR; RRtimes1_filt; RRtimes2_filt;RRtimes3_filt];
    ExcitRR_BL.(an{a}).AllBins_prethresh=combinedData;

end

%% Now filter through ,  ADAPT stim 
for a = 1:length(an)
    RRtimes1= T0_divide_AD.(an{a}).StimRR1;
    RRtimes2= T0_divide_AD.(an{a}).StimRR2;
    RRtimes3= T0_divide_AD.(an{a}).StimRR3;

    % Keep only values >= -0.08 and <= -0.02
    RRtimes1_filt = RRtimes1(RRtimes1 >= -0.08 & RRtimes1 <= -0.02);
    RRtimes2_filt = RRtimes2(RRtimes2 >= -0.08 & RRtimes2 <= -0.02);
    RRtimes3_filt = RRtimes3(RRtimes3 >= -0.08 & RRtimes3 <= -0.02);

    % (Optional) Save back the filtered values, if desired:
    ExcitRR_AD.(an{a}).StimRR1_filt=  RRtimes1_filt;
    ExcitRR_AD.(an{a}).StimRR2_filt=  RRtimes2_filt;
    ExcitRR_AD.(an{a}).StimRR3_filt=  RRtimes3_filt;

    AllbinsRR=ExcitRR_AD.(an{a}).AllBins;

    combinedData=[];
    combinedData = [AllbinsRR; RRtimes1_filt; RRtimes2_filt;RRtimes3_filt];
    ExcitRR_AD.(an{a}).AllBins_prethresh=combinedData;

end

%% Washout all
for a = 1:length(an)
    RRtimes1= T0_divide_WO.(an{a}).StimRR1;
    RRtimes2= T0_divide_WO.(an{a}).StimRR2;
    RRtimes3= T0_divide_WO.(an{a}).StimRR3;

    % Keep only values >= -0.08 and <= -0.02
    RRtimes1_filt = RRtimes1(RRtimes1 >= -0.08 & RRtimes1 <= -0.02);
    RRtimes2_filt = RRtimes2(RRtimes2 >= -0.08 & RRtimes2 <= -0.02);
    RRtimes3_filt = RRtimes3(RRtimes3 >= -0.08 & RRtimes3 <= -0.02);

    % (Optional) Save back the filtered values, if desired:
    ExcitRR_WO.(an{a}).StimRR1_filt=  RRtimes1_filt;
    ExcitRR_WO.(an{a}).StimRR2_filt=  RRtimes2_filt;
    ExcitRR_WO.(an{a}).StimRR3_filt=  RRtimes3_filt;

    AllbinsRR=ExcitRR_WO.(an{a}).AllBins;

    combinedData=[];
    combinedData = [AllbinsRR; RRtimes1_filt; RRtimes2_filt;RRtimes3_filt];
    ExcitRR_WO.(an{a}).AllBins_prethresh=combinedData;

end

%% Washout early 
for a = 1:length(an)
    RRtimes1= T0_divide_WOE.(an{a}).StimRR1;
    RRtimes2= T0_divide_WOE.(an{a}).StimRR2;
    RRtimes3= T0_divide_WOE.(an{a}).StimRR3;

    % Keep only values >= -0.08 and <= -0.02
    RRtimes1_filt = RRtimes1(RRtimes1 >= -0.08 & RRtimes1 <= -0.02);
    RRtimes2_filt = RRtimes2(RRtimes2 >= -0.08 & RRtimes2 <= -0.02);
    RRtimes3_filt = RRtimes3(RRtimes3 >= -0.08 & RRtimes3 <= -0.02);

    % (Optional) Save back the filtered values, if desired:
    ExcitRR_WO_E.(an{a}).StimRR1_filt=  RRtimes1_filt;
    ExcitRR_WO_E.(an{a}).StimRR2_filt=  RRtimes2_filt;
    ExcitRR_WO_E.(an{a}).StimRR3_filt=  RRtimes3_filt;

    AllbinsRR=ExcitRR_WO_E.(an{a}).AllBins;

    combinedData=[];
    combinedData = [AllbinsRR; RRtimes1_filt; RRtimes2_filt;RRtimes3_filt];
    ExcitRR_WO_E.(an{a}).AllBins_prethresh=combinedData;

end

%% Washout late 
for a = 1:length(an)
    RRtimes1= T0_divide_WOL.(an{a}).StimRR1;
    RRtimes2= T0_divide_WOL.(an{a}).StimRR2;
    RRtimes3= T0_divide_WOL.(an{a}).StimRR3;

    % Keep only values >= -0.08 and <= -0.02
    RRtimes1_filt = RRtimes1(RRtimes1 >= -0.08 & RRtimes1 <= -0.02);
    RRtimes2_filt = RRtimes2(RRtimes2 >= -0.08 & RRtimes2 <= -0.02);
    RRtimes3_filt = RRtimes3(RRtimes3 >= -0.08 & RRtimes3 <= -0.02);

    % (Optional) Save back the filtered values, if desired:
    ExcitRR_WO_L.(an{a}).StimRR1_filt=  RRtimes1_filt;
    ExcitRR_WO_L.(an{a}).StimRR2_filt=  RRtimes2_filt;
    ExcitRR_WO_L.(an{a}).StimRR3_filt=  RRtimes3_filt;

    AllbinsRR=ExcitRR_WO_L.(an{a}).AllBins;

    combinedData=[];
    combinedData = [AllbinsRR; RRtimes1_filt; RRtimes2_filt;RRtimes3_filt];
    ExcitRR_WO_L.(an{a}).AllBins_prethresh=combinedData;

end
%% Plot All Groups Together with KS Test

% Define colors
stimcolor = rgb('DodgerBlue'); % ExcitRR
ctrlcolor = rgb('Gray'); % ExcitRR_Ctrl
shamcolor = rgb('Black'); % ExcitRR_Sham (Dashed line)
washecolor=rgb('violet');
washlcolor=rgb('Magenta');

% Initialize storage for KDEs
num_points = 100; % Number of points for KDE interpolation
xi_common = linspace(-0.1, 0.2, num_points); % Common x-axis


% -----------------------------
% Compute KDE for ExcitRR_BL Group
% -----------------------------
all_f_excitrr_ctrl = [];
excitam_excitrr_ctrl = fieldnames(ExcitRR_BL);
num_excitrr_ctrl = length(excitam_excitrr_ctrl);
all_data_excitrr_ctrl = [];

for m = 1:num_excitrr_ctrl
    allbins_excitrr_ctrl = ExcitRR_BL.(excitam_excitrr_ctrl{m}).AllBins_prethresh;
    all_data_excitrr_ctrl = [all_data_excitrr_ctrl; allbins_excitrr_ctrl];  % Collect for KS test
    %[f_excitrr_ctrl, xi_excitrr_ctrl] = ksdensity(allbins_excitrr_ctrl, 'Bandwidth', 0.008);
    [f_excitrr_ctrl, xi_excitrr_ctrl] = ksdensity(allbins_excitrr_ctrl, 'Bandwidth', 0.004);
    f_interp_excitrr_ctrl = interp1(xi_excitrr_ctrl, f_excitrr_ctrl, xi_common, 'linear', 0);
    all_f_excitrr_ctrl = [all_f_excitrr_ctrl; f_interp_excitrr_ctrl];
end


% -----------------------------
% Compute KDE for ExcitRR_AD Group
% -----------------------------
all_f_excitrr = [];
excitam_excitrr = fieldnames(ExcitRR_AD);
num_excitrr = length(excitam_excitrr);
all_data_excitrr = [];

for m = 1:num_excitrr
    allbins_excitrr = ExcitRR_AD.(excitam_excitrr{m}).AllBins_prethresh;
    all_data_excitrr = [all_data_excitrr; allbins_excitrr];  % Collect for KS test
    %[f_excitrr, xi_excitrr] = ksdensity(allbins_excitrr, 'Bandwidth', 0.008);
     [f_excitrr, xi_excitrr] = ksdensity(allbins_excitrr, 'Bandwidth', 0.004);
    f_interp_excitrr = interp1(xi_excitrr, f_excitrr, xi_common, 'linear', 0);
    all_f_excitrr = [all_f_excitrr; f_interp_excitrr];
end


% % -----------------------------
% % Compute KDE for ExcitRR_WO_ EGroup
% % -----------------------------
% all_f_excitrr_woe = [];
% excitam_excitrr_woe = fieldnames(ExcitRR_WO_E);
% num_excitrr_woe = length(excitam_excitrr_woe);
% all_data_excitrr_woe = [];
% 
% for m = 1:num_excitrr_woe
%     allbins_excitrr_woe = ExcitRR_WO_E.(excitam_excitrr_woe{m}).AllBins_prethresh;
%     all_data_excitrr_woe = [all_data_excitrr_woe; allbins_excitrr_woe];  % Collect for KS test
%     [f_excitrr_woe, xi_excitrr_woe] = ksdensity(allbins_excitrr_woe, 'Bandwidth', 0.008);
%     f_interp_excitrr_woe = interp1(xi_excitrr_woe, f_excitrr_woe, xi_common, 'linear', 0);
%     all_f_excitrr_woe = [all_f_excitrr_woe; f_interp_excitrr_woe];
% end

% -----------------------------
% Compute KDE for ExcitRR_WO_ EGroup
% -----------------------------
all_f_excitrr_woe = [];
excitam_excitrr_woe = fieldnames(ExcitRR_WO);
num_excitrr_woe = length(excitam_excitrr_woe);
all_data_excitrr_woe = [];

for m = 1:num_excitrr_woe
    allbins_excitrr_woe = ExcitRR_WO.(excitam_excitrr_woe{m}).AllBins_prethresh;
    all_data_excitrr_woe = [all_data_excitrr_woe; allbins_excitrr_woe];  % Collect for KS test
    %[f_excitrr_woe, xi_excitrr_woe] = ksdensity(allbins_excitrr_woe, 'Bandwidth', 0.008);
    [f_excitrr_woe, xi_excitrr_woe] = ksdensity(allbins_excitrr_woe, 'Bandwidth', 0.004);
    f_interp_excitrr_woe = interp1(xi_excitrr_woe, f_excitrr_woe, xi_common, 'linear', 0);
    all_f_excitrr_woe = [all_f_excitrr_woe; f_interp_excitrr_woe];
end

% % -----------------------------
% % Compute KDE for ExcitRR_WO_ LGroup
% % -----------------------------
% all_f_excitrr_wol = [];
% excitam_excitrr_wol = fieldnames(ExcitRR_WO_L);
% num_excitrr_wol = length(excitam_excitrr_wol);
% all_data_excitrr_wol = [];
% 
% for m = 1:num_excitrr_wol
%     allbins_excitrr_wol = ExcitRR_WO_L.(excitam_excitrr_wol{m}).AllBins_prethresh;
%     all_data_excitrr_wol = [all_data_excitrr_wol; allbins_excitrr_wol];  % Collect for KS test
%     [f_excitrr_wol, xi_excitrr_wol] = ksdensity(allbins_excitrr_wol, 'Bandwidth', 0.008);
%     f_interp_excitrr_wol = interp1(xi_excitrr_wol, f_excitrr_wol, xi_common, 'linear', 0);
%     all_f_excitrr_wol = [all_f_excitrr_wol; f_interp_excitrr_wol];
% end

% -----------------------------
% Compute KDE for ExcitRR_WO_ LGroup
% -----------------------------
all_f_excitrr_wol = [];
excitam_excitrr_wol = fieldnames(ExcitRR_WO_E);
num_excitrr_wol = length(excitam_excitrr_wol);
all_data_excitrr_wol = [];

for m = 1:num_excitrr_wol
    allbins_excitrr_wol = ExcitRR_WO_E.(excitam_excitrr_wol{m}).AllBins_prethresh;
    all_data_excitrr_wol = [all_data_excitrr_wol; allbins_excitrr_wol];  % Collect for KS test
    %[f_excitrr_wol, xi_excitrr_wol] = ksdensity(allbins_excitrr_wol, 'Bandwidth', 0.008);
    [f_excitrr_wol, xi_excitrr_wol] = ksdensity(allbins_excitrr_wol, 'Bandwidth', 0.004);
    f_interp_excitrr_wol = interp1(xi_excitrr_wol, f_excitrr_wol, xi_common, 'linear', 0);
    all_f_excitrr_wol = [all_f_excitrr_wol; f_interp_excitrr_wol];
end




% -----------------------------
% Compute KDE for ExcitRR_Sham Group (Dashed Line Only, No Fill)
% -----------------------------
all_f_sham = [];
excitam_sham = fieldnames(ExcitRR_Sham);
num_sham = length(excitam_sham);
all_data_sham = [];

for m = 1:num_sham
    allbins_sham = ExcitRR_Sham.(excitam_sham{m}).AllBins_prethresh;
    all_data_sham = [all_data_sham; allbins_sham];  % Collect for KS test
    %[f_sham, xi_sham] = ksdensity(allbins_sham, 'Bandwidth', 0.008);
    [f_sham, xi_sham] = ksdensity(allbins_sham, 'Bandwidth', 0.004);
    f_interp_sham = interp1(xi_sham, f_sham, xi_common, 'linear', 0);
    all_f_sham = [all_f_sham; f_interp_sham];
end

%% Run KS Test for comparisons
[h_excitrr_ctrl, p_excitrr_ctrl, ks_stat_excitrr_ctrl] = kstest2(all_data_excitrr, all_data_excitrr_ctrl);
[h_excitrr_sham, p_excitrr_sham, ks_stat_excitrr_sham] = kstest2(all_data_excitrr, all_data_sham);
[h_excitrr_ctrl_sham, p_excitrr_ctrl_sham, ks_stat_excitrr_ctrl_sham] = kstest2(all_data_excitrr_ctrl, all_data_sham);

[h_excitrr_ctrl_early, p_excitrr_ctrl_early, ks_stat_excitrr_ctrl_early] = kstest2(all_data_excitrr_ctrl, all_data_excitrr_woe);
[h_excitrr_ctrl_late, p_excitrr_ctrl_late, ks_stat_excitrr_ctrl_late] = kstest2(all_data_excitrr_ctrl, all_data_excitrr_wol);
[h_excitrr_wash, p_excitrr_wash, ks_stat_excitrr_wash] = kstest2(all_data_excitrr_woe, all_data_excitrr_wol);

% Display KS Test results in command window
fprintf('KS test between ExcitRR and ExcitRR_Ctrl:\n');
fprintf('  h = %d (1 = reject null hypothesis)\n', h_excitrr_ctrl);
fprintf('  p = %.4f\n', p_excitrr_ctrl);
fprintf('  KS statistic = %.4f\n\n', ks_stat_excitrr_ctrl);

fprintf('KS test between ExcitRR and ExcitRR_Sham:\n');
fprintf('  h = %d (1 = reject null hypothesis)\n', h_excitrr_sham);
fprintf('  p = %.4f\n', p_excitrr_sham);
fprintf('  KS statistic = %.4f\n\n', ks_stat_excitrr_sham);

fprintf('KS test between ExcitRR_Ctrl and ExcitRR_Sham:\n');
fprintf('  h = %d (1 = reject null hypothesis)\n', h_excitrr_ctrl_sham);
fprintf('  p = %.4f\n', p_excitrr_ctrl_sham);
fprintf('  KS statistic = %.4f\n\n', ks_stat_excitrr_ctrl_sham);

fprintf('KS test between ExcitRR_Ctrl and ExcitRR_woe\n');
fprintf('  h = %d (1 = reject null hypothesis)\n', h_excitrr_ctrl_early);
fprintf('  p = %.4f\n', p_excitrr_ctrl_early);
fprintf('  KS statistic = %.4f\n\n', ks_stat_excitrr_ctrl_early);

fprintf('KS test between ExcitRR_Ctrl and ExcitRR_wol\n');
fprintf('  h = %d (1 = reject null hypothesis)\n', h_excitrr_ctrl_late);
fprintf('  p = %.4f\n', p_excitrr_ctrl_late);
fprintf('  KS statistic = %.4f\n\n', ks_stat_excitrr_ctrl_late);

fprintf('KS test between ExcitRR_wol and ExcitRR_woe\n');
fprintf('  h = %d (1 = reject null hypothesis)\n', h_excitrr_wash);
fprintf('  p = %.4f\n', p_excitrr_wash);
fprintf('  KS statistic = %.4f\n\n', ks_stat_excitrr_wash);


%% Plot All Groups Together
mean_positive_excitrr = mean(all_data_excitrr(all_data_excitrr > 0));
%figure; hold on;
figure('Position', [100, 100, 600, 300]);  hold on

% Fill areas for each group
fill([xi_common, fliplr(xi_common)], [mean(all_f_excitrr, 1), fliplr(zeros(1, num_points))], ...
    stimcolor, 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % ExcitRR (Teal)

fill([xi_common, fliplr(xi_common)], [mean(all_f_excitrr_ctrl, 1), fliplr(zeros(1, num_points))], ...
    ctrlcolor, 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % ExcitRR_Ctrl (LightSlateGray)

fill([xi_common, fliplr(xi_common)], [mean(all_f_excitrr_woe, 1), fliplr(zeros(1, num_points))], ...
    washecolor, 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % ExcitRR_Ctrl (LightSlateGray)

fill([xi_common, fliplr(xi_common)], [mean(all_f_excitrr_wol, 1), fliplr(zeros(1, num_points))], ...
    washlcolor, 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % ExcitRR_Ctrl (LightSlateGray)

% Plot mean KDE lines
%plot(xi_common, mean(all_f_excitrr, 1), 'Color', stimcolor, 'LineWidth', 3); % ExcitRR (Teal)
plot(xi_common, mean(all_f_excitrr_ctrl, 1), 'Color', ctrlcolor, 'LineWidth', 3); % ExcitRR_Ctrl (LightSlateGray)
plot(xi_common, mean(all_f_excitrr_woe, 1), 'Color', washecolor, 'LineWidth', 3); % ExcitRR_Ctrl (LightSlateGray)
plot(xi_common, mean(all_f_excitrr_wol, 1), 'Color', washlcolor, 'LineWidth', 3); % ExcitRR_Ctrl (LightSlateGray)


% Plot Sham group with dashed line, no fill
plot(xi_common, mean(all_f_sham, 1), '--', 'Color', shamcolor, 'LineWidth', 3); % Sham (Blue, Dashed)

% Set axis limits
xlim([-0.1 0.2]);
xline(mean_positive_excitrr,'--')

% Customize x-ticks
xticks([-0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2]);

% Labels and title
xlabel('Time from threshold crossing (s)');
ylabel('Density');
title('Average KDE across Animals (ExcitRR, ExcitRR_Ctrl, ExcitRR_Sham)');

% Annotate KS test results on plot
ymax = max([mean(all_f_excitrr, 1), mean(all_f_excitrr_ctrl, 1), mean(all_f_sham, 1)], [], 'all');
ks_text_1 = sprintf('KS (ExcitRR_Ctrlend vs ExcitRR_Ctrl1): p = %.4f, D = %.3f', p_excitrr_ctrl, ks_stat_excitrr_ctrl);
ks_text_2 = sprintf('KS (ExcitRR_Ctrlend vs ExcitRR_Sham): p = %.4f, D = %.3f', p_excitrr_sham, ks_stat_excitrr_sham);
ks_text_3 = sprintf('KS (ExcitRR_Ctrl1 vs Sham): p = %.4f, D = %.3f', p_excitrr_ctrl_sham, ks_stat_excitrr_ctrl_sham);

% Display the KS test text on the plot
text(0.05, ymax * 0.9, ks_text_1, 'FontSize', 12, 'FontWeight', 'bold');
text(0.05, ymax * 0.85, ks_text_2, 'FontSize', 12, 'FontWeight', 'bold');
text(0.05, ymax * 0.8, ks_text_3, 'FontSize', 12, 'FontWeight', 'bold');

% Formatting
set(gca, 'FontSize', 12, 'LineWidth', 1);
%ylim([0 10])

hold off;
