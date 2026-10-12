function [ReReach,StimRR,MaxVel]=ReReachFilter_accel_ctrlreaches_edit(Reaches)
%Rereach is structure only including reaches that had a re-reach
%StimRR has the time of the rereach time and position 
    
%Input StimReaches 
    Sreaches=Reaches;
    reach_ind = fieldnames(Sreaches);
    ReReach=[];
    
    ReReachInfo={};
    StimInfo={};
    MaxVs={};
 
    thresh=0.16; % For headfixed reaches 

    for reachnum = 1:length(reach_ind)
        reach = Sreaches.(reach_ind{reachnum}).kin;
        stim_inx = find(reach(:,2)>=thresh,1); 
        accel_reach=reach(stim_inx:end,:);
        
        % Window of time for corrective movement
        windowstart=reach(stim_inx,1); 
        windowend=windowstart+0.20; 
        
        differences = abs(accel_reach(:,1) - windowend); 
        [minDifference, index] = min(differences);
        accel_reachclipped = accel_reach(1:index,:);
        
        accel_col = accel_reachclipped(:,11); % Outward accel column
        
        % Find indices where acceleration is positive
        pos_accel_idx = find(accel_col > 0);
        
        % Ensure acceleration is sustained for at least 25 ms
        sampling_interval = mean(diff(accel_reach(:,1))); % Get time step size
        min_samples = ceil(0.025 / sampling_interval); % Convert 25 ms to number of samples

        %min_samples = ceil(0.05 / sampling_interval); % Convert 50 ms to number of samples

        startidx = NaN; % Default value if no sustained positive acceleration is found
        for i = 1:length(pos_accel_idx) - min_samples + 1
            if all(accel_col(pos_accel_idx(i:i+min_samples-1)) > 0)
                startidx = pos_accel_idx(i); 
                break;
            end
        end

        if isnan(startidx)
            continue;
        else
            ReReach.(reach_ind{reachnum}).kin=reach;    
            ReReach.(reach_ind{reachnum}).inx=Sreaches.(reach_ind{reachnum}).inx;
            ReReach.(reach_ind{reachnum}).success=Sreaches.(reach_ind{reachnum}).success;
            ReReach.(reach_ind{reachnum}).stim=Sreaches.(reach_ind{reachnum}).stim;

            % Get time and position of rereach
            rrtime=accel_reach(startidx,1);
            rrpos=accel_reach(startidx,2);
            
            % Find time and position of stim
            stiminx = find(accel_reach(:,2)>=thresh,1); 
            stimtime=accel_reach(stiminx,1);
            stimpos=accel_reach(stiminx,2);

            % Store rereach and stimulus info
            RR={rrtime,rrpos};
            Stiminfo={stimtime,stimpos};
            ReReachInfo=[ReReachInfo;RR];
            StimInfo=[StimInfo;Stiminfo];

            % Find max velocity before and after rereach
            reachbefore=accel_reach(1:startidx,7); % Velocity before re-reach
            maxvbefore=max(reachbefore);
            reachafter=accel_reach(startidx:end,7); % Velocity after re-reach
            maxvafter=max(reachafter);
            Maxvs={maxvbefore,maxvafter};
            MaxVs=[MaxVs;Maxvs];
        end
    end 
    
    ReReachInfo=cell2mat(ReReachInfo);
    StimInfo=cell2mat(StimInfo);
    StimRR=ReReachInfo-StimInfo;
    MaxVel=cell2mat(MaxVs);
end
