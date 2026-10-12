function [Events] = GetEvents_NP(Reaches)

%get event times (eg point of max speed)
rnames = fieldnames(Reaches);

for rnum = 1:length(rnames)
    kin = Reaches.(rnames{rnum}).kin;
    
    [~,maxs] = max(kin(:,6));
    Events.max_speed(rnum) = (kin(maxs,1));
    
    [~,maxvx] = max(kin(:,7)); % align to max outward velicity
    Events.max_vx(rnum) = (kin(maxvx,1));
    
    [~,maxvy] = max(kin(:,8));
    Events.max_vy(rnum) = (kin(maxvy,1));
    
    [~,maxvz] = max(kin(:,9));
    Events.max_vz(rnum) = (kin(maxvz,1));
    
    [~,endp] = max(kin(:,2));
    Events.end_point(rnum) = (kin(endp,1));
    
    Events.reach_start(rnum) = (kin(1,1));
    
    pthr = find(kin(:,2) > .1,1);
    Events.pos_thresh(rnum) = (kin(pthr,1));



    
end