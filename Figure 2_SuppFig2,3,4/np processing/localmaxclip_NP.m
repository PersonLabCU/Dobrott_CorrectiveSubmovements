
function [Reaches] = localmaxclip_NP(Reaches)


reach_ind = fieldnames(Reaches);

for reachnum = 1:length(reach_ind)
    
    reach = Reaches.(reach_ind{reachnum}).kin;
%     scatter(1:length(reach),reach(:,2)); hold on
    inx_start = find(reach(:,2)>=.1,1);
    
    if length(reach(inx_start:end,2)) < 3  % if too short
        Reaches.(reach_ind{reachnum}).kin = reach(1:inx_start,:);
    else
    
        [~,inx_peak] = findpeaks(reach(inx_start:end,2),'NPeaks',1);
    %     scatter(1:(inx_start+inx_peak),reach(1:(inx_start+inx_peak),2));
        Reaches.(reach_ind{reachnum}).kin = reach(1:(inx_start+inx_peak-1),:);
    %    close all
    end
    
end