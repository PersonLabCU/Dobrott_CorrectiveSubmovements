%GetReaches.  Extracts reaches from trial kinematic data structure.  no
%overlap. Doesn't count reaches less than 200 ms (can modify in line 68).
%Doesn't count reaches shorter than .16 in x direction (line 73).
function [Expt, Reaches] = GetReachesHF (Cut,thresh)
  


%creates logical array for x values above/equal or below threshold.
above_thresh = Cut(:,2) >= thresh;

%find diff values; values of 1 are upward crossing, -1 downward
% concat a '0' to align to original index
thresh_diffs = [0;diff(above_thresh)];

%get upward crossings (where x(t+1)>x(t))
cross_up = find(thresh_diffs==1);

%get downward crossings (where x(t+1)<x(t))
cross_down = find(thresh_diffs==-1);

if (isempty(cross_up) || isempty(cross_down))
    Expt.reachinds = {};
    Reaches = struct;
else
    
if length(cross_up) < length(cross_down)
    crosses = [cross_up cross_down(2:end)];
else
    if cross_up(1) > cross_down(1)
        crosses = [cross_up(1:(length(cross_up)-1)) cross_down(2:end)];
    else
        crosses = [cross_up(1:length(cross_down)) cross_down];
    end
end

Expt.reachinds = crosses;

%get ready to loop through and cut out each reach.
%set up reach naming convention
reach_names = 'R%d';
reach_count = 1;

    for reachnum = 1:length(crosses(:,1));

        % start and stop at any direction reversal
       x = Cut(:,2);
       pre_x = x(1:crosses(reachnum,1));
       diffprex = diff(pre_x);
       start_inx  = find(diffprex<0,1,'last') + 1;
       if isempty(start_inx)
           start_inx = 1;
       end
       
       post_x = x(crosses(reachnum,2):end);
       diffpostx = diff(post_x);
       end_inx = find(diffpostx>0,1) + crosses(reachnum,2) - 1;
       if isempty(end_inx)
           start_inx = length(x);
       end

       %pull out kin data at crossed indices
       reach = Cut(start_inx:end_inx,:);
       
       %don't count reach if empty
       if isempty(reach)
           continue
       end
%        %don't count reach if less than 200 msec
%        if reach(end,1) - reach(1,1) <= .2
%            continue
%        end
       %don't count reach if doesn't get beyond .1
       if all (reach(:,2)<.1)
           continue
       end
       
       %establish reach name
       reach_name = sprintf(reach_names,reach_count);
       reach_count = reach_count+1;

       %Put reach data into structure
       Reaches.(reach_name).kin = reach;
       Reaches.(reach_name).inx = crosses(reachnum,:);

    end
    
    if ~exist('Reaches','var')
        Reaches  = struct;
    end
    
end










