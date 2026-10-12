%establish directory
fname = uigetdir;

direc = dir(strcat(fname,'\*.mat'));

%%
%get indices of files for desired animal. stored in anim_inx


a_names = {direc(:).name}';

get_name = 'SAM009';

anim_inx = [];

for i = 1:length(a_names)
    
    %extract full name
    full_name = a_names{i};
    
    %separate date and animal ID
    piece_name = strsplit(full_name,'_');
    
    %if no underscore, data file is not in proper format--> skip it
    if length(piece_name) < 2
        continue
    end
    
    %remove the '.mat' extension
    inx = find(contains(piece_name,'.mat'));
    mouse_name = regexprep(piece_name{inx},'.mat','');
    
    
    %if a recording animal, test for that here
    if contains(get_name, 'njc')
        
        %do the thing you want to do to keep track of IDs and indices
        if strcmp(get_name,mouse_name)
            anim_inx = [anim_inx i];
        end
        
        
    elseif contains(get_name, 'okd')
        
        %do the thing you want to do to keep track of IDs and indices
        if strcmp(get_name,mouse_name)
            anim_inx = [anim_inx i];
        end
        
    else
%         for everything other than recording animals ('njc' label)
        
        %remove the MIB label if it exists
        final_name = regexprep(mouse_name,'MIB','','ignorecase');
        
        
%         %activate this setion if just looking for MIBs. Otherwise,
%         %inactivate
%         %if any other letters, we don't want it
%         if any(isstrprop(final_name,'alpha'))
%             continue
%         end
%         
%         %add the MIB label back
%         final_name = ['MIB' final_name];
        
%         compare to desired name
        if strcmp(get_name,final_name)
            anim_inx = [anim_inx i];
        end
        
    end
    
    

end


%%
%get all data out from files, organized by days. will be stored in
%'anim_data'

%set up struct to store data
anim_data = struct;


for dnum = 1:length(anim_inx)
    
    %initialize naming convention first
    dnames = 'd%d';
    dname = sprintf(dnames,dnum);
    
    
    try
        load(fullfile(fname,direc(anim_inx(dnum)).name));
    catch
        warning('data prolly bad');
    end
    
    
    %if dta doesnt exist, data is formatted wrong for analysis
    if ~exist('dta','var')
        anim_data.(dname).Raw = [];
        anim_data.(dname).Reaches = struct;
        continue
    end
    if isempty(dta)
        anim_data.(dname).Raw = [];
        anim_data.(dname).Reaches = struct;
        continue
    end
    
    
    %analyze
    mark = MocapOrg(dta);
    ided = assign_identity3(mark);
    %Cut out noise (nonreach data per DJC).
    x = .3; y = .25; %these are the thresholds dylan uses (MIB 12/2021)
    nohead = ided.m1(~(ided.m1(:,2)>x | ided.m1(:,3)>y),:);
    Filt = MasterFilt2(nohead);
    [Expt,tempReaches] = GetReachesHF(Filt,.1);
    gapped = endgap(tempReaches,.05);

    %suc doesnt exist sometimes
    if exist('suc','var')
        Succed = sucassign(gapped,suc);
        Stimmed = StimAssign(Succed);
        Blocked = blockassign(Stimmed,blk_track);%add in block to allow for easy org 
        if exist('r_act_time','var')
            %[Key,NoKey] = Keypress(Stimmed,r_act_time);
            [Key,NoKey] = Keypress(Blocked,r_act_time);
        end
        
    else
        Succed = sucassign_nan(gapped);
        Stimmed = StimAssign(Succed);
       %Blocked = blockassign(Stimmed,blk_track);%add in block to allow for easy org 
        if exist('r_act_time','var')
            %[Key,NoKey] = Keypress(Stimmed,r_act_time);
            [Key,NoKey] = Keypress(Blocked,r_act_time);
        end
    end
    
    %blk_track doesnt exist sometimes
    if exist('blk_track','var')
         anim_data.(dname).blk_track = blk_track;
    end
    
    
    %organize data, raw and reaches
    anim_data.(dname).Raw = ided.m1;
    anim_data.(dname).Reaches = Stimmed;
    if exist('r_act_time','var')
        anim_data.(dname).Key = Key;
        anim_data.(dname).NoKey = NoKey;
    end

    %         master_data.(anames{anum}).(dname).thresh = .16;
    %         master_data.(anames{anum}).(dname).dir = 1;
    
    clearvars dta suc blk_track r_act_time
end