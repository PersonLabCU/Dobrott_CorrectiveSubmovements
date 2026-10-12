%%behav analysis per animal


%loop through sessions, extract blocks
anim = anim_data;
day=fieldnames(anim); 

for d = 1:length(anim)

    %check to see if keypress- unsure how to get keypress logical here 
%     if anim(d).keypress_logical == 0
%         continue
%     end

    %Reaches = anim(d).Reaches;
    Reaches = anim.(day{d}).Key;
   

    %check to see whether any stim blocks
    rn = fieldnames(Reaches);
    stim_inx = [];
    for r = 1:length(rn)
        stim_inx(r,1) = Reaches.(rn{r}).stim;
    end
    if sum(stim_inx) < 5
        continue
    end


    %pull out blocks
    [R_Adapt,R_Adapt_inx,BL,ST,WA] = block_adapt_org_HF(Reaches);

    anim(d).R_Adapt = R_Adapt;
    anim(d).R_Adapt_inx = R_Adapt_inx;
    anim(d).BL = BL;
    anim(d).ST = ST;
    anim(d).WA = WA;
end

%% Behav analysis for Random Stim 
%loop through sessions, extract blocks
anim = anim_data;
day=fieldnames(anim);


for d = 1:length(day)

    Reaches = anim.(day{d}).Key;
    %Reaches = anim.(day{d}).NoKey;

    [Stim, Control]=StimFilter(Reaches);
   
    anim.(day{d}).Stim = Stim;
    anim.(day{d}).Control=Control;
end

%% PLot random stim 

%analyze variables of interest and plot

%find sessions with significant stim effect

thresh = .1;
direction = 1;
dn=fieldnames(anim);

stim_effect_inx = false(length(anim),1);
folder_path="C:\Users\Person Lab\OneDrive - The University of Colorado Denver\Documents\Summary Data\HF_BlockStim_Matt Data\individual days\SAM010_random";
for d = 1:length(dn)


    C = AlignThresh(anim.(dn{d}).Control,thresh,direction,200);
    S = AlignThresh(Extract(anim.(dn{d}).Stim),thresh,direction,200);
    
    x1=0 ;x2=0.05;
    y_upper=2; y_lower=-1;

    dim = 5;
    StimvelPlot(C,dim,[.5 .5 .5],1)
    StimvelPlot(S,dim,[.2 .6 .9],0)
    hold on 
    % Fill the area between the vertical lines
    fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
    yline(0,'--k')    
    title(num2str(d))
   % ylim([0 0.3])
    xlim([-0.2 0.2])

    filename=sprintf('SAM010_outward_vel_random%d.svg',d);
    %filename=sprintf('vmax_raster%d.svg',c);
    full_path=fullfile(folder_path,filename);
    saveas(gcf,full_path);

    close;
    Ctest = C;
    Stest = S;

    parray = NaN(size(Ctest.stimvel,1),2);

    for pnum = 1: length(parray)
        a4b = squeeze(Ctest.stimvel(pnum,5,:));
        b4b = squeeze(Stest.stimvel(pnum,5,:));
        
        if all(isnan(a4b)) || all(isnan(b4b))
            continue
        end

        [parray(pnum,1),b,pstat] = ranksum(a4b,b4b,'method','approximate');
        parray(pnum,2) = pstat.zval;
    end

    if any(parray(23:31,1) < .05)
        stim_effect_inx(d) = true;
    end

end

%manually change stim_effect_inx if graphs looks weird


%%
%analyze variables of interest and plot

%find sessions with significant stim effect

thresh = .1;
direction = 1;

stim_effect_inx = false(length(anim),1);

for d = 1:length(anim)

    if isempty(anim(d).BL)
        continue
    elseif length(fieldnames(anim(d).BL)) < 5
        continue
    end


    C = AlignThresh(anim(d).BL,thresh,direction,200);
    S = AlignThresh(Extract(anim(d).ST,1:10),thresh,direction,200);

    dim = 5;
    StimvelPlot(C,dim,[.5 .5 .5],1)
    StimvelPlot(S,dim,[.2 .6 .9],0)
    title(num2str(d))
    ylim([-1 2.5])

    Ctest = C;
    Stest = S;

    parray = NaN(size(Ctest.stimvel,1),2);

    for pnum = 1: length(parray)
        a4b = squeeze(Ctest.stimvel(pnum,5,:));
        b4b = squeeze(Stest.stimvel(pnum,5,:));
        
        if all(isnan(a4b)) || all(isnan(b4b))
            continue
        end

        [parray(pnum,1),b,pstat] = ranksum(a4b,b4b,'method','approximate');
        parray(pnum,2) = pstat.zval;
    end

    if any(parray(23:31,1) < .05)
        stim_effect_inx(d) = true;
    end

end

%manually change stim_effect_inx if graphs looks weird


%%

%collate data per stim_effect_inx
good_days = find(stim_effect_inx == 1);

BASE = struct;
ADAP_E = struct;
ADAP_L = struct;
WASH = struct;

chunk_size = 1;
chnk = chunk_size - 1;

for d = 8:11

    BL_len = length(fieldnames(anim(good_days(d)).BL));
    BL_temp = Extract(anim(good_days(d)).BL,BL_len-chnk:BL_len);
    BASE = Combine(BASE,BL_temp);

    ST_len = length(fieldnames(anim(good_days(d)).ST));
    ADAPE_temp = Extract(anim(good_days(d)).ST,1:chnk+1);
    ADAP_E = Combine(ADAP_E,ADAPE_temp);

    ADAPL_temp = Extract(anim(good_days(d)).ST,ST_len-chnk:ST_len);
    ADAP_L = Combine(ADAP_L,ADAPL_temp);

    WA_len = length(fieldnames(anim(good_days(d)).WA));
    WA_temp = Extract(anim(good_days(d)).WA,1:chnk+1);
    WASH = Combine(WASH,WA_temp);

end

thresh = 0.1;
BLN = AlignThresh(BASE,thresh,1,200);
STN_E = AlignThresh(ADAP_E,thresh,1,200);
STN_L = AlignThresh(ADAP_L,thresh,1,200);
WAN = AlignThresh(WASH,thresh,1,200);

%plotting
dim = 5;
%figure with basline (all or late), adapt early, adapt late, and washout
%early
er_ln = 1;
figure; hold on
StimvelPlot2(BLN,dim,[.7 .5 .5],0,er_ln);
StimvelPlot2(STN_E,dim,[.2 .4 .8],0,er_ln);
StimvelPlot2(STN_L,dim,[.1 .3 .5],0,er_ln);
StimvelPlot2(WAN,dim,[.5 .7 .5],0,er_ln);
xlim([-.2 .2])
ylim auto

anim(1).BLN = BLN;
anim(1).STN_E = STN_E;
anim(1).STN_L = STN_L;
anim(1).WAN = WAN;


%%

SAM010 = anim;
