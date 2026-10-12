%HF recording organization script

%%
%load NP data % For some reason won't run in Matlab 2020- run in Matlabb
%2018 
%establish directory: navigate to 'recording' subfolder
pathy = uigetdir;
pathx = string(pathy);
data=np_holder(pathx); % view 'np_holder_usage' for more info

data.get_event_times(.05);

%% plot event times 
% can get a general sense of when the optotaggin bbegan-2seconds bbetween
% laser pulses
diffEvent=diff(data.event_times);
figure()
plot(diffEvent,'o','MarkerSize', 8)

figure()
histogram(diffEvent)
 figure()
plot(data.event_times,'o','MarkerSize',0.5)

figure()
scatter(ided.m1(:,1),ided.m1(:,5))
%%
%load behavior data
[fz,pz] = uigetfile;
load(fullfile(pz,fz));

%process behav file into reaches
[Reaches,ided,Kins] = HFKinsWrap(dta,suc,r_act_time,blk_track); % adding blk_track so that blocks can bbe added to reaches 

%%
%unify timebase
NP_events = data.event_times;
time_dif = TimeOffset_NP(NP_events,ided);

%correct reach timebase
ided.m1(:,1) = ided.m1(:,1) + time_dif;
Kins(:,1) = Kins(:,1) + time_dif;
rn = fieldnames(Reaches);
for r = 1:length(rn)
    Reaches.(rn{r}).kin(:,1) = Reaches.(rn{r}).kin(:,1) + time_dif;
end

%%
%process reaches and identify kinematic event times for alignment

%process reaches
Local = localmaxclip_NP(Reaches);

%get event times
Events = GetEvents_NP(Local);

%%
%make rasters
% define baseline duration and event duration (time pre and post 'event')
% units are msec
% this is used for creating referenced spike times
baseline_epoch = 1000;
event_epoch = 1000;

%choose event to align to. relies on assigning exactly correct event_name
%above
event_name = ["end_point","pos_thresh"]; 
enames = fieldnames(Events);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000;
for enum = 1:length(evnt)
    event_times = Events.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    
    data.units(c).(enames{evnt(enum)}).Referenced = Ref;
    data.units(c).(enames{evnt(enum)}).EventTrial = EvTr;
    data.units(c).(enames{evnt(enum)}).Raster = BinMast;
    
    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).(enames{evnt(enum)}).PSTH_raw = psth_raw;
    data.units(c).(enames{evnt(enum)}).PSTH_smooth = psth_smooth;
    
    
end
end

%%

%plotting

epo = 500;

evnt = "end_point";
%evnt="pos_thresh";

for c = 1:22
    
    
    figure; hold on
    
    subplot(2,1,1); hold on
    h = title(['cluster' ' ' num2str(data.units(c).cluster_ID) ' ' 'depth' ' ' num2str(data.units(c).depth)]) ;
    set(h,'interpreter','none')
    
    % Plots Raster from binary matrix
    tick_length = 1; % length of tick marks for spikes
    spike_tick_width = 0.5; % width of lines marking spikes
    event_tick_width = 2; % width of lines marking secondary events
    
    
    rastr = data.units(c).(evnt).Raster(:,(1000-epo):(2000-(1000-epo)));
    [y_co,x_co] = find(rastr);
    x_co = x_co - epo;
    line([x_co';x_co'],[y_co'-(tick_length/2);y_co'+(tick_length/2)],'Color', [0 0 0], 'LineWidth', spike_tick_width);
    line([zeros(1,size(rastr,1));zeros(1,size(rastr,1))],[(1:1:size(rastr,1))-.5;(1:1:size(rastr,1))+.5],'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);
    axis ij;
    ylim([0,size(rastr,1)])
    xlim([-epo,epo])
    
    
    
    subplot(2,1,2); hold on
    title('peth')
    
    psth_smooth = data.units(c).(evnt).PSTH_smooth;
    psth_raw = data.units(c).(evnt).PSTH_raw;
    td1 = length(psth_raw)-length(psth_smooth)+1;
    td2 = td1/2;
    plot((-epo:epo),psth_smooth(((1000-epo)-td2):((2000-(1000-epo))-td2)),'k','LineWidth',3)
    xlim([-epo,epo])
%     ylim([0,100])

if mod(c,10) == 0
    drawnow
    input("press enter when done looking at figures")
    close all
end

end


%%
%heatmap

%choose event
%evnt = "end_point";
evnt="pos_thresh";

%optional indexing step to order units by depth or maxfr

% % % %depth
% dep = [data(:).units.depth]';
% [dep_sort,inx] = sort(dep);
% inx = flip(inx);

% %time of maxfr
maxfr_inx = NaN(length(data.units),1);
for c = 1:length(data.units)
    psth = data.units(c).(evnt).PSTH_smooth;
    
    [~,maxfr_inx(c)] = max(psth);
end
[maxfr_t,inx] = sort(maxfr_inx);
% 
% %order in array
inx = 1:length(data.units);


%organizing into master matrix
psth_master = NaN(length(data.units),length(data.units(inx(1)).(evnt).PSTH_smooth));
for c = 1:length(data.units)
    
    psth_master(c,:) = data.units(inx(c)).(evnt).PSTH_smooth;
    
    %normalize to max fr version
    psth_max = max(data.units(inx(c)).(evnt).PSTH_smooth);
    psth_master(c,:) = data.units(inx(c)).(evnt).PSTH_smooth/psth_max;
end



td1 = 2000-length(psth_master)+1;
td2 = td1/2;

figure;
imagesc(-1000+td2,1,psth_master)
h = colorbar;
cmap = brewermap(128,'*RdPu');
colormap(cmap)
% caxis([20 150])
xlim([-1000 1000])  
%%
% optotagging analysis - get event times for opto 
NP_opto=data.event_times;
diffNP_opto=diff(NP_opto);
idx=find(diffNP_opto>=2.95 & diffNP_opto<=3.1);%find first stim given during 
%eventopto=NP_opto(87:end); %
eventopto=NP_opto(idx(1,:):end)';
Events.eventopto=eventopto; %add opto event times to structure 

%% Get rasters for opto
%think about changing this - shorten to 500
baseline_epoch = 500;
event_epoch = 500;

event_name = ["eventopto"]; 
enames = fieldnames(Events);
evnt= find(ismember(enames, event_name));

%choose optoevents to align to 

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000;% putting into seconds)
for enum = 1:length(evnt)
    event_times = Events.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    
    data.units(c).eventopto(enum).Referenced = Ref;
    data.units(c).eventopto(enum).EventTrial = EvTr;
    data.units(c).eventopto(enum).Raster = BinMast;
    
    [psth_raw,psth_smooth] = Psth1(BinMast,5);
    data.units(c).eventopto(enum).PSTH_raw = psth_raw;
    data.units(c).eventopto(enum).PSTH_smooth = psth_smooth;
    
    
end
end



% need totake the difference but refer back to NP_opto data to get the
% actual times 
    
%%

%plotting

epo = 250;

%evnt = "end_point";
evnt="eventopto";

for c = 38

    
    
    figure; hold on
    
    subplot(2,1,1); hold on
    h = title(['cluster' ' ' num2str(data.units(c).cluster_ID) ' ' 'depth' ' ' num2str(data.units(c).depth)]) ;
    set(h,'interpreter','none')
    
    % Plots Raster from binary matrix
    tick_length = 1; % length of tick marks for spikes
    spike_tick_width = 0.5; % width of lines marking spikes
    event_tick_width = 2; % width of lines marking secondary events
    
    
    rastr = data.units(c).(evnt).Raster(:,(500-epo):(1000-(500-epo)));
    [y_co,x_co] = find(rastr);
    x_co = x_co - epo;
    line([x_co';x_co'],[y_co'-(tick_length/2);y_co'+(tick_length/2)],'Color', [0 0 0], 'LineWidth', spike_tick_width);
    line([zeros(1,size(rastr,1));zeros(1,size(rastr,1))],[(1:1:size(rastr,1))-.5;(1:1:size(rastr,1))+.5],'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);
    axis ij;
    ylim([0,size(rastr,1)])
    xlim([-epo,epo])
    
    
    
    subplot(2,1,2); hold on
    title('peth')
    
    psth_smooth = data.units(c).(evnt).PSTH_smooth;
    psth_raw = data.units(c).(evnt).PSTH_raw;
    td1 = length(psth_raw)-length(psth_smooth)+1;
    td2 = td1/2;
    plot((-epo:epo),psth_smooth(((500-epo)-td2):((1000-(500-epo))-td2)),'Color',[.2 .4 .8],'LineWidth',3)
    xlim([-epo,epo])
%     ylim([0,100])

if mod(c,10) == 0
    drawnow
    input("press enter when done looking at figures")
    close all
end

end


%% heatmap of unit activity from optotagged cells 
%heatmap

%choose event
%evnt = "end_point";
evnt="eventopto";

%optional indexing step to order units by depth or maxfr

% % %depth
dep = [data(:).units.depth]';
[dep_sort,inx] = sort(dep);
inx = flip(inx);

% %time of maxfr
% maxfr_inx = NaN(length(data.units),1);
% for c = 1:length(data.units)
%     psth = data.units(c).(evnt).PSTH_smooth;
%     
%     [~,maxfr_inx(c)] = max(psth);
% end
% [maxfr_t,inx] = sort(maxfr_inx);
% 
% %order in array
% inx = 1:length(data.units);


%organizing into master matrix
psth_master = NaN(length(data.units),length(data.units(inx(1)).(evnt).PSTH_smooth));
for c = 1:length(data.units)
    
    psth_master(c,:) = data.units(inx(c)).(evnt).PSTH_smooth;
    
    %normalize to max fr version
    psth_max = max(data.units(inx(c)).(evnt).PSTH_smooth);
    psth_master(c,:) = data.units(inx(c)).(evnt).PSTH_smooth/psth_max;
end



td1 = 2000-length(psth_master)+1;
td2 = td1/2;

figure;
imagesc(-1000+td2,1,psth_master)
h = colorbar;
cmap = brewermap(128,'*RdBu');
colormap(cmap)
% caxis([20 150])
xlim([-500 500]) 
xline(0)
%%
%adaptation ephys analysis
% Extract blocks and find the index so that it can be used 
rn = fieldnames(Reaches);
blk_array = [];
for r = 1:length(rn)
    blk_array(r) = Reaches.(rn{r}).block;
end
blk_array=blk_array';
%assume at least three blocks (baseline, adaptation, washout). separate
%blocks out 
bl_inx = find(blk_array == 1);
BASE = Extract(Reaches,bl_inx);
ad_inx = find(blk_array == 2);
ADAPT = Extract(Reaches,ad_inx);
wo_inx = find(blk_array == 3);
WASH = Extract(Reaches,wo_inx);


event_name = ["end_point","pos_thresh"]; 
enames = fieldnames(Events);
evnt= find(ismember(enames, event_name));

for c = 1:length(data.units)
    spike_times = data.units(c).spike_times*1000;
for enum = 1:length(evnt)
    event_times = Events.(enames{evnt(enum)})*1000;
    window = [baseline_epoch;event_epoch];
    
    [Ref, EvTr, BinMast] = Raster1(spike_times, event_times, window);
    
    data.units(c).(enames{evnt(enum)}).Referenced = Ref;
    data.units(c).(enames{evnt(enum)}).EventTrial = EvTr;
    data.units(c).(enames{evnt(enum)}).Raster = BinMast;
    % seaprate out psth 
    %baseline 
    BinMastBL=BinMast(bl_inx,:);
    [psth_rawBL,psth_smoothBL] = Psth1(BinMastBL,5);
    data.units(c).(enames{evnt(enum)}).PSTH_rawBL = psth_rawBL;
    data.units(c).(enames{evnt(enum)}).PSTH_smoothBL = psth_smoothBL;
    %adapt
     BinMastAD=BinMast(ad_inx,:);
    [psth_rawAD,psth_smoothAD] = Psth1(BinMastAD,5);
    data.units(c).(enames{evnt(enum)}).PSTH_rawAD = psth_rawAD;
    data.units(c).(enames{evnt(enum)}).PSTH_smoothAD = psth_smoothAD;
    % washout 
    wo_indx=wo_inx(1:5); % can change this but keeping consistent with Dylan
    BinMastWO=BinMast(wo_indx,:); % consider shortening washout 
    [psth_rawWO,psth_smoothWO] = Psth1(BinMastWO,5);
    data.units(c).(enames{evnt(enum)}).PSTH_rawWO = psth_rawWO;
    data.units(c).(enames{evnt(enum)}).PSTH_smoothWO = psth_smoothWO;
    
    
end
end



%% Plotting 

epo = 500;

evnt = "end_point";

%evnt="pos_thresh";

for c = 1:5
    
    
    figure; hold on
    
    subplot(4,1,1); hold on
    h = title(['cluster' ' ' num2str(data.units(c).cluster_ID) ' ' 'depth' ' ' num2str(data.units(c).depth)]) ;
    set(h,'interpreter','none')
    
    % Plots Raster from binary matrix
    tick_length = 1; % length of tick marks for spikes
    spike_tick_width = 0.5; % width of lines marking spikes
    event_tick_width = 2; % width of lines marking secondary events
    
    raster=data.units(c).(evnt).Raster;
    Start=length(bl_inx)-5;
    End=length(bl_inx);
    rastr=raster(Start:End,:); % extract indexfor trials from baseline 
    rastr = rastr(:,(1000-epo):(2000-(1000-epo)));
    [y_co,x_co] = find(rastr);
    x_co = x_co - epo;
    line([x_co';x_co'],[y_co'-(tick_length/2);y_co'+(tick_length/2)],'Color', [0 0 0], 'LineWidth', spike_tick_width);
    line([zeros(1,size(rastr,1));zeros(1,size(rastr,1))],[(1:1:size(rastr,1))-.5;(1:1:size(rastr,1))+.5],'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);
    axis ij;
    ylim([0,size(rastr,1)])
    xlim([-epo,epo])
    
    subplot(4,1,2); hold on
    h = title(['cluster' ' ' num2str(data.units(c).cluster_ID) ' ' 'depth' ' ' num2str(data.units(c).depth)]) ;
    set(h,'interpreter','none')
    
    % Plots Raster from binary matrix
    tick_length = 1; % length of tick marks for spikes
    spike_tick_width = 0.5; % width of lines marking spikes
    event_tick_width = 2; % width of lines marking secondary events
    
    raster=data.units(c).(evnt).Raster;
    rastr=raster(ad_inx,:); % extract indexfor trials from baseline 
    rastr = rastr(:,(1000-epo):(2000-(1000-epo)));
    [y_co,x_co] = find(rastr);
    x_co = x_co - epo;
    line([x_co';x_co'],[y_co'-(tick_length/2);y_co'+(tick_length/2)],'Color', [0 0 1], 'LineWidth', spike_tick_width);
    line([zeros(1,size(rastr,1));zeros(1,size(rastr,1))],[(1:1:size(rastr,1))-.5;(1:1:size(rastr,1))+.5],'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);
    axis ij;
    ylim([0,size(rastr,1)])
    xlim([-epo,epo])
    
    subplot(4,1,3); hold on
    h = title(['cluster' ' ' num2str(data.units(c).cluster_ID) ' ' 'depth' ' ' num2str(data.units(c).depth)]) ;
    set(h,'interpreter','none')
    
    % Plots Raster from binary matrix
    tick_length = 1; % length of tick marks for spikes
    spike_tick_width = 0.5; % width of lines marking spikes
    event_tick_width = 2; % width of lines marking secondary events
    
    raster=data.units(c).(evnt).Raster;
    wo_indx=wo_inx(1:5);
    rastr=raster(wo_indx,:); % extract indexfor trials from baseline 
    rastr = rastr(:,(1000-epo):(2000-(1000-epo)));
    [y_co,x_co] = find(rastr);
    x_co = x_co - epo;
    line([x_co';x_co'],[y_co'-(tick_length/2);y_co'+(tick_length/2)],'Color', [1 0 1], 'LineWidth', spike_tick_width);
    line([zeros(1,size(rastr,1));zeros(1,size(rastr,1))],[(1:1:size(rastr,1))-.5;(1:1:size(rastr,1))+.5],'Color', [0.7 0.7 0.7], 'LineWidth', event_tick_width);
    axis ij;
    ylim([0,size(rastr,1)])
    xlim([-epo,epo])
    
    % plot difference between baseline and washout ( washout - bbaseline ) 
%     subplot(4,1,4); hold on
%     title('peth (Washout-Baseline)')
%     
%     psth_smoothBL = data.units(c).(evnt).PSTH_smoothBL;
%     psth_raw = data.units(c).(evnt).PSTH_rawBL;
%     psth_smoothWO=data.units(c).(evnt).PSTH_smoothWO;
%     psth_smooth=psth_smoothWO-psth_smoothBL;
%     td1 = length(psth_raw)-length(psth_smoothBL)+1;
%     td2 = td1/2;
%     plot((-epo:epo),psth_smooth(((1000-epo)-td2):((2000-(1000-epo))-td2)),'k','LineWidth',3)
%     xlim([-epo,epo])
%     ylim([0,100])

    % plot PSTH for baseline and washout on top of eachother when looking
    % at positional threshold 
    subplot(4,1,4); hold on
    title('peth (Washout and Baseline)')
    
    psth_smoothBL = data.units(c).(evnt).PSTH_smoothBL;
    psth_raw = data.units(c).(evnt).PSTH_rawBL;
    psth_smoothWO=data.units(c).(evnt).PSTH_smoothWO;
    td1 = length(psth_raw)-length(psth_smoothBL)+1;
    td2 = td1/2;
    plot((-epo:epo),psth_smoothBL(((1000-epo)-td2):((2000-(1000-epo))-td2)),'k','LineWidth',1.5)
    hold on 
    plot((-epo:epo),psth_smoothWO(((1000-epo)-td2):((2000-(1000-epo))-td2)),'m','LineWidth',1.5)
    xlim([-epo,epo])
%     ylim([0,100])

if mod(c,10) == 0
    drawnow
    input("press enter when done looking at figures")
    close all
end

end

%% %%
%heatmap

%choose event
%evnt = "end_point";
evnt="pos_thresh";

%optional indexing step to order units by depth or maxfr

% % %depth
% dep = [data(:).units.depth]';
% [dep_sort,inx] = sort(dep);
% inx = flip(inx);

% %time of maxfr
maxfr_inx = NaN(length(data.units),1);
for c = 1:length(data.units)
    psth = data.units(c).(evnt).PSTH_smoothBL;
    
    [~,maxfr_inx(c)] = max(psth);
end
[maxfr_t,inx] = sort(maxfr_inx);
% 
% %order in array
inx = 1:length(data.units);


%organizing into master matrix
psth_masterBL = NaN(length(data.units),length(data.units(inx(1)).(evnt).PSTH_smoothBL));
for c = 1:length(data.units)
    
    psth_masterBL(c,:) = data.units(inx(c)).(evnt).PSTH_smoothBL;
    
    %normalize to max fr version
    psth_max = max(data.units(inx(c)).(evnt).PSTH_smoothBL);
    psth_masterBL(c,:) = data.units(inx(c)).(evnt).PSTH_smoothBL/psth_max;
end

psth_masterWO = NaN(length(data.units),length(data.units(inx(1)).(evnt).PSTH_smoothWO));
for c = 1:length(data.units)
    
    psth_masterWO(c,:) = data.units(inx(c)).(evnt).PSTH_smoothWO;
    
    %normalize to max fr version
    psth_max = max(data.units(inx(c)).(evnt).PSTH_smoothWO);
    psth_masterWO(c,:) = data.units(inx(c)).(evnt).PSTH_smoothWO/psth_max;
end



td1 = 2000-length(psth_masterBL)+1;
td2 = td1/2;

figure;
subplot(2,1,1);hold on
imagesc(-1000+td2,1,psth_masterBL)
h = colorbar;
cmap = brewermap(128,'*RdPu');
colormap(cmap)
% caxis([20 150])
xlim([-1000 1000])  
title('Baseline')

subplot(2,1,2); hold on
imagesc(-1000+td2,1,psth_masterWO)
h = colorbar;
cmap = brewermap(128,'*RdPu');
colormap(cmap)
% caxis([20 150])
xlim([-1000 1000])  
title('Washout')

%%
% find the difference in firing rates 
psth_masterdiff= psth_masterWO-psth_masterBL;
figure;
imagesc(-1000+td2,1,psth_masterdiff)
h = colorbar;
cmap = brewermap(128,'*RdBu');
colormap(cmap)
% caxis([20 150])
%xlim([-1000 1000])
xlim([-250 250])
title('Washout-Baseline(sorted by depth)')



%%
%get all reach names from original reach struct
%get subset of reach names from R_Adapt (blockadaptation_HF)
%get indices of R_Adapt reach names in original struct
%extract subset of rows from raster
%utilize stim_inx to group reaches into categories (BL, ST_E, ST_L, WA)
%Quality control metrics ? Should I be implementing some of these ? 
% need to make paramters for figuring out what cells are PC or NC also
% quality control metrics 

