Reaches=anim_data.d2.Key;
R_Adapt = Reaches;

% %remove reaches that mess up the block
% R_Adapt = Extract(R_Adapt,[1:12 14:18 21 22 24:45]);

rn = fieldnames(R_Adapt);
stim_inx = [];
for r = 1:length(rn)
    stim_inx(r,1) = R_Adapt.(rn{r}).stim;
end



stim_block_start  = find(stim_inx == 1,1);
wash_block_start = find(stim_inx(stim_block_start:end) == 0,1)+stim_block_start-1;

if any(stim_inx(wash_block_start:end) == 1)
    second_stim_start = find(stim_inx(wash_block_start:end) == 1,1)+wash_block_start-1;
    R_Adapt = Extract(R_Adapt,1:second_stim_start-1);
end

BL = Extract(R_Adapt,1:(stim_block_start-1));
ST_E = Extract(R_Adapt,stim_block_start:(stim_block_start+4));
ST_L = Extract(R_Adapt,(wash_block_start-5):(wash_block_start-1));
% WA = Extract(R_Adapt,wash_block_start:length(fieldnames(R_Adapt)));
WA = Extract(R_Adapt,wash_block_start:(wash_block_start+4));


%%
thresh = .1;
BLN = AlignThresh(BL,thresh,1,200);
STN_E = AlignThresh(ST_E,thresh,1,200);
STN_L = AlignThresh(ST_L,thresh,1,200);
WAN = AlignThresh(WA,thresh,1,200);

%plotting
dim = 2;
%figure with basline (all or late), adapt early, adapt late, and washout
%early
er_ln = 1;
figure; hold on
StimvelPlot2(BLN,dim,[.5 .5 .5],0);
StimvelPlot2(STN_E,dim,[0 1 1],0);
StimvelPlot2(STN_L,dim,[0 0 1],0);
StimvelPlot2(WAN,dim,[1 0 1],0);
xlim([-.2 .2])
ylim auto



