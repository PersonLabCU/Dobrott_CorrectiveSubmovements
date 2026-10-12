
% Add alignments to data 
an=fieldnames(ExcitRR_BL);
EndPtMean = struct();

for a=1:length(an)

    %BLR=ExcitRR_BL.(an{a}).BL;
    % ADR=ExcitRR_AD.(an{a}).AD;
    % ADER=ExcitRR_AD_E.(an{a}).AD_E;
    % ADLR=ExcitRR_AD_L.(an{a}).AD_L;
    % WOR=ExcitRR_WO.(an{a}).WO;
    % WOER=ExcitRR_WO_E.(an{a}).WO_E;
     WOLR=ExcitRR_WO_L.(an{a}).WO;

     %BLR=ExcitRR_BL.(an{a}).ReReaches;
    ADR=ExcitRR_AD.(an{a}).ReReaches;
    ADER=ExcitRR_AD_E.(an{a}).ReReaches;
    ADLR=ExcitRR_AD_L.(an{a}).ReReaches;
    %  WOR=ExcitRR_WO.(an{a}).ReReaches;
    % WOER=ExcitRR_WO_E.(an{a}).ReReaches;
    %WOLR=ExcitRR_WO_L.(an{a}).ReReaches;

     BLR=ExcitRR_BL.(an{a}).RRT2; % use for analysis
    % % ADR=ExcitRR_AD.(an{a}).RRT2;
    % % ADER=ExcitRR_AD_E.(an{a}).RRT2;
    % % ADLR=ExcitRR_AD_L.(an{a}).RRT2;
      WOR=ExcitRR_WO.(an{a}).RRT2; % use for analysis 
      WOER=ExcitRR_WO_E.(an{a}).RRT2; % use for analysis 
    %WOLR=ExcitRR_WO_L.(an{a}).RRT2;
    
      % ---- endpoints + mean endpoints ----
    BL_end  = endpoint(BLR);   BL_mean  = mean(BL_end,  1, 'omitnan');
    AD_end  = endpoint(ADR);   AD_mean  = mean(AD_end,  1, 'omitnan');
    ADE_end = endpoint(ADER);  ADE_mean = mean(ADE_end, 1, 'omitnan');
    ADL_end = endpoint(ADLR);  ADL_mean = mean(ADL_end, 1, 'omitnan');
    WO_end  = endpoint(WOR);   WO_mean  = mean(WO_end,  1, 'omitnan');
    WOE_end = endpoint(WOER);  WOE_mean = mean(WOE_end, 1, 'omitnan');

    % ---- save into struct (by animal) ----
    EndPtMean.(an{a}).BL  = BL_mean;
    EndPtMean.(an{a}).AD  = AD_mean;
    EndPtMean.(an{a}).AD_E = ADE_mean;
    EndPtMean.(an{a}).AD_L = ADL_mean;
    EndPtMean.(an{a}).WO  = WO_mean;
    EndPtMean.(an{a}).WO_E = WOE_mean;

    % (optional but useful) save N reaches too
    EndPtMean.(an{a}).N.BL   = size(BL_end,1);
    EndPtMean.(an{a}).N.AD   = size(AD_end,1);
    EndPtMean.(an{a}).N.AD_E = size(ADE_end,1);
    EndPtMean.(an{a}).N.AD_L = size(ADL_end,1);
    EndPtMean.(an{a}).N.WO   = size(WO_end,1);
    EndPtMean.(an{a}).N.WO_E = size(WOE_end,1);
    
    thresh=0.1; dirdir=1;
    % % Align all reaches to CSM- average time of induced CSM 
    BL_CSM=AlignThreshCSMHF(BLR,thresh,dirdir,250);
    BL_CSM=stimvel_accel(BL_CSM);
    ExcitRR_BL.(an{a}).BL_CSM=BL_CSM;

    % % Stim all
    AD_CSM=AlignThreshCSMHF(ADR,thresh,dirdir,250);
    AD_CSM=stimvel_accel(AD_CSM);
    ExcitRR_AD.(an{a}).AD_CSM=AD_CSM;
    %Stim early
    ADE_CSM=AlignThreshCSMHF(ADER,thresh,dirdir,250);
    ADE_CSM=stimvel_accel(ADE_CSM);
    ExcitRR_AD_E.(an{a}).ADE_CSM=ADE_CSM;
    %Stim late
    ADL_CSM=AlignThreshCSMHF(ADLR,thresh,dirdir,250);
    ADL_CSM=stimvel_accel(ADL_CSM);
    ExcitRR_AD_L.(an{a}).ADL_CSM=ADL_CSM;

    % % Washout all
    WO_CSM=AlignThreshCSMHF(WOR,thresh,dirdir,250);
    WO_CSM=stimvel_accel(WO_CSM);
    ExcitRR_WO.(an{a}).WO_CSM=WO_CSM;
    %Washout early
    WOE_CSM=AlignThreshCSMHF(WOER,thresh,dirdir,250);
    WOE_CSM=stimvel_accel(WOE_CSM);
    ExcitRR_WO_E.(an{a}).WOE_CSM=WOE_CSM;
    %Washout early
    WOL_CSM=AlignThreshCSMHF(WOLR,thresh,dirdir,250);
    WOL_CSM=stimvel_accel(WOL_CSM);
    ExcitRR_WO_L.(an{a}).WOL_CSM=WOL_CSM;

    % % Align all reaches to Max V 
    dim=7; dir=1;
    BL_peakV= AlignPeakV(BLR,dim,250);
    BL_peakV=stimvel_accel(BL_peakV);
    ExcitRR_BL.(an{a}).BL_peakV=BL_peakV;

    % % Stimulated trials all
    AD_peakV= AlignPeakV(ADR,dim,250);
    AD_peakV=stimvel_accel(AD_peakV);
    ExcitRR_AD.(an{a}).AD_peakV=AD_peakV;
    % Adapt early
    ADE_peakV= AlignPeakV(ADER,dim,250);
    ADE_peakV=stimvel_accel(ADE_peakV);
    ExcitRR_AD_E.(an{a}).ADE_peakV=ADE_peakV;
     % Adapt late
    ADL_peakV= AlignPeakV(ADLR,dim,250);
    ADL_peakV=stimvel_accel(ADL_peakV);
    ExcitRR_AD_L.(an{a}).ADL_peakV=ADL_peakV;
    
    % % All washout trials 

    WO_peakV= AlignPeakV(WOR,dim,250);
    WO_peakV=stimvel_accel(WO_peakV);
    ExcitRR_WO.(an{a}).WO_peakV=WO_peakV;

    WOE_peakV= AlignPeakV(WOER,dim,250);
    WOE_peakV=stimvel_accel(WOE_peakV);
    ExcitRR_WO_E.(an{a}).WOE_peakV=WOE_peakV;

    WOL_peakV= AlignPeakV(WOLR,dim,250);
    WOL_peakV=stimvel_accel(WOL_peakV);
    ExcitRR_WO_L.(an{a}).WOL_peakV=WOL_peakV;

    % % % Align all trials to thresh

    BL_thresh=AlignThreshHF(BLR,thresh,dirdir,250);
    BL_thresh=stimvel_accel(BL_thresh);
    ExcitRR_BL.(an{a}).BLN=BL_thresh;

    % % Stim all
    AD_thresh=AlignStim(ADR,250);
    AD_thresh=stimvel_accel(AD_thresh);
    ExcitRR_AD.(an{a}).ADP=AD_thresh;
    %Stim early
    ADE_thresh=AlignStim(ADER,250);
    ADE_thresh=stimvel_accel(ADE_thresh);
    ExcitRR_AD_E.(an{a}).ADP_E=ADE_thresh;
    %Stim late
    ADL_thresh=AlignStim(ADLR,250);
    ADL_thresh=stimvel_accel(ADL_thresh);
    ExcitRR_AD_L.(an{a}).ADP_L=ADL_thresh;

    % % Washout all
    WO_thresh=AlignThreshHF(WOR,thresh,dirdir,250);
    WO_thresh=stimvel_accel(WO_thresh);
    ExcitRR_WO.(an{a}).WSH =WO_thresh;
    %Washout early
    WOE_thresh=AlignThreshHF(WOER,thresh,dirdir,250);
    WOE_thresh=stimvel_accel(WOE_thresh);
    ExcitRR_WO_E.(an{a}).WSH_E=WOE_thresh;
    %Washout early
    WOL_thresh=AlignThreshHF(WOLR,thresh,dirdir,250);
    WOL_thresh=stimvel_accel(WOL_thresh);
    ExcitRR_WO_L.(an{a}).WSH_L=WOL_thresh;


  
end
%% SHAM
am=fieldnames(ExcitRR_Sham);

for b=1:length(am)

    %BLR=ExcitRR_Sham.(am{b}).CtrlReaches;
    BLR=ExcitRR_Sham.(am{b}).ReReaches;
    BLRAlign=AlignThreshHF(BLR,thresh,dirdir,250);
    BLRAlign=stimvel_accel(BLRAlign);
    ExcitRR_Sham.(am{b}).BLRAlign=BLRAlign;

end
%% ===================== PLOT MEAN ENDPOINTS BY CATEGORY =====================
% Requires: EndPtMean struct from your previous loop
% Each EndPtMean.(animal).<cat> is a 1x3 mean endpoint [outward, y, z] (or [x,y,z])

an = fieldnames(EndPtMean);

cats = {'BL','AD','AD_E','AD_L','WO','WO_E'};  % edit if you changed names

% ---------- collect per-animal mean endpoints into matrices ----------
E = struct();
for ci = 1:numel(cats)
    cat = cats{ci};
    E.(cat) = nan(numel(an), 3);
end

for a = 1:numel(an)
    animal = an{a};
    for ci = 1:numel(cats)
        cat = cats{ci};
        if isfield(EndPtMean.(animal), cat) && ~isempty(EndPtMean.(animal).(cat))
            E.(cat)(a,:) = EndPtMean.(animal).(cat); % 1x3
        end
    end
end

% ---------- 1) 3D plot (grand mean across animals for each category) ----------
Mean3D = nan(numel(cats),3);
for ci = 1:numel(cats)
    Mean3D(ci,:) = mean(E.(cats{ci}), 1, 'omitnan'); % 1x3
end

figure('Color','w','Position',[200 200 700 550]); hold on;
for ci = 1:numel(cats)
    plot3(Mean3D(ci,1), Mean3D(ci,2), Mean3D(ci,3), 'o', 'MarkerSize', 10, 'LineWidth', 2);
    text(Mean3D(ci,1), Mean3D(ci,2), Mean3D(ci,3), ['  ' cats{ci}], 'FontSize', 11);
end
grid on; axis equal;
xlabel('Outward endpoint (col 1)');
ylabel('Endpoint dim 2');
zlabel('Endpoint dim 3');
title('Mean endpoint by reach category (across animals)');

% (optional) connect in a sequence if you want
% order = [1 2 3 4 5 6]; % BL->AD->AD_E->AD_L->WO->WO_E
% plot3(Mean3D(order,1), Mean3D(order,2), Mean3D(order,3), '-', 'LineWidth', 1.5);

% ---------- 2) Outward only (first column) ----------
MeanOut = nan(numel(cats),1);
SemOut  = nan(numel(cats),1);

for ci = 1:numel(cats)
    x = E.(cats{ci})(:,1);              % outward endpoint per animal
    MeanOut(ci) = mean(x,'omitnan');
    SemOut(ci)  = std(x,'omitnan') ./ sqrt(sum(~isnan(x)));
end

figure('Color','w','Position',[250 250 700 320]); hold on;
xpos = 1:numel(cats);
errorbar(xpos, MeanOut, SemOut, 'o', 'LineWidth', 2, 'MarkerSize', 8);
xlim([0.5 numel(cats)+0.5]);
xticks(xpos);
xticklabels(cats);
xtickangle(30);
ylabel('Outward endpoint (col 1)');
title('Outward endpoint by reach category (mean ± SEM across animals)');
grid on;

%% ===================== 3D ENDPOINT PLOT WITH SEM (COLORED) =====================

an   = fieldnames(EndPtMean);
%cats = {'BL','AD','AD_E','AD_L','WO','WO_E'};   % order matters for colors
cats = {'BL','AD_E','AD_L','WO_E'};   % order matters for colors
stimLcolor = rgb('RoyalBlue');
stimEcolor = rgb('LightSteelBlue');

% ---------- collect endpoints across animals ----------
E = struct();
for ci = 1:numel(cats)
    E.(cats{ci}) = [];
end

for a = 1:numel(an)
    animal = an{a};
    for ci = 1:numel(cats)
        cat = cats{ci};
        if isfield(EndPtMean.(animal), cat)
            E.(cat) = [E.(cat); EndPtMean.(animal).(cat)];
        end
    end
end

% ---------- compute mean + SEM ----------
Mean3D = nan(numel(cats),3);
SEM3D  = nan(numel(cats),3);

for ci = 1:numel(cats)
    X = E.(cats{ci});
    Mean3D(ci,:) = mean(X,1,'omitnan');
    SEM3D(ci,:)  = std(X,[],1,'omitnan') ./ sqrt(size(X,1));
end

% ---------- plot ----------
figure('Color','w','Position',[200 200 650 550]); hold on;

for ci = 1:numel(cats)

    m = Mean3D(ci,:);
    s = SEM3D(ci,:);

    % -------- color logic --------
    switch cats{ci}
        case 'BL'
            col = 'k';
        case 'WO_E'
            col = rgb('Magenta');
        case 'AD_L'
            col = stimLcolor;
        case 'AD_E'
            col = stimEcolor;
        otherwise
            col = [0.6 0.6 0.6]; % neutral gray for AD + WO
    end

    % % mean point
    % plot3(m(1), m(2), m(3), 'o', ...
    %     'MarkerSize',10, 'LineWidth',2, ...
    %     'MarkerFaceColor',col, 'MarkerEdgeColor',col);
    % 
    % % x error
    % plot3([m(1)-s(1) m(1)+s(1)], [m(2) m(2)], [m(3) m(3)], ...
    %     'Color',col, 'LineWidth',2)
    % 
    % % y error
    % plot3([m(1) m(1)], [m(2)-s(2) m(2)+s(2)], [m(3) m(3)], ...
    %     'Color',col, 'LineWidth',2)
    % 
    % % z error
    % plot3([m(1) m(1)], [m(2) m(2)], [m(3)-s(3) m(3)+s(3)], ...
    %     'Color',col, 'LineWidth',2)

    % ---------- axis mapping to match your example ----------
    % X = outward (dim1)
    % Y = lateral (dim3)
    % Z = upward (dim2)

    % mean point
    plot3(m(1), m(3), m(2), 'o', ...
        'MarkerSize',10, 'LineWidth',2, ...
        'MarkerFaceColor',col, 'MarkerEdgeColor',col);

    % x error (OUTWARD = dim1)
    plot3([m(1)-s(1) m(1)+s(1)], [m(3) m(3)], [m(2) m(2)], ...
        'Color',col, 'LineWidth',2);

    % y error (LATERAL = dim3)
    plot3([m(1) m(1)], [m(3)-s(3) m(3)+s(3)], [m(2) m(2)], ...
        'Color',col, 'LineWidth',2);

    % z error (UPWARD = dim2)
    plot3([m(1) m(1)], [m(3) m(3)], [m(2)-s(2) m(2)+s(2)], ...
        'Color',col, 'LineWidth',2);



end

% grid on
% axis equal
% xlabel('Outward endpoint (col 1)')
% ylabel('Dim 2')
% zlabel('Dim 3')
% title('Mean 3D endpoints ± SEM across animals')
%view(135,25)
%view(140,20)

% ===================== AXES STYLE (CUBE WITH TICKS) =====================

axis vis3d
axis equal
box on
grid off   % keep off to match clean cube look

% axis labels (matching your mapping)
xlabel('Outward Position')     % dim1
ylabel('Lateral Position')     % dim3
zlabel('Upward Position')      % dim2

% keep ticks but make them subtle / clean
set(gca,'LineWidth',1.5,'FontSize',11)
set(gca,'TickDir','out')   % ticks outside like paper figures

% THIS is what usually makes it match the cube orientation
set(gca,'YDir','reverse')   % lateral axis runs toward the back

% camera angle similar to figure
view(-35,18)

%% ===================== 3D ENDPOINT PLOT WITH SEM (COLORED + FITS WINDOW) =====================

an   = fieldnames(EndPtMean);
cats = {'BL','AD_E','AD_L','WO_E'};   % ONLY plot these
stimLcolor = rgb('RoyalBlue');
stimEcolor = rgb('LightSteelBlue');
woeColor   = rgb('Magenta');

% ---------- collect endpoints across animals ----------
E = struct();
for ci = 1:numel(cats)
    E.(cats{ci}) = [];
end

for a = 1:numel(an)
    animal = an{a};
    for ci = 1:numel(cats)
        cat = cats{ci};
        if isfield(EndPtMean.(animal), cat) && ~isempty(EndPtMean.(animal).(cat))
            E.(cat) = [E.(cat); EndPtMean.(animal).(cat)];
        end
    end
end

% ---------- compute mean + SEM ----------
Mean3D = nan(numel(cats),3);
SEM3D  = nan(numel(cats),3);

for ci = 1:numel(cats)
    X = E.(cats{ci});
    Mean3D(ci,:) = mean(X,1,'omitnan');
    SEM3D(ci,:)  = std(X,[],1,'omitnan') ./ sqrt(size(X,1));
end

% ---------- plot ----------
figure('Color','w'); hold on;
%set(gcf,'Units','normalized','Position',[0.08 0.10 0.70 0.75]);  % bigger window
set(gcf,'Units','normalized','Position',[0.12 0.15 0.55 0.60]);


for ci = 1:numel(cats)

    m = Mean3D(ci,:);
    s = SEM3D(ci,:);

    % -------- color logic --------
    switch cats{ci}
        case 'BL'
            col = 'k';
        case 'WO_E'
            col = woeColor;
        case 'AD_L'
            col = stimLcolor;
        case 'AD_E'
            col = stimEcolor;
    end

    % ---------- axis mapping to match your example ----------
    % X = outward (dim1)
    % Y = lateral (dim3)
    % Z = upward (dim2)

    % mean point
    plot3(m(1), m(3), m(2), 'o', ...
        'MarkerSize',15, 'LineWidth',2, ...
        'MarkerFaceColor',col, 'MarkerEdgeColor',col);

    % x error (OUTWARD = dim1)
    plot3([m(1)-s(1) m(1)+s(1)], [m(3) m(3)], [m(2) m(2)], ...
        'Color',col, 'LineWidth',2);

    % y error (LATERAL = dim3)
    plot3([m(1) m(1)], [m(3)-s(3) m(3)+s(3)], [m(2) m(2)], ...
        'Color',col, 'LineWidth',2);

    % z error (UPWARD = dim2)
    plot3([m(1) m(1)], [m(3) m(3)], [m(2)-s(2) m(2)+s(2)], ...
        'Color',col, 'LineWidth',2);

end

% ---------- expand axis limits so nothing gets clipped ----------
allPts = [];
for ci = 1:numel(cats)
    % convert each mean±sem into plotted coordinate system [X Y Z] = [1 3 2]
    m = Mean3D(ci,:);
    s = SEM3D(ci,:);
    Pplus  = [m(1) m(3) m(2)] + [s(1) s(3) s(2)];
    Pminus = [m(1) m(3) m(2)] - [s(1) s(3) s(2)];
    allPts = [allPts; Pplus; Pminus];
end

pad = 0.08 * range(allPts,1);                 % 8% padding
if any(pad==0), pad(pad==0) = 1; end          % fallback if range is 0

xlim([min(allPts(:,1))-pad(1)  max(allPts(:,1))+pad(1)]);
ylim([min(allPts(:,2))-pad(2)  max(allPts(:,2))+pad(2)]);
zlim([min(allPts(:,3))-pad(3)  max(allPts(:,3))+pad(3)]);

% ===================== AXES STYLE (CUBE WITH TICKS) =====================
axis vis3d
axis equal

% ===================== SHRINK ACTUAL PLOT AREA =====================
ax = gca;
ax.Units = 'normalized';
ax.Position = [0.20 0.20 0.60 0.60];  
% [left bottom width height]

box on
grid on

xlabel('Outward Position')   % dim1
ylabel('Lateral Position')   % dim3 (mapped to Y)
zlabel('Upward Position')    % dim2 (mapped to Z)

set(gca,'LineWidth',1.5,'FontSize',11);
set(gca,'TickDir','out');
set(gca,'YDir','reverse');   % usually needed for the "paper cube" orientation

view(-35,18);
camzoom(0.92);               % slight zoom-out so all bars stay visible

set(gcf,'Renderer','painters');   % ensures vector output
exportgraphics(gcf,'Endpoint3D.pdf','ContentType','vector');
%% Align Stimvel Zscore - aligned to average time of induced CSM 
BL_CSM_big = struct;
AD_CSM_big = struct;
ADE_CSM_big = struct;
ADL_CSM_big = struct;
WO_CSM_big = struct;
WOE_CSM_big = struct;
WOL_CSM_big = struct;

an = fieldnames(ExcitRR_BL);

for a = 1:length(an)
    BL_CSM_big.stimvel(:,:,a)  = ExcitRR_BL.(an{a}).BL_CSM.Average;
    AD_CSM_big.stimvel(:,:,a)  = ExcitRR_AD.(an{a}).AD_CSM.Average;
    ADE_CSM_big.stimvel(:,:,a) = ExcitRR_AD_E.(an{a}).ADE_CSM.Average;
    ADL_CSM_big.stimvel(:,:,a) = ExcitRR_AD_L.(an{a}).ADL_CSM.Average;
    WO_CSM_big.stimvel(:,:,a)  = ExcitRR_WO.(an{a}).WO_CSM.Average;
    WOE_CSM_big.stimvel(:,:,a) = ExcitRR_WO_E.(an{a}).WOE_CSM.Average;
    WOL_CSM_big.stimvel(:,:,a) = ExcitRR_WO_L.(an{a}).WOL_CSM.Average;
end

% Compute nanmean, nanstd, and sample count for each time point
groupList = {'BL_CSM_big', 'AD_CSM_big', 'ADE_CSM_big', ...
             'ADL_CSM_big', 'WO_CSM_big', 'WOE_CSM_big', 'WOL_CSM_big'};

for g = 1:length(groupList)
    grp = eval(groupList{g});
    grp.Average = nanmean(grp.stimvel, 3);
    grp.Std     = nanstd(grp.stimvel, 0, 3);

    % Count number of contributing samples at each time point
    N = zeros(size(grp.Average,1), 1);
    for nnum = 1:length(N)
        samp = squeeze(grp.stimvel(nnum,2,:));  % using dim 2 = outward vel
        N(nnum) = sum(~isnan(samp));
    end
    grp.N = N;

    % Run stimvel_accel
    %grp = stimvel_accel(grp);

    % Save back to workspace
    assignin('base', groupList{g}, grp);
end

%% Align Stimvel to Peak Velocity
BL_peakV_big = struct;
AD_peakV_big = struct;
ADE_peakV_big = struct;
ADL_peakV_big = struct;
WO_peakV_big = struct;
WOE_peakV_big = struct;
WOL_peakV_big = struct;

an = fieldnames(ExcitRR_BL);

for a = 1:length(an)
    BL_peakV_big.stimvel(:,:,a)  = ExcitRR_BL.(an{a}).BL_peakV.Average;
    AD_peakV_big.stimvel(:,:,a)  = ExcitRR_AD.(an{a}).AD_peakV.Average;
    ADE_peakV_big.stimvel(:,:,a) = ExcitRR_AD_E.(an{a}).ADE_peakV.Average;
    ADL_peakV_big.stimvel(:,:,a) = ExcitRR_AD_L.(an{a}).ADL_peakV.Average;
    WO_peakV_big.stimvel(:,:,a)  = ExcitRR_WO.(an{a}).WO_peakV.Average;
    WOE_peakV_big.stimvel(:,:,a) = ExcitRR_WO_E.(an{a}).WOE_peakV.Average;
    WOL_peakV_big.stimvel(:,:,a) = ExcitRR_WO_L.(an{a}).WOL_peakV.Average;
end

% Compute nanmean, nanstd, and sample count for each time point
groupList = {'BL_peakV_big', 'AD_peakV_big', 'ADE_peakV_big', ...
             'ADL_peakV_big', 'WO_peakV_big', 'WOE_peakV_big', 'WOL_peakV_big'};

for g = 1:length(groupList)
    grp = eval(groupList{g});
    grp.Average = nanmean(grp.stimvel, 3);
    grp.Std     = nanstd(grp.stimvel, 0, 3);

    % Count number of contributing samples at each time point
    N = zeros(size(grp.Average,1), 1);
    for nnum = 1:length(N)
        samp = squeeze(grp.stimvel(nnum,2,:));  % using dim 2 = outward velocity
        N(nnum) = sum(~isnan(samp));
    end
    grp.N = N;

    % Run stimvel_accel
    grp = stimvel_accel(grp);

    % Save back to workspace
    assignin('base', groupList{g}, grp);
end

%% Align Stimvel to Peak Velocity
BL_big = struct;
AD_big = struct;
ADE_big = struct;
ADL_big = struct;
WO_big = struct;
WOE_big = struct;
WOL_big = struct;

an = fieldnames(ExcitRR_BL);

for a = 1:length(an)
    BL_big.stimvel(:,:,a)  = ExcitRR_BL.(an{a}).BLN.Average;
    AD_big.stimvel(:,:,a)  = ExcitRR_AD.(an{a}).ADP.Average;
    ADE_big.stimvel(:,:,a) = ExcitRR_AD_E.(an{a}).ADP_E.Average;
    ADL_big.stimvel(:,:,a) = ExcitRR_AD_L.(an{a}).ADP_L.Average;
    WO_big.stimvel(:,:,a)  = ExcitRR_WO.(an{a}).WSH.Average;
    WOE_big.stimvel(:,:,a) = ExcitRR_WO_E.(an{a}).WSH_E.Average;
    %WOL_big.stimvel(:,:,a) = ExcitRR_WO_L.(an{a}).WSH_L.Average;
    WOL_big.stimvel(:,:,a) = ExcitRR_WO_L.(an{a}).WSH.Average;
end

% Compute nanmean, nanstd, and sample count for each time point
groupList = {'BL_big', 'AD_big', 'ADE_big', ...
             'ADL_big', 'WO_big', 'WOE_big', 'WOL_big'};

for g = 1:length(groupList)
    grp = eval(groupList{g});
    grp.Average = nanmean(grp.stimvel, 3);
    grp.Std     = nanstd(grp.stimvel, 0, 3);

    % Count number of contributing samples at each time point
    N = zeros(size(grp.Average,1), 1);
    for nnum = 1:length(N)
        samp = squeeze(grp.stimvel(nnum,2,:));  % using dim 2 = outward velocity
        N(nnum) = sum(~isnan(samp));
    end
    grp.N = N;

    % Run stimvel_accel
    grp = stimvel_accel(grp);

    % Save back to workspace
    assignin('base', groupList{g}, grp);
end

 %% Align Sham
Sham_big=struct;


for a = 1:length(am) 
  Sham_big.stimvel(:,:,a) = ExcitRR_Sham.(am{a}).BLRAlign.Average;
   
end

Sham_big.Average = nanmean(Sham_big.stimvel,3);
Sham_big.Std = nanstd(Sham_big.stimvel,0,3);

%loop through to count the number of reaches contributing to each time pt
N = zeros(length(Sham_big.Average),1); 
for nnum = 1:length(Sham_big.Average)
    samp = squeeze(Sham_big.stimvel(nnum,2,:));
    N(nnum) = length(samp(~isnan(samp)));

end
Sham_big.N = N;

%% plot stimvel - Aligned to Stim
baselinecolor=rgb('Black');
stimLcolor= rgb('RoyalBlue');
stimEcolor= rgb('LightSteelBlue');
washEcolor=rgb('Violet');
washLcolor=rgb('Magenta');
stimcolor=rgb('DodgerBlue');
Shamcolor=rgb('gray');

dim = 5;

x1=0 ;x2=0.05;
y_upper=40; y_lower=-60;
figure('Position', [100, 100, 600, 300]); hold on
%figure('Position', [100, 100, 400, 300]); hold on


StimvelPlot(BL_big,dim,baselinecolor,0);
%StimvelPlot(Sham_big,dim,Shamcolor,0);
StimvelPlot(ADE_big,dim,stimEcolor,0);
StimvelPlot(ADL_big,dim,stimLcolor,0);
% StimvelPlot(AD_big,dim,stimcolor,0);
% StimvelPlot(WOE_big,dim,washEcolor,0);
% StimvelPlot(WOL_big,dim,washLcolor,0);

xline(0.0815,'--')

hold on 
% Fill the area between the vertical lines
%fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
%xlim([-.15 .25])
xlim([0 .05])

xline(0.14)
%xlim([-0.2 0.2])
ylim('auto')
%ylim([-.5,1.5])
%ylim([0 0.27])
%ylim([-60 30])
yline(0,'--')
xlabel('Time from Stim.(s)')
ylabel('Outward Velocity(cm)')

%%
dim=2;

x1=0 ;x2=0.05;
y_upper=0.27; y_lower=0;
subplot(1,4,1); hold on 
StimvelPlot(BL_big,dim,baselinecolor,0);
StimvelPlot(ADE_big,dim,stimEcolor,0);
hold on 
% Fill the area between the vertical lines
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');

xlim([-.25 .25])
ylim([0 0.2])

xlabel('Time from Stim.(s)')
ylabel('Outward Position (cm)')
pbaspect([1 2 1]) % for position plots 

subplot(1,4,2); hold on 
y_upper=0.27; y_lower=0;
StimvelPlot(BL_big,dim,baselinecolor,0);
StimvelPlot(ADL_big,dim,stimLcolor,0);
hold on 
% Fill the area between the vertical lines
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');

xlim([-.25 .25])
ylim([0 0.2])
xlabel('Time from Stim.(s)')
ylabel('Outward Position (cm)')
pbaspect([1 2 1]) % for position plots 

subplot(1,4,3); hold on
StimvelPlot(BL_big,dim,baselinecolor,0);
StimvelPlot(WOE_big,dim,washEcolor,0);

xlim([-.25 .25])
ylim([0 0.2])
xlabel('Time from Stim.(s)')
ylabel('Outward Position (cm)')
pbaspect([1 2 1]) % for position plots 

subplot(1,4,4); hold on
StimvelPlot(BL_big,dim,baselinecolor,0);
StimvelPlot(WOL_big,dim,washLcolor,0);

xlim([-.25 .25])
ylim([0 0.2])
xlabel('Time from Stim.(s)')
ylabel('Outward Position (cm)')
pbaspect([1 2 1]) % for position plots 
%% Individual animal
% plot single stimvel curve for animal
an=fieldnames(ExcitRR_BL);
a =1;


BLN= ExcitRR_BL.(an{a}).BLN;
ADP= ExcitRR_AD.(an{a}).ADP;
ADP_E= ExcitRR_AD_E.(an{a}).ADP_E;
ADP_L= ExcitRR_AD_L.(an{a}).ADP_L;
WSH= ExcitRR_WO_E.(an{a}).WSH_E;
WSH_L=ExcitRR_WO_L.(an{a}).WSH_L;

stimLcolor= rgb('RoyalBlue');
stimEcolor= rgb('LightSteelBlue');
adaptEcolor=rgb('Violet');
adaptLcolor=rgb('Magenta');

dim = 2;

subplot(1,4,1); hold on 
StimvelPlot(BLN,dim,[.5 .5 .5],0);
StimvelPlot(ADP_E,dim,stimEcolor,0);
xlim([-.2 .2])
ylim([0 0.27])
hold on 
% Fill the area between the vertical lines
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');

xlabel('Time from Stim.(s)')
ylabel('Outward Position (cm)')
pbaspect([1 2 1]) % for position plots 

subplot(1,4,2); hold on 
StimvelPlot(BLN,dim,[.5 .5 .5],0);
StimvelPlot(ADP_L,dim,stimLcolor,0);
hold on 
% Fill the area between the vertical lines
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
xlim([-.2 .2])
ylim([0 0.27])
xlabel('Time from Stim.(s)')
pbaspect([1 2 1]) % for position plots 

subplot(1,4,3); hold on
StimvelPlot(BLN,dim,[.5 .5 .5],0);
StimvelPlot(WSH,dim,adaptEcolor,0);
xlim([-.2 .2])
ylim([0 0.27])
xlabel('Time from Stim.(s)')
pbaspect([1 2 1]) 

subplot(1,4,4); hold on
StimvelPlot(BLN,dim,[.5 .5 .5],0);
StimvelPlot(WSH_L,dim,adaptLcolor,0);
xlim([-.2 .2])
ylim([0 0.27])
xlabel('Time from Stim.(s)')
pbaspect([1 2 1]) 

%% 
dim = 5;


figure('Position', [100, 100, 600, 300]); hold on


StimvelPlot(BLN,dim,baselinecolor,0);
%StimvelPlot(Sham_big,dim,Shamcolor,0);
StimvelPlot(ADP_E,dim,stimEcolor,0);
StimvelPlot(ADP_L,dim,stimLcolor,0);
StimvelPlot(ADP,dim,stimcolor,0);
StimvelPlot(WSH,dim,washEcolor,0);
% StimvelPlot(WOL_big,dim,washLcolor,0);

xline(0.0815,'--')

hold on 
% Fill the area between the vertical lines
%fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
xlim([-.15 .25])


%xlim([-0.2 0.2])
ylim('auto')
%ylim([-.5,1.5])
%ylim([0 0.27])
%ylim([-60 30])
yline(0,'--')
xlabel('Time from Stim.(s)')
ylabel('Outward Velocity(cm)')
%%

dim = 5;
figure('Position', [100, 100, 600, 300]); hold on

StimvelPlot(BL_peakV_big,dim,baselinecolor,0);
StimvelPlot(ADE_peakV_big,dim,stimEcolor,0);
StimvelPlot(ADL_peakV_big,dim,stimLcolor,0);
StimvelPlot(WOE_peakV_big,dim,washEcolor,0);
StimvelPlot(WOL_peakV_big,dim,washLcolor,0);

hold on 

%xlim([-.1 .25])
xlim([-0.25 0.25])
ylim('auto')
%ylim([-.5,1.5])
%ylim([0 0.27])
%ylim([-20 40])
yline(0,'--')
xlabel('Time from MaxV')
ylabel('Outward Velocity(cm)')
%pbaspect([1 2 1]) % for position plots 
%%

dim = 5;
figure('Position', [100, 100, 600, 300]); hold on

StimvelPlot(BL_CSM_big,dim,baselinecolor,0);
StimvelPlot(AD_CSM_big,dim,stimcolor,0);
% StimvelPlot(ADE_CSM_big,dim,stimEcolor,0);
% StimvelPlot(ADL_CSM_big,dim,stimLcolor,0);
StimvelPlot(WOE_CSM_big,dim,washEcolor,0);
%StimvelPlot(WOL_CSM_big,dim,washLcolor,0);

hold on 

xlim([-.25 .25])
%xlim([-0.2 0.2])
ylim('auto')
%ylim([-.5,1.5])
%ylim([0 0.27])
%ylim([-20 20])
xline(0)
yline(0,'--')

xline(-0.0815,'r')
xline(-0.1365,'b')
xlabel('Time from CSM')
ylabel('Outward Velocity(cm)')
%$pbaspect([1 2 1]) % for position plots 
%% Difference Plots for stim aligned
%
% For example, just plot the y-dimension (2nd column)
% diff_avg_y = WOE_big.Average(:,8) - BL_big.Average(:,8);
% diff_avg_y=diff_avg_y*10;
% SE_diff_y = sqrt((WOE_big.Std(:,8).^2 ./ WOE_big.N) + (BL_big.Std(:,8).^2 ./ BL_big.N));

diff_avg_y = WOE_big.Average(:,5) - BL_big.Average(:,5);
diff_avg_y=diff_avg_y*10;
SE_diff_y = sqrt((WOE_big.Std(:,5).^2 ./ WOE_big.N) + (BL_big.Std(:,5).^2 ./ BL_big.N));

time = (BL_big.Average(:,1))'; % ensure it's a column vector

% Now plot
x1=0 ;x2=0.05;
y_upper=8; y_lower=-8;

figure('Position', [100, 100, 600, 300]);  
shadedErrorBar(time, diff_avg_y, SE_diff_y, 'lineProps', {'k', 'LineWidth', 2});
hold on 
% Fill the area between the vertical lines
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
csmmean_sec=0.0815;
xline(csmmean_sec,':')
xlabel('Time from Stim (ms)');
xlim([-0.25 0.25])
%ylim([-8 8])
ylabel('Delta- WOE-BL');
title('WOE-BL Kin Diff');
yline(0, '--');

%%
% For example, just plot the y-dimension (2nd column)
diff_avg_y = AD_big.Average(:,5) - BL_big.Average(:,5);
diff_avg_y=diff_avg_y*10;
SE_diff_y = sqrt((AD_big.Std(:,5).^2 ./ AD_big.N) + (BL_big.Std(:,5).^2 ./ BL_big.N));
% 
% diff_avg_y = AD_big.Average(:,8) - BL_big.Average(:,8);
% diff_avg_y=diff_avg_y*10;
% SE_diff_y = sqrt((AD_big.Std(:,8).^2 ./ AD_big.N) + (BL_big.Std(:,8).^2 ./ BL_big.N));

time = (BL_big.Average(:,1))'; % ensure it's a column vector
x1=0 ;x2=0.05;
y_upper=8; y_lower=-8;
% Now plot
figure('Position', [100, 100, 600, 300]);  
shadedErrorBar(time, diff_avg_y, SE_diff_y, 'lineProps', {'k', 'LineWidth', 2});
hold on 
% Fill the area between the vertical lines
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
csmmean_sec=0.0815;
xline(csmmean_sec,':')
xlabel('Time from Stim (ms)');
%ylim([-200 200])
ylim([-8 8])
xlim([-0.25 0.25])
ylabel('Y-Velocity Difference (AD - BL)');
title('AD-BL Kin Diff');
yline(0, '--');

%% WOE and Adapt on the same plot 

% For example, just plot the y-dimension (2nd column)
diff_avg_y = WOE_big.Average(:,5) - BL_big.Average(:,5);
diff_avg_y=diff_avg_y*10;
SE_diff_y = sqrt((WOE_big.Std(:,5).^2 ./ WOE_big.N) + (BL_big.Std(:,5).^2 ./ BL_big.N));




diff_avg_y_AD = AD_big.Average(:,5) - BL_big.Average(:,5);
diff_avg_y_AD=diff_avg_y_AD*10;
SE_diff_y_AD = sqrt((AD_big.Std(:,5).^2 ./ AD_big.N) + (BL_big.Std(:,5).^2 ./ BL_big.N));
% 
% % For example, just plot the y-dimension (2nd column)
% diff_avg_y = WOE_big.Average(:,8) - BL_big.Average(:,8);
% diff_avg_y=diff_avg_y*10;
% SE_diff_y = sqrt((WOE_big.Std(:,8).^2 ./ WOE_big.N) + (BL_big.Std(:,8).^2 ./ BL_big.N));
% 
% 
% diff_avg_y_AD = AD_big.Average(:,8) - BL_big.Average(:,8);
% diff_avg_y_AD=diff_avg_y_AD*10;
% SE_diff_y_AD = sqrt((AD_big.Std(:,8).^2 ./ AD_big.N) + (BL_big.Std(:,8).^2 ./ BL_big.N));

time = (BL_big.Average(:,1))'; % ensure it's a column vector

x1=0 ;x2=0.05;
y_upper=8; y_lower=-8;
% Now plot
figure('Position', [100, 100, 600, 300]);  

% Plot BL (or other group) in Violet
shadedErrorBar(time, diff_avg_y, SE_diff_y, 'lineProps', {'Color', [.9297    0.5078    0.9297], 'LineWidth', 2});  % Violet RGB

hold on 

% Plot AD in Dodger Blue
shadedErrorBar(time, diff_avg_y_AD, SE_diff_y_AD, 'lineProps', {'Color', [0.1172    0.5625    1.0000], 'LineWidth', 2});  % Dodger Blue RGB

% Fill the area between the vertical lines
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], ...
    [0.1172, 0.5625, 1.0000], 'FaceAlpha', 0.2, 'EdgeColor', 'none');

csmmean_sec = 0.0815;
xline(csmmean_sec, ':');

xlabel('Time from Stim (ms)');
ylabel('Y-Velocity Difference (AD - BL)');
title('AD-BL Kin Diff');

xlim([-0.15 0.25]);
%ylim([-200 200]);
%ylim('auto')
ylim([-6 8])
yline(0, '--');

%% ============================================================
% WOE + AD delta traces with p-value heatmap above
% Testing whether each delta trace is different from 0
% ============================================================

% ---------------- DELTAS ----------------
diff_WOE = (WOE_big.Average(:,5) - BL_big.Average(:,5)) * 10;
SE_WOE   = sqrt((WOE_big.Std(:,5).^2 ./ WOE_big.N) + (BL_big.Std(:,5).^2 ./ BL_big.N));

diff_AD = (AD_big.Average(:,5) - BL_big.Average(:,5)) * 10;
SE_AD   = sqrt((AD_big.Std(:,5).^2 ./ AD_big.N) + (BL_big.Std(:,5).^2 ./ BL_big.N));

time = BL_big.Average(:,1)';

% ---------------- APPROXIMATE TEST VS 0 ----------------
% t = mean / SE, then two-sided p-value
% df approximated conservatively from condition N and BL N

df_WOE = max(min(WOE_big.N, BL_big.N) - 1, 1);
df_AD  = max(min(AD_big.N,  BL_big.N) - 1, 1);

t_WOE = diff_WOE ./ SE_WOE;
t_AD  = diff_AD  ./ SE_AD;

p_WOE = 2 * (1 - tcdf(abs(t_WOE), df_WOE));
p_AD  = 2 * (1 - tcdf(abs(t_AD),  df_AD));

% avoid zeros for plotting on log scale
p_WOE(p_WOE < 1e-10) = 1e-10;
p_AD(p_AD < 1e-10)   = 1e-10;

% convert to -log10(p) so stronger significance = hotter color
logp_WOE = -log10(p_WOE);
logp_AD  = -log10(p_AD);

Pmap = [logp_WOE'; logp_AD'];   % row 1 = WOE, row 2 = AD
% ---------------- FIGURE LAYOUT ----------------
fig = figure('Color','w','Position',[100 100 700 380]);

ax_top  = axes('Position',[0.10 0.76 0.84 0.12]); hold(ax_top,'on')
ax_main = axes('Position',[0.10 0.14 0.84 0.56]); hold(ax_main,'on')

% ---------------- TOP PANEL: P-VALUE HEATMAP ----------------
% imagesc(ax_top, time, [1 2], Pmap);
% set(ax_top,'YDir','normal')
% 
% % choose color scaling
% % 1.3 = p<0.05
% % 2   = p<0.01
% % 3   = p<0.001
% caxis(ax_top, [0 3])
% colormap(ax_top, hot)
% 
% set(ax_top, 'YTick', [1 2], ...
%             'YTickLabel', {'WOE','AD'}, ...
%             'XTick', [], ...
%             'TickLength', [0 0], ...
%             'Box', 'off')
% 
% ylabel(ax_top, '-log10(p)')
% 
% xline(ax_top, 0.0815, ':', 'LineWidth', 1.2, 'Color', 'k')
% 
% % optional threshold markers on colorbar meaning:
% cb = colorbar(ax_top, 'eastoutside');
% cb.Label.String = '-log10(p)';
% cb.Ticks = [0 1.3 2 3];
% cb.TickLabels = {'1','0.05','0.01','0.001'};
% 
% xlim(ax_top, [-0.25 0.25])
% ylim(ax_top, [0.5 2.5])

% ---------------- TOP PANEL: P-VALUE HEATMAP ----------------
imagesc(ax_top, time, [1 2], Pmap);
set(ax_top,'YDir','normal')

% ----- threshold (in -log10 space) -----
pthresh = 0.05;
thresh_log = -log10(pthresh);   % ~1.3

% ----- color scaling -----
caxis(ax_top, [0 3])   % adjust upper bound if needed

% ----- custom colormap -----
n = 256;

% fraction of map below threshold
frac_gray = thresh_log / 3;   % matches caxis upper bound (3)

n_gray = round(n * frac_gray);
n_red  = n - n_gray;

% gray part (non-significant)
gray_map = repmat(linspace(0.85, 0.6, n_gray)', 1, 3);

% red gradient (significant)
red_map = [linspace(1, 0.6, n_red)', ...   % R
           linspace(0.8, 0, n_red)', ...   % G
           linspace(0.8, 0, n_red)'];     % B

cmap = [gray_map; red_map];
colormap(ax_top, cmap)

% ----- labels -----
set(ax_top, 'YTick', [1 2], ...
            'YTickLabel', {'WOE','AD'}, ...
            'XTick', [], ...
            'TickLength', [0 0], ...
            'Box', 'off')

ylabel(ax_top, '-log10(p)')

xline(ax_top, 0.0815, ':', 'LineWidth', 1.2, 'Color', 'k')

xlim(ax_top, [-0.25 0.25])
ylim(ax_top, [0.5 2.5])

% ----- colorbar -----
cb = colorbar(ax_top, 'eastoutside');
cb.Label.String = 'p-value';

cb.Ticks = [0 thresh_log 2 3];
cb.TickLabels = {'1', '0.05', '0.01', '0.001'};

% ---------------- MAIN PANEL: DELTA TRACES ----------------
shadedErrorBar(time, diff_WOE, SE_WOE, ...
    'lineProps', {'Color', [.9297 0.5078 0.9297], 'LineWidth', 2});

hold(ax_main,'on')

shadedErrorBar(time, diff_AD, SE_AD, ...
    'lineProps', {'Color', [0.1172 0.5625 1.0000], 'LineWidth', 2});

% stimulation patch
x1 = 0; 
x2 = 0.05;
y_lower = -6; 
y_upper = 8;

fill(ax_main, [x1 x2 x2 x1], [y_lower y_lower y_upper y_upper], ...
    [0.1172 0.5625 1.0000], 'FaceAlpha', 0.2, 'EdgeColor', 'none');

xline(ax_main, 0.0815, ':', 'LineWidth', 1.2);
yline(ax_main, 0, '--', 'LineWidth', 1);

xlabel(ax_main, 'Time from Stim (s)');
ylabel(ax_main, 'Y-Velocity Difference from BL');
title(ax_main, 'WOE-BL and AD-BL with p-value heatmap');

xlim(ax_main, [-0.25 0.25]);
ylim(ax_main, [-6 8]);
box(ax_main, 'off')

legend(ax_main, {'WOE-BL','AD-BL'}, 'Location', 'northwest');

% ---------------- LINK X AXES ----------------
linkaxes([ax_top ax_main], 'x')

%% CSM aligned
% WOE and Adapt on the same plot 

% For example, just plot the y-dimension (2nd column)
diff_avg_y = WOE_CSM_big.Average(:,5) - BL_CSM_big.Average(:,5);
diff_avg_y=diff_avg_y*10;
SE_diff_y = sqrt((WOE_CSM_big.Std(:,5).^2 ./ WOE_CSM_big.N) + (BL_CSM_big.Std(:,5).^2 ./ BL_CSM_big.N));


% For example, just plot the y-dimension (2nd column)
% diff_avg_y = WOL_CSM_big.Average(:,5) - BL_CSM_big.Average(:,5);
% diff_avg_y=diff_avg_y*10;
% SE_diff_y = sqrt((WOL_CSM_big.Std(:,5).^2 ./ WOL_CSM_big.N) + (BL_CSM_big.Std(:,5).^2 ./ BL_CSM_big.N));


diff_avg_y_AD = AD_CSM_big.Average(:,5) - BL_CSM_big.Average(:,5);
diff_avg_y_AD=diff_avg_y_AD*10;
SE_diff_y_AD = sqrt((AD_CSM_big.Std(:,5).^2 ./ AD_CSM_big.N) + (BL_CSM_big.Std(:,5).^2 ./ BL_CSM_big.N));

% % For example, just plot the y-dimension (2nd column)
% diff_avg_y = WOE_CSM_big.Average(:,8) - BL_CSM_big.Average(:,8);
% diff_avg_y=diff_avg_y*10;
% SE_diff_y = sqrt((WOE_CSM_big.Std(:,8).^2 ./ WOE_CSM_big.N) + (BL_CSM_big.Std(:,8).^2 ./ BL_CSM_big.N));
% 
% 
% diff_avg_y_AD = AD_CSM_big.Average(:,8) - BL_CSM_big.Average(:,8);
% diff_avg_y_AD=diff_avg_y_AD*10;
% SE_diff_y_AD = sqrt((AD_CSM_big.Std(:,8).^2 ./ AD_CSM_big.N) + (BL_CSM_big.Std(:,8).^2 ./ BL_CSM_big.N));

time = (BL_CSM_big.Average(:,1))'; % ensure it's a column vector

x1=-0.05 ;x2=0.05;
%y_upper=210; y_lower=-200;
y_upper=8; y_lower=-8;
% Now plot
figure('Position', [100, 100, 600, 300]);  

% Plot BL (or other group) in Violet
shadedErrorBar(time, diff_avg_y, SE_diff_y, 'lineProps', {'Color', [.9297    0.5078    0.9297], 'LineWidth', 2});  % Violet RGB

hold on 

% Plot AD in Dodger Blue
shadedErrorBar(time, diff_avg_y_AD, SE_diff_y_AD, 'lineProps', {'Color', [0.1172    0.5625    1.0000], 'LineWidth', 2});  % Dodger Blue RGB

% Fill the area between the vertical lines
% fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], ...
%     [0.1172, 0.5625, 1.0000], 'FaceAlpha', 0.2, 'EdgeColor', 'none');


% xline(-0.05, ':');
% xline(0.05, ':');
xline(0,'--')
xline(-0.0815,'r')
xline(-0.1365,'b')

xlabel('Time from CSM (s)');
ylabel('Y-Velocity Difference (AD - BL)');
title('AD-BL Kin Diff');

xlim([-0.25 0.25]);
%ylim([-200 210]);
%ylim('auto')
ylim([-6 8])
yline(0, '--');
%% Similar plot as above but with heatmap below 

% Z-score = (difference) / (standard error of the difference)

z_AD  = diff_avg_y_AD ./ SE_diff_y_AD;   % AD vs BL
z_WOE = diff_avg_y    ./ SE_diff_y;      % WOE vs BL

% significance threshold (z > ~1.96 means p<0.05 two-tailed)
z_thresh = 1.96;

% =====================  MAIN FIGURE  ==================================
figure('Position', [100, 100, 650, 450]);

% --- Create vertical layout: trace on top, heatmap on bottom
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% =====================  TOP PANEL: Kinematic Differences  =============
nexttile; hold on;

% WOE (violet)
shadedErrorBar(time, diff_avg_y, SE_diff_y,...
    'lineProps', {'Color', [.9297 0.5078 0.9297], 'LineWidth', 2});

% AD (dodger blue)
shadedErrorBar(time, diff_avg_y_AD, SE_diff_y_AD,...
    'lineProps', {'Color', [0.1172 0.5625 1.0000], 'LineWidth', 2});

% Stim window fill
% fill([x1 x2 x2 x1], [y_lower y_lower y_upper y_upper],...
%     [0.1172 0.5625 1.0000], 'FaceAlpha', 0.15, 'EdgeColor', 'none');

xline(-0.05, ':'); xline(0.05, ':');
yline(0,'--');

xlabel('Time from CSM (s)');
ylabel('Δ Y-Vel (cm/s)');
title('Kinematic Difference Relative to Baseline');
xlim([-0.25 0.25]);
%ylim([y_lower y_upper]);

% ======================== EFFECT SIZE HEATMAP BELOW ==========================
% If you haven't already defined the tile layout, do it once:
% tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile;  % Move to lower tile

% ---- Compute Cohen's d ----

% Pooled standard deviations
SD_pooled_AD  = sqrt((AD_CSM_big.Std(:,5).^2 + BL_CSM_big.Std(:,5).^2) / 2);
SD_pooled_WOE = sqrt((WOE_CSM_big.Std(:,5).^2 + BL_CSM_big.Std(:,5).^2) / 2);

% Cohen's d effect size
d_AD  = (AD_CSM_big.Average(:,5)  - BL_CSM_big.Average(:,5))  ./ SD_pooled_AD;
d_WOE = (WOE_CSM_big.Average(:,5) - BL_CSM_big.Average(:,5)) ./ SD_pooled_WOE;

% Effect size threshold
%d_thresh = 0.5;   % medium effect size
 d_thresh = 0.8; % switch to larger effect if preferred

% ---- Build heatmap matrix ----
Dmat = [d_AD(:)'; d_WOE(:)'];

imagesc(time, [1 2], Dmat);
colormap(brewermap([], 'PrGn'));      % diverging colormap recommended
colorbar;
caxis([-2 2]);                % keep symmetric around zero

yticks([1 2]);
yticklabels({'AD–BL', 'WOE–BL'});
xlabel('Time from CSM (s)');
title('Effect Size (Cohen''s d)');

hold on;

% Mark meaningful effect sizes
idx_AD  = find(abs(d_AD)  > d_thresh);
idx_WOE = find(abs(d_WOE) > d_thresh);

plot(time(idx_AD),   ones(size(idx_AD))*1,  'k.', 'MarkerSize', 10);
plot(time(idx_WOE),  ones(size(idx_WOE))*2, 'k.', 'MarkerSize', 10);




%% Similar plot as above but with heatmap below (Cohen's d with auto-scaled caxis)

% =====================  MAIN FIGURE  ==================================
figure('Position', [100, 100, 600, 400]);

% --- Create vertical layout: trace on top, heatmap on bottom
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% =====================  TOP PANEL: Kinematic Differences  =============
nexttile; hold on;



% WOE (violet)
shadedErrorBar(time, diff_avg_y, SE_diff_y,...
    'lineProps', {'Color', [.9297 0.5078 0.9297], 'LineWidth', 2});

% AD (dodger blue)
shadedErrorBar(time, diff_avg_y_AD, SE_diff_y_AD,...
    'lineProps', {'Color', [0.1172 0.5625 1.0000], 'LineWidth', 2});

% Stim window fill
fill([x1 x2 x2 x1], [y_lower y_lower y_upper y_upper],...
    [0.1172 0.5625 1.0000], 'FaceAlpha', 0.15, 'EdgeColor', 'none');

xline(-0.05, ':'); xline(0.05, ':');
yline(0,'--');


xlabel('Time from CSM (s)');
ylabel('Δ Y-Vel (cm/s)');
title('Kinematic Difference Relative to Baseline');
xlim([-0.25 0.25]);
%ylim([y_lower y_upper]);
ylim('auto')

% ======================== EFFECT SIZE HEATMAP BELOW ==========================
nexttile;  % Move to lower tile

% ---- Compute Cohen's d ----
SD_pooled_AD  = sqrt((AD_CSM_big.Std(:,5).^2 + BL_CSM_big.Std(:,5).^2) / 2);
SD_pooled_WOE = sqrt((WOE_CSM_big.Std(:,5).^2 + BL_CSM_big.Std(:,5).^2) / 2);

d_AD  = (AD_CSM_big.Average(:,5)  - BL_CSM_big.Average(:,5))  ./ SD_pooled_AD;
d_WOE = (WOE_CSM_big.Average(:,5) - BL_CSM_big.Average(:,5)) ./ SD_pooled_WOE;

% SD_pooled_AD  = sqrt((AD_CSM_big.Std(:,8).^2 + BL_CSM_big.Std(:,8).^2) / 2);
% SD_pooled_WOE = sqrt((WOE_CSM_big.Std(:,8).^2 + BL_CSM_big.Std(:,8).^2) / 2);
% 
% d_AD  = (AD_CSM_big.Average(:,8)  - BL_CSM_big.Average(:,8))  ./ SD_pooled_AD;
% d_WOE = (WOE_CSM_big.Average(:,8) - BL_CSM_big.Average(:,8)) ./ SD_pooled_WOE;

% Effect size threshold
d_thresh = 0.8; % large effect

% ---- Build heatmap matrix ----
Dmat = [d_AD(:)'; d_WOE(:)'];

imagesc(time, [1 2], Dmat);
%colormap(brewermap([], 'PrGn'));      % diverging colormap
colormap(slanCM('PrGn'))

% ---- Auto-scale caxis based on actual max absolute effect size ----
all_d = [d_AD(:); d_WOE(:)];
cmax = max(abs(all_d));
%caxis([-cmax cmax]);
%caxis([-3 3])
caxis([-2 2])

colorbar;

yticks([1 2]);
yticklabels({'AD–BL', 'WOE–BL'});
xlabel('Time from CSM (s)');
title('Effect Size (Cohen''s d)');

hold on;

% ---- Mark meaningful effect sizes ----
idx_AD  = find(abs(d_AD)  > d_thresh);
idx_WOE = find(abs(d_WOE) > d_thresh);

plot(time(idx_AD),   ones(size(idx_AD))*1,  'k.', 'MarkerSize', 10);
plot(time(idx_WOE),  ones(size(idx_WOE))*2, 'k.', 'MarkerSize', 10);

%% Similar plot as above but with heatmap below (Z-score vs Baseline)

% =====================  MAIN FIGURE  ==================================
figure('Position', [100, 100, 600, 400]);

tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% =====================  TOP PANEL: Kinematic Differences  =============
nexttile; hold on;

% WOE (violet)
shadedErrorBar(time, diff_avg_y, SE_diff_y,...
    'lineProps', {'Color', [.9297 0.5078 0.9297], 'LineWidth', 2});

% AD (dodger blue)
shadedErrorBar(time, diff_avg_y_AD, SE_diff_y_AD,...
    'lineProps', {'Color', [0.1172 0.5625 1.0000], 'LineWidth', 2});

% Stim window
fill([x1 x2 x2 x1], [y_lower y_lower y_upper y_upper],...
    [0.1172 0.5625 1.0000], 'FaceAlpha', 0.15, 'EdgeColor', 'none');

xline(-0.05, ':'); 
xline(0.05, ':');
yline(0,'--');

xlabel('Time from CSM (s)');
ylabel('Δ Y-Vel (cm/s)');
title('Kinematic Difference Relative to Baseline');
xlim([-0.25 0.25]);
ylim([y_lower y_upper]);

% =====================  BOTTOM PANEL: Z-SCORE HEATMAP  =================
nexttile; hold on;

% ---- Compute Z-scores relative to baseline ----
% Mean difference normalized by baseline variability
z_AD  = (AD_CSM_big.Average(:,5)  - BL_CSM_big.Average(:,5)) ./ BL_CSM_big.Std(:,5);
z_WOE = (WOE_CSM_big.Average(:,5) - BL_CSM_big.Average(:,5)) ./ BL_CSM_big.Std(:,5);

% ---- Build heatmap matrix ----
Zmat = [z_AD(:)'; z_WOE(:)'];

imagesc(time, [1 2], Zmat);
colormap(slanCM('RdBu'));   % diverging, centered at 0

% ---- Symmetric color limits ----
zmax = max(abs(Zmat(:)));
caxis([-3 3]);  % or fixed [-3 3] if you prefer

colorbar;

yticks([1 2]);
yticklabels({'AD–BL', 'WOE–BL'});
xlabel('Time from CSM (s)');
title('Z-score vs Baseline');

% ---- Statistical threshold ----
z_thresh = 2;   % ~p < 0.05 two-sided
idx_AD  = find(abs(z_AD)  > z_thresh);
idx_WOE = find(abs(z_WOE) > z_thresh);

plot(time(idx_AD),  ones(size(idx_AD))*1, 'k.', 'MarkerSize', 10);
plot(time(idx_WOE), ones(size(idx_WOE))*2, 'k.', 'MarkerSize', 10);

%% Similar plot as above but with heatmap below (BOOTSTRAP Z vs Baseline)

% =====================  PARAMETERS  ====================================
nBoot = 2000;
alpha = 0.05;

% define baseline window (adjust if needed)
baseline_idx = time >= -0.25 & time <= -0.15;

% =====================  MAIN FIGURE  ===================================
figure('Position', [100, 100, 600, 400]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% =====================  TOP PANEL: Kinematic Differences  ===============
nexttile; hold on;

% WOE (violet)
shadedErrorBar(time, diff_avg_y, SE_diff_y,...
    'lineProps', {'Color', [.9297 0.5078 0.9297], 'LineWidth', 2});

% AD (dodger blue)
shadedErrorBar(time, diff_avg_y_AD, SE_diff_y_AD,...
    'lineProps', {'Color', [0.1172 0.5625 1.0000], 'LineWidth', 2});

% Stim window
fill([x1 x2 x2 x1], [y_lower y_lower y_upper y_upper],...
    [0.1172 0.5625 1.0000], 'FaceAlpha', 0.15, 'EdgeColor', 'none');

xline(-0.05, ':');
xline(0.05, ':');
yline(0,'--');

xlabel('Time from CSM (s)');
ylabel('Δ Y-Vel (cm/s)');
title('Kinematic Difference Relative to Baseline');
xlim([-0.25 0.25]);
ylim([y_lower y_upper]);

% =====================  BOTTOM PANEL: BOOTSTRAP Z HEATMAP  ===============
nexttile; hold on;

% -----------------------------------------------------------------------
% Data assumed shape: [nTrials x nTime]
AD_vals  = squeeze(AD_CSM_big.Average(:,5,:));
WOE_vals = squeeze(WOE_CSM_big.Average(:,5,:));
BL_vals  = squeeze(BL_CSM_big.Average(:,5,:));

nT = size(BL_vals,2);

% -----------------------------------------------------------------------
% Bootstrap baseline distribution
boot_baseline = zeros(nBoot,1);

for b = 1:nBoot
    samp = BL_vals(randi(size(BL_vals,1),[],1),baseline_idx);
    boot_baseline(b) = mean(samp(:));
end

mu0  = mean(boot_baseline);
sig0 = std(boot_baseline);

% -----------------------------------------------------------------------
% Bootstrap time-resolved means
boot_AD  = zeros(nBoot,nT);
boot_WOE = zeros(nBoot,nT);

for b = 1:nBoot
    idxAD  = randi(size(AD_vals,1),[],1);
    idxWOE = randi(size(WOE_vals,1),[],1);

    boot_AD(b,:)  = mean(AD_vals(idxAD,:),1);
    boot_WOE(b,:) = mean(WOE_vals(idxWOE,:),1);
end

% -----------------------------------------------------------------------
% Bootstrap Z-scores
z_AD  = (mean(boot_AD,1)  - mu0) ./ sig0;
z_WOE = (mean(boot_WOE,1) - mu0) ./ sig0;

Zmat = [z_AD; z_WOE];

imagesc(time, [1 2], Zmat);
colormap(slanCM('RdBu'));
caxis([-3 3]);
colorbar;

yticks([1 2]);
yticklabels({'AD–BL', 'WOE–BL'});
xlabel('Time from CSM (s)');
title('Bootstrap Z-score vs Baseline');

% -----------------------------------------------------------------------
% Bootstrap p-values
p_AD  = mean(abs(boot_AD  - mu0) >= abs(mean(boot_AD,1)  - mu0),1);
p_WOE = mean(abs(boot_WOE - mu0) >= abs(mean(boot_WOE,1) - mu0),1);

sig_AD  = p_AD  < alpha;
sig_WOE = p_WOE < alpha;

% -----------------------------------------------------------------------
% Overlay significance dots
plot(time(sig_AD),  ones(sum(sig_AD),1)*1, 'k.', 'MarkerSize', 10);
plot(time(sig_WOE), ones(sum(sig_WOE),1)*2, 'k.', 'MarkerSize', 10);

%% Kinematics + Heatmap with exact x-axis match and custom heights
figure('Position', [100, 100, 650, 450]);

% Define axes heights (as fraction of figure height)
topHeight = 0.65;
bottomHeight = 0.30;
leftMargin = 0.12;
bottomMargin = 0.1;
axWidth = 0.8;

% =====================  TOP PANEL: Kinematic Differences  =============
ax1 = axes('Position', [leftMargin, bottomMargin + bottomHeight + 0.05, axWidth, topHeight]);
axes(ax1); hold on;   % Set current axes

% WOE (violet)
shadedErrorBar(time, diff_avg_y, SE_diff_y,...
    'lineProps', {'Color', [.9297 0.5078 0.9297], 'LineWidth', 2});

% AD (dodger blue)
shadedErrorBar(time, diff_avg_y_AD, SE_diff_y_AD,...
    'lineProps', {'Color', [0.1172 0.5625 1.0000], 'LineWidth', 2});

% Stim window fill
fill([x1 x2 x2 x1], [y_lower y_lower y_upper y_upper],...
    [0.1172 0.5625 1.0000], 'FaceAlpha', 0.15, 'EdgeColor', 'none');

xline(-0.05, ':'); xline(0.05, ':');
yline(0,'--');

xlabel('Time from CSM (s)');
ylabel('Δ Y-Vel (cm/s)');
title('Kinematic Difference Relative to Baseline');
xlim([-0.25 0.25]);
ylim([y_lower y_upper]);

% ======================== EFFECT SIZE HEATMAP BELOW ==========================
ax2 = axes('Position', [leftMargin, bottomMargin, axWidth, bottomHeight]);
axes(ax2); hold on;  % Set current axes

% Compute Cohen's d
SD_pooled_AD  = sqrt((AD_CSM_big.Std(:,5).^2 + BL_CSM_big.Std(:,5).^2) / 2);
SD_pooled_WOE = sqrt((WOE_CSM_big.Std(:,5).^2 + BL_CSM_big.Std(:,5).^2) / 2);

d_AD  = (AD_CSM_big.Average(:,5)  - BL_CSM_big.Average(:,5))  ./ SD_pooled_AD;
d_WOE = (WOE_CSM_big.Average(:,5) - BL_CSM_big.Average(:,5)) ./ SD_pooled_WOE;

Dmat = [d_AD(:)'; d_WOE(:)'];

imagesc(time, [1 2], Dmat);
colormap(slanCM('PrGn'))
caxis([-3 3]);

colorbar;

yticks([1 2]);
yticklabels({'AD–BL', 'WOE–BL'});
xlabel('Time from CSM (s)');
title('Effect Size (Cohen''s d)');

% Mark meaningful effect sizes
d_thresh = 0.8;
idx_AD  = find(abs(d_AD)  > d_thresh);
idx_WOE = find(abs(d_WOE) > d_thresh);

plot(time(idx_AD),   ones(size(idx_AD))*1,  'k.', 'MarkerSize', 10);
plot(time(idx_WOE),  ones(size(idx_WOE))*2, 'k.', 'MarkerSize', 10);

% ======================== LINK X-AXES ==========================
linkaxes([ax1, ax2],'x');
xlim(ax1, [-0.25 0.25]);
xlim(ax2, [-0.25 0.25]);

%
 %=====================  STATISTICAL TEST INSTEAD OF COHEN'S D =====================

% extract time-series distributions
% each column = timepoint, each row = trial/unit (match your data structure)
BL_vals  = BL_CSM_big.stimvel(:,5,:);    % [nTrials x 1 x nTime]
AD_vals  = AD_CSM_big.stimvel(:,5,:);
WOE_vals = WOE_CSM_big.stimvel(:,5,:);

BL_vals  = squeeze(BL_vals);   % now [nTrials x nTime]
AD_vals  = squeeze(AD_vals);
WOE_vals = squeeze(WOE_vals);

nTime = size(BL_vals,2);

p_AD  = nan(1,nTime);
p_WOE = nan(1,nTime);

% --- perform per-timepoint t-test ---
for t = 1:nTime
    p_AD(t)  = ttest2(AD_vals(:,t),  BL_vals(:,t),  'Vartype','unequal');
    p_WOE(t) = ttest2(WOE_vals(:,t), BL_vals(:,t), 'Vartype','unequal');
end

% FDR correction
[~, crit_AD]  = fdr_bh(p_AD, 0.05);
[~, crit_WOE] = fdr_bh(p_WOE, 0.05);

sig_AD  = p_AD  < crit_AD;
sig_WOE = p_WOE < crit_WOE;

% ---- convert to heatmap matrix: use -log10(p) ----
logp_AD  = -log10(p_AD);
logp_WOE = -log10(p_WOE);

Pmat = [logp_AD; logp_WOE];

% plotting replaces the Cohen’s d block:
imagesc(time, [1 2], Pmat);
colormap(slanCM('PrGn'));
caxis([0 5]);                     % 5 = p = 1e-5

colorbar;
yticks([1 2]);
yticklabels({'AD–BL', 'WOE–BL'});
xlabel('Time from CSM (s)');
title('-log_{10}(p) From Two-sample t-Test');

hold on;

% ---- Mark statistically significant timepoints ----
plot(time(sig_AD),  ones(1,sum(sig_AD))*1,  'k.', 'MarkerSize', 10);
plot(time(sig_WOE), ones(1,sum(sig_WOE))*2, 'k.', 'MarkerSize', 10);

%% Just heatmap
% Now plot
figure('Position', [100, 100, 600, 100]);  
imagesc(time, [1 2], Dmat);
%colormap(brewermap([], 'PrGn'));      % diverging colormap
colormap(slanCM('PrGn'))

% ---- Auto-scale caxis based on actual max absolute effect size ----
all_d = [d_AD(:); d_WOE(:)];
cmax = max(abs(all_d));
%caxis([-cmax cmax]);
caxis([-2 2])

colorbar;

yticks([1 2]);
yticklabels({'AD–BL', 'WOE–BL'});
xlabel('Time from CSM (s)');
title('Effect Size (Cohen''s d)');

hold on;

% ---- Mark meaningful effect sizes ----
idx_AD  = find(abs(d_AD)  > d_thresh);
idx_WOE = find(abs(d_WOE) > d_thresh);

plot(time(idx_AD),   ones(size(idx_AD))*1,  'k.', 'MarkerSize', 10);
plot(time(idx_WOE),  ones(size(idx_WOE))*2, 'k.', 'MarkerSize', 10);

%%
% % For example, just plot the y-dimension (2nd column)
% diff_avg_y = AD_CSM_big.Average(:,5) - BL_CSM_big.Average(:,5);
% diff_avg_y=diff_avg_y*10;
% SE_diff_y = sqrt((AD_CSM_big.Std(:,5).^2 ./ AD_CSM_big.N) + (BL_CSM_big.Std(:,5).^2 ./ BL_CSM_big.N));

diff_avg_y = AD_CSM_big.Average(:,8) - BL_CSM_big.Average(:,8);
diff_avg_y=diff_avg_y*10;
SE_diff_y = sqrt((AD_CSM_big.Std(:,8).^2 ./ AD_CSM_big.N) + (BL_CSM_big.Std(:,8).^2 ./ BL_CSM_big.N));


time = (BL_CSM_big.Average(:,1))'; % ensure it's a column vector
x1=0 ;x2=0.05;
y_upper=20; y_lower=-20;
% Now plot
figure('Position', [100, 100, 600, 300]);  
shadedErrorBar(time, diff_avg_y, SE_diff_y, 'lineProps', {'k', 'LineWidth', 2});

xline(0,'--')
xlabel('Time from CSM (ms)');
%ylim([-200 200])
%ylim([-8 8])
ylim('auto')
xlim([-0.25 0.25])
ylabel('Y-Velocity Difference (AD - BL)');
title('AD-BL Kin Diff');
yline(0, '--');

%% Difference Plots for Max Velocity aligned
% For example, just plot the y-dimension (2nd column)
% diff_avg_y = WOE_big.Average(:,8) - BL_big.Average(:,8);
% diff_avg_y=diff_avg_y*10;
% SE_diff_y = sqrt((WOE_big.Std(:,8).^2 ./ WOE_big.N) + (BL_big.Std(:,8).^2 ./ BL_big.N));

diff_avg_y = WOE_big.Average(:,5) - BL_big.Average(:,5);
diff_avg_y=diff_avg_y*10;
SE_diff_y = sqrt((WOE_big.Std(:,5).^2 ./ WOE_big.N) + (BL_big.Std(:,5).^2 ./ BL_big.N));

time = (BL_big.Average(:,1))'; % ensure it's a column vector

% Now plot
x1=0 ;x2=0.05;
y_upper=8; y_lower=-8;

figure('Position', [100, 100, 600, 300]);  
shadedErrorBar(time, diff_avg_y, SE_diff_y, 'lineProps', {'k', 'LineWidth', 2});
hold on 
% Fill the area between the vertical lines
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
csmmean_sec=0.0815;
xline(csmmean_sec,':')
xlabel('Time from Stim (ms)');
xlim([-0.25 0.25])
ylim([-8 8])
ylabel('Delta- WOE-BL');
title('WOE-BL Kin Diff');
yline(0, '--');

%%
% For example, just plot the y-dimension (2nd column)
diff_avg_y = AD_big.Average(:,5) - BL_big.Average(:,5);
diff_avg_y=diff_avg_y*10;
SE_diff_y = sqrt((AD_big.Std(:,5).^2 ./ AD_big.N) + (BL_big.Std(:,5).^2 ./ BL_big.N));
% 
% diff_avg_y = AD_big.Average(:,8) - BL_big.Average(:,8);
% diff_avg_y=diff_avg_y*10;
% SE_diff_y = sqrt((AD_big.Std(:,8).^2 ./ AD_big.N) + (BL_big.Std(:,8).^2 ./ BL_big.N));

time = (BL_big.Average(:,1))'; % ensure it's a column vector
x1=0 ;x2=0.05;
y_upper=8; y_lower=-8;
% Now plot
figure('Position', [100, 100, 600, 300]);  
shadedErrorBar(time, diff_avg_y, SE_diff_y, 'lineProps', {'k', 'LineWidth', 2});
hold on 
% Fill the area between the vertical lines
fill([x1, x2, x2, x1], [y_lower, y_lower, y_upper, y_upper], [ 0.1172    0.5625    1.0000],'FaceAlpha',0.2,'EdgeColor','none');
csmmean_sec=0.0815;
xline(csmmean_sec,':')
xlabel('Time from Stim (ms)');
%ylim([-200 200])
ylim([-8 8])
xlim([-0.25 0.25])
ylabel('Y-Velocity Difference (AD - BL)');
title('AD-BL Kin Diff');
yline(0, '--');

%% New method for quantification 
%% ============================================================
% TIME-RESOLVED STATS (AD, WOE vs BL) USING stimvel(:,:,animal)
% Plots delta traces + significance above traces
% ============================================================



% ---------------- SETTINGS ----------------
var_col = 5;        % velocity column
scale_factor = 10;
alpha = 0.05;
min_cluster = 5;

col_AD  = [0.1172 0.5625 1.0000];      % blue
col_WOE = [.9297 0.5078 0.9297];       % violet
col_cmp = [0 0 0];

xlim_vals = [-0.25 0.25];
ylim_vals = [-6 8];

% ---------------- GET TIME ----------------
time = BL_big.stimvel(:,1,1)';   % assumes col 1 = time

% ---------------- EXTRACT PER-ANIMAL ----------------
nAnimals = size(BL_big.stimvel,3);
nTime    = size(BL_big.stimvel,1);

BL = nan(nTime, nAnimals);
AD = nan(nTime, nAnimals);
WOE = nan(nTime, nAnimals);

for a = 1:nAnimals
    BL(:,a)  = BL_big.stimvel(:,var_col,a);
    AD(:,a)  = AD_big.stimvel(:,var_col,a);
    WOE(:,a) = WOE_big.stimvel(:,var_col,a);
end

% ---------------- DELTAS ----------------
delta_AD  = (AD  - BL) * scale_factor;
delta_WOE = (WOE - BL) * scale_factor;

% ---------------- MEAN + SEM ----------------
mean_AD  = mean(delta_AD,2,'omitnan')';
mean_WOE = mean(delta_WOE,2,'omitnan')';

sem_AD  = std(delta_AD,0,2,'omitnan') ./ sqrt(sum(~isnan(delta_AD),2));
sem_WOE = std(delta_WOE,0,2,'omitnan') ./ sqrt(sum(~isnan(delta_WOE),2));

sem_AD  = sem_AD';
sem_WOE = sem_WOE';

% ---------------- STATS ----------------
p_AD      = nan(nTime,1);
p_WOE     = nan(nTime,1);
p_ADvsWOE = nan(nTime,1);

for t = 1:nTime
    
    x_ad  = delta_AD(t,:);
    x_woe = delta_WOE(t,:);
    
    x_ad  = x_ad(~isnan(x_ad));
    x_woe = x_woe(~isnan(x_woe));
    
    % AD vs baseline (delta != 0)
    if ~isempty(x_ad)
        p_AD(t) = signrank(x_ad,0);
    end
    
    % WOE vs baseline
    if ~isempty(x_woe)
        p_WOE(t) = signrank(x_woe,0);
    end
    
    % AD vs WOE
    both = ~isnan(delta_AD(t,:)) & ~isnan(delta_WOE(t,:));
    if any(both)
        p_ADvsWOE(t) = signrank(delta_AD(t,both), delta_WOE(t,both));
    end
end

% ---------------- FDR ----------------
p_AD_corr      = fdr_bh(p_AD);
p_WOE_corr     = fdr_bh(p_WOE);
p_ADvsWOE_corr = fdr_bh(p_ADvsWOE);

sig_AD      = p_AD_corr      < alpha;
sig_WOE     = p_WOE_corr     < alpha;
sig_ADvsWOE = p_ADvsWOE_corr < alpha;

% ---------------- CLUSTER FILTER ----------------
sig_AD      = keep_clusters(sig_AD, min_cluster);
sig_WOE     = keep_clusters(sig_WOE, min_cluster);
sig_ADvsWOE = keep_clusters(sig_ADvsWOE, min_cluster);

% ---------------- PLOT ----------------
figure('Position',[100 100 700 340]); hold on

% WOE
shadedErrorBar(time, mean_WOE, sem_WOE, ...
    'lineProps', {'Color', col_WOE, 'LineWidth', 2});

% AD
shadedErrorBar(time, mean_AD, sem_AD, ...
    'lineProps', {'Color', col_AD, 'LineWidth', 2});

% refs
yline(0,'--');
xline(0.0815,':');

xlim(xlim_vals)
ylim(ylim_vals)

xlabel('Time (s)')
ylabel('Velocity Δ from BL')
title('AD & WOE vs BL with stats')

% ---------------- SIGNIFICANCE ABOVE ----------------
y1 = ylim_vals(2) - 0.5;
y2 = ylim_vals(2) - 1.0;
y3 = ylim_vals(2) - 1.5;

plot(time(sig_AD),      y1*ones(1,sum(sig_AD)),      '.', 'Color', col_AD,  'MarkerSize',10)
plot(time(sig_WOE),     y2*ones(1,sum(sig_WOE)),     '.', 'Color', col_WOE, 'MarkerSize',10)
plot(time(sig_ADvsWOE), y3*ones(1,sum(sig_ADvsWOE)), '.', 'Color', col_cmp, 'MarkerSize',10)

text(xlim_vals(1)+0.01, y1, 'AD≠BL','Color',col_AD)
text(xlim_vals(1)+0.01, y2, 'WOE≠BL','Color',col_WOE)
text(xlim_vals(1)+0.01, y3, 'AD≠WOE','Color',col_cmp)

legend({'WOE-BL','AD-BL'}, 'Location','northwest')
box off

% ---------------- SUMMARY ----------------
fprintf('\nSignificant bins:\n')
fprintf('AD vs BL   : %d\n', sum(sig_AD))
fprintf('WOE vs BL  : %d\n', sum(sig_WOE))
fprintf('AD vs WOE  : %d\n', sum(sig_ADvsWOE))
%% Slighly different method 

%% ============================================================
% TEST WHETHER DELTA TRACES DIFFER FROM 0
% Uses:
%   BL_big.stimvel(:,:,animal)
%   AD_big.stimvel(:,:,animal)
%   WOE_big.stimvel(:,:,animal)
% ============================================================

% ---------------- SETTINGS ----------------
var_col = 5;        % velocity column
scale_factor = 10;
alpha = 0.05;

col_AD  = [0.1172 0.5625 1.0000];      % blue
col_WOE = [.9297 0.5078 0.9297];       % violet

xlim_vals = [-0.25 0.25];
ylim_vals = [-6 8];

% ---------------- GET TIME ----------------
time = BL_big.stimvel(:,1,1)';   % assumes col 1 = time

% ---------------- EXTRACT PER-ANIMAL DATA ----------------
nAnimals = size(BL_big.stimvel,3);
nTime    = size(BL_big.stimvel,1);

BL  = nan(nTime,nAnimals);
AD  = nan(nTime,nAnimals);
WOE = nan(nTime,nAnimals);

for animal = 1:nAnimals
    BL(:,animal)  = BL_big.stimvel(:,var_col,animal);
    AD(:,animal)  = AD_big.stimvel(:,var_col,animal);
    WOE(:,animal) = WOE_big.stimvel(:,var_col,animal);
end

% ---------------- COMPUTE DELTAS ----------------
delta_AD  = (AD  - BL) * scale_factor;    % time x animal
delta_WOE = (WOE - BL) * scale_factor;    % time x animal

% ---------------- MEAN + SEM FOR PLOTTING ----------------
mean_AD  = mean(delta_AD,2,'omitnan')';
mean_WOE = mean(delta_WOE,2,'omitnan')';

sem_AD  = std(delta_AD,0,2,'omitnan') ./ sqrt(sum(~isnan(delta_AD),2));
sem_WOE = std(delta_WOE,0,2,'omitnan') ./ sqrt(sum(~isnan(delta_WOE),2));

sem_AD  = sem_AD';
sem_WOE = sem_WOE';

% ---------------- TEST DELTA AGAINST 0 ----------------
p_AD  = nan(nTime,1);
p_WOE = nan(nTime,1);

for t = 1:nTime
    
    x_ad = delta_AD(t,:);
    x_ad = x_ad(~isnan(x_ad));
    
    x_woe = delta_WOE(t,:);
    x_woe = x_woe(~isnan(x_woe));
    
    if ~isempty(x_ad)
        p_AD(t) = signrank(x_ad,0);   % is AD-BL different from 0?
    end
    
    if ~isempty(x_woe)
        p_WOE(t) = signrank(x_woe,0); % is WOE-BL different from 0?
    end
end

sig_AD  = p_AD  < alpha;
sig_WOE = p_WOE < alpha;

% ---------------- PLOT ----------------
figure('Position',[100 100 700 340]); 
hold on

% WOE delta
shadedErrorBar(time, mean_WOE, sem_WOE, ...
    'lineProps', {'Color', col_WOE, 'LineWidth', 2});

% AD delta
shadedErrorBar(time, mean_AD, sem_AD, ...
    'lineProps', {'Color', col_AD, 'LineWidth', 2});

% reference lines
yline(0,'--','LineWidth',1);
xline(0.0815,':','LineWidth',1.2);

xlim(xlim_vals)
ylim(ylim_vals)
box off

xlabel('Time from stim (s)')
ylabel('\Delta velocity from BL')
title('Delta velocity relative to baseline')

% ---------------- SIGNIFICANCE ABOVE TRACES ----------------
y1 = ylim_vals(2) - 0.4;   % AD row
y2 = ylim_vals(2) - 0.9;   % WOE row

plot(time(sig_AD),  y1*ones(1,sum(sig_AD)),  '.', 'Color', col_AD,  'MarkerSize', 10)
plot(time(sig_WOE), y2*ones(1,sum(sig_WOE)), '.', 'Color', col_WOE, 'MarkerSize', 10)

text(xlim_vals(1)+0.01, y1, 'AD \neq 0',  'Color', col_AD)
text(xlim_vals(1)+0.01, y2, 'WOE \neq 0', 'Color', col_WOE)

legend({'WOE-BL','AD-BL'}, 'Location','northwest')

% ---------------- OPTIONAL: PRINT COUNTS ----------------
fprintf('\nSignificant bins vs 0:\n')
fprintf('AD  : %d / %d bins\n', sum(sig_AD),  nTime)
fprintf('WOE : %d / %d bins\n', sum(sig_WOE), nTime)

%% With zscore 
%% ============================================================
% WOE + AD delta traces with Z-SCORE heatmap above
% ============================================================

% ---------------- DELTAS ----------------
diff_WOE = (WOE_big.Average(:,5) - BL_big.Average(:,5)) * 10;
SE_WOE   = sqrt((WOE_big.Std(:,5).^2 ./ WOE_big.N) + (BL_big.Std(:,5).^2 ./ BL_big.N));

diff_AD = (AD_big.Average(:,5) - BL_big.Average(:,5)) * 10;
SE_AD   = sqrt((AD_big.Std(:,5).^2 ./ AD_big.N) + (BL_big.Std(:,5).^2 ./ BL_big.N));

time = BL_big.Average(:,1)';

% ---------------- Z-SCORE RELATIVE TO BL ----------------
BL_std = BL_big.Std(:,5);

% avoid divide-by-zero
BL_std(BL_std == 0) = NaN;

z_WOE = (WOE_big.Average(:,5) - BL_big.Average(:,5)) ./ BL_std;
z_AD  = (AD_big.Average(:,5)  - BL_big.Average(:,5)) ./ BL_std;

% row 1 = WOE, row 2 = AD
Zmap = [z_WOE'; z_AD'];

% ---------------- FIGURE LAYOUT ----------------
fig = figure('Color','w','Position',[100 100 700 380]);

ax_top  = axes('Position',[0.10 0.76 0.84 0.12]); hold(ax_top,'on')
ax_main = axes('Position',[0.10 0.14 0.84 0.56]); hold(ax_main,'on')

% ---------------- TOP PANEL: Z-SCORE HEATMAP ----------------
imagesc(ax_top, time, [1 2], Zmap);
set(ax_top,'YDir','normal')

% color scale for z-score
caxis(ax_top, [-3 3])

% diverging colormap: blue -> white/gray -> red
n = 256;
n_half = n/2;

blue_map = [linspace(0,0.85,n_half)', ...
            linspace(0.2,0.85,n_half)', ...
            linspace(0.8,0.85,n_half)'];

red_map  = [linspace(0.85,0.8,n_half)', ...
            linspace(0.85,0,n_half)', ...
            linspace(0.85,0,n_half)'];


colormap(slanCM('PrGn'))


set(ax_top, 'YTick', [1 2], ...
            'YTickLabel', {'WOE','AD'}, ...
            'XTick', [], ...
            'TickLength', [0 0], ...
            'Box', 'off')

ylabel(ax_top, 'z')

xline(ax_top, 0.0815, ':', 'LineWidth', 1.2, 'Color', 'k')

xlim(ax_top, [-0.25 0.25])
ylim(ax_top, [0.5 2.5])

cb = colorbar(ax_top, 'eastoutside');
cb.Label.String = 'z-score';
cb.Ticks = [-3 -2 -1 0 1 2 3];

% ---------------- MAIN PANEL: DELTA TRACES ----------------
axes(ax_main)
hold(ax_main,'on')

shadedErrorBar(time, diff_WOE, SE_WOE, ...
    'lineProps', {'Color', [.9297 0.5078 0.9297], 'LineWidth', 2});

shadedErrorBar(time, diff_AD, SE_AD, ...
    'lineProps', {'Color', [0.1172 0.5625 1.0000], 'LineWidth', 2});

% stimulation patch
x1 = 0; 
x2 = 0.05;
y_lower = -6; 
y_upper = 8;

fill(ax_main, [x1 x2 x2 x1], [y_lower y_lower y_upper y_upper], ...
    [0.1172 0.5625 1.0000], 'FaceAlpha', 0.2, 'EdgeColor', 'none');

xline(ax_main, 0.0815, ':', 'LineWidth', 1.2);
yline(ax_main, 0, '--', 'LineWidth', 1);

xlabel(ax_main, 'Time from Stim (s)');
ylabel(ax_main, 'Y-Velocity Difference from BL');
title(ax_main, 'WOE-BL and AD-BL with z-score heatmap');

xlim(ax_main, [-0.25 0.25]);
ylim(ax_main, [-6 8]);
box(ax_main, 'off')

legend(ax_main, {'WOE-BL','AD-BL'}, 'Location', 'northwest');

% ---------------- LINK X AXES ----------------
linkaxes([ax_top ax_main], 'x')



%% Other methods 

%% ============================================================
% DELTA VELOCITY TRACE + TOP HEATMAP OF MEAN Z-SCORE
%
% Main panel:
%   mean delta velocity (AD-BL and WOE-BL), no SEM
%
% Top panel:
%   2-row heatmap of mean z-score relative to BL
%   z = (mean(cond) - mean(BL)) ./ std(BL)
%
% Uses:
%   BL_big.stimvel(:,:,animal)
%   AD_big.stimvel(:,:,animal)
%   WOE_big.stimvel(:,:,animal)
% ============================================================


% ---------------- SETTINGS ----------------
var_col      = 5;      % velocity column
scale_factor = 10;     % convert to desired units

col_AD  = [0.1172 0.5625 1.0000];
col_WOE = [.9297 0.5078 0.9297];

xlim_vals = [-0.15 0.25];
ylim_main = [-6 8];
csmmean_sec = 0.0815;

% heatmap color limits (adjust if needed)
zlim_vals = [-3 3];

% ---------------- GET TIME ----------------
time = BL_big.stimvel(:,1,1)';   % assumes col 1 is time
nAnimals = size(BL_big.stimvel,3);
nTime    = size(BL_big.stimvel,1);

% ---------------- EXTRACT DATA ----------------
BL  = nan(nTime,nAnimals);
AD  = nan(nTime,nAnimals);
WOE = nan(nTime,nAnimals);

for a = 1:nAnimals
    BL(:,a)  = BL_big.stimvel(:,var_col,a);
    AD(:,a)  = AD_big.stimvel(:,var_col,a);
    WOE(:,a) = WOE_big.stimvel(:,var_col,a);
end

% ---------------- MEAN TRACES ----------------
BL_mean  = mean(BL,  2, 'omitnan');
AD_mean  = mean(AD,  2, 'omitnan');
WOE_mean = mean(WOE, 2, 'omitnan');

% ---------------- DELTA VELOCITY (MAIN PANEL) ----------------
delta_AD_mean  = (AD_mean  - BL_mean)  * scale_factor;
delta_WOE_mean = (WOE_mean - BL_mean)  * scale_factor;

% ---------------- MEAN Z-SCORE RELATIVE TO BL ----------------
BL_std = std(BL, 0, 2, 'omitnan');
BL_std(BL_std == 0) = NaN;

z_AD_mean  = (AD_mean  - BL_mean) ./ BL_std;
z_WOE_mean = (WOE_mean - BL_mean) ./ BL_std;

% 2-row heatmap matrix
% Row 1 = WOE
% Row 2 = AD
Zmap = [z_WOE_mean'; z_AD_mean'];

% ---------------- FIGURE ----------------
fig = figure('Color','w','Position',[100 100 720 400]);

% --- top heatmap axis
ax_top = axes('Position',[0.10 0.77 0.84 0.12]);
hold(ax_top,'on')

% --- main trace axis
ax_main = axes('Position',[0.10 0.14 0.84 0.56]);
hold(ax_main,'on')

% ---------------- TOP PANEL: HEATMAP ----------------
axes(ax_top)

% Build edges so image fills bins nicely
if numel(time) > 1
    %dt = median(diff(time));
    dt = median(builtin('diff', time));
else
    dt = 0.001;
end

x_edges = [time - dt/2, time(end) + dt/2];
y_edges = [0.5 1.5 2.5];

imagesc(ax_top, time, [1 2], Zmap);
set(ax_top,'YDir','normal')


colormap(slanCM('PrGn'))
caxis(ax_top, zlim_vals)

xline(ax_top, csmmean_sec, ':', 'LineWidth', 1.2, 'Color', 'k')

xlim(ax_top, xlim_vals)
ylim(ax_top, [0.5 2.5])

set(ax_top, 'XTick', [], ...
            'YTick', [1 2], ...
            'YTickLabel', {'WOE','AD'}, ...
            'Box', 'off', ...
            'TickLength', [0 0])

ylabel(ax_top, 'z')

% ---------------- MAIN PANEL: MEAN DELTA TRACES ----------------
axes(ax_main)

plot(ax_main, time, delta_WOE_mean, 'Color', col_WOE, 'LineWidth', 2.5);
plot(ax_main, time, delta_AD_mean,  'Color', col_AD,  'LineWidth', 2.5);

yline(ax_main, 0, '--', 'LineWidth', 1);
xline(ax_main, csmmean_sec, ':', 'LineWidth', 1.2);

xlim(ax_main, xlim_vals)
ylim(ax_main, ylim_main)

xlabel(ax_main, 'Time from stim (s)')
ylabel(ax_main, '\Delta velocity from BL')
title(ax_main, 'Mean delta velocity with mean z-score heatmap')
box(ax_main, 'off')

legend(ax_main, {'WOE-BL','AD-BL'}, 'Location', 'northwest')

% ---------------- COLORBAR ----------------
% cb = colorbar(ax_top, 'eastoutside');
% cb.Label.String = 'Mean z-score';

% ---------------- LINK X AXES ----------------
linkaxes([ax_top ax_main], 'x')


%% ---------------- Z-SCORE HEATMAP FIGURE ----------------

figure('Color','w','Position',[100 100 700 150]);
ax = axes; hold(ax,'on')

imagesc(ax, time, [1 2], Zmap);
set(ax,'YDir','normal')

colormap(slanCM('PrGn'))   % change if you want
caxis(ax, [-3 3])      % adjust range as needed

% labels
set(ax,'YTick',[1 2], ...
       'YTickLabel',{'WOE','AD'}, ...
       'TickLength',[0 0])

xlabel('Time (s)')
ylabel('')

% reference line
xline(ax, 0.0815, ':', 'LineWidth', 1.2, 'Color', 'k')

xlim([min(time) max(time)])
box off

% colorbar
cb = colorbar;
cb.Label.String = 'Mean z-score';
%% ============================================================
function sig_out = keep_clusters(sig_in, min_len)

sig_in = logical(sig_in(:)');
sig_out = false(size(sig_in));

d = diff([0 sig_in 0]);
start_idx = find(d==1);
end_idx   = find(d==-1)-1;

for i = 1:length(start_idx)
    if (end_idx(i)-start_idx(i)+1) >= min_len
        sig_out(start_idx(i):end_idx(i)) = true;
    end
end
end

function p_fdr = fdr_bh(p)

    p = p(:);
    n = length(p);
    
    % handle NaNs
    nan_mask = isnan(p);
    p_valid = p(~nan_mask);
    
    [p_sorted, sort_idx] = sort(p_valid);
    m = length(p_sorted);
    
    % BH correction
    q = (1:m)'/m * 0.05;  % alpha = 0.05 baked in
    
    p_adj = p_sorted .* m ./ (1:m)';
    
    % enforce monotonicity
    p_adj = flipud(cummin(flipud(p_adj)));
    
    % cap at 1
    p_adj(p_adj > 1) = 1;
    
    % reinsert into original shape
    p_fdr = nan(size(p));
    p_fdr(~nan_mask) = p_adj(inv_perm(sort_idx));
    
    p_fdr = reshape(p_fdr, size(p));
end

function idx = inv_perm(p)
    idx = zeros(size(p));
    idx(p) = 1:length(p);
end
%% Quanrification of difference in acceleration during 50ms period aligned to stim 
BL_diff=[];
for a=1:3

    bl=BL_big.stimvel(:,:,a);
    % stimon= bl(26,8);
    % stimoff=bl(31,8);

    stimon= bl(26,5);
    stimoff=bl(31,5);
    diff=stimoff-stimon;diff=diff*10;
    BL_diff=[BL_diff;diff];

end

BL_diffCSM=[];
for a=1:3

    bl=BL_CSM_big.stimvel(:,:,a);
    % stimon= bl(21,8);
    % stimoff=bl(31,8);

    stimon= bl(21,5);
    stimoff=bl(31,5);
    diff=stimoff-stimon;diff=diff*10;
    BL_diffCSM=[BL_diffCSM;diff];

end
%%
AD_diff=[];
for a=1:3

    ad=AD_big.stimvel(:,:,a);
    %  stimon= ad(26,8);
    % stimoff=ad(31,8);
    stimon= ad(26,5);
    stimoff=ad(31,5);
    diff=stimoff-stimon; diff=diff*10;
    AD_diff=[AD_diff;diff];

end

AD_diffCSM=[];
for a=1:3

    ad=AD_CSM_big.stimvel(:,:,a);
    %  stimon= ad(21,8);
    % stimoff=ad(31,8);
    stimon= ad(21,5);
    stimoff=ad(31,5);
    diff=stimoff-stimon; diff=diff*10;
    AD_diffCSM=[AD_diffCSM;diff];

end

%%
WOE_diff=[];
for a=1:3

    woe=WOE_big.stimvel(:,:,a);
    %  stimon= woe(26,8);
    % stimoff=woe(31,8);
     stimon= woe(26,5);
    stimoff=woe(31,5);
    diff=stimoff-stimon; diff=diff*10;
    WOE_diff=[WOE_diff;diff];

end

WOE_diffCSM=[];
for a=1:3

    woe=WOE_CSM_big.stimvel(:,:,a);
    %  stimon= woe(21,8);
    % stimoff=woe(31,8);
     stimon= woe(21,5);
    stimoff=woe(31,5);
    diff=stimoff-stimon; diff=diff*10;
    WOE_diffCSM=[WOE_diffCSM;diff];

end
%%
AD_differences=AD_diff-BL_diff;

WOE_differences=WOE_diff-BL_diff;

AD_differencesCSM=AD_diffCSM-BL_diffCSM;

WOE_differencesCSM=WOE_diffCSM-BL_diffCSM;
%% Generate stimvel zscore
% Initialize arrays for p-values and Z-scores
parray0 = NaN(size(BL_big.stimvel, 1), 1);  % Store p-values
zarray0 =   NaN(size(BL_big.stimvel, 1), 1)   ;    % Store Z-scores (same shape as stimvel)
Zscores0= NaN(size(BL_big.stimvel, 1), 5);

for pnum = 1:length(parray0)
    % Extract trial data for C_big (control) and T0_big (experimental)
    % a4b = squeeze(WOE_CSM_big.stimvel(pnum, 5, :));  % Experimental data (T0_big)
    % b4b = squeeze(BL_CSM_big.stimvel(pnum, 5 , :));   % Control data (C_big), matching trial
    % 
    % a4b = squeeze(AD_big.stimvel(pnum, 8, :));  % Experimental data (T0_big)
    % b4b = squeeze(BL_big.stimvel(pnum, 8 , :));   % Control data (C_big), matching trial

    a4b = squeeze(WOE_big.stimvel(pnum, 5, :));  % Experimental data (T0_big)
    b4b = squeeze(BL_big.stimvel(pnum, 5 , :));   % Control data (C_big), matching trial

    % Calculate p-value using Wilcoxon signed-rank test (trial-wise)
    [parray0(pnum, 1), ~, pstat] = signrank(a4b, b4b, 'method', 'approximate');
    parray0(pnum, 2) = pstat.zval;  % Store Z-value from signrank

    % Calculate mean and standard deviation for the control trial (b4b)
    mu_b4b = mean(b4b);             % Mean of the current trial in C_big
    sigma_b4b = std(b4b);           % Standard deviation of the current trial in C_big

    % Calculate Z-scores for each time point in a4b relative to b4b stats
    zscores = (a4b - mu_b4b) ./ sigma_b4b;  % Point-wise Z-score
    Zscores(pnum,:)=zscores;
    mean_zscore=mean(zscores);
    zarray0(pnum,:)=mean_zscore;
end

%%
BLcolor=rgb('Gray');
CSMcolor=rgb('Orchid');
stimcolor=rgb('DodgerBlue');
%zscore plot 
% Adjust the figure to provide more vertical space
figure('Position', [100, 100, 600, 300]);  hold on


dim = 8;
StimvelPlot(BL_CSM_big, dim,BLcolor, 0);
StimvelPlot(WOE_CSM_big, dim,CSMcolor, 0);
StimvelPlot(AD_CSM_big,dim,stimcolor,0);
% 
% StimvelPlot(BL_big, dim,BLcolor, 1);
% StimvelPlot(WOE_big, dim,CSMcolor, 0);
% StimvelPlot(AD_big,dim,stimcolor,0);

% Plot reference line
yline(0, '--k');

% Axis limits and labels
ylim('auto')
xlim([-0.25 0.25])
%ylim([-20 20])
%xlabel('Time from MaxV (s)');
%xlabel('Time from Stim(s)');
xlabel('Time from CSM')

ylabel('Outward Acceleration');


% Adjust the main plot's position to create space below
ax1 = gca;
ax1.Position(2) = ax1.Position(2) + 0.1;
ax1.Position(4) = ax1.Position(4) - 0.1;

% Add Z-score color bar below
ax2 = axes('Position', [ax1.Position(1), ax1.Position(2) - 0.2, ax1.Position(3), 0.05]);

% Plot Z-scores as a heatmap
imagesc(zarray0', 'Parent', ax2);  % Transpose to match time x trial layout
colormap(ax2, brewermap([], 'PrGn'));  % Diverging colormap
caxis([-3, 3]);  % Adjust the Z-score color limits
colorbar(ax2, 'Position', [0.92, ax2.Position(2), 0.02, ax2.Position(4)]);

% Hide ticks and labels for the Z-score axis
set(ax2, 'XTick', [], 'YTick', [], 'YColor', 'none', 'XColor', 'none');

%% New method for plotting 
BLcolor=rgb('Gray');
CSMcolor=rgb('Orchid');
stimcolor=rgb('DodgerBlue');

%zscore plot 
figure('Position', [100, 100, 600, 300]);  hold on

dim = 5;
StimvelPlot(BL_CSM_big, dim,BLcolor, 0);
StimvelPlot(WOE_CSM_big, dim,CSMcolor, 0);
StimvelPlot(AD_CSM_big,dim,stimcolor,0);

yline(0, '--k');

ylim('auto')
xlim([-0.25 0.25])
xlabel('Time from Stim(s)');
ylabel('Outward Acceleration');

ax1 = gca;
ax1.Position(2) = ax1.Position(2) + 0.1;
ax1.Position(4) = ax1.Position(4) - 0.1;

% ======== NEW: compute and plot arrows for Z-score peaks ========
% assume zarray0 is 1 x timeBins
[~, maxIdx] = max(zarray0);
[~, minIdx] = min(zarray0);

% create time vector matching PSTH range
timevec = linspace(-0.25, 0.25, length(zarray0));

% arrow coordinates
xMax = timevec(maxIdx);
xMin = timevec(minIdx);

% arrow y positions
yTop = max(ylim) * 0.95;    % near top
yBot = min(ylim) * 0.95;    % near bottom

% plot arrows on main axis
annotation('textarrow', ...
    ax1.Position(1:2) + [ (xMax - ax1.XLim(1))/diff(ax1.XLim)*ax1.Position(3), 0 ], ...
    ax1.Position(2:2) + [ 0, (yTop - ax1.YLim(1))/diff(ax1.YLim)*ax1.Position(4)], ...
    'String','max Z', 'Color','k');

annotation('textarrow', ...
    ax1.Position(1:2) + [ (xMin - ax1.XLim(1))/diff(ax1.XLim)*ax1.Position(3), 0 ], ...
    ax1.Position(2:2) + [ 0, (yBot - ax1.YLim(1))/diff(ax1.YLim)*ax1.Position(4)], ...
    'String','min Z', 'Color','k');

% ================================================================


% Add Z-score color bar below
ax2 = axes('Position', [ax1.Position(1), ax1.Position(2) - 0.2, ax1.Position(3), 0.05]);

imagesc(zarray0', 'Parent', ax2);
colormap(ax2, brewermap([], 'PrGn'));
caxis([-3, 3]);
colorbar(ax2, 'Position', [0.92, ax2.Position(2), 0.02, ax2.Position(4)]);

set(ax2, 'XTick', [], 'YTick', [], 'YColor', 'none', 'XColor', 'none');


%%
CSMcolor = rgb('Orchid');

% zscore plot 
set(gcf, 'Position', [100, 100, 800, 600]); % [left, bottom, width, height]

dim = 5;

% % Plot both datasets on the same axes
% StimvelPlot(BL_CSM_big, dim,BLcolor, 1);
% hold on
% StimvelPlot(WO_CSM_big, dim,CSMcolor, 0);
% hold on
% StimvelPlot(AD_CSM_big,dim,stimcolor,0);
% hold off;
% Plot both datasets on the same axes
StimvelPlot(BL_big, dim,BLcolor, 1);
hold on
StimvelPlot(WOE_big, dim,CSMcolor, 0);
hold on
StimvelPlot(AD_big,dim,stimcolor,0);
hold off;

% Plot reference line
yline(0, '--k');
xlim([-0.25 0.25])
ylim('auto')
% Axis limits and labels
ylim('auto')
xlabel('Time from MaxV (s)');
ylabel('Outward Acceleration');

% Adjust the main plot's position to make it shorter
ax1 = gca;
ax1.Position(2) = ax1.Position(2) + 0.15;  % Raise the plot up
ax1.Position(4) = ax1.Position(4) * 0.6;   % Shrink the height (e.g., 60% of original)

% Add Z-score color bar below
ax2 = axes('Position', [ax1.Position(1), ax1.Position(2) - 0.18, ax1.Position(3), 0.05]);

% Plot Z-scores as a heatmap
imagesc(zarray0', 'Parent', ax2);  % Transpose to match time x trial layout
colormap(ax2, brewermap([], 'PrGn'));  % Diverging colormap
caxis([-2, 2]);  % Adjust the Z-score color limits
colorbar(ax2, 'Position', [0.92, ax2.Position(2), 0.02, ax2.Position(4)]);

% Hide ticks and labels for the Z-score axis
set(ax2, 'XTick', [], 'YTick', [], 'YColor', 'none', 'XColor', 'none');


%%
% Initialize
an = fieldnames(ExcitRR_BL);
stimveldiffX = [];
stimveldiffY = [];
maxStimVel = zeros(length(an),1);
maxStimVelIdx = zeros(length(an),1);
maxStimVelTime = zeros(length(an),1);  % Optional, if you have a time vector

% Assuming a time vector exists and is the same length as the rows of Average
% E.g., time_vector = -500:10:1000; (adjust accordingly)
% and it's stored in ChR2.(an{a}).S.time or given externally

for a = 1:length(an)
    % Get averages at 63rd row (110 ms)
    CavgX = ExcitRR_BL.(an{a}).BLN.Average(31,5);
    SavgX = ExcitRR_AD.(an{a}).ADP.Average(31,5);

    CavgY = ExcitRR_BL.(an{a}).BLN.Average(31,6);
    SavgY =ExcitRR_AD.(an{a}).ADP.Average(31,6);

    % Compute difference at 110 ms
    DiffX = SavgX - CavgX;
    DiffY = SavgY - CavgY;

    % Store
    stimveldiffX(a,1) = DiffX;
    stimveldiffY(a,1) = DiffY;

    % === Get max stim velocity and index between rows 56 to 72 ===
    stim_segment = ExcitRR_AD.(an{a}).ADP.Average(31:51,5);
    [maxStimVel(a), rel_idx] = max(stim_segment);
    maxStimVelIdx(a) = rel_idx + 30;  % Convert to absolute row index (56 + rel_idx - 1)

    % === If you have a time vector ===
    % Uncomment and adjust if you have time_vector
     % maxStimVelTime(a) = time_vector(maxStimVelIdx(a));
end

% Optional: Combine into one output matrix
ExcitAnimalsDiff_CSM = [stimveldiffX, stimveldiffY, maxStimVel, maxStimVelIdx];

%% Look at poisiton and velocity 
% Preallocate output array: 5 rows, 2 columns (for col 2 and col 5)
BL_values0 = zeros(3, 2);
BL_values50 = zeros(3, 2);
BL_values100 = zeros(3, 2);

zeroidx=26;
endidx=31;
csmidx=21;
csmidxbl=24;
%csmidx=38;
% csmidx=40;
%csmidx=34:46;

% Loop through each slice in 3rd dimension
for i = 1:3
    % Extract values at row 61, columns 2 and 5 for each slice
    % stim on
    BL_values0(i, 1) = BL_big.stimvel(zeroidx, 2, i);
    BL_values0(i, 2) = BL_big.stimvel(zeroidx, 5, i);

    BL_values50(i, 1) = BL_big.stimvel(endidx, 2, i);
    BL_values50(i, 2) = BL_big.stimvel(endidx, 5, i);
    
    % % 50 ms after
    % BL_values100(i, 1) = BL_big.stimvel(csmidx, 2, i);
    % BL_values100(i, 2) = BL_big.stimvel(csmidx, 5, i);
    BL_values100(i, 1) = BL_CSM_big.stimvel(csmidxbl, 2, i);
    BL_values100(i, 2) = BL_CSM_big.stimvel(csmidxbl, 5, i);

    % % CSM window: take max across range for col 2 and col 5
    % BL_values100(i, 1) = max(BL_big.stimvel(csmidx, 2, i));
    % BL_values100(i, 2) = max(BL_big.stimvel(csmidx, 5, i));
   
end


AD_values0 = zeros(3, 2);
AD_values50 = zeros(3, 2);
AD_values100 = zeros(3, 2);

% Loop through each slice in 3rd dimension
for i = 1:3
    % Extract values at row 61, columns 2 and 5 for each slice
    AD_values0(i, 1) = AD_big.stimvel(zeroidx, 2, i);
    AD_values0(i, 2) = AD_big.stimvel(zeroidx, 5, i);

    AD_values50(i, 1) = AD_big.stimvel(endidx, 2, i);
    AD_values50(i, 2) = AD_big.stimvel(endidx, 5, i);
    % 
    % ADE_values100(i, 1) = ADE_big.stimvel(csmidx, 2, i);
    % ADE_values100(i, 2) = ADE_big.stimvel(csmidx, 5, i);

    AD_values100(i, 1) = AD_CSM_big.stimvel(csmidx, 2, i);
    AD_values100(i, 2) = AD_CSM_big.stimvel(csmidx, 5, i);

        % CSM window: take max across range for col 2 and col 5
    % AD_values100(i, 1) = max(AD_big.stimvel(csmidx, 2, i));
    % AD_values100(i, 2) = max(AD_big.stimvel(csmidx, 5, i));
end

ADE_values0 = zeros(3, 2);
ADE_values50 = zeros(3, 2);
ADE_values100 = zeros(3, 2);

% Loop through each slice in 3rd dimension
for i = 1:3
    % Extract values at row 61, columns 2 and 5 for each slice
    ADE_values0(i, 1) = ADE_big.stimvel(zeroidx, 2, i);
    ADE_values0(i, 2) = ADE_big.stimvel(zeroidx, 5, i);

    ADE_values50(i, 1) = ADE_big.stimvel(endidx, 2, i);
    ADE_values50(i, 2) = ADE_big.stimvel(endidx, 5, i);
    % 
    % ADE_values100(i, 1) = ADE_big.stimvel(csmidx, 2, i);
    % ADE_values100(i, 2) = ADE_big.stimvel(csmidx, 5, i);

    ADE_values100(i, 1) = ADE_CSM_big.stimvel(csmidx, 2, i);
    ADE_values100(i, 2) = ADE_CSM_big.stimvel(csmidx, 5, i);

        % CSM window: take max across range for col 2 and col 5
    % ADE_values100(i, 1) = max(ADE_big.stimvel(csmidx, 2, i));
    % ADE_values100(i, 2) = max(ADE_big.stimvel(csmidx, 5, i));
end

ADL_values0 = zeros(3, 2);
ADL_values50 = zeros(3, 2);
ADL_values100 = zeros(3, 2);

% Loop through each slice in 3rd dimension
for i = 1:3
    % Extract values at row 61, columns 2 and 5 for each slice
    ADL_values0(i, 1) = ADL_big.stimvel(zeroidx, 2, i);
    ADL_values0(i, 2) = ADL_big.stimvel(zeroidx, 5, i);

    ADL_values50(i, 1) = ADL_big.stimvel(endidx, 2, i);
    ADL_values50(i, 2) = ADL_big.stimvel(endidx, 5, i);

    % ADL_values100(i, 1) = ADL_big.stimvel(csmidx, 2, i);
    % ADL_values100(i, 2) = ADL_big.stimvel(csmidx, 5, i);

     ADL_values100(i, 1) = ADL_CSM_big.stimvel(csmidx, 2, i);
    ADL_values100(i, 2) = ADL_CSM_big.stimvel(csmidx, 5, i);

            % CSM window: take max across range for col 2 and col 5
    % ADL_values100(i, 1) = max(ADL_big.stimvel(csmidx, 2, i));
    % ADL_values100(i, 2) = max(ADL_big.stimvel(csmidx, 5, i));
end

WOE_values0 = zeros(3, 2);
WOE_values50 = zeros(3, 2);
WOE_values100 = zeros(3, 2);

% Loop through each slice in 3rd dimension
for i = 1:3
    % Extract values at row 61, columns 2 and 5 for each slice
    WOE_values0(i, 1) = WOE_big.stimvel(zeroidx, 2, i);
    WOE_values0(i, 2) = WOE_big.stimvel(zeroidx, 5, i);

    WOE_values50(i, 1) = WOE_big.stimvel(endidx, 2, i);
    WOE_values50(i, 2) = WOE_big.stimvel(endidx, 5, i);

    % WOE_values100(i, 1) = WOE_big.stimvel(csmidx, 2, i);
    % WOE_values100(i, 2) = WOE_big.stimvel(csmidx, 5, i);

     WOE_values100(i, 1) = WOE_CSM_big.stimvel(csmidxbl, 2, i);
    WOE_values100(i, 2) = WOE_CSM_big.stimvel(csmidxbl, 5, i);

    %         % CSM window: take max across range for col 2 and col 5
    % WOE_values100(i, 1) = max(WOE_big.stimvel(csmidx, 2, i));
    % WOE_values100(i, 2) = max(WOE_big.stimvel(csmidx, 5, i));

end

WOL_values0 = zeros(3, 2);
WOL_values50 = zeros(3, 2);
WOL_values100 = zeros(3, 2);

% Loop through each slice in 3rd dimension
for i = 1:3
    % Extract values at row 61, columns 2 and 5 for each slice
    WOL_values0(i, 1) = WOL_big.stimvel(zeroidx, 2, i);
    WOL_values0(i, 2) = WOL_big.stimvel(zeroidx, 5, i);

    WOL_values50(i, 1) = WOL_big.stimvel(endidx, 2, i);
    WOL_values50(i, 2) = WOL_big.stimvel(endidx, 5, i);

    % WOL_values100(i, 1) = WOL_big.stimvel(csmidx, 2, i);
    % WOL_values100(i, 2) = WOL_big.stimvel(csmidx, 5, i);


    WOL_values100(i, 1) = WOL_CSM_big.stimvel(csmidx, 2, i);
    WOL_values100(i, 2) = WOL_CSM_big.stimvel(csmidx, 5, i);

              % CSM window: take max across range for col 2 and col 5
    % WOL_values100(i, 1) = max(WOL_big.stimvel(csmidx, 2, i));
    % WOL_values100(i, 2) = max(WOL_big.stimvel(csmidx, 5, i));

end
%%
ADStimdiff=AD_values50-BL_values50;
ADPostdiff=AD_values100-BL_values100;


ADEStimdiff=ADE_values50-BL_values50;
ADEPostdiff=ADE_values100-BL_values100;

ADLStimdiff=ADL_values50-BL_values50;
ADLPostdiff=ADL_values100-BL_values100;

WOEStimdiff=WOE_values50-BL_values50;
WOEPostdiff=WOE_values100-BL_values100;

WOLStimdiff=WOL_values50-BL_values50;
WOLPostdiff=WOL_values100-BL_values100;
%% More dynamic process for looking at peaks 

% Define ranges for each time window
range0 = 26;    % Example range for time window 0 ms
range50 = 27:33;   % Example range for ~50 ms
range100 = 34:40;  % Example range for ~100 ms

% Define groups and their data
groups = {'BL_big', 'ADE_big', 'ADL_big', 'WOE_big', 'WOL_big'};
results = struct();

for g = 1:length(groups)
    groupName = groups{g};
    data = eval(groupName);  % Load the 3D matrix for current group

    % Preallocate result arrays: 3 slices × 2 columns
    values0 = zeros(3, 2);
    values50 = zeros(3, 2);
    values100 = zeros(3, 2);

    % Loop through slices (3rd dimension)
    for i = 1:3
        slice = data.stimvel(:, :, i);  % Extract 2D slice

        % Find max in column 2 and 5 for each range
        values0(i, 1) = max(slice(range0, 2));
        values0(i, 2) = max(slice(range0, 5));

        values50(i, 1) = max(slice(range50, 2));
        values50(i, 2) = max(slice(range50, 5));

        values100(i, 1) = max(slice(range100, 2));
        values100(i, 2) = max(slice(range100, 5));
    end

    % Store results in a struct
    results.(groupName).values0 = values0;
    results.(groupName).values50 = values50;
    results.(groupName).values100 = values100;
end

%% 

BLStimdiff=BL_values50-BL_values0;
BLPostdiff=BL_values100-BL_values50;

ADEStimdiff=ADE_values50-ADE_values0;
ADEPostdiff=ADE_values100-ADE_values50;

ADLStimdiff=ADL_values50-ADL_values0;
ADLPostdiff=ADL_values100-ADL_values50;

WOEStimdiff=WOE_values50-WOE_values0;
WOEPostdiff=WOE_values100-WOE_values50;

WOLStimdiff=WOL_values50-WOL_values0;
WOLPostdiff=WOL_values100-WOL_values50;

%%
% Stim phase relative to BL_values0
ADEStim_rel2BL = (ADE_values50 - BL_values0) ./ BL_values0;
ADLStim_rel2BL = (ADL_values50 - BL_values0) ./ BL_values0;
WOEStim_rel2BL = (WOE_values50 - BL_values0) ./ BL_values0;
WOLStim_rel2BL = (WOL_values50 - BL_values0) ./ BL_values0;

% Post phase relative to BL_values50
ADEPost_rel2BL = (ADE_values100 - BL_values50) ./ BL_values50;
ADLPost_rel2BL = (ADL_values100 - BL_values50) ./ BL_values50;
WOEPost_rel2BL = (WOE_values100 - BL_values50) ./ BL_values50;
WOLPost_rel2BL = (WOL_values100 - BL_values50) ./ BL_values50;

%%
% Organize labels and colors
labels = {'ADE', 'ADL', 'WOE', 'WOL'};
colors = lines(4);

% For each velocity dimension (col 2 and 5)
for dim = 1:2
    figure('Name', sprintf('Relative to Baseline - Velocity Dim %d', dim));
    hold on;

    % Collect stim and post data for plotting
    stim_data = [ADEStim_rel2BL(:,dim), ADLStim_rel2BL(:,dim), ...
                 WOEStim_rel2BL(:,dim), WOLStim_rel2BL(:,dim)];

    post_data = [ADEPost_rel2BL(:,dim), ADLPost_rel2BL(:,dim), ...
                 WOEPost_rel2BL(:,dim), WOLPost_rel2BL(:,dim)];

    % Plot stim data
    subplot(1,2,1);
    boxplot(stim_data, 'Labels', labels);
    title(sprintf('Stim Relative to BL (Col %d)', dim+1));
    ylabel('Fractional Change from BL');
    hold on;
    for i = 1:4
        scatter(repmat(i,5,1), stim_data(:,i), 60, ...
                'filled', 'MarkerFaceColor', colors(i,:), ...
                'MarkerEdgeColor', 'k');
    end
    yline(0, '--k');
    ylim padded
    grid on

    % Plot post data
    subplot(1,2,2);
    boxplot(post_data, 'Labels', labels);
    title(sprintf('Post Relative to BL (Col %d)', dim+1));
    ylabel('Fractional Change from BL');
    hold on;
    for i = 1:4
        scatter(repmat(i,5,1), post_data(:,i), 60, ...
                'filled', 'MarkerFaceColor', colors(i,:), ...
                'MarkerEdgeColor', 'k');
    end
    yline(0, '--k');
    ylim padded
    grid on
end

%% %% ============================================================
%  CSM-ALIGNED WINDOW SUMMARY: -50 to -30 ms (idx 21:23)
%  - Extract per-animal mean position/velocity (cols 2 and 5)
%  - Plot grouped bar plots (mean ± SEM) with individual animal points
%  - Run significance tests vs BL (paired t-test + Wilcoxon as backup)
%
%  Your mapping:
%   -50 ms -> idx 21
%   -40 ms -> idx 22
%   -30 ms -> idx 23
% ============================================================

% ---------------- USER-DEFINED INDICES ----------------
idxWin   = 21:23;   % -50:-30 ms window (CSM-aligned)
col_pos  = 2;       % position column
col_vel  = 5;       % velocity column

% If some conditions use different CSM-aligned structures/indices, keep these:
% (You had csmidx=21; csmidxbl=24 previously; NOW we use idxWin directly.)
% Make sure these exist and are aligned to CSM:
%   *_CSM_big.stimvel(timeIdx, col, animalIdx)

% ---------------- HELPER: extract window mean per animal ----------------
% Returns [nAnimals x 2] matrix: [:,1]=pos, [:,2]=vel
getWinMean = @(S) squeeze( mean( S.stimvel(idxWin, [col_pos col_vel], :), 1 ) )';

% ---------------- EXTRACT PER-ANIMAL WINDOW MEANS ----------------
% Each is nAnimals x 2 (pos, vel). Your code assumes 3 animals.
BL_win  = getWinMean(BL_CSM_big);

AD_win  = getWinMean(AD_CSM_big);
ADE_win = getWinMean(ADE_CSM_big);
ADL_win = getWinMean(ADL_CSM_big);

WOE_win = getWinMean(WOE_CSM_big);
WOL_win = getWinMean(WOL_CSM_big);

% ---------------- PICK METRIC TO PLOT/TEST ----------------
% Choose which column to use:
metric_name = 'Velocity';   % 'Position' or 'Velocity'
if strcmpi(metric_name,'Position')
    mcol = 1; ylabel_txt = 'Position (col 2 units)';
else
    mcol = 2; ylabel_txt = 'Velocity (col 5 units)';
end

% Pack into matrix for plotting/testing
conds   = {'BL','AD','ADE','ADL','WOE','WOL'};
WinAll  = [BL_win(:,mcol), AD_win(:,mcol), ADE_win(:,mcol), ADL_win(:,mcol), WOE_win(:,mcol), WOL_win(:,mcol)];
nA      = size(WinAll,1);

% ---------------- BAR PLOT (MEAN ± SEM) + INDIVIDUAL POINTS ----------------
mu  = mean(WinAll,1,'omitnan');
sem = std(WinAll,0,1,'omitnan') ./ sqrt(sum(~isnan(WinAll),1));

figure('Color','w','Position',[200 200 900 350]); hold on
b = bar(1:numel(conds), mu, 0.7); %#ok<NASGU>
errorbar(1:numel(conds), mu, sem, 'k', 'LineStyle','none', 'LineWidth', 1.5);

% overlay individual animal points
for ci = 1:numel(conds)
    xj = ci + (rand(nA,1)-0.5)*0.12;
    scatter(xj, WinAll(:,ci), 60, 'k', 'filled', 'MarkerFaceAlpha', 0.65);
end

xticks(1:numel(conds)); xticklabels(conds);
ylabel(ylabel_txt);
title(sprintf('%s mean in CSM window (-50 to -30 ms; idx 21:23)', metric_name));
yline(0,'--'); box off

% ---------------- SIGNIFICANCE TESTING (vs BL) ----------------
% Paired tests are appropriate because each row is the same animal across conditions.
% We report both:
%   - paired t-test (parametric)
%   - Wilcoxon signed-rank (nonparametric)
alpha = 0.05;

p_t   = nan(1,numel(conds));
p_wsr = nan(1,numel(conds));

for ci = 2:numel(conds)   % compare each to BL
    xBL = WinAll(:,1);
    x   = WinAll(:,ci);

    % keep only animals that exist in BOTH
    ok = ~isnan(xBL) & ~isnan(x);
    if sum(ok) < 2
        fprintf('%s vs BL: not enough paired samples.\n', conds{ci});
        continue
    end

    [~, p_t(ci)] = ttest(x(ok), xBL(ok));              % paired t-test
    p_wsr(ci)    = signrank(x(ok), xBL(ok));           % Wilcoxon signed-rank
end

% FDR correction across the 5 comparisons (AD, ADE, ADL, WOE, WOL)
pvec = p_t(2:end);
[~, ~, p_t_fdr] = fdr_bh(pvec, alpha, 'pdep', 'yes');  % requires fdr_bh on path
% If you don't have fdr_bh, comment out above and use raw p-values.

% Print results
fprintf('\n=== Window stats: %s (-50 to -30 ms; idx 21:23) ===\n', metric_name);
for ci = 2:numel(conds)
    fprintf('%s vs BL: p_t=%.4g, p_wilcoxon=%.4g\n', conds{ci}, p_t(ci), p_wsr(ci));
end

% ---------------- OPTIONAL: annotate significant bars (using paired t-test raw) ----------------
% stars: * <0.05, ** <0.01, *** <0.001
yl = ylim;
yStar = yl(2) - 0.06*(yl(2)-yl(1));

for ci = 2:numel(conds)
    if isnan(p_t(ci)), continue; end
    if p_t(ci) < 0.001
        s='***';
    elseif p_t(ci) < 0.01
        s='**';
    elseif p_t(ci) < 0.05
        s='*';
    else
        s='';
    end
    if ~isempty(s)
        text(ci, yStar, s, 'HorizontalAlignment','center', 'FontSize', 14, 'FontWeight','bold');
    end
end
ylim(yl);


%% ============================================================
% SUMMARY (3 animals) for BL, ADE, ADL, WOE
% Metric: MAX within CSM window (-50 to -30 ms; idx 21:23)
% Source: *_CSM_big.stimvel(timeIdx, col, animalIdx)
% Plots: (1) bars of raw summary (mean ± SEM + animal points)
%        (2) bars of Δ vs BL (mean ± SEM + animal points)
% Stats : paired t-test + Wilcoxon signed-rank vs BL


% ---------------- SETTINGS ----------------
idxWin  = 21:23;   % -50:-30 ms relative to CSM
col_vel = 5;       % velocity column (set to 2 for position)

conds   = {'BL','ADE','ADL','WOE'};

% ---------------- EXTRACT: max in window per animal ----------------
% Robust extractor that returns [nAnimals x 1]
getWinMax = @(S) reshape( max(S.stimvel(idxWin, col_vel, :), [], 1), [], 1 );

BL_vals  = getWinMax(BL_CSM_big);
ADE_vals = getWinMax(ADE_CSM_big);
ADL_vals = getWinMax(ADL_CSM_big);
WOE_vals = getWinMax(WOE_CSM_big);

SummaryMat = [BL_vals ADE_vals ADL_vals WOE_vals];   % [nAnimals x 4]
nA = size(SummaryMat,1);

fprintf('\nSummaryMat (rows=animals, cols=BL ADE ADL WOE):\n');
disp(SummaryMat);
fprintf('size(SummaryMat) = [%d %d]\n\n', size(SummaryMat,1), size(SummaryMat,2));

%% ---------------- PLOT 1: RAW SUMMARY BARS ----------------
mu  = mean(SummaryMat,1,'omitnan');
sem = std(SummaryMat,0,1,'omitnan') ./ sqrt(sum(~isnan(SummaryMat),1));

figure('Color','w','Position',[250 250 750 320]); hold on
bar(1:4, mu, 0.65);
errorbar(1:4, mu, sem, 'k', 'LineStyle','none', 'LineWidth', 1.5);

% overlay animal points
for ci = 1:4
    xj = ci + (rand(nA,1)-0.5)*0.10;
    scatter(xj, SummaryMat(:,ci), 60, 'k', 'filled', 'MarkerFaceAlpha', 0.7);
end

xticks(1:4); xticklabels(conds);
ylabel('Peak value (col 5)');
title('Max within CSM window (-50 to -30 ms; idx 21:23)');
yline(0,'--'); box off;

% ---------------- STATS (paired vs BL) on RAW VALUES ----------------
alpha = 0.05;
p_t   = nan(1,4);
p_w   = nan(1,4);

for ci = 2:4
    xBL = SummaryMat(:,1);
    x   = SummaryMat(:,ci);
    ok  = ~isnan(xBL) & ~isnan(x);
    if sum(ok) < 2
        fprintf('%s vs BL: not enough paired samples.\n', conds{ci});
        continue
    end
    [~, p_t(ci)] = ttest(x(ok), xBL(ok));      % paired t-test
    p_w(ci)      = signrank(x(ok), xBL(ok));   % Wilcoxon signed-rank
end

fprintf('\nRAW (max in window) paired tests vs BL:\n');
for ci = 2:4
    fprintf('%s vs BL: ttest p = %.4g | signrank p = %.4g\n', conds{ci}, p_t(ci), p_w(ci));
end

%% ---------------- PLOT 2: DELTA VS BL BARS (recommended) ----------------
DeltaMat = SummaryMat(:,2:4) - SummaryMat(:,1);   % [nAnimals x 3] for ADE/ADL/WOE
condsD   = {'ADE-BL','ADL-BL','WOE-BL'};

muD  = mean(DeltaMat,1,'omitnan');
semD = std(DeltaMat,0,1,'omitnan') ./ sqrt(sum(~isnan(DeltaMat),1));

figure('Color','w','Position',[250 620 650 320]); hold on
bar(1:3, muD, 0.65);
errorbar(1:3, muD, semD, 'k', 'LineStyle','none', 'LineWidth', 1.5);

for ci = 1:3
    xj = ci + (rand(nA,1)-0.5)*0.10;
    scatter(xj, DeltaMat(:,ci), 60, 'k', 'filled', 'MarkerFaceAlpha', 0.7);
end

xticks(1:3); xticklabels(condsD);
ylabel('Δ Peak value (vs BL)');
title('Δ Max within CSM window (-50 to -30 ms; idx 21:23)');
yline(0,'--'); box off;

%% ---------------- STATS (paired vs BL) on DELTAS ----------------
% Here we test whether Δ differs from 0 (equivalent to paired test vs BL)
p_tD = nan(1,3);
p_wD = nan(1,3);

for ci = 1:3
    x = DeltaMat(:,ci);
    ok = ~isnan(x);
    if sum(ok) < 2
        fprintf('%s: not enough samples.\n', condsD{ci});
        continue
    end
    [~, p_tD(ci)] = ttest(x(ok));      % test vs 0
    p_wD(ci)      = signrank(x(ok));   % test vs 0
end

fprintf('\nDELTA (vs BL) tests vs 0:\n');
for ci = 1:3
    fprintf('%s: ttest p = %.4g | signrank p = %.4g\n', condsD{ci}, p_tD(ci), p_wD(ci));
end

% ---------------- OPTIONAL: annotate stars on DELTA plot (ttest) ----------------
figure(gcf); % ensure we're on delta figure
yl = ylim;
yStar = yl(2) - 0.07*(yl(2)-yl(1));

for ci = 1:3
    if isnan(p_tD(ci)), continue; end
    if p_tD(ci) < 0.001
        s='***';
    elseif p_tD(ci) < 0.01
        s='**';
    elseif p_tD(ci) < 0.05
        s='*';
    else
        s='';
    end
    if ~isempty(s)
        text(ci, yStar, s, 'HorizontalAlignment','center', 'FontSize', 14, 'FontWeight','bold');
    end
end
ylim(yl);


%% ============================================================
% DECELERATION GAIN: slope of velocity vs time
% Windows: -100->0 ms and -50->0 ms (CSM-aligned)
% Conditions: BL, ADE, ADL, WOE
%
% OUTPUTS:
%   slopes.(cond).w100_0  : [nAnimals x 1] slope (vel units / ms)
%   slopes.(cond).w50_0   : [nAnimals x 1] slope (vel units / ms)
%
% PLOTS:
%   (1) raw slopes by condition (two figures: one per window)
%   (2) optional "decel gain" = -slope (positive = stronger decel)
%
% STATS:
%   paired vs BL for each window (ttest + signrank)


% ---------------- SETTINGS ----------------
col_vel = 5;    % velocity column (change if needed)
zeroidx = 26;   % 0 ms index (relative to CSM), from your previous code

dt_ms = 10;     % IMPORTANT: set to your bin width in ms (10ms if true)

% Window endpoints (indices)
idx_0      = zeroidx;
idx_m50    = 21;
%idx_m50    = 22; % 40
idx_m100   =16;
%idx_m100   =19;

% Sanity check with your known mapping (you said -50ms is idx 21)
% If dt_ms=10 and zeroidx=26 => -50ms = 21 (matches)
fprintf('Computed idx for -100ms: %d | -50ms: %d | 0ms: %d\n', idx_m100, idx_m50, idx_0);

conds = {'BL','ADE','ADL','WOE'};

% ---------------- Helper: slope via linear fit ----------------
% S.stimvel is [time x dim x animal]
% returns slopes [nAnimals x 1] in (velocity units / ms)
getSlope = @(S, idxStart, idxEnd) localSlopePerAnimal(S.stimvel(:,col_vel,:), idxStart, idxEnd, zeroidx, dt_ms);

%% ---------------- Compute slopes for each condition ----------------
slopes = struct();

slopes.BL.w100_0  = getSlope(BL_CSM_big,  idx_m100, idx_0);
slopes.BL.w50_0   = getSlope(BL_CSM_big,  idx_m50,  idx_0);

slopes.ADE.w100_0 = getSlope(ADE_CSM_big, idx_m100, idx_0);
slopes.ADE.w50_0  = getSlope(ADE_CSM_big, idx_m50,  idx_0);

slopes.ADL.w100_0 = getSlope(ADL_CSM_big, idx_m100, idx_0);
slopes.ADL.w50_0  = getSlope(ADL_CSM_big, idx_m50,  idx_0);

slopes.WOE.w100_0 = getSlope(WOE_CSM_big, idx_m100, idx_0);
slopes.WOE.w50_0  = getSlope(WOE_CSM_big, idx_m50,  idx_0);

% Pack matrices [nAnimals x nConds]
S100 = [slopes.BL.w100_0, slopes.ADE.w100_0, slopes.ADL.w100_0, slopes.WOE.w100_0];
S50  = [slopes.BL.w50_0,  slopes.ADE.w50_0,  slopes.ADL.w50_0,  slopes.WOE.w50_0];

nA = size(S100,1);
fprintf('\nSlope matrices (rows=animals, cols=BL ADE ADL WOE):\n');
disp('S100 (−100→0):'); disp(S100);
disp('S50  (−50→0):');  disp(S50);

% ---------------- Optional: define "Decel gain" as -slope ----------------
% If velocity is decreasing toward 0, slope will be negative.
% DecelGain = -slope makes "more decel" positive.
useDecelGain = true;
if useDecelGain
    S100_plot = -S100;
    S50_plot  = -S50;
    ylab = 'Decel gain (−slope; vel units / ms)';
else
    S100_plot = S100;
    S50_plot  = S50;
    ylab = 'Slope (vel units / ms)';
end

% ---------------- Plot function (mean ± SEM + animal points) ----------------
plotBars(S100_plot, conds, ylab, 'Deceleration gain from −100 to 0 ms');
plotBars(S50_plot,  conds, ylab, 'Deceleration gain from −50 to 0 ms');

% ---------------- Stats: paired vs BL for each window ----------------
% Compare each condition to BL (paired across animals)
fprintf('\n=== Paired tests vs BL (window −100→0 ms) ===\n');
pairedVsBL(S100_plot, conds);

fprintf('\n=== Paired tests vs BL (window −50→0 ms) ===\n');
pairedVsBL(S50_plot, conds);

%% ============================================================
% Local functions
% ============================================================

function slopes = localSlopePerAnimal(stimvel_2d3d, idxStart, idxEnd, zeroidx, dt_ms)
    % stimvel_2d3d: [time x 1 x animal] OR [time x animal] depending on squeeze
    V = stimvel_2d3d;
    if ndims(V) == 3
        V = squeeze(V); % -> [time x animal]
    end

    nA = size(V,2);
    idx = idxStart:idxEnd;

    % time vector in ms relative to 0 at zeroidx
    t_ms = (idx - zeroidx) * dt_ms;   % column vector
    t_ms = t_ms(:);

    slopes = nan(nA,1);
    for a = 1:nA
        y = V(idx,a);
        ok = ~isnan(y(:)) & ~isnan(t_ms);
        if sum(ok) < 2
            slopes(a) = NaN;
            continue
        end

        % linear fit: y = m*t + b; slope m in vel units / ms
        p = polyfit(t_ms(ok), y(ok), 1);
        slopes(a) = p(1);
    end
end

function plotBars(M, conds, ylab, ttl)
    mu  = mean(M,1,'omitnan');
    sem = std(M,0,1,'omitnan') ./ sqrt(sum(~isnan(M),1));
    nA = size(M,1);

    figure('Color','w','Position',[250 250 720 320]); hold on
    bar(1:numel(conds), mu, 0.65);
    errorbar(1:numel(conds), mu, sem, 'k', 'LineStyle','none', 'LineWidth', 1.5);

    for ci = 1:numel(conds)
        xj = ci + (rand(nA,1)-0.5)*0.10;
        scatter(xj, M(:,ci), 60, 'k', 'filled', 'MarkerFaceAlpha', 0.7);
    end

    xticks(1:numel(conds)); xticklabels(conds);
    ylabel(ylab);
    title(ttl);
    yline(0,'--');
    box off
end

function pairedVsBL(M, conds)
    % M: [nAnimals x nConds] with col1 = BL
    for ci = 2:size(M,2)
        xBL = M(:,1);
        x   = M(:,ci);
        ok  = ~isnan(xBL) & ~isnan(x);

        if sum(ok) < 2
            fprintf('%s vs BL: not enough paired samples.\n', conds{ci});
            continue
        end

        [~, p_t] = ttest(x(ok), xBL(ok));      % paired t-test
        p_w      = signrank(x(ok), xBL(ok));   % Wilcoxon signed-rank
        fprintf('%s vs BL: ttest p=%.4g | signrank p=%.4g\n', conds{ci}, p_t, p_w);
    end
end

%% DIFFERENT DECEL GAIN METHOD
%% ============================================================
% DECELERATION GAIN: slope of velocity curve pre-CSM
% Windows: [-100, 0] ms and [-50, 0] ms
% Conditions: BL, ADE, ADL, WOE
% Metric: slope from linear fit (vel vs time)
%   - more negative slope = stronger deceleration (higher gain)
%
% Requires: *_CSM_big.stimvel(timeIdx, dim, animalIdx)


% ---------------- SETTINGS ----------------
dimVel   = 5;   % velocity column
zeroidx  = 26;  % 0 ms index in your aligned structure
idx_m50  = 21;  % -50 ms index
idx_m100 = 16;  % -100 ms index (assumes 10 ms/bin given your mapping)

% Build index ranges
idxWin_100_0 = idx_m100:zeroidx;
idxWin_50_0  = idx_m50:zeroidx;

% Time vectors (ms) for regression (only relative spacing matters)
t_100_0 = (-numel(idxWin_100_0)+1:0)' * 10;  % e.g., -100:-90:...:0
t_50_0  = (-numel(idxWin_50_0 )+1:0)' * 10;  % e.g., -50:-40:...:0

conds   = {'BL','ADE','ADL','WOE'};

% ---------------- GRAB PER-ANIMAL VELOCITY TRACES ----------------
% Each returns [nTime x nAnimals]
getVel = @(S) squeeze(S.stimvel(:, dimVel, :));

V_BL  = getVel(BL_CSM_big);
V_ADE = getVel(ADE_CSM_big);
V_ADL = getVel(ADL_CSM_big);
V_WOE = getVel(WOE_CSM_big);

% Stack for looping: [nTime x nAnimals x nConds]
V = cat(3, V_BL, V_ADE, V_ADL, V_WOE);

nA = size(V,2);
nC = size(V,3);

%% ---------------- COMPUTE SLOPES (per animal, per condition) ----------------
% slope units: (velocity units) / ms
slope_100_0 = nan(nA, nC);
slope_50_0  = nan(nA, nC);

for c = 1:nC
    for a = 1:nA
        y1 = V(idxWin_100_0, a, c);
        y2 = V(idxWin_50_0,  a, c);

        ok1 = ~isnan(y1);
        ok2 = ~isnan(y2);

        if sum(ok1) >= 2
            p1 = polyfit(t_100_0(ok1), y1(ok1), 1);
            slope_100_0(a,c) = p1(1);
        end
        if sum(ok2) >= 2
            p2 = polyfit(t_50_0(ok2), y2(ok2), 1);
            slope_50_0(a,c) = p2(1);
        end
    end
end

% “Gain” version (positive number = stronger decel), optional:
gain_100_0 = -slope_100_0;
gain_50_0  = -slope_50_0;

%% ---------------- PLOT: grouped bars + animal points ----------------
plotMetric = 'gain';  % 'slope' or 'gain'
switch lower(plotMetric)
    case 'slope'
        M1 = slope_100_0;  M2 = slope_50_0;
        ylab = 'Slope (vel / ms)  (more negative = more decel)';
    otherwise
        M1 = gain_100_0;   M2 = gain_50_0;
        ylab = 'Deceleration gain = -slope (vel / ms)  (bigger = more decel)';
end

mu1  = mean(M1,1,'omitnan');  sem1 = std(M1,0,1,'omitnan') ./ sqrt(sum(~isnan(M1),1));
mu2  = mean(M2,1,'omitnan');  sem2 = std(M2,0,1,'omitnan') ./ sqrt(sum(~isnan(M2),1));

figure('Color','w','Position',[200 200 900 350]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% --- Panel A: -100 to 0
nexttile; hold on
bar(1:nC, mu1, 0.65);
errorbar(1:nC, mu1, sem1, 'k','LineStyle','none','LineWidth',1.5);
for c = 1:nC
    xj = c + (rand(nA,1)-0.5)*0.10;
    scatter(xj, M1(:,c), 60, 'k','filled','MarkerFaceAlpha',0.7);
end
xticks(1:nC); xticklabels(conds);
ylabel(ylab);
title('Pre-CSM decel: -100 to 0 ms');
yline(0,'--'); box off

% --- Panel B: -50 to 0
nexttile; hold on
bar(1:nC, mu2, 0.65);
errorbar(1:nC, mu2, sem2, 'k','LineStyle','none','LineWidth',1.5);
for c = 1:nC
    xj = c + (rand(nA,1)-0.5)*0.10;
    scatter(xj, M2(:,c), 60, 'k','filled','MarkerFaceAlpha',0.7);
end
xticks(1:nC); xticklabels(conds);
title('Pre-CSM decel: -50 to 0 ms');
yline(0,'--'); box off

% ---------------- STATS: paired vs BL (per window) ----------------
% Pairing is by animal (row)
alpha = 0.05;

fprintf('\n=== Paired tests vs BL (%s) ===\n', plotMetric);

% helper to run paired tests
run_paired = @(X, name) deal( ...
    arrayfun(@(c) ttest_paired(X(:,1), X(:,c)), 2:nC), ...
    arrayfun(@(c) signrank_paired(X(:,1), X(:,c)), 2:nC), ...
    name);

% local inline helpers (keep at bottom of script if you prefer)
function p = ttest_paired(xBL, x)
    ok = ~isnan(xBL) & ~isnan(x);
    if sum(ok) < 2, p = nan; return; end
    [~,p] = ttest(x(ok), xBL(ok));
end
function p = signrank_paired(xBL, x)
    ok = ~isnan(xBL) & ~isnan(x);
    if sum(ok) < 2, p = nan; return; end
    p = signrank(x(ok), xBL(ok));
end

% window -100:0
p_t_100 = nan(1,nC); p_w_100 = nan(1,nC);
for c = 2:nC
    p_t_100(c) = ttest_paired(M1(:,1), M1(:,c));
    p_w_100(c) = signrank_paired(M1(:,1), M1(:,c));
end

% window -50:0
p_t_50 = nan(1,nC); p_w_50 = nan(1,nC);
for c = 2:nC
    p_t_50(c) = ttest_paired(M2(:,1), M2(:,c));
    p_w_50(c) = signrank_paired(M2(:,1), M2(:,c));
end

fprintf('\nWindow -100:0 ms\n');
for c = 2:nC
    fprintf('%s vs BL: ttest p=%.4g | signrank p=%.4g\n', conds{c}, p_t_100(c), p_w_100(c));
end

fprintf('\nWindow -50:0 ms\n');
for c = 2:nC
    fprintf('%s vs BL: ttest p=%.4g | signrank p=%.4g\n', conds{c}, p_t_50(c), p_w_50(c));
end


%% ============================================================
% DECEL GAIN across trial types using pooled reach structs
% Uses: ExcitRR_*.(animal).RRT2 or ReReaches (your choice)
% ============================================================



an = fieldnames(ExcitRR_BL);

% parameters (match your originals)
thresh    = 0.1;
minpoints = 5;
cov_pct   = 0.8;

% Choose which reach field to use per condition:
% For washout you said you want RRT2. For ADE/ADL, set whichever is correct.
USE_FIELD_BL  = 'RRT2';
USE_FIELD_ADE = 'ReReaches';  % change to 'RRT2' if you have it & want apples-to-apples
USE_FIELD_ADL = 'ReReaches';  % change to 'RRT2' if you have it
USE_FIELD_WO  = 'RRT2';
USE_FIELD_WOE = 'RRT2';
USE_FIELD_WOL = 'RRT2';

conds = {'BL','ADE','ADL','WO','WOE','WOL'};
nA = numel(an);

Gain   = nan(nA, numel(conds));
Rval   = nan(nA, numel(conds));
MedRG  = nan(nA, numel(conds));
Nreach = nan(nA, numel(conds));

DGout = cell(nA, numel(conds));

for a = 1:nA

    % -------- pull reaches using your structure --------
    BLR  = ExcitRR_BL.(an{a}).(USE_FIELD_BL);

    ADER = ExcitRR_AD_E.(an{a}).(USE_FIELD_ADE);
    ADLR = ExcitRR_AD_L.(an{a}).(USE_FIELD_ADL);

    WOR  = ExcitRR_WO.(an{a}).(USE_FIELD_WO);
    WOER = ExcitRR_WO_E.(an{a}).(USE_FIELD_WOE);
    WOLR = ExcitRR_WO_L.(an{a}).(USE_FIELD_WOL);

    ReachSets = {BLR, ADER, ADLR, WOR, WOER, WOLR};

    for ci = 1:numel(conds)

        DG = decel_gain_pooled_reaches(ReachSets{ci}, thresh, minpoints, cov_pct);

        DGout{a,ci} = DG;

        Gain(a,ci)   = DG.gain;
        Rval(a,ci)   = DG.R;
        Nreach(a,ci) = DG.nReaches;

        if ~isempty(DG.reach_gains)
            MedRG(a,ci) = median(DG.reach_gains(:,1), 'omitnan'); % x component
        end
    end
end

disp('Decel gain (slope) rows=animals cols=BL ADE ADL WO WOE WOL'); disp(Gain)
disp('N reaches rows=animals cols=BL ADE ADL WO WOE WOL'); disp(Nreach)

%% ---------------- PLOT: GAIN (mean ± SEM + animal points) ----------------
mu  = mean(Gain,1,'omitnan');
sem = std(Gain,0,1,'omitnan') ./ sqrt(sum(~isnan(Gain),1));

figure('Color','w','Position',[300 300 950 340]); hold on
bar(1:numel(conds), mu, 0.65);
errorbar(1:numel(conds), mu, sem, 'k','LineStyle','none','LineWidth',1.5);

for ci = 1:numel(conds)
    xj = ci + (rand(nA,1)-0.5)*0.10;
    scatter(xj, Gain(:,ci), 70, 'k','filled','MarkerFaceAlpha',0.7);
end

xticks(1:numel(conds)); xticklabels(conds);
ylabel('Decel gain (slope of decel vs vmax, x component)');
title('Decel gain by trial type (pooled reaches)');
yline(0,'--'); box off

% ---------------- STATS vs BL (paired across animals) ----------------
BL = Gain(:,1);

fprintf('\n=== Paired tests vs BL (across animals) ===\n');
for ci = 2:numel(conds)
    X = Gain(:,ci);
    ok = ~isnan(BL) & ~isnan(X);
    if sum(ok) >= 2
        [~,pT] = ttest(X(ok), BL(ok));
        pSR = signrank(X(ok), BL(ok));
        fprintf('%s vs BL: ttest p=%.4g | signrank p=%.4g\n', conds{ci}, pT, pSR);
    else
        fprintf('%s vs BL: not enough paired samples\n', conds{ci});
    end
end

%% ---------------- OPTIONAL: within-animal reach-level tests (like your old script) ----------------
% This tests whether reach_gains distributions differ within each animal.
% (Not independent across reaches, but useful descriptively.)
fprintf('\n=== Within-animal ranksum on reach_gains(:,1) vs BL ===\n');
for a = 1:nA
    BLrg = DGout{a,1}.reach_gains(:,1);
    for ci = 2:numel(conds)
        rg = DGout{a,ci}.reach_gains(:,1);
        if isempty(BLrg) || isempty(rg), continue; end
        p = ranksum(BLrg, rg);
        fprintf('Animal %s: %s vs BL ranksum p=%.4g\n', an{a}, conds{ci}, p);
    end
end


%% ============================================================
% PLOTTING SUITE for Decel Gain analyses (BL, ADE, ADL, WO, WOE, WOL)
%
% ASSUMES you already ran your decel-gain loop and have:
%   - an (cellstr)              : animal names (fieldnames)
%   - conds (cellstr)           : {'BL','ADE','ADL','WO','WOE','WOL'}
%   - DGout{a,ci} (struct) with fields:
%       .vmax        (nReaches x 1)   [x component]
%       .avgDecel    (nReaches x 1)
%       .gain        (scalar slope)
%       .lm          (fitlm object) (optional but nice)
%       .reach_gains (nReaches x 1) OR (nReaches x 3; we use col1)
%   - Gain (nAnimals x nConds)  : slopes (optional; can be derived)
%
% This script makes:
%   (1) Per-animal scatter: avgDecel vs vmax + regression line (small multiples)
%   (2) Spaghetti plot: slope (gain) per animal across conditions
%   (3) Spaghetti plot: intercept per animal across conditions (if lm exists)
%   (4) ECDF of reach_gains pooled across animals (one panel)
%   (5) Violin-like (box + swarm) of reach_gains pooled, plus per-animal medians
%   (6) Predicted decel at fixed vmax* (standardized braking) spaghetti
%   (7) Residual decel relative to BL model (per-animal distribution + summary)
% ============================================================



%---------------- SANITY CHECKS ----------------
assert(exist('DGout','var')==1, 'DGout not found. Run decel-gain loop first.');
assert(exist('conds','var')==1, 'conds not found.');
assert(exist('an','var')==1, 'an not found.');

nA = size(DGout,1);
nC = size(DGout,2);
assert(nC == numel(conds), 'DGout columns must match conds.');

% Build Gain/Intercept arrays if not present
GainPlot = nan(nA,nC);
IntPlot  = nan(nA,nC);
RPlot    = nan(nA,nC);
NPlot    = nan(nA,nC);

for a = 1:nA
    for ci = 1:nC
        if isempty(DGout{a,ci}), continue; end
        if isfield(DGout{a,ci},'gain'), GainPlot(a,ci) = DGout{a,ci}.gain; end
        if isfield(DGout{a,ci},'R'),    RPlot(a,ci)    = DGout{a,ci}.R;    end
        if isfield(DGout{a,ci},'nReaches'), NPlot(a,ci)= DGout{a,ci}.nReaches; end

        if isfield(DGout{a,ci},'lm') && ~isempty(DGout{a,ci}.lm)
            lm = DGout{a,ci}.lm;
            try
                IntPlot(a,ci) = lm.Coefficients.Estimate(1);
            catch
                IntPlot(a,ci) = NaN;
            end
        end
    end
end

%% ============================================================
% (1) PER-ANIMAL SCATTER: avgDecel vs vmax + regression line
% ============================================================
figure('Color','w','Position',[100 80 1400 750]);
tiledlayout(nA, nC, 'TileSpacing','compact','Padding','compact');

for a = 1:nA
    for ci = 1:nC
        nexttile; hold on

        DG = DGout{a,ci};
        if isempty(DG) || ~isfield(DG,'vmax') || isempty(DG.vmax)
            title(sprintf('%s %s (empty)', an{a}, conds{ci}), 'Interpreter','none');
            axis off; continue
        end

        x = DG.vmax(:);
        y = DG.avgDecel(:);

        ok = ~isnan(x) & ~isnan(y);
        x = x(ok); y = y(ok);

        scatter(x, y, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.35);

        % regression line
        if isfield(DG,'lm') && ~isempty(DG.lm)
            b0 = DG.lm.Coefficients.Estimate(1);
            b1 = DG.lm.Coefficients.Estimate(2);

            xx = linspace(min(x), max(x), 50);
            yy = b0 + b1*xx;
            plot(xx, yy, 'k-', 'LineWidth', 1.5);

            txt = sprintf('slope=%.3g\nR=%.2f\nN=%d', b1, DG.R, DG.nReaches);
        else
            txt = sprintf('R=%.2f\nN=%d', DG.R, DG.nReaches);
        end

        title(sprintf('%s | %s', an{a}, conds{ci}), 'Interpreter','none', 'FontSize', 9);
        xlabel('v_{max}');
        ylabel('avg decel');

        % small annotation
        xlim auto; ylim auto;
        xl = xlim; yl = ylim;
        text(xl(1)+0.02*range(xl), yl(2)-0.05*range(yl), txt, 'FontSize', 8, ...
            'VerticalAlignment','top');

        box off
    end
end

%% ============================================================
% (2) SPAGHETTI: SLOPE (GAIN) per animal across conditions
% ============================================================
figure('Color','w','Position',[150 150 900 320]); hold on
for a = 1:nA
    plot(1:nC, GainPlot(a,:), '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
end
xticks(1:nC); xticklabels(conds);
ylabel('Decel gain (slope)');
title('Decel gain per animal (spaghetti)');
yline(0,'--'); box off

%% ============================================================
% (3) SPAGHETTI: INTERCEPT per animal across conditions (if available)
% ============================================================
figure('Color','w','Position',[150 520 900 320]); hold on
for a = 1:nA
    plot(1:nC, IntPlot(a,:), '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
end
xticks(1:nC); xticklabels(conds);
ylabel('Intercept (b0)');
title('Intercept per animal (from decel ~ vmax)');
yline(0,'--'); box off



%% ============================================================
% (5) BOX + SWARM of reach_gains (pooled) + per-animal medians overlay
% ============================================================
% Build pooled vectors for boxchart
allRG = [];
allGrp = [];
allAnimal = [];

for ci = 1:nC
    for a = 1:nA
        DG = DGout{a,ci};
        if isempty(DG) || ~isfield(DG,'reach_gains') || isempty(DG.reach_gains), continue; end
        rg = DG.reach_gains;
        if size(rg,2) > 1, rg = rg(:,1); end
        rg = rg(:);
        rg = rg(~isnan(rg));
        if isempty(rg), continue; end

        allRG = [allRG; rg];
        allGrp = [allGrp; repmat(ci, numel(rg), 1)];
        allAnimal = [allAnimal; repmat(a, numel(rg), 1)];
    end
end

figure('Color','w','Position',[150 160 1050 420]); hold on
boxchart(allGrp, allRG);
xlim([0.5 nC+0.5]);
xticks(1:nC); xticklabels(conds);
ylabel('reach gain = decel / vmax');
title('Reach-gain distributions (pooled) + animal medians');

% overlay a light swarm (subsample if huge)
rng(0);
maxPts = 3000;
if numel(allRG) > maxPts
    idx = randperm(numel(allRG), maxPts);
else
    idx = 1:numel(allRG);
end
xj = allGrp(idx) + (rand(numel(idx),1)-0.5)*0.18;
scatter(xj, allRG(idx), 10, 'k', 'filled', 'MarkerFaceAlpha', 0.15);

% per-animal medians as connected points
MedRG_an = nan(nA,nC);
for a = 1:nA
    for ci = 1:nC
        DG = DGout{a,ci};
        if isempty(DG) || ~isfield(DG,'reach_gains') || isempty(DG.reach_gains), continue; end
        rg = DG.reach_gains; if size(rg,2)>1, rg=rg(:,1); end
        MedRG_an(a,ci) = median(rg,'omitnan');
    end
end
for a = 1:nA
    plot(1:nC, MedRG_an(a,:), '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
end
box off


%% OTHER ANALYSES


% Pick which two conditions to compare
BL = BL_CSM_big;         % baseline struct with .stimvel
WO = WOE_CSM_big;        % washout struct with .stimvel  (change to WO_CSM_big or WOL_CSM_big)

% Which dimension in stimvel is outward velocity / acceleration?
VEL_DIM = 5;             % you said outward velocity is col 7
ACC_DIM = 8;             % CHANGE THIS to your outward acceleration dim in stimvel

% Time vector (seconds) corresponding to stimvel time axis:
% If you already have it as 'time', comment this out and use yours.
time = BL.Average(:,1);  % or: time = your_time_vector(:);

% Windows (seconds) – adjust if needed
win_prepeak  = [-0.20 -0.08];
win_trough   = [-0.08  0.00];
win_rebound  = [ 0.00  0.12];
win_postosc  = [ 0.00  0.20];

% Cluster permutation settings
nPerm = 5000;
alpha = 0.05;

%% ========================= FIG 1: VELOCITY =========================
[muBL, semBL, YBL] = mean_sem_from_stimvel(BL.stimvel, VEL_DIM);
[muWO, semWO, YWO] = mean_sem_from_stimvel(WO.stimvel, VEL_DIM);

% Cluster permutation on paired per-animal difference (WO - BL)
sigMask_vel = cluster_perm_paired(YWO - YBL, nPerm, alpha);

figure('Color','w','Position',[200 200 900 380]); hold on
plot_shaded(time, muBL, semBL, [0.2 0.2 0.2]);   % BL gray/black
plot_shaded(time, muWO, semWO, [0.75 0.3 0.85]); % WO purple

xline(0,'--k'); yline(0,'--k');
xlabel('Time from CSM (s)');
ylabel('Outward Velocity');
title('Velocity: BL vs Washout');

% Significance bar
yl = ylim;
ySig = yl(1) + 0.05*range(yl);
plot_sig_bar(time, sigMask_vel, ySig);

box off

%% ========================= FIG 2: ACCELERATION =========================
[muBL_a, semBL_a, ABL] = mean_sem_from_stimvel(BL.stimvel, ACC_DIM);
[muWO_a, semWO_a, AWO] = mean_sem_from_stimvel(WO.stimvel, ACC_DIM);

sigMask_acc = cluster_perm_paired(AWO - ABL, nPerm, alpha);

figure('Color','w','Position',[200 620 900 380]); hold on
plot_shaded(time, muBL_a, semBL_a, [0.2 0.2 0.2]);
plot_shaded(time, muWO_a, semWO_a, [0.75 0.3 0.85]);

xline(0,'--k'); yline(0,'--k');
xlabel('Time from CSM (s)');
ylabel('Outward Acceleration');
title('Acceleration: BL vs Washout');

yl = ylim;
ySig = yl(1) + 0.05*range(yl);
plot_sig_bar(time, sigMask_acc, ySig);

box off
%%
fprintf('\nVEL sigMask: %d of %d timepoints TRUE\n', nnz(sigMask_vel), numel(sigMask_vel));
fprintf('ACC sigMask: %d of %d timepoints TRUE\n', nnz(sigMask_acc), numel(sigMask_acc));

% Quick check: biggest effect time
dVel = mean(YWO - YBL, 2, 'omitnan');
[~, imax] = max(abs(dVel));
fprintf('Max |mean diff| at t=%.3f s: diff=%.3g\n', time(imax), dVel(imax));


%% ========================= FIG 3: FEATURE METRICS =========================
M = compute_metrics(time, YBL, YWO, win_prepeak, win_trough, win_rebound, win_postosc);

metricNames = {'prePeak','trough','rebound','A_decel','A_rebound','t_rebound','postOsc'};
metricLabels = {'Pre-peak max','Trough min','Rebound max','Decel area','Rebound area','Rebound time','Post osc (std)'};

figure('Color','w','Position',[150 120 1100 520]);
tiledlayout(2,4,'TileSpacing','compact','Padding','compact');

for k = 1:numel(metricNames)
    nexttile; hold on
    xBLm = M.BL.(metricNames{k});
    xWOm = M.WO.(metricNames{k});

    % paired spaghetti
    for a = 1:numel(xBLm)
        plot([1 2],[xBLm(a) xWOm(a)],'-o','LineWidth',1.5,'MarkerSize',5);
    end

    % paired stats
    ok = ~isnan(xBLm) & ~isnan(xWOm);
    if sum(ok) >= 2
        [~,pT] = ttest(xWOm(ok), xBLm(ok));
        pSR = signrank(xWOm(ok), xBLm(ok));
        ttl = sprintf('%s\np(t)=%.3g  p(sr)=%.3g', metricLabels{k}, pT, pSR);
    else
        ttl = sprintf('%s\nn too small', metricLabels{k});
    end

    xlim([0.7 2.3]); xticks([1 2]); xticklabels({'BL','WO'});
    title(ttl);
    box off
end

%% Helper functions 

function [mu, sem, Y] = mean_sem_from_stimvel(stimvel, dim)
% stimvel: [nTime x nDim x nAnimals]
Y = squeeze(stimvel(:, dim, :)); % [nTime x nAnimals]
mu = mean(Y, 2, 'omitnan');
sem = std(Y, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(Y),2));
end

function plot_shaded(t, mu, sem, col)
% simple shaded error without toolboxes
fill([t(:); flipud(t(:))], [mu(:)+sem(:); flipud(mu(:)-sem(:))], col, ...
    'FaceAlpha',0.18,'EdgeColor','none');
plot(t, mu, 'Color', col, 'LineWidth', 3);
end

function plot_sig_bar(t, sigMask, y)
% draw thick horizontal segments where sigMask==true
sigMask = sigMask(:) > 0;
if ~any(sigMask), return; end
d = diff([false; sigMask; false]);
starts = find(d==1);
ends   = find(d==-1)-1;

for i = 1:numel(starts)
    plot([t(starts(i)) t(ends(i))], [y y], 'k-', 'LineWidth', 6);
end
end

function sigMask = cluster_perm_paired(D, nPerm, alpha)
% D: [nTime x nAnimals] paired difference (WO - BL)
% Cluster-based sign-flip permutation (two-sided) across animals

[tN, nA] = size(D);
sigMask = false(tN,1);

% observed t-stat vs 0 at each time
tObs = nan(tN,1);
for ti = 1:tN
    x = D(ti,:);
    ok = ~isnan(x);
    if sum(ok) >= 2
        tObs(ti) = mean(x(ok)) ./ (std(x(ok))/sqrt(sum(ok)));
    end
end

% threshold for forming clusters (uncorrected)
df = max(1, nA-1);
tThresh = tinv(1 - alpha/2, df);  % two-sided

% get observed cluster masses
obsMasses = cluster_masses(tObs, tThresh);

% permutation max cluster mass distribution
maxMass = zeros(nPerm,1);

for p = 1:nPerm
    flips = (rand(1,nA) > 0.5)*2 - 1;   % +/-1
    Dp = D .* flips;                    % sign flip by animal

    tP = nan(tN,1);
    for ti = 1:tN
        x = Dp(ti,:);
        ok = ~isnan(x);
        if sum(ok) >= 2
            tP(ti) = mean(x(ok)) ./ (std(x(ok))/sqrt(sum(ok)));
        end
    end

    masses = cluster_masses(tP, tThresh);
    if isempty(masses)
        maxMass(p) = 0;
    else
        maxMass(p) = max(masses);
    end
end

% cluster-mass significance threshold
crit = prctile(maxMass, 100*(1-alpha));

% mark significant observed clusters
if ~isempty(obsMasses)
    sigMask = mark_sig_clusters(tObs, tThresh, crit);
end
end

function masses = cluster_masses(tvec, tThresh)
% return absolute cluster mass = sum(|t|) for contiguous supra-threshold points
sup = abs(tvec) > tThresh;
masses = [];
if ~any(sup), return; end
d = diff([false; sup; false]);
starts = find(d==1);
ends   = find(d==-1)-1;
for i = 1:numel(starts)
    masses(end+1) = sum(abs(tvec(starts(i):ends(i))), 'omitnan'); %#ok<AGROW>
end
end

function sigMask = mark_sig_clusters(tvec, tThresh, crit)
sup = abs(tvec) > tThresh;
sigMask = false(size(tvec));
d = diff([false; sup; false]);
starts = find(d==1);
ends   = find(d==-1)-1;
for i = 1:numel(starts)
    mass = sum(abs(tvec(starts(i):ends(i))), 'omitnan');
    if mass >= crit
        sigMask(starts(i):ends(i)) = true;
    end
end
end

function M = compute_metrics(t, YBL, YWO, win_prepeak, win_trough, win_rebound, win_postosc)
% YBL/YWO: [nTime x nAnimals]
idx_pre = t>=win_prepeak(1) & t<=win_prepeak(2);
idx_tr  = t>=win_trough(1)  & t<=win_trough(2);
idx_re  = t>=win_rebound(1) & t<=win_rebound(2);
idx_po  = t>=win_postosc(1) & t<=win_postosc(2);

nA = size(YBL,2);

M.BL.prePeak  = nan(nA,1);  M.WO.prePeak  = nan(nA,1);
M.BL.trough   = nan(nA,1);  M.WO.trough   = nan(nA,1);
M.BL.rebound  = nan(nA,1);  M.WO.rebound  = nan(nA,1);
M.BL.A_decel  = nan(nA,1);  M.WO.A_decel  = nan(nA,1);
M.BL.A_rebound= nan(nA,1);  M.WO.A_rebound= nan(nA,1);
M.BL.t_rebound= nan(nA,1);  M.WO.t_rebound= nan(nA,1);
M.BL.postOsc  = nan(nA,1);  M.WO.postOsc  = nan(nA,1);

for a = 1:nA
    v1 = YBL(:,a); v2 = YWO(:,a);

    % pre-peak max
    M.BL.prePeak(a) = max(v1(idx_pre), [], 'omitnan');
    M.WO.prePeak(a) = max(v2(idx_pre), [], 'omitnan');

    % trough min
    M.BL.trough(a) = min(v1(idx_tr), [], 'omitnan');
    M.WO.trough(a) = min(v2(idx_tr), [], 'omitnan');

    % rebound max + time
    [M.BL.rebound(a), i1] = max(v1(idx_re), [], 'omitnan');
    [M.WO.rebound(a), i2] = max(v2(idx_re), [], 'omitnan');

    tre = t(idx_re);
    if ~isempty(tre)
        M.BL.t_rebound(a) = tre(i1);
        M.WO.t_rebound(a) = tre(i2);
    end

    % decel area: integrate negative velocity magnitude during trough window
    v1t = v1(idx_tr); v2t = v2(idx_tr);
    tt  = t(idx_tr);
    M.BL.A_decel(a) = trapz(tt, abs(v1t).* (v1t<0));
    M.WO.A_decel(a) = trapz(tt, abs(v2t).* (v2t<0));

    % rebound area: integrate positive velocity during rebound window
    v1r = v1(idx_re); v2r = v2(idx_re);
    tr  = t(idx_re);
    M.BL.A_rebound(a) = trapz(tr, v1r.* (v1r>0));
    M.WO.A_rebound(a) = trapz(tr, v2r.* (v2r>0));

    % post-oscillation variability
    M.BL.postOsc(a) = std(v1(idx_po), 'omitnan');
    M.WO.postOsc(a) = std(v2(idx_po), 'omitnan');
end
end
