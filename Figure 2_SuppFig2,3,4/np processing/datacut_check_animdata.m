

cn = fieldnames(anim_data);
for c = [1:3]
    
    
    
figure; hold on

scatter3(anim_data.(cn{c}).Raw(:,1),anim_data.(cn{c}).Raw(:,2),anim_data.(cn{c}).Raw(:,3),30,'filled')


%Cut out noise (nonreach data per DJC). 
x = .3; y = .25; %these are the thresholds dylan uses (MIB 12/2021)
nohead = anim_data.(cn{c}).Raw(~(anim_data.(cn{c}).Raw(:,2)>x | anim_data.(cn{c}).Raw(:,3)>y),:);

scatter3(nohead(:,1),nohead(:,2),nohead(:,3),15,'filled');



end


