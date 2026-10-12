function  accg_3d= get_3d_ACCG(mat,n,win)
time=mat(:,1);
spk=mat(:,2);
ix= find(spk==1);
fr=zeros(1,length(spk));
bin=1;
Tax=-win:win;
for i=    2:length(ix)-1
    T1=ix(i-1);
    T=ix(i);
    T2=ix(i+1);
    for t=round(T1+.5*(T-T1)):round(T+.5*(T2-T))
        fr(t)=(2*bin)/(T2-T1);
    end
end
sfr= movmean(fr, 250);
ifr=fr(ix);
all_index=repmat(ix,1,length(Tax))+repmat(Tax,length(ix),1);
tmp=find(all_index(:,1)>0 & all_index(:,end)<length(fr));
ix= ix(tmp); all_index=all_index(tmp,:);
deciles= quantile(sfr(ix),0:1/n:1);
all_fr_real=spk(all_index);
for d=1:length(deciles)-1
    index= find(sfr(ix)>=deciles(d) & sfr(ix)<=deciles(d+1));
    this_fr=all_fr_real(index,:);
    accg_3d(d,:)=nanmean(this_fr,1);
end
accg_3d(:,find(Tax==0))=0;
end