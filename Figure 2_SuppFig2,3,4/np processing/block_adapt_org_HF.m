function [R_Adapt,R_Adapt_inx,BL,ST,WA] = block_adapt_org_HF(Reaches)


R_Adapt = Reaches;

% %remove reaches that mess up the blocks (spurious stims or unstims)
% R_Adapt_inx = [1:31 33:100];
% R_Adapt = Extract(R_Adapt,R_Adapt_inx);

rn = fieldnames(R_Adapt);
stim_inx = [];
for r = 1:length(rn)
    stim_inx(r,1) = R_Adapt.(rn{r}).stim;
end

%check sum of abs(diff). if greater than 2, that means there were extra
%'blocks' or spurious stims, which need to be inspected and corrected
%manually above
if sum(abs(diff(stim_inx))) > 2
    warning('too many blocks');
    [];
end


stim_block_start  = find(stim_inx == 1,1);
wash_block_start = find(stim_inx(stim_block_start:end) == 0,1)+stim_block_start-1;

BL = Extract(R_Adapt,1:(stim_block_start-1));
ST = Extract(R_Adapt,stim_block_start:(wash_block_start-1));
WA = Extract(R_Adapt,wash_block_start:length(fieldnames(R_Adapt)));

if ~exist('R_Adapt_inx','var')
    R_Adapt_inx = 1:length(rn);
end

