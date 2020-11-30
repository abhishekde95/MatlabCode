function [OCidx, LMidx, LUMidx, SOidx, hardtoclassifyidx] = classifycells(S1LMS,S2LMS)
% Classifying cells into Orange-cyan, lime-magenta DO, luminance cells, hardtoclassifycells
idx = find(sum(sign(S1LMS).*sign(S2LMS),1)==-3);
OCidx = idx(find(sign(S1LMS(2,idx)).*sign(S1LMS(3,idx))==1 & sign(S1LMS(1,idx)).*sign(S1LMS(3,idx))==-1));
LMidx = idx(find(sign(S1LMS(1,idx)).*sign(S1LMS(3,idx))==1 & sign(S1LMS(2,idx)).*sign(S1LMS(3,idx))==-1));


LUMidx = idx(find(sign(S1LMS(1,idx)).*sign(S1LMS(2,idx))==1));
hardtoclassifyidx = find(sum(sign(S1LMS).*sign(S2LMS),1)>-3);
moreLUMidxs = find(sum(sign(S1LMS(1:2,:)).*sign(S2LMS(1:2,:)))==-2 & sign(S1LMS(1,:)).*sign(S1LMS(2,:))==1 & sum(sign(S1LMS).*sign(S2LMS),1)>-3 & abs(S1LMS(3,:))+abs(S2LMS(3,:))<0.2*2);
LUMidx = [LUMidx moreLUMidxs];


SOidx = find(sign(S1LMS(1,:)).*sign(S2LMS(1,:))==1 & sign(S1LMS(2,:)).*sign(S2LMS(2,:))==1 & sign(S1LMS(3,:)).*sign(S2LMS(3,:))==1);
hardtoclassifyidx(ismember(hardtoclassifyidx,LUMidx)) = [];
hardtoclassifyidx(ismember(hardtoclassifyidx,SOidx)) = [];
end

