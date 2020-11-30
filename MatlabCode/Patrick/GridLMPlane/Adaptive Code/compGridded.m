function [GInit, finit, minneglogev] = compGridded(InitialValues, datastruct)

HowManyInit = size(InitialValues, 1);
neglogev = zeros(HowManyInit,1);
fmapRmat = cell(HowManyInit,1);

for i=1:HowManyInit 
    prs = InitialValues(i,:)';
    [neglogev(i), fmapR]  = updateFmap(prs, datastruct);
    fmapRmat{i} = fmapR;
end

minneglogev =(min(neglogev)); % find a minimum neg log evid and its index
indx = (neglogev== minneglogev);
GInit = InitialValues(indx,:)'; % choose the best initial value for optimization
GInit = GInit(:,1); % choose first one if there are more than one minimium value
finit = fmapRmat{indx};