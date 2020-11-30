%function calculate z score 

function [zscoresG] = zscorefunc(meanG, uniqueors)
mnsG = mean(meanG,2);
stdsG = std(meanG,0,2);
zscoresG = (meanG-repmat(mnsG,[1,length(uniqueors),1]))./repmat(stdsG,[1,length(uniqueors),1]);
end
