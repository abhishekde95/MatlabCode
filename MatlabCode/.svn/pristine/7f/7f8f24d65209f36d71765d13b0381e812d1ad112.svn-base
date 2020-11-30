function [oog]=isItOog(currentLMS,stimIntensity)
global M
global bkgndLMS

if nargin==1
    stimIntensity=1;
end

currentLMS=currentLMS*stimIntensity;
currentLMS=bkgndLMS.*(currentLMS+1);
rgb=currentLMS*inv(M);

if min(rgb)>0 && max(rgb)<1
    oog=0;
else
    oog=1;
end