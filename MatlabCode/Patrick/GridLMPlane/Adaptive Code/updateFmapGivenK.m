function [neglogev, fmapR, a, W, sqrtL]  = updateFmapGivenK(datastruct)

% Given K, find fmap
f0 = datastruct.finit;

thresh = 0.05;
diffInF = 1;
maxIter = 10;
count = 1;
while (diffInF>thresh)&&(count<=maxIter)
    [diffInF, fmapR, a, W, sqrtL, obj] = objFun_NewtonsMethod_GivenK(f0, datastruct);
    if count ==1
        datastruct.obj = obj;
    end
%     diffInF
    % check if the obj is increasing
    if count >1
        if obj < datastruct.obj
            f0 = (fmapR+f0)/2; % if obj is decreasing, choose a smaller step
        else
            datastruct.obj = obj;
            f0 = fmapR;
        end
    end
    count = count+1;
end

% 3. given fmap, compute neglogev
neglogev = computeLogevidence(fmapR, W, datastruct.muf, a, datastruct);