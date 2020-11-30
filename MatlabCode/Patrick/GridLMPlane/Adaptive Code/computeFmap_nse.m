function [neglogev, hmapR, a, W, sqrtL]  = computeFmap_nse(datastruct)

% given H, find hmap
h0 = datastruct.hinit;
thresh = 0.05;
diffIn_h = 1;
maxIter = 10; 
count = 1;
while (diffIn_h>thresh)&&(count<=maxIter)
    
    [diffIn_h, hmapR, a, W, sqrtL, obj] = objFun_NewtonsMethod_nse(h0, datastruct);
    if count ==1
        datastruct.obj = obj;
    end
    %     diffInF
    % check if the obj is increasing
    if count >1
        if obj < datastruct.obj
            h0 = (hmapR+h0)/2; % if obj is decreasing, choose a smaller step
        else
            datastruct.obj = obj;
            h0 = hmapR;
        end
    end
    count = count+1;
    
    
end
    
% 3. given fmap, compute neglogev
neglogev = computeLogevidence_nse(hmapR, W, datastruct.muf, a, datastruct.r);