function [neglogev, fmapR, a, W, sqrtL]  = updateFmap(p, datastruct)

%% linear scale for hyperparameters:

% 1. given theta, form K
mufFinal = abs(p(1));
datastruct.muf = mufFinal;
param = abs(p(2:end));
norm_mat = datastruct.norm_mat; 
datastruct.K = param(1)*exp(-.5/param(2).*norm_mat);

% 2. given K, find fmap
f0 = datastruct.finit;
% options = optimset('Display','iter','TolFun',1e-6, 'algorithm',{'levenberg-marquardt',.005}, 'maxFunEvals', 1000, 'maxIter', 5);
% fmap = fsolve(@(f)objFun_NewtonsMethod(f, datastruct), f0, options);
% [differece_in_f, fmapR, a, W, sqrtL] = objFun_NewtonsMethod(fmap, datastruct);

thresh = 1e-3;
diffInF = 1;
maxIter = 10;
count = 1;
while (diffInF>thresh)&&(count<=maxIter)
    [diffInF, fmapR, a, W, sqrtL, obj] = objFun_NewtonsMethod(f0, datastruct);
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