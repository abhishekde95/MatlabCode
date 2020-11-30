function [flag, gamutViolation] = Check_ifcompletely_probed(stro,basisvec,bkgnd_monitor,idxs1,weights)
%Check if the a direction was succesfully/completely probed and if there were any gamut violations in particular directions
% flag = 0 - incompletely probed
% flag = 1 - completely probed
% gamutViolation = 1 - gun intensity exceeded the gamut
global reversalflagidx stepsizescale stepsize nreversals
% keyboard;
gamutViolation = 0;
flag = 0;
num_reversals = stro.trial(idxs1(end),reversalflagidx);
img = weights(1)*(basisvec{1}-bkgnd_monitor) + weights(2)*(basisvec{2}-bkgnd_monitor);
tmp_stepsize = stepsize*(stepsizescale^(num_reversals-1));
img_decrease = img*(1-tmp_stepsize) + bkgnd_monitor;
img_increase = img*(1+tmp_stepsize) + bkgnd_monitor;
if num_reversals == nreversals-1
    flag = 1; % completely probed
else
    if (any(img_increase(:)>1) | any(img_increase(:)<0)) | (any(img_decrease(:)>1) | any(img_decrease(:)<0))
        % checking if the next image would be satisfy out of gamut condition
        flag = 1;
        gamutViolation = 1;
    end
end
end

