function response = getresponse(inpimage,STA,PC,subunit1,subunit2)
% Assuming that each filter does a linear spatial summation but combines
% the output of the filters linearly/ quadratically
max_FR = 100;
mode = 1;
if mode == 1
    netdrive = (inpimage(:)'*STA(:))^2 + (inpimage(:)'*PC(:))^2; % squaring and adding outputs of STA and PC filters
elseif mode == 2
    netdrive = (inpimage(:)'*STA(:))*(inpimage(:)'*PC(:)); % multiplying the outputs of the STA and PC filters
elseif mode == 3
    netdrive = (inpimage(:)'*subunit1(:))^2 + (inpimage(:)'*subunit2(:))^2; % squaring the drives from the 2 subunits
elseif mode == 4
    netdrive = (inpimage(:)'*subunit1(:)) + (inpimage(:)'*subunit2(:)); % linearly adding the drives from the 2 subunits
elseif mode == 5
    netdrive = (inpimage(:)'*STA(:)) + 1*(inpimage(:)'*PC(:))^2; % linear drive from STA filter + quadratic drive from PC filter
elseif mode ==6 
    netdrive = (inpimage(:)'*STA(:)); % only STA matters 
elseif mode == 7
    netdrive = max(0,inpimage(:)'*STA(:))^2 + 10*(inpimage(:)'*PC(:))^2; % rectified quadratui drive from STA filter + quadratic drive from PC filter
end
response = max_FR * 1./(1+exp(-1*(netdrive-1))); % Spike generation mechanism: sigmoidal non-linearity 

% response = poissrnd(response); % adding poisson variability
end

