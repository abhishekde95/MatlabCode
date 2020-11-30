% This script will run gabors for Rieke Lab Rotation Presentation

function []=presentation(rgb,stro,Stim)

contrasts(1,:) = rgb(15,:)*1.5
contrasts(2,:) = rgb(15,:)
contrasts(3,:) = rgb(2,:)

for s = 1:size(contrasts,1)
        [Stim.Power,Stim.Pixel.Gaussian] = stimDynamics(contrasts(s,:),stro,Stim); % Retuns a matrix at each frame refresh
        keyboard
end